"""
Auto-calibration loop for flood models.
AutoResearch-style iterative calibration with lessons learned.
"""

import numpy as np
import json
import os
from typing import Dict, List, Any, Tuple
from datetime import datetime

from .metrics import CalibrationMetrics
from ..models.adapter import ModelAdapter


class CalibrationResult:
    """Container for calibration iteration results"""
    
    def __init__(self, iteration: int, params: Dict[str, float], 
                 metrics: Dict[str, float], score: float):
        self.iteration = iteration
        self.params = params.copy()
        self.metrics = metrics.copy()
        self.score = score
        self.timestamp = datetime.now().isoformat()


class AutoCalibrationLoop:
    """
    AutoResearch-style automatic calibration loop.
    
    Iteratively proposes parameters, runs simulations, and evaluates metrics.
    Maintains Pareto front and learns from past runs.
    """
    
    def __init__(self, 
                 model_adapter: ModelAdapter,
                 storm_event: Dict[str, Any],
                 target_metrics: Dict[str, float],
                 weights: Dict[str, float]):
        """
        Initialize calibration loop.
        
        Args:
            model_adapter: Model adapter (LISFLOOD-FP, etc.)
            storm_event: Storm event configuration
            target_metrics: Target values for NSE, KGE, IoU
            weights: Objective weights for composite score
        """
        self.model = model_adapter
        self.storm_event = storm_event
        self.target_metrics = target_metrics
        self.weights = weights
        
        self.metrics = CalibrationMetrics()
        
        # State
        self.iteration = 0
        self.best_score = -np.inf
        self.best_params = {}
        self.best_metrics = {}
        self.pareto_front: List[CalibrationResult] = []
        self.history: List[CalibrationResult] = []
        self.lessons: List[Dict[str, Any]] = []
        
        # Convergence tracking
        self.converged = False
        self.convergence_window = 10  # Check last N iterations
        self.tolerance = 0.001  # Score change threshold
    
    def propose_parameters(self) -> Dict[str, float]:
        """
        Propose new parameter set (agent-driven).
        
        Strategy evolves over iterations:
        - Iterations 0-9: Latin Hypercube sampling (exploration)
        - Iterations 10-29: Bayesian optimization (refinement)
        - Iterations 30+: Local search around best (exploitation)
        
        Returns:
            Dict of parameter names to values
        """
        param_bounds, param_names = self.model.get_parameter_bounds()
        
        if self.iteration < 10:
            # Initial exploration: Latin Hypercube sampling
            params = self._latin_hypercube_sample(param_bounds, param_names)
        
        elif self.iteration < 30:
            # Refinement: Around best (Gaussian perturbation)
            params = self._local_search(self.best_params, param_bounds, param_names)
        
        else:
            # Exploitation: Smaller perturbations
            params = self._local_search(self.best_params, param_bounds, param_names, scale=0.1)
        
        return params
    
    def evaluate(self, params: Dict[str, float]) -> Tuple[Dict[str, float], bool, str]:
        """
        Run simulation and compute metrics.
        
        Args:
            params: Parameter set to evaluate
        
        Returns:
            Tuple of (metrics_dict, success, error_message)
        """
        # 1. Validate parameters
        if not self.model.validate_parameters(params):
            metrics = {"nse": -np.inf, "kge": -np.inf, "iou": 0.0}
            return metrics, False, "Invalid parameters"
        
        # 2. Run simulation
        result = self.model.run_simulation(self.storm_event)
        
        if not result.success:
            metrics = {"nse": -np.inf, "kge": -np.inf, "iou": 0.0}
            return metrics, False, result.error
        
        # 3. Extract outputs
        discharge_obs = self.storm_event["observations"]["gauge"]["discharge"]
        discharge_sim = result.discharge["values"]
        inundation_obs = self.storm_event["observations"]["satellite"]["mask"]
        inundation_sim = result.inundation["grid"]
        
        # Align discharge time series (if lengths differ)
        if len(discharge_obs) != len(discharge_sim):
            min_len = min(len(discharge_obs), len(discharge_sim))
            discharge_obs = discharge_obs[:min_len]
            discharge_sim = discharge_sim[:min_len]
        
        # 4. Compute metrics
        metrics = {
            "nse": self.metrics.nash_sutcliffe(discharge_obs, discharge_sim),
            "kge": self.metrics.kling_gupta(discharge_obs, discharge_sim),
            "iou": self.metrics.iou(inundation_obs, inundation_sim),
        }
        
        # 5. Validate metrics
        is_valid, msg = self.metrics.validate_metrics(metrics)
        if not is_valid:
            return metrics, False, msg
        
        # 6. Composite score
        score = self.metrics.composite_score(metrics, self.weights)
        metrics["score"] = score
        
        return metrics, True, ""
    
    def update_pareto_front(self, result: CalibrationResult):
        """
        Update Pareto front with new result.
        
        A result is non-dominated if no other result is better in all objectives.
        """
        dominated = False
        
        # Check if dominated by any existing result
        for existing in self.pareto_front:
            if self._dominates(existing, result):
                dominated = True
                break
            elif self._dominates(result, existing):
                # Remove dominated result
                self.pareto_front.remove(existing)
        
        if not dominated:
            self.pareto_front.append(result)
    
    def _dominates(self, a: CalibrationResult, b: CalibrationResult) -> bool:
        """
        Check if solution a dominates solution b.
        
        a dominates b if:
        - a is better or equal in all objectives
        - a is strictly better in at least one objective
        """
        nse_better = a.metrics["nse"] >= b.metrics["nse"]
        kge_better = a.metrics["kge"] >= b.metrics["kge"]
        iou_better = a.metrics["iou"] >= b.metrics["iou"]
        
        strictly_better = (
            a.metrics["nse"] > b.metrics["nse"] or
            a.metrics["kge"] > b.metrics["kge"] or
            a.metrics["iou"] > b.metrics["iou"]
        )
        
        return nse_better and kge_better and iou_better and strictly_better
    
    def check_convergence(self) -> bool:
        """
        Check if calibration has converged.
        
        Convergence criteria:
        - Score improvement < tolerance over last N iterations
        - Target metrics achieved
        """
        if len(self.history) < self.convergence_window:
            return False
        
        # Check score improvement
        recent_scores = [r.score for r in self.history[-self.convergence_window:]]
        score_range = max(recent_scores) - min(recent_scores)
        
        if score_range < self.tolerance:
            return True
        
        # Check target metrics
        latest = self.history[-1]
        targets_met = all(
            latest.metrics[k] >= self.target_metrics[k]
            for k in self.target_metrics
        )
        
        return targets_met
    
    def run_iteration(self) -> CalibrationResult:
        """
        Run single calibration iteration.
        
        Returns:
            CalibrationResult with metrics
        """
        # 1. Propose parameters
        params = self.propose_parameters()
        
        # 2. Evaluate
        metrics, success, error = self.evaluate(params)
        
        # 3. Create result
        score = metrics.get("score", -np.inf)
        result = CalibrationResult(self.iteration, params, metrics, score)
        
        # 4. Update state
        self.history.append(result)
        
        if success and score > self.best_score:
            self.best_score = score
            self.best_params = params.copy()
            self.best_metrics = metrics.copy()
        
        # 5. Update Pareto front
        if success:
            self.update_pareto_front(result)
        
        # 6. Record lessons
        if self.iteration > 0:
            self._record_lessons(result)
        
        # 7. Check convergence
        self.converged = self.check_convergence()
        
        # 8. Increment iteration
        self.iteration += 1
        
        return result
    
    def _record_lessons(self, result: CalibrationResult):
        """
        Record lessons from this iteration.
        
        Lessons include:
        - Performance drops (why did it get worse?)
        - Numerical instabilities (invalid parameters)
        - Parameter correlations (patterns in good solutions)
        """
        previous = self.history[-2]
        
        # Lesson: Performance dropped
        if result.score < previous.score - 0.05:  # Significant drop
            diff_params = {
                k: result.params[k] - previous.params[k]
                for k in result.params
                if abs(result.params[k] - previous.params[k]) > 0.01
            }
            
            self.lessons.append({
                "type": "performance_drop",
                "iteration": self.iteration,
                "previous_score": previous.score,
                "new_score": result.score,
                "param_changes": diff_params,
                "suggestion": f"Avoid increasing {[k for k in diff_params.values() if diff_params[k] > 0]}"
            })
        
        # Lesson: Metric-specific trade-offs
        if (result.metrics["nse"] > previous.metrics["nse"] and 
            result.metrics["iou"] < previous.metrics["iou"]):
            self.lessons.append({
                "type": "metric_tradeoff",
                "iteration": self.iteration,
                "tradeoff": "nse vs iou",
                "suggestion": "NSE improved but IoU worsened - check spatial accuracy"
            })
    
    def run(self, max_iterations: int = 100, 
            verbose: bool = True) -> Dict[str, Any]:
        """
        Run full calibration loop.
        
        Args:
            max_iterations: Maximum number of iterations
            verbose: Print progress
        
        Returns:
            Dict with calibration results
        """
        print(f"\n🌊 Starting automatic calibration (max {max_iterations} iterations)...")
        print(f"   Target: NSE≥{self.target_metrics.get('nse', 0.85)}, "
              f"KGE≥{self.target_metrics.get('kge', 0.80)}, "
              f"IoU≥{self.target_metrics.get('iou', 0.70)}")
        
        while self.iteration < max_iterations and not self.converged:
            result = self.run_iteration()
            
            if verbose:
                print(f"Iter {self.iteration:3d}: "
                      f"NSE={result.metrics['nse']:.3f}, "
                      f"KGE={result.metrics['kge']:.3f}, "
                      f"IoU={result.metrics['iou']:.3f}, "
                      f"Score={result.score:.3f}")
        
        # Summary
        print(f"\n✅ Calibration complete after {self.iteration} iterations")
        print(f"   Best: NSE={self.best_metrics['nse']:.3f}, "
              f"KGE={self.best_metrics['kge']:.3f}, "
              f"IoU={self.best_metrics['iou']:.3f}, "
              f"Score={self.best_score:.3f}")
        print(f"   Pareto front: {len(self.pareto_front)} non-dominated solutions")
        print(f"   Lessons learned: {len(self.lessons)}")
        
        return {
            "best_params": self.best_params,
            "best_metrics": self.best_metrics,
            "best_score": self.best_score,
            "pareto_front": [self._serialize_result(r) for r in self.pareto_front],
            "history": [self._serialize_result(r) for r in self.history],
            "lessons": self.lessons,
            "converged": self.converged,
            "iterations": self.iteration
        }
    
    def _latin_hypercube_sample(self, bounds: List[Tuple[float, float]], 
                                  names: List[str]) -> Dict[str, float]:
        """Latin Hypercube sampling for initial exploration"""
        np.random.seed(self.iteration * 42)  # Reproducible
        
        n_params = len(bounds)
        samples = np.zeros(n_params)
        
        for i in range(n_params):
            # LHS: each dimension is stratified into n_params intervals
            perm = np.random.permutation(n_params)
            sample = (perm[i] + np.random.random()) / n_params
            samples[i] = bounds[i][0] + sample * (bounds[i][1] - bounds[i][0])
        
        return {name: val for name, val in zip(names, samples)}
    
    def _local_search(self, best_params: Dict[str, float],
                     bounds: List[Tuple[float, float]],
                     names: List[str],
                     scale: float = 0.2) -> Dict[str, float]:
        """Local search around best parameters with Gaussian perturbation"""
        np.random.seed(self.iteration * 42)
        
        new_params = {}
        for i, name in enumerate(names):
            center = best_params.get(name, (bounds[i][0] + bounds[i][1]) / 2)
            range_size = (bounds[i][1] - bounds[i][0]) * scale
            
            perturbation = np.random.normal(0, range_size)
            new_val = center + perturbation
            
            # Clamp to bounds
            new_val = max(bounds[i][0], min(bounds[i][1], new_val))
            new_params[name] = new_val
        
        return new_params
    
    def _serialize_result(self, result: CalibrationResult) -> Dict[str, Any]:
        """Convert CalibrationResult to JSON-serializable dict"""
        return {
            "iteration": result.iteration,
            "params": result.params,
            "metrics": result.metrics,
            "score": result.score,
            "timestamp": result.timestamp
        }
    
    def save_results(self, output_path: str):
        """Save calibration results to JSON"""
        results = self.run(max_iterations=0)  # Get current state
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\n📁 Calibration results saved to: {output_path}")
