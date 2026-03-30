"""
Calibration metrics for flood models.
Computes NSE, KGE (temporal) and IoU (spatial) metrics.
"""

import numpy as np
from typing import Dict, Tuple
from scipy import stats


class CalibrationMetrics:
    """Compute multi-objective calibration metrics"""
    
    @staticmethod
    def nash_sutcliffe(observed: np.ndarray, simulated: np.ndarray) -> float:
        """
        Nash-Sutcliffe Efficiency (NSE).
        
        Range: -∞ to 1, higher is better.
        NSE = 1 - (sum((obs - sim)^2) / sum((obs - mean(obs))^2))
        
        Args:
            observed: Observed discharge time series
            simulated: Simulated discharge time series
        
        Returns:
            NSE value
        """
        # Remove NaN values
        mask = ~np.isnan(observed) & ~np.isnan(simulated)
        obs = observed[mask]
        sim = simulated[mask]
        
        if len(obs) == 0:
            return -np.inf  # No valid data
        
        obs_mean = np.mean(obs)
        
        numerator = np.sum((obs - sim) ** 2)
        denominator = np.sum((obs - obs_mean) ** 2)
        
        if denominator == 0:
            return 1.0 if numerator == 0 else -np.inf
        
        nse = 1 - (numerator / denominator)
        return nse
    
    @staticmethod
    def kling_gupta(observed: np.ndarray, simulated: np.ndarray) -> float:
        """
        Kling-Gupta Efficiency (KGE).
        
        Range: -∞ to 1, higher is better.
        KGE = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
        
        Where:
        - r: Correlation coefficient
        - α: Variability ratio (std(sim) / std(obs))
        - β: Bias ratio (mean(sim) / mean(obs))
        
        Args:
            observed: Observed discharge time series
            simulated: Simulated discharge time series
        
        Returns:
            KGE value
        """
        # Remove NaN values
        mask = ~np.isnan(observed) & ~np.isnan(simulated)
        obs = observed[mask]
        sim = simulated[mask]
        
        if len(obs) == 0:
            return -np.inf
        
        # Compute components
        r = np.corrcoef(obs, sim)[0, 1]
        
        if np.isnan(r):
            return -np.inf
        
        α = np.std(sim) / np.std(obs) if np.std(obs) > 0 else 1.0
        β = np.mean(sim) / np.mean(obs) if np.mean(obs) > 0 else 1.0
        
        kge = 1 - np.sqrt((r - 1)**2 + (α - 1)**2 + (β - 1)**2)
        return kge
    
    @staticmethod
    def iou(observed_mask: np.ndarray, 
            simulated_mask: np.ndarray,
            depth_threshold: float = 0.1) -> float:
        """
        Intersection over Union (IoU) for flood extent.
        
        Range: 0 to 1, higher is better.
        IoU = |A ∩ B| / |A ∪ B|
        
        Args:
            observed_mask: Observed flood depth grid (meters) or binary mask
            simulated_mask: Simulated flood depth grid (meters) or binary mask
            depth_threshold: Threshold to consider a cell flooded (meters)
        
        Returns:
            IoU value
        """
        # Convert depth grids to binary masks
        if observed_mask.dtype in [np.float32, np.float64]:
            obs_binary = (observed_mask > depth_threshold).astype(np.uint8)
        else:
            obs_binary = observed_mask.astype(np.uint8)
        
        if simulated_mask.dtype in [np.float32, np.float64]:
            sim_binary = (simulated_mask > depth_threshold).astype(np.uint8)
        else:
            sim_binary = simulated_mask.astype(np.uint8)
        
        # Ensure same shape
        if obs_binary.shape != sim_binary.shape:
            print(f"Warning: Mask shapes differ: {obs_binary.shape} vs {sim_binary.shape}")
            # Resize to smaller shape (or use intersection of valid areas)
            min_rows = min(obs_binary.shape[0], sim_binary.shape[0])
            min_cols = min(obs_binary.shape[1], sim_binary.shape[1])
            obs_binary = obs_binary[:min_rows, :min_cols]
            sim_binary = sim_binary[:min_rows, :min_cols]
        
        # Compute IoU
        intersection = np.logical_and(obs_binary, sim_binary).sum()
        union = np.logical_or(obs_binary, sim_binary).sum()
        
        if union == 0:
            return 0.0  # Both masks are empty
        
        iou = intersection / union
        return iou
    
    @staticmethod
    def f1_score(observed_mask: np.ndarray,
                 simulated_mask: np.ndarray,
                 depth_threshold: float = 0.1) -> float:
        """
        F1 score for flood extent binary classification.
        
        F1 = 2 * (precision * recall) / (precision + recall)
        
        Args:
            observed_mask: Observed flood depth grid or binary mask
            simulated_mask: Simulated flood depth grid or binary mask
            depth_threshold: Threshold to consider a cell flooded (meters)
        
        Returns:
            F1 score
        """
        # Convert to binary masks
        if observed_mask.dtype in [np.float32, np.float64]:
            obs_binary = (observed_mask > depth_threshold).astype(np.uint8)
        else:
            obs_binary = observed_mask.astype(np.uint8)
        
        if simulated_mask.dtype in [np.float32, np.float64]:
            sim_binary = (simulated_mask > depth_threshold).astype(np.uint8)
        else:
            sim_binary = simulated_mask.astype(np.uint8)
        
        # Ensure same shape
        if obs_binary.shape != sim_binary.shape:
            min_rows = min(obs_binary.shape[0], sim_binary.shape[0])
            min_cols = min(obs_binary.shape[1], sim_binary.shape[1])
            obs_binary = obs_binary[:min_rows, :min_cols]
            sim_binary = sim_binary[:min_rows, :min_cols]
        
        # Compute TP, FP, FN
        tp = np.logical_and(obs_binary == 1, sim_binary == 1).sum()
        fp = np.logical_and(obs_binary == 0, sim_binary == 1).sum()
        fn = np.logical_and(obs_binary == 1, sim_binary == 0).sum()
        
        # Compute precision and recall
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        
        # Compute F1 score
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        return f1
    
    @staticmethod
    def peak_flow_error(observed: np.ndarray, 
                       simulated: np.ndarray) -> float:
        """
        Absolute error in peak flow.
        
        Args:
            observed: Observed discharge time series
            simulated: Simulated discharge time series
        
        Returns:
            Absolute peak flow error (m³/s)
        """
        obs_max = np.max(observed[~np.isnan(observed)])
        sim_max = np.max(simulated[~np.isnan(simulated)])
        
        return abs(obs_max - sim_max)
    
    @staticmethod
    def timing_error(observed: np.ndarray,
                    simulated: np.ndarray,
                    timestamps: np.ndarray) -> float:
        """
        Timing error in peak flow (hours).
        
        Args:
            observed: Observed discharge time series
            simulated: Simulated discharge time series
            timestamps: Timestamps (hours or seconds)
        
        Returns:
            Timing error (peak lag in same unit as timestamps)
        """
        obs_max_idx = np.argmax(observed[~np.isnan(observed)])
        sim_max_idx = np.argmax(simulated[~np.isnan(simulated)])
        
        # Convert to actual timestamps
        time_diff = abs(timestamps[sim_max_idx] - timestamps[obs_max_idx])
        
        return time_diff
    
    @staticmethod
    def composite_score(metrics: Dict[str, float], 
                      weights: Dict[str, float]) -> float:
        """
        Weighted composite calibration score.
        
        Normalizes each metric to [0, 1] range (higher is better).
        
        Args:
            metrics: Dict with "nse", "kge", "iou" (or other metrics)
            weights: Dict with same keys as metrics
        
        Returns:
            Composite score (0 to 1)
        """
        # Normalize metrics to [0, 1]
        normalized = {}
        
        for key in metrics:
            if key == "nse":
                normalized[key] = max(0.0, metrics[key])  # NSE: -∞ to 1
            elif key == "kge":
                normalized[key] = max(0.0, metrics[key])  # KGE: -∞ to 1
            elif key == "iou":
                normalized[key] = metrics[key]  # IoU: 0 to 1
            elif key == "f1":
                normalized[key] = metrics[key]  # F1: 0 to 1
            else:
                # Assume already normalized or handle specially
                normalized[key] = metrics[key]
        
        # Compute weighted sum
        score = sum(weights.get(k, 0) * normalized[k] for k in normalized)
        
        return score
    
    @staticmethod
    def validate_metrics(metrics: Dict[str, float]) -> Tuple[bool, str]:
        """
        Validate calibration metrics.
        
        Args:
            metrics: Dict with calibration metrics
        
        Returns:
            Tuple of (is_valid, message)
        """
        # Check for NaN or infinite values
        for key, value in metrics.items():
            if np.isnan(value) or np.isinf(value):
                return False, f"Invalid {key}: {value}"
        
        # Check metric ranges
        if "nse" in metrics and metrics["nse"] > 1.0:
            return False, f"NSE > 1.0: {metrics['nse']}"
        
        if "kge" in metrics and metrics["kge"] > 1.0:
            return False, f"KGE > 1.0: {metrics['kge']}"
        
        if "iou" in metrics and (metrics["iou"] < 0.0 or metrics["iou"] > 1.0):
            return False, f"IoU out of [0,1]: {metrics['iou']}"
        
        return True, "Metrics are valid"
