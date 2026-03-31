"""
Tests for calibration loop (parameter proposal, Pareto front, convergence).
"""

import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from calibration.loop import AutoCalibrationLoop, CalibrationResult
from models.adapter import ModelAdapter, SimulationResult


class FakeAdapter(ModelAdapter):
    """Fake model adapter for testing"""

    def __init__(self):
        super().__init__({})
        self.param_ranges = {
            "friction": {"value": 0.035, "min": 0.01, "max": 0.1},
            "theta": {"value": 0.7, "min": 0.5, "max": 1.0},
        }
        self.current_params = {}
        self._call_count = 0

    def get_parameters(self):
        return self.param_ranges.copy()

    def set_parameters(self, params):
        self.current_params = params.copy()

    def run_simulation(self, storm_event):
        self._call_count += 1
        # Return simple synthetic output
        return SimulationResult(
            discharge={
                "timestamps": np.arange(10, dtype=float),
                "values": np.random.random(10) + 1.0
            },
            inundation={"grid": np.random.random((5, 5))},
            success=True
        )

    def get_outputs(self):
        return {
            "discharge": {"timestamps": np.arange(10), "values": np.ones(10)},
            "inundation": {"grid": np.ones((5, 5))}
        }


def make_loop():
    adapter = FakeAdapter()
    storm = {
        "observations": {
            "gauge": {"discharge": np.ones(10)},
            "satellite": {"mask": np.ones((5, 5))}
        }
    }
    targets = {"nse": 0.5, "kge": 0.5, "iou": 0.5}
    weights = {"nse": 0.4, "kge": 0.3, "iou": 0.3}

    return AutoCalibrationLoop(adapter, storm, targets, weights)


class TestParameterProposal:
    """Tests for parameter proposal strategies"""

    def test_latin_hypercube_within_bounds(self):
        loop = make_loop()
        loop.iteration = 0
        params = loop.propose_parameters()

        adapter = FakeAdapter()
        bounds = adapter.get_parameter_bounds()

        for name, val in params.items():
            param_info = adapter.get_parameters()[name]
            assert param_info["min"] <= val <= param_info["max"], \
                f"{name}={val} out of range [{param_info['min']}, {param_info['max']}]"

    def test_local_search_within_bounds(self):
        loop = make_loop()
        loop.iteration = 20
        loop.best_params = {"friction": 0.05, "theta": 0.75}

        params = loop.propose_parameters()

        adapter = FakeAdapter()
        for name, val in params.items():
            param_info = adapter.get_parameters()[name]
            assert param_info["min"] <= val <= param_info["max"]

    def test_lhs_diversity_across_iterations(self):
        """Different iterations should produce different parameters"""
        loop = make_loop()
        all_params = []

        for i in range(5):
            loop.iteration = i
            all_params.append(loop.propose_parameters())

        # Check that not all parameter sets are identical
        friction_values = [p["friction"] for p in all_params]
        assert len(set(round(v, 6) for v in friction_values)) > 1


class TestParetoFront:
    """Tests for Pareto front maintenance"""

    def test_non_dominated_added(self):
        loop = make_loop()

        # First result always added
        r1 = CalibrationResult(0, {}, {"nse": 0.8, "kge": 0.7, "iou": 0.6}, 0.7)
        loop.update_pareto_front(r1)
        assert len(loop.pareto_front) == 1

    def test_dominating_replaces(self):
        loop = make_loop()

        r1 = CalibrationResult(0, {}, {"nse": 0.8, "kge": 0.7, "iou": 0.6}, 0.7)
        loop.update_pareto_front(r1)

        # r2 dominates r1 (better in all metrics)
        r2 = CalibrationResult(1, {}, {"nse": 0.9, "kge": 0.8, "iou": 0.7}, 0.8)
        loop.update_pareto_front(r2)

        assert len(loop.pareto_front) == 1
        assert loop.pareto_front[0].iteration == 1

    def test_non_comparable_both_kept(self):
        loop = make_loop()

        r1 = CalibrationResult(0, {}, {"nse": 0.9, "kge": 0.6, "iou": 0.6}, 0.7)
        loop.update_pareto_front(r1)

        # r2 better in KGE/IoU but worse in NSE → non-dominated
        r2 = CalibrationResult(1, {}, {"nse": 0.7, "kge": 0.9, "iou": 0.9}, 0.83)
        loop.update_pareto_front(r2)

        assert len(loop.pareto_front) == 2

    def test_dominated_not_added(self):
        loop = make_loop()

        r1 = CalibrationResult(0, {}, {"nse": 0.9, "kge": 0.8, "iou": 0.7}, 0.8)
        loop.update_pareto_front(r1)

        # r2 dominated by r1 (worse in all metrics)
        r2 = CalibrationResult(1, {}, {"nse": 0.5, "kge": 0.4, "iou": 0.3}, 0.4)
        loop.update_pareto_front(r2)

        assert len(loop.pareto_front) == 1
        assert loop.pareto_front[0].iteration == 0


class TestConvergence:
    """Tests for convergence detection"""

    def test_not_converged_early(self):
        loop = make_loop()
        loop.convergence_window = 5
        # Add only 3 results
        for i in range(3):
            loop.history.append(CalibrationResult(i, {}, {"nse": 0.5, "kge": 0.5, "iou": 0.5}, 0.5))
        assert loop.check_convergence() is False

    def test_converged_on_targets_met(self):
        loop = make_loop()
        loop.convergence_window = 3
        loop.target_metrics = {"nse": 0.5, "kge": 0.5, "iou": 0.5}

        # Fill history past the window
        for i in range(5):
            loop.history.append(CalibrationResult(i, {}, {"nse": 0.6, "kge": 0.6, "iou": 0.6}, 0.6))

        assert loop.check_convergence() is True

    def test_not_converged_targets_not_met(self):
        loop = make_loop()
        loop.convergence_window = 3
        loop.target_metrics = {"nse": 0.9, "kge": 0.9, "iou": 0.9}

        for i in range(5):
            loop.history.append(CalibrationResult(i, {}, {"nse": 0.5, "kge": 0.5, "iou": 0.5}, 0.5))

        # Targets not met and not enough history for plateau check (2x window)
        assert loop.check_convergence() is False


class TestEvaluateSetsParams:
    """Tests that evaluate() applies parameters to model (issue #1 fix)"""

    def test_params_set_before_simulation(self):
        adapter = FakeAdapter()
        storm = {
            "observations": {
                "gauge": {"discharge": np.ones(10)},
                "satellite": {"mask": np.ones((5, 5))}
            }
        }
        loop = AutoCalibrationLoop(
            adapter, storm,
            {"nse": 0.5, "kge": 0.5, "iou": 0.5},
            {"nse": 0.4, "kge": 0.3, "iou": 0.3}
        )

        test_params = {"friction": 0.05, "theta": 0.8}
        loop.evaluate(test_params)

        # Adapter should have the proposed parameters
        assert adapter.current_params == test_params
