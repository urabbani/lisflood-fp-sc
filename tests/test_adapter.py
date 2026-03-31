"""
Tests for model adapter base class and parameter validation.
"""

import numpy as np
import pytest
from models.adapter import ModelAdapter, SimulationResult


class ConcreteAdapter(ModelAdapter):
    """Concrete implementation for testing"""

    def __init__(self):
        super().__init__({})
        self.param_ranges = {
            "friction": {"value": 0.035, "min": 0.01, "max": 0.1},
            "theta": {"value": 0.7, "min": 0.5, "max": 1.0},
        }

    def get_parameters(self):
        return self.param_ranges.copy()

    def set_parameters(self, params):
        self.current_params = params.copy()

    def run_simulation(self, storm_event):
        return SimulationResult(success=True)

    def get_outputs(self):
        return {}


class TestSimulationResult:
    """Tests for SimulationResult defaults (issue #21 fix)"""

    def test_none_defaults_to_empty_dict(self):
        result = SimulationResult()
        assert result.discharge == {}
        assert result.inundation == {}

    def test_explicit_none_same_as_default(self):
        result = SimulationResult(discharge=None, inundation=None)
        assert result.discharge == {}
        assert result.inundation == {}

    def test_passed_values_preserved(self):
        d = {"values": np.array([1, 2, 3])}
        result = SimulationResult(discharge=d)
        assert result.discharge is d


class TestValidateParameters:
    """Tests for parameter validation (issue #17 fix)"""

    def setup_method(self):
        self.adapter = ConcreteAdapter()

    def test_valid_params(self):
        assert self.adapter.validate_parameters({"friction": 0.05, "theta": 0.8}) is True

    def test_unknown_param_rejected(self):
        # Issue #17: unknown params should cause validation failure
        assert self.adapter.validate_parameters({"nonexistent": 0.5}) is False

    def test_out_of_range_rejected(self):
        assert self.adapter.validate_parameters({"friction": 0.5}) is False

    def test_boundary_values_accepted(self):
        assert self.adapter.validate_parameters({"friction": 0.01}) is True
        assert self.adapter.validate_parameters({"friction": 0.1}) is True

    def test_empty_params_accepted(self):
        assert self.adapter.validate_parameters({}) is True


class TestParameterBounds:
    """Tests for get_parameter_bounds"""

    def setup_method(self):
        self.adapter = ConcreteAdapter()

    def test_returns_sorted_bounds(self):
        bounds, names = self.adapter.get_parameter_bounds()
        assert len(bounds) == 2
        assert len(names) == 2
        # Sorted alphabetically
        assert names == ["friction", "theta"]
        assert bounds[0] == (0.01, 0.1)
        assert bounds[1] == (0.5, 1.0)
