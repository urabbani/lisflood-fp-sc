"""
Base adapter interface for flood model integration.
All model adapters must inherit from this class.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any
import numpy as np

class SimulationResult:
    """Container for model simulation results"""
    
    def __init__(self, discharge: Dict[str, np.ndarray] = None, 
                 inundation: Dict[str, Any] = None,
                 success: bool = True,
                 error: str = None):
        self.discharge = discharge or {}  # {"timestamps": [], "values": []}
        self.inundation = inundation or {}  # {"grid": 2D_array, "resolution": ..., "crs": ...}
        self.success = success
        self.error = error


class ModelAdapter(ABC):
    """Base class for all flood model adapters"""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.param_ranges = {}
        self.current_params = {}
    
    @abstractmethod
    def get_parameters(self) -> Dict[str, Dict[str, float]]:
        """
        Return current parameters with ranges.
        
        Returns:
            Dict mapping parameter names to {"value": float, "min": float, "max": float}
        """
        pass
    
    @abstractmethod
    def set_parameters(self, params: Dict[str, float]):
        """
        Update model parameters.
        
        Args:
            params: Dict mapping parameter names to values
        """
        pass
    
    @abstractmethod
    def run_simulation(self, storm_event: Dict[str, Any]) -> SimulationResult:
        """
        Run model simulation and return outputs.
        
        Args:
            storm_event: Dict with start_date, end_date, rainfall_data, etc.
        
        Returns:
            SimulationResult with discharge time series and inundation grid
        """
        pass
    
    @abstractmethod
    def get_outputs(self) -> Dict[str, Any]:
        """
        Get model outputs for calibration.
        
        Returns:
            Dict with "discharge" and "inundation" data
        """
        pass
    
    def validate_parameters(self, params: Dict[str, float]) -> bool:
        """
        Check if parameters are within valid ranges.
        
        Args:
            params: Dict mapping parameter names to values
        
        Returns:
            True if valid, False otherwise
        """
        param_ranges = self.get_parameters()
        for name, value in params.items():
            if name not in param_ranges:
                print(f"Warning: Unknown parameter {name}")
                continue
            
            min_val = param_ranges[name]["min"]
            max_val = param_ranges[name]["max"]
            
            if not (min_val <= value <= max_val):
                print(f"Warning: {name} = {value} out of range [{min_val}, {max_val}]")
                return False
        
        return True
    
    def get_parameter_bounds(self) -> tuple:
        """
        Get parameter bounds for optimization algorithms.
        
        Returns:
            Tuple of (bounds_list, param_names)
            bounds_list: List of (min, max) tuples
            param_names: List of parameter names
        """
        params = self.get_parameters()
        bounds = []
        names = []
        
        for name in sorted(params.keys()):
            bounds.append((params[name]["min"], params[name]["max"]))
            names.append(name)
        
        return bounds, names
