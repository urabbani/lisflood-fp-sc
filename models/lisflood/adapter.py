"""
LISFLOOD-FP model adapter.
Parses .par files, runs simulations, extracts outputs.
"""

import os
import re
import shutil
import subprocess
import logging
from pathlib import Path
import numpy as np
import rasterio
from rasterio.transform import from_origin
from typing import Dict, Any
from models.adapter import ModelAdapter, SimulationResult

logger = logging.getLogger(__name__)


class LISFLOODAdapter(ModelAdapter):
    """Adapter for LISFLOOD-FP hydrodynamic model"""

    @staticmethod
    def _resolve_executable(executable: str) -> str:
        """Resolve 'auto' to the vendored LISFLOOD-FP binary."""
        if executable != "auto":
            return executable

        from models.lisflood82.adapter import LISFLOOD82Adapter
        return LISFLOOD82Adapter._resolve_executable("auto")
    
    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)

        self.executable = LISFLOODAdapter._resolve_executable(config["executable"])
        self.param_template = config["param_template"]
        self.output_dir = config["output_dir"]
        
        # LISFLOOD-FP parameters (with physical ranges)
        self.param_ranges = {
            "friction": {"value": 0.035, "min": 0.01, "max": 0.1},
            "infiltration": {"value": 0.002, "min": 0.0, "max": 0.01},
            "manning_n": {"value": 0.035, "min": 0.02, "max": 0.06},
            "evaporation": {"value": 0.0, "min": 0.0, "max": 0.001},
            "theta": {"value": 0.7, "min": 0.5, "max": 1.0},  # Diffusion coefficient
        }
        
        # Load template .par file
        self.template_lines = self._load_template()
    
    def _load_template(self) -> list:
        """Load template .par file lines"""
        with open(self.param_template, 'r') as f:
            return f.readlines()
    
    def get_parameters(self) -> Dict[str, Dict[str, float]]:
        """Return LISFLOOD-FP parameters with ranges"""
        return self.param_ranges.copy()
    
    def set_parameters(self, params: Dict[str, float]):
        """
        Update LISFLOOD-FP parameters in .par file.
        
        Parameters:
        - friction: Channel/friction coefficient (Manning's n)
        - infiltration: Infiltration rate (m/s)
        - manning_n: Alternative Manning's n parameter
        - evaporation: Evaporation rate (m/s)
        - theta: Diffusion coefficient (0.5-1.0, 1.0 = full diffusion)
        """
        self.current_params = params.copy()
        
        # Create new .par file with updated parameters
        new_lines = []
        for line in self.template_lines:
            updated_line = line
            for param_name, param_value in params.items():
                # Pattern: PARAM_NAME: VALUE
                pattern = re.compile(rf'^{param_name.upper()}:\s*\d+\.?\d*', re.IGNORECASE)
                if pattern.search(line):
                    updated_line = pattern.sub(f'{param_name.upper()}: {param_value:.6f}', line)
            new_lines.append(updated_line)
        
        # Write to output directory
        os.makedirs(self.output_dir, exist_ok=True)
        output_par = os.path.join(self.output_dir, "sim.par")
        with open(output_par, 'w') as f:
            f.writelines(new_lines)
        
        return output_par
    
    def run_simulation(self, storm_event: Dict[str, Any]) -> SimulationResult:
        """
        Run LISFLOOD-FP simulation.
        
        Args:
            storm_event: Dict with:
                - start_date: str (ISO format)
                - end_date: str (ISO format)
                - rainfall_data: str (path to rainfall file)
                - boundary_discharge: str (path to boundary condition file)
        
        Returns:
            SimulationResult with discharge and inundation outputs
        """
        try:
            # 1. Prepare simulation
            par_file = self.set_parameters(self.current_params)
            
            # 2. Create input file (if not using .par file directly)
            input_file = self._create_input_file(storm_event, par_file)
            
            # 3. Run LISFLOOD-FP
            cmd = [self.executable, input_file]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                return SimulationResult(
                    success=False,
                    error=f"LISFLOOD-FP failed: {result.stderr}"
                )
            
            # 4. Extract outputs
            discharge = self._extract_discharge()
            inundation = self._extract_inundation()
            
            return SimulationResult(
                discharge=discharge,
                inundation=inundation,
                success=True
            )
        
        except subprocess.TimeoutExpired:
            return SimulationResult(
                success=False,
                error="LISFLOOD-FP simulation timeout (1 hour)"
            )
        except Exception as e:
            return SimulationResult(
                success=False,
                error=f"LISFLOOD-FP error: {str(e)}"
            )
    
    def _create_input_file(self, storm_event: Dict[str, Any], par_file: str) -> str:
        """
        Create LISFLOOD-FP input file.
        
        LISFLOOD-FP typically uses a single input file that references the .par file.
        """
        input_file = os.path.join(self.output_dir, "input.txt")
        
        with open(input_file, 'w') as f:
            f.write(f"PARFILE: {par_file}\n")
            f.write(f"STARTDATE: {storm_event['start_date']}\n")
            f.write(f"ENDDATE: {storm_event['end_date']}\n")
            f.write(f"RAINDATA: {storm_event['rainfall_data']}\n")
            f.write(f"BOUNDARY: {storm_event['boundary_discharge']}\n")
        
        return input_file
    
    def _extract_discharge(self) -> Dict[str, np.ndarray]:
        """
        Extract discharge time series from LISFLOOD-FP output.
        
        LISFLOOD-FP typically outputs discharge in a .dat or .out file.
        Format: timestamp (seconds), discharge (m³/s)
        """
        discharge_file = os.path.join(self.output_dir, "discharge.out")
        
        if not os.path.exists(discharge_file):
            logger.warning("Discharge file not found: %s", discharge_file)
            return {"timestamps": np.array([]), "values": np.array([])}
        
        # Parse discharge file (format may vary)
        timestamps = []
        values = []
        
        with open(discharge_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        timestamps.append(float(parts[0]))  # Time in seconds
                        values.append(float(parts[1]))      # Discharge in m³/s
                    except ValueError:
                        continue
        
        return {
            "timestamps": np.array(timestamps),
            "values": np.array(values)
        }
    
    def _extract_inundation(self) -> Dict[str, Any]:
        """
        Extract inundation depth grid from LISFLOOD-FP output.
        
        LISFLOOD-FP outputs spatial results as .asc or .tif files.
        """
        # Try multiple output formats
        possible_files = [
            os.path.join(self.output_dir, "inundation.tif"),
            os.path.join(self.output_dir, "depth.asc"),
            os.path.join(self.output_dir, "max_depth.tif"),
        ]
        
        inundation_file = None
        for f in possible_files:
            if os.path.exists(f):
                inundation_file = f
                break
        
        if inundation_file is None:
            logger.warning("Inundation file not found in %s", self.output_dir)
            return {"grid": np.array([]), "resolution": None, "crs": None}
        
        # Read with rasterio
        with rasterio.open(inundation_file) as src:
            grid = src.read(1)  # 2D array of depths (meters)
            resolution = src.res[0]  # Grid resolution (meters)
            crs = src.crs  # Coordinate reference system
            transform = src.transform  # Geotransform
        
        return {
            "grid": grid,
            "resolution": resolution,
            "crs": crs,
            "transform": transform
        }
    
    def get_outputs(self) -> Dict[str, Any]:
        """Get model outputs for calibration"""
        discharge = self._extract_discharge()
        inundation = self._extract_inundation()
        
        return {
            "discharge": discharge,
            "inundation": inundation
        }
    
    def save_config(self, output_path: str):
        """Save calibrated parameters to .par file"""
        par_file = self.set_parameters(self.current_params)
        
        # Copy to output path
        shutil.copy2(par_file, output_path)
        
        print(f"Calibrated parameters saved to: {output_path}")
