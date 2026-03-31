"""
LISFLOOD-FP 8.2 model adapter.
Specifically adapted for urabbani/lisflood-fp_8.2_update with GPU acceleration,
multiple solvers, and enhanced features.
"""

import os
import re
import shutil
import subprocess
import logging
import numpy as np
import rasterio
from pathlib import Path
from typing import Dict, Any, Tuple
import pandas as pd

from models.adapter import ModelAdapter, SimulationResult

logger = logging.getLogger(__name__)


class LISFLOOD82Adapter(ModelAdapter):
    """Adapter for LISFLOOD-FP 8.2 with GPU acceleration and advanced solvers"""

    def __init__(self, config: Dict[str, Any]):
        super().__init__(config)

        self.executable = self._resolve_executable(config["executable"])
        self.param_template = config["param_template"]
        self.output_dir = config["output_dir"]
        
        # LISFLOOD-FP 8.2 calibratable parameters
        self.param_ranges = {
            # Hydraulic parameters
            "fpfric": {
                "value": 0.035,
                "min": 0.01,
                "max": 0.1,
                "description": "Floodplain Manning's n",
                "unit": "s/m^(1/3)"
            },
            "SGCn": {
                "value": 0.035,
                "min": 0.02,
                "max": 0.06,
                "description": "Sub-grid channel Manning's n",
                "unit": "s/m^(1/3)"
            },
            
            # Infiltration
            "infinity": {
                "value": 0.0001,
                "min": 0.0,
                "max": 0.001,
                "description": "Infiltration rate",
                "unit": "m/s"
            },
            
            # Numerical parameters
            "courant": {
                "value": 0.7,
                "min": 0.3,
                "max": 0.9,
                "description": "CFL number (time step control)",
                "unit": "-"
            },
            "theta": {
                "value": 1.0,
                "min": 0.5,
                "max": 1.0,
                "description": "Time-stepping parameter",
                "unit": "-"
            },
            "depth_thresh": {
                "value": 0.001,
                "min": 0.0005,
                "max": 0.01,
                "description": "Minimum depth threshold",
                "unit": "m"
            },
            
            # Physical parameters
            "gravity": {
                "value": 9.8,
                "min": 9.7,
                "max": 9.9,
                "description": "Gravitational acceleration",
                "unit": "m/s²"
            },
            "max_Froude": {
                "value": 10000.0,
                "min": 1.0,
                "max": 10000.0,
                "description": "Maximum Froude number",
                "unit": "-"
            },
            
            # Sub-grid channel parameters
            "SGCr": {
                "value": 0.15,
                "min": 0.1,
                "max": 0.2,
                "description": "Sub-grid channel multiplier",
                "unit": "-"
            },
            "SGCp": {
                "value": 0.78,
                "min": 0.7,
                "max": 0.85,
                "description": "Sub-grid channel exponent",
                "unit": "-"
            },
        }
        
        # Solver options
        self.solver_options = {
            "ACC": "Acceleration solver (simplified momentum)",
            "FV1": "First-order finite volume",
            "FV2": "Second-order finite volume",
            "DG2": "Discontinuous Galerkin (higher order)",
            "ACC_NUGRID": "Non-uniform grid acceleration"
        }
        
        # Load template .par file
        self.template_lines = self._load_template()
        
        # Simulation options
        self.use_cuda = config.get("use_cuda", True)
        self.verbose = config.get("verbose", True)
        self.solver = config.get("solver", "ACC")  # Default: ACC (fastest)
    
    @staticmethod
    def _resolve_executable(executable: str) -> str:
        """Resolve 'auto' to the integrated LISFLOOD-FP binary, compiling if needed."""
        if executable != "auto":
            return executable

        # Look for integrated LISFLOOD-FP source
        project_root = Path(__file__).resolve().parent.parent.parent
        source_dir = project_root / "src" / "lisflood-fp"
        build_dir = source_dir / "build"

        if not source_dir.exists():
            raise FileNotFoundError(
                "LISFLOOD-FP source not found in src/lisflood-fp/"
            )

        # Check for existing binary
        for name in ["lisflood", "lisflood.exe"]:
            binary = build_dir / name
            if binary.exists():
                logger.info("Using integrated LISFLOOD-FP binary: %s", binary)
                return str(binary)

        # Also check Release subdirectory (Windows)
        release_binary = build_dir / "Release" / "lisflood.exe"
        if release_binary.exists():
            logger.info("Using integrated LISFLOOD-FP binary: %s", release_binary)
            return str(release_binary)

        # Binary not found — attempt to compile
        logger.info("LISFLOOD-FP binary not found, compiling from integrated source...")
        # Import build script from the new location
        import sys
        sys.path.append(str(project_root / "scripts"))
        from build_lisflood import build
        binary = build(clean=False, cuda=False)
        logger.info("Compiled LISFLOOD-FP: %s", binary)
        return str(binary)

    def _load_template(self) -> list:
        """Load template .par file lines"""
        if not os.path.exists(self.param_template):
            raise FileNotFoundError(f"Template file not found: {self.param_template}")
        
        with open(self.param_template, 'r') as f:
            return f.readlines()
    
    def get_parameters(self) -> Dict[str, Dict[str, float]]:
        """Return LISFLOOD-FP 8.2 parameters with ranges"""
        return self.param_ranges.copy()
    
    def set_parameters(self, params: Dict[str, float]) -> str:
        """
        Update LISFLOOD-FP 8.2 parameters in .par file.
        
        Parameters:
        - fpfric: Floodplain Manning's n
        - SGCn: Sub-grid channel Manning's n
        - infinity: Infiltration rate (m/s)
        - courant: CFL number
        - theta: Time-stepping parameter
        - depth_thresh: Minimum depth threshold
        - gravity: Gravitational acceleration
        - max_Froude: Maximum Froude number
        - SGCr: Sub-grid channel multiplier
        - SGCp: Sub-grid channel exponent
        
        Returns:
            Path to updated .par file
        """
        self.current_params = params.copy()
        
        # Create new .par file with updated parameters
        new_lines = []
        for line in self.template_lines:
            updated_line = line
            for param_name, param_value in params.items():
                # Pattern: PARAM_NAME VALUE (case-insensitive, escaped param name)
                pattern = rf'^{re.escape(param_name)}\s+\d+\.?\d*e?[+-]?\d*'
                if re.search(pattern, line, re.IGNORECASE):
                    # Preserve comments
                    comment = ""
                    if "#" in line:
                        comment = "#" + line.split("#")[1]
                    updated_line = f"{param_name} {param_value:.6f} {comment}"
            new_lines.append(updated_line)
        
        # Write to output directory
        os.makedirs(self.output_dir, exist_ok=True)
        output_par = os.path.join(self.output_dir, "simulation.par")
        with open(output_par, 'w') as f:
            f.writelines(new_lines)
        
        return output_par
    
    def run_simulation(self, storm_event: Dict[str, Any]) -> SimulationResult:
        """
        Run LISFLOOD-FP 8.2 simulation.
        
        Args:
            storm_event: Dict with:
                - rainfall_path: Path to rainfall data
                - boundary_path: Path to boundary condition
                - duration: Simulation duration (seconds)
        
        Returns:
            SimulationResult with discharge and inundation outputs
        """
        try:
            # 1. Update parameters
            par_file = self.set_parameters(self.current_params)
            
            # 2. Build command
            cmd = [self.executable]
            
            if self.verbose:
                cmd.append("-v")
            
            if self.use_cuda:
                cmd.append("-cuda")
            
            cmd.append(par_file)
            
            # 3. Run simulation
            print(f"   Running LISFLOOD-FP 8.2...")
            print(f"   Command: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1 hour timeout
            )
            
            if result.returncode != 0:
                error_msg = result.stderr if result.stderr else "Unknown error"
                return SimulationResult(
                    success=False,
                    error=f"LISFLOOD-FP failed (code {result.returncode}): {error_msg}"
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
    
    def _extract_discharge(self) -> Dict[str, np.ndarray]:
        """
        Extract discharge time series from LISFLOOD-FP 8.2 output.
        
        LISFLOOD-FP 8.2 outputs discharge to .discharge file at gauge points.
        Format: timestamp (s), discharge (m³/s) per gauge
        """
        discharge_file = os.path.join(self.output_dir, "simulation.discharge")
        
        if not os.path.exists(discharge_file):
            logger.warning("Primary discharge file not found: %s", discharge_file)
            # Try alternative names
            for prefix in ["sim1", "output"]:
                for ext in [".discharge", ".stage"]:
                    test_file = os.path.join(self.output_dir, f"{prefix}{ext}")
                    if os.path.exists(test_file):
                        discharge_file = test_file
                        break
                if os.path.exists(discharge_file):
                    break
        
        if not os.path.exists(discharge_file):
            logger.warning("Discharge file not found in %s", self.output_dir)
            return {"timestamps": np.array([]), "values": np.array([])}
        
        # Parse discharge file
        timestamps = []
        values = []
        
        try:
            with open(discharge_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('#') or not line:
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            timestamps.append(float(parts[0]))  # Time (seconds)
                            values.append(float(parts[1]))      # Discharge (m³/s)
                        except ValueError:
                            continue
        except Exception as e:
            logger.warning("Error reading discharge file: %s", str(e))
            return {"timestamps": np.array([]), "values": np.array([])}
        
        print(f"   Extracted {len(values)} discharge data points")
        
        return {
            "timestamps": np.array(timestamps),
            "values": np.array(values)
        }
    
    def _extract_inundation(self) -> Dict[str, Any]:
        """
        Extract inundation depth grid from LISFLOOD-FP 8.2 output.
        
        LISFLOOD-FP 8.2 outputs spatial results as:
        - .max: Maximum water depths
        - .wd: Water depths at save intervals
        - .elev: Water surface elevations
        
        Formats: ASCII (.asc) or Binary (.bin), optionally NetCDF
        """
        # Try multiple output files (priority: .max, .wd)
        possible_files = [
            os.path.join(self.output_dir, "simulation.max"),
            os.path.join(self.output_dir, "sim1.max"),
            os.path.join(self.output_dir, "simulation.wd"),
            os.path.join(self.output_dir, "sim1.wd"),
            os.path.join(self.output_dir, "simulation.elev"),
            os.path.join(self.output_dir, "sim1.elev"),
        ]
        
        inundation_file = None
        for f in possible_files:
            if os.path.exists(f):
                inundation_file = f
                break
        
        if inundation_file is None:
            logger.warning("Inundation file not found in %s", self.output_dir)
            return {"grid": np.array([]), "resolution": None, "crs": None}
        
        print(f"   Found inundation file: {inundation_file}")
        
        # Read with rasterio (handles both ASCII and binary)
        try:
            with rasterio.open(inundation_file) as src:
                grid = src.read(1)  # 2D array of depths (meters)
                resolution = src.res[0]  # Grid resolution (meters)
                crs = src.crs  # Coordinate reference system
                transform = src.transform  # Geotransform
                
                print(f"   Grid shape: {grid.shape}")
                print(f"   Resolution: {resolution} meters")
                print(f"   CRS: {crs}")
                print(f"   Depth range: {grid.min():.3f} - {grid.max():.3f} m")
                print(f"   Flooded cells: {np.sum(grid > 0.1)} / {grid.size}")
            
            return {
                "grid": grid,
                "resolution": resolution,
                "crs": crs,
                "transform": transform
            }
        
        except Exception as e:
            logger.warning("Error reading inundation file: %s", str(e))
            
            # Try reading as ASCII text if rasterio fails
            try:
                if inundation_file.endswith('.asc'):
                    grid = self._read_ascii_raster(inundation_file)
                    return {
                        "grid": grid,
                        "resolution": None,
                        "crs": None,
                        "transform": None
                    }
            except Exception:
                return {"grid": np.array([]), "resolution": None, "crs": None}
    
    def _read_ascii_raster(self, file_path: str) -> np.ndarray:
        """Read ASCII raster file manually (ESRI format)"""
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        ncols = None
        nrows = None
        xllcorner = None
        yllcorner = None
        cellsize = None
        nodata_value = None
        
        header_lines = []
        for i, line in enumerate(lines):
            line_lower = line.strip().lower()
            if line_lower.startswith('ncols'):
                ncols = int(line.split()[1])
                header_lines.append(i)
            elif line_lower.startswith('nrows'):
                nrows = int(line.split()[1])
                header_lines.append(i)
            elif line_lower.startswith('xllcorner'):
                xllcorner = float(line.split()[1])
                header_lines.append(i)
            elif line_lower.startswith('yllcorner'):
                yllcorner = float(line.split()[1])
                header_lines.append(i)
            elif line_lower.startswith('cellsize'):
                cellsize = float(line.split()[1])
                header_lines.append(i)
            elif line_lower.startswith('nodata_value'):
                nodata_value = float(line.split()[1])
                header_lines.append(i)
        
        # Read data (skip header)
        data_lines = [line for i, line in enumerate(lines) if i not in header_lines]
        data = np.array([list(map(float, line.split())) for line in data_lines])
        
        print(f"   Read ASCII raster: {nrows}x{ncols}, cellsize={cellsize}m")
        
        return data
    
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
    
    def set_solver(self, solver: str):
        """Set numerical solver (ACC, FV1, FV2, DG2, ACC_NUGRID)"""
        if solver not in self.solver_options:
            raise ValueError(f"Unknown solver: {solver}. Options: {list(self.solver_options.keys())}")
        
        self.solver = solver
        print(f"Solver set to: {solver} ({self.solver_options[solver]})")
    
    def set_gpu(self, enable: bool):
        """Enable/disable GPU acceleration"""
        self.use_cuda = enable
        print(f"GPU acceleration: {'enabled' if enable else 'disabled'}")
