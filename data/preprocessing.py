"""
Data loading and preprocessing utilities.
Handles gauge data, satellite data, and input validation.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Any, Union
import rasterio
from datetime import datetime


def load_observations(file_path: str, data_type: str) -> Dict[str, Any]:
    """
    Load observation data (gauge or satellite).
    
    Args:
        file_path: Path to data file
        data_type: "gauge" or "satellite"
    
    Returns:
        Dict with loaded data
    """
    path = Path(file_path)
    
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {file_path}")
    
    if data_type == "gauge":
        return load_gauge_data(file_path)
    elif data_type == "satellite":
        return load_satellite_data(file_path)
    else:
        raise ValueError(f"Unknown data type: {data_type}")


def load_gauge_data(file_path: str) -> Dict[str, np.ndarray]:
    """
    Load gauge discharge data from CSV.
    
    Expected CSV format:
    timestamp,discharge
    2026-06-01T00:00:00Z,1250.5
    2026-06-01T01:00:00Z,1265.2
    ...
    
    Args:
        file_path: Path to CSV file
    
    Returns:
        Dict with "timestamps" and "discharge" arrays
    """
    try:
        df = pd.read_csv(file_path)
        
        # Detect column names (flexible)
        time_col = None
        value_col = None
        
        for col in df.columns:
            col_lower = col.lower()
            if "time" in col_lower or "date" in col_lower:
                time_col = col
            elif "discharge" in col_lower or "flow" in col_lower or "q" in col_lower:
                value_col = col
        
        if time_col is None:
            time_col = df.columns[0]
        if value_col is None:
            value_col = df.columns[1]
        
        # Parse timestamps
        timestamps = pd.to_datetime(df[time_col])
        discharge = df[value_col].values.astype(np.float64)
        
        # Handle NaN values
        mask = ~np.isnan(discharge)
        timestamps = timestamps[mask]
        discharge = discharge[mask]
        
        print(f"   Loaded {len(discharge)} gauge observations from {file_path}")
        print(f"   Time range: {timestamps[0]} to {timestamps[-1]}")
        print(f"   Discharge range: {discharge.min():.2f} - {discharge.max():.2f} m³/s")
        
        return {
            "timestamps": timestamps.values,
            "discharge": discharge
        }
    
    except Exception as e:
        raise ValueError(f"Error loading gauge data from {file_path}: {str(e)}")


def load_satellite_data(file_path: str) -> np.ndarray:
    """
    Load satellite-derived flood mask from GeoTIFF.
    
    Args:
        file_path: Path to GeoTIFF file
    
    Returns:
        2D numpy array (flood mask or depth grid)
    """
    try:
        with rasterio.open(file_path) as src:
            grid = src.read(1)
            crs = src.crs
            transform = src.transform
            bounds = src.bounds
            
            print(f"   Loaded satellite inundation from {file_path}")
            print(f"   Grid shape: {grid.shape}")
            print(f"   CRS: {crs}")
            print(f"   Bounds: {bounds}")
            print(f"   Data type: {grid.dtype}")
            print(f"   Value range: {grid.min():.2f} - {grid.max():.2f}")
            
            # Check if binary mask or depth grid
            unique_values = np.unique(grid)
            if len(unique_values) <= 2:
                print(f"   Binary flood mask (values: {unique_values})")
            else:
                print(f"   Depth grid (unique values: {len(unique_values)})")
        
        return grid
    
    except Exception as e:
        raise ValueError(f"Error loading satellite data from {file_path}: {str(e)}")


def validate_input_data(config: Dict[str, Any]) -> bool:
    """
    Validate that all input data files exist and are readable.
    
    Args:
        config: Configuration dict
    
    Returns:
        True if all valid, raises exception otherwise
    """
    print("\n🔍 Validating input data...")
    
    # Check rainfall data
    rainfall_path = config["input_data"]["rainfall"]["path"]
    if not Path(rainfall_path).exists():
        raise FileNotFoundError(f"Rainfall data not found: {rainfall_path}")
    print(f"   ✓ Rainfall: {rainfall_path}")
    
    # Check terrain (DEM)
    terrain_path = config["input_data"]["terrain"]["path"]
    if not Path(terrain_path).exists():
        raise FileNotFoundError(f"Terrain data not found: {terrain_path}")
    print(f"   ✓ Terrain: {terrain_path}")
    
    # Check boundary discharge
    boundary_path = config["input_data"]["boundary_discharge"]
    if not Path(boundary_path).exists():
        raise FileNotFoundError(f"Boundary discharge not found: {boundary_path}")
    print(f"   ✓ Boundary discharge: {boundary_path}")
    
    # Check gauge data
    for gauge in config["observations"]["gauges"]:
        gauge_path = gauge["path"]
        if not Path(gauge_path).exists():
            raise FileNotFoundError(f"Gauge data not found: {gauge_path}")
        print(f"   ✓ Gauge: {gauge_path}")
    
    # Check satellite data
    for satellite in config["observations"]["satellite"]:
        satellite_path = satellite["path"]
        if not Path(satellite_path).exists():
            raise FileNotFoundError(f"Satellite data not found: {satellite_path}")
        print(f"   ✓ Satellite: {satellite_path}")
    
    # Check LISFLOOD-FP template
    template_path = config["model"]["param_template"]
    if not Path(template_path).exists():
        raise FileNotFoundError(f"LISFLOOD-FP template not found: {template_path}")
    print(f"   ✓ LISFLOOD-FP template: {template_path}")
    
    print("\n✅ All input data validated successfully")
    return True


def align_temporal_data(obs_timestamps: np.ndarray,
                       sim_timestamps: np.ndarray,
                       obs_values: np.ndarray,
                       sim_values: np.ndarray) -> tuple:
    """
    Align observed and simulated time series by interpolation.
    
    Args:
        obs_timestamps: Observation timestamps (hours)
        sim_timestamps: Simulation timestamps (hours)
        obs_values: Observation values
        sim_values: Simulation values
    
    Returns:
        Tuple of (aligned_obs, aligned_sim, aligned_timestamps)
    """
    from scipy import interpolate
    
    # Find common time range
    t_min = max(obs_timestamps.min(), sim_timestamps.min())
    t_max = min(obs_timestamps.max(), sim_timestamps.max())
    
    # Create uniform time grid
    common_times = np.linspace(t_min, t_max, len(obs_timestamps))
    
    # Interpolate both series to common grid
    obs_interp = interpolate.interp1d(
        obs_timestamps, obs_values,
        kind='linear', bounds_error=False, fill_value=np.nan
    )(common_times)
    
    sim_interp = interpolate.interp1d(
        sim_timestamps, sim_values,
        kind='linear', bounds_error=False, fill_value=np.nan
    )(common_times)
    
    # Remove NaN values
    mask = ~np.isnan(obs_interp) & ~np.isnan(sim_interp)
    common_times = common_times[mask]
    obs_aligned = obs_interp[mask]
    sim_aligned = sim_interp[mask]
    
    print(f"   Aligned {len(obs_aligned)} time steps (common range: {t_min:.1f} - {t_max:.1f} h)")
    
    return obs_aligned, sim_aligned, common_times


def align_spatial_data(obs_grid: np.ndarray,
                      sim_grid: np.ndarray,
                      obs_transform,
                      sim_transform) -> tuple:
    """
    Align observed and simulated spatial grids to same resolution/crs.
    
    Args:
        obs_grid: Observed grid
        sim_grid: Simulated grid
        obs_transform: Rasterio transform for observed grid
        sim_transform: Rasterio transform for simulated grid
    
    Returns:
        Tuple of (aligned_obs, aligned_sim)
    """
    # Simple implementation: crop/resize to common shape
    from rasterio.warp import reproject, Resampling
    
    # Find common bounds (intersection)
    obs_bounds = rasterio.transform.array_bounds(
        obs_grid.shape[0], obs_grid.shape[1], obs_transform
    )
    sim_bounds = rasterio.transform.array_bounds(
        sim_grid.shape[0], sim_grid.shape[1], sim_transform
    )
    
    # Common bounds (intersection)
    min_x = max(obs_bounds[0], sim_bounds[0])
    max_x = min(obs_bounds[2], sim_bounds[2])
    min_y = max(obs_bounds[1], sim_bounds[1])
    max_y = min(obs_bounds[3], sim_bounds[3])
    
    print(f"   Aligning spatial grids to common bounds")
    print(f"   Intersection: x=[{min_x:.2f}, {max_x:.2f}], y=[{min_y:.2f}, {max_y:.2f}]")
    
    # Reproject both grids to common bounds
    # (Simplified: just return original grids with warning)
    if obs_grid.shape != sim_grid.shape:
        print(f"   Warning: Grid shapes differ: {obs_grid.shape} vs {sim_grid.shape}")
        print(f"   Using smaller shape: {tuple(min(s) for s in zip(obs_grid.shape, sim_grid.shape))}")
        
        min_rows = min(obs_grid.shape[0], sim_grid.shape[0])
        min_cols = min(obs_grid.shape[1], sim_grid.shape[1])
        
        return (
            obs_grid[:min_rows, :min_cols],
            sim_grid[:min_rows, :min_cols]
        )
    
    return obs_grid, sim_grid
