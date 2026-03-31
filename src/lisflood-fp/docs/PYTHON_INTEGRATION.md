# LISFLOOD-FP Python Integration Guide

This guide demonstrates how to integrate LISFLOOD-FP with Python workflows for pre-processing, job management, and post-processing tasks.

## 1. Prerequisites

Install the required Python packages:

```bash
pip install numpy pandas xarray matplotlib rasterio geopandas scipy netCDF4 h5py
```

## 2. Pre-processing with Python

### 2.1 DEM Preparation

```python
import rasterio
import numpy as np
from scipy import ndimage

def prepare_dem(dem_path, output_path, resolution=10, fill_nodata=True, smooth=True):
    """Prepare a DEM for LISFLOOD-FP simulation."""
    # Read the DEM
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        meta = src.meta.copy()
        
        # Resample to target resolution if needed
        if meta['transform'][0] != resolution:
            # Implement resampling here
            pass
            
        # Fill no-data values
        if fill_nodata:
            mask = dem == src.nodata
            dem = ndimage.gaussian_filter(dem, sigma=1)
            dem[~mask] = src.read(1)[~mask]  # Keep original values for non-nodata cells
            
        # Apply smoothing if requested
        if smooth:
            dem = ndimage.gaussian_filter(dem, sigma=0.5)
            
        # Write the processed DEM
        meta.update({
            'dtype': 'float32',
            'nodata': -9999
        })
        
        with rasterio.open(output_path, 'w', **meta) as dst:
            dst.write(dem.astype(np.float32), 1)
            
    return output_path
```

### 2.2 Generating Boundary Conditions

```python
def create_boundary_file(dem_path, output_path, inflow_points=None, stage_points=None):
    """Create a boundary condition file."""
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        meta = src.meta
        
    with open(output_path, 'w') as f:
        # Write header
        f.write("# LISFLOOD-FP boundary condition file\n")
        f.write("# Created with Python script\n\n")
        
        # Add inflow points
        if inflow_points:
            for point in inflow_points:
                row, col = point['row'], point['col']
                hydrograph = point['hydrograph']
                f.write(f"QPOINT {col} {row} {hydrograph}\n")
                
        # Add stage boundaries
        if stage_points:
            for point in stage_points:
                edge = point['edge']  # N, E, S, W
                start = point['start']
                end = point['end']
                timeseries = point['timeseries']
                f.write(f"HVAR {edge} {start} {end} {timeseries}\n")
    
    return output_path
```

### 2.3 Generating Manning's Roughness Map

```python
import geopandas as gpd

def create_roughness_map(dem_path, landuse_path, roughness_values, output_path):
    """Create a Manning's roughness map based on land use data."""
    # Read DEM for template
    with rasterio.open(dem_path) as src:
        dem = src.read(1)
        meta = src.meta.copy()
        
    # Read land use vector data
    landuse = gpd.read_file(landuse_path)
    
    # Create roughness raster
    roughness = np.ones_like(dem) * 0.035  # Default value
    
    # Assign roughness values based on land use
    for idx, row in landuse.iterrows():
        ltype = row['type']
        if ltype in roughness_values:
            # Rasterize this polygon to the DEM grid and assign roughness
            # Implementation depends on specific data format
            pass
    
    # Write the roughness raster
    meta.update({
        'dtype': 'float32',
        'nodata': -9999
    })
    
    with rasterio.open(output_path, 'w', **meta) as dst:
        dst.write(roughness.astype(np.float32), 1)
    
    return output_path
```

## 3. Running LISFLOOD-FP from Python

### 3.1 Single Simulation

```python
import subprocess
import os

def run_lisflood(parameter_file, output_dir, use_gpu=True):
    """Run a LISFLOOD-FP simulation."""
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = ['lisflood']
    if use_gpu:
        cmd.append('-cuda')
    cmd.append(parameter_file)
    
    # Run the simulation
    result = subprocess.run(cmd, 
                           capture_output=True, 
                           text=True, 
                           check=True)
    
    return result.stdout
```

### 3.2 Parameter Sensitivity Analysis

```python
import itertools
import pandas as pd

def sensitivity_analysis(base_parameter_file, output_dir, parameters):
    """Run sensitivity analysis by varying parameters."""
    # Create parameter combinations
    param_names = parameters.keys()
    param_values = [parameters[name] for name in param_names]
    combinations = list(itertools.product(*param_values))
    
    results = []
    
    for i, combination in enumerate(combinations):
        # Create parameter file for this combination
        param_dict = {name: value for name, value in zip(param_names, combination)}
        param_file = create_parameter_file(base_parameter_file, 
                                          f"{output_dir}/param_set_{i}.par",
                                          param_dict)
        
        # Run the simulation
        run_dir = f"{output_dir}/run_{i}"
        stdout = run_lisflood(param_file, run_dir)
        
        # Extract results
        max_depth, inundation_area = extract_results(run_dir)
        
        # Store results
        results.append({
            **param_dict,
            'max_depth': max_depth,
            'inundation_area': inundation_area
        })
    
    # Create summary dataframe
    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{output_dir}/sensitivity_results.csv", index=False)
    
    return results_df
```

## 4. Post-processing with Python

### 4.1 Reading Model Output

```python
import xarray as xr

def read_results(output_dir, timestep=None):
    """Read LISFLOOD-FP results."""
    # For NetCDF output
    if os.path.exists(f"{output_dir}/results.nc"):
        ds = xr.open_dataset(f"{output_dir}/results.nc")
        return ds
        
    # For traditional ASCII output
    depth_files = sorted(glob.glob(f"{output_dir}/*.wd"))
    
    if timestep is not None:
        # Read specific timestep
        depth_file = [f for f in depth_files if f.endswith(f"{timestep:05d}.wd")]
        if depth_file:
            with rasterio.open(depth_file[0]) as src:
                depth = src.read(1)
                meta = src.meta
            return depth, meta
    else:
        # Read all timesteps
        depths = []
        for file in depth_files:
            with rasterio.open(file) as src:
                depth = src.read(1)
                depths.append(depth)
                if len(depths) == 1:
                    meta = src.meta
        
        return np.array(depths), meta
```

### 4.2 Calculating Flood Statistics

```python
def calculate_flood_statistics(depth_array, cell_area, threshold=0.1):
    """Calculate flood statistics from depth array."""
    stats = {}
    
    # Maximum water depth
    stats['max_depth'] = np.max(depth_array)
    
    # Inundation area
    flood_cells = np.sum(depth_array > threshold)
    stats['inundation_area'] = flood_cells * cell_area
    
    # Flood volume
    stats['flood_volume'] = np.sum(depth_array[depth_array > threshold]) * cell_area
    
    # Inundation duration (if time dimension exists)
    if len(depth_array.shape) > 2:
        duration = np.sum(depth_array > threshold, axis=0)
        stats['max_duration'] = np.max(duration)
        stats['mean_duration'] = np.mean(duration[duration > 0])
    
    return stats
```

### 4.3 Visualization

```python
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def plot_flood_map(depth, dem, output_path, title="Flood Depth"):
    """Create a flood depth visualization."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create a custom colormap for flood depths
    colors = [(0.8, 0.8, 1), (0, 0, 0.8), (0, 0, 0.5)]
    cmap = LinearSegmentedColormap.from_list('flood_cmap', colors)
    
    # Plot hillshade of DEM for context
    ls = plt.matplotlib.colors.LightSource(azdeg=315, altdeg=45)
    hillshade = ls.hillshade(dem, vert_exag=2)
    ax.imshow(hillshade, cmap='gray', alpha=0.5)
    
    # Plot flood depths
    masked_depth = np.ma.masked_where(depth < 0.01, depth)
    im = ax.imshow(masked_depth, cmap=cmap, vmin=0, vmax=2)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Flood Depth (m)')
    
    # Set title and save
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path
```

## 5. Complete Workflow Example

```python
def complete_workflow(dem_path, landuse_path, rainfall_path, output_dir):
    """Run a complete LISFLOOD-FP workflow."""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Preprocess data
    processed_dem = prepare_dem(dem_path, f"{output_dir}/processed_dem.asc")
    
    # Create roughness map
    roughness_values = {
        'urban': 0.05,
        'forest': 0.1,
        'grass': 0.035,
        'water': 0.03
    }
    roughness_map = create_roughness_map(processed_dem, landuse_path, 
                                         roughness_values, 
                                         f"{output_dir}/roughness.asc")
    
    # 2. Create parameter file
    param_file = f"{output_dir}/simulation.par"
    with open(param_file, 'w') as f:
        f.write("# LISFLOOD-FP parameter file\n")
        f.write(f"DEMfile         {processed_dem}\n")
        f.write(f"manningsn       {roughness_map}\n")
        f.write(f"rain            {rainfall_path}\n")
        f.write(f"resroot         {output_dir}/results\n")
        f.write("saveint         3600\n")
        f.write("sim_time        86400\n")
        f.write("solver          FV1\n")
        f.write("cuda            1\n")
    
    # 3. Run simulation
    run_output = run_lisflood(param_file, output_dir)
    print(run_output)
    
    # 4. Process results
    depth_array, meta = read_results(f"{output_dir}/results")
    
    cell_area = abs(meta['transform'][0] * meta['transform'][4])
    stats = calculate_flood_statistics(depth_array[-1], cell_area)
    print(stats)
    
    # 5. Visualize results
    with rasterio.open(processed_dem) as src:
        dem = src.read(1)
    
    plot_flood_map(depth_array[-1], dem, f"{output_dir}/flood_map.png")
    
    return stats
```

## 6. Tips for Working with Large Datasets

1. Use memory-mapped arrays (with numpy.memmap) for large raster files
2. Process large domains in tiles
3. Use dask for parallel processing of large datasets
4. With NetCDF outputs, use xarray with dask to process data without loading everything into memory
5. Consider using the LISFLOOD-FP Docker container for consistent environments
6. For multiple simulations, take advantage of multiprocessing to run them in parallel

## 7. Advanced Integration: Real-time Forecasting

For real-time flood forecasting systems, you can integrate LISFLOOD-FP with weather forecast data:

```python
def forecast_workflow(dem_path, forecast_data, previous_state=None):
    """Run a forecast simulation using weather forecast data."""
    # Extract rainfall forecast
    rainfall = extract_rainfall_from_forecast(forecast_data)
    
    # Create rainfall file
    rainfall_file = create_rainfall_file(rainfall, "forecast_rain.asc")
    
    # Set up simulation
    if previous_state:
        # Continue from previous state
        start_option = f"startfile {previous_state}"
    else:
        # Start from dry conditions
        start_option = ""
    
    # Create parameter file
    param_file = "forecast.par"
    with open(param_file, 'w') as f:
        f.write("# Forecast simulation\n")
        f.write(f"DEMfile         {dem_path}\n")
        f.write(f"rain            {rainfall_file}\n")
        f.write(f"{start_option}\n")
        f.write("resroot         forecast_results\n")
        f.write("saveint         3600\n")
        f.write("sim_time        86400\n")
        f.write("solver          FV1\n")
        f.write("cuda            1\n")
    
    # Run forecast
    run_lisflood(param_file, "forecast_results")
    
    # Save final state for next forecast
    latest_state = find_latest_state("forecast_results")
    
    return "forecast_results", latest_state
```