# LISFLOOD-FP Rainfall Input Guide

This document provides comprehensive information about the various rainfall input options available in LISFLOOD-FP, their formats, and usage examples.

## Table of Contents

1. [Introduction](#introduction)
2. [Rainfall Input Methods](#rainfall-input-methods)
3. [File Formats](#file-formats)
4. [Unit Conversions](#unit-conversions)
5. [Rainfall Masks](#rainfall-masks)
6. [Dynamic Rainfall](#dynamic-rainfall)
7. [Rainfall with Routing](#rainfall-with-routing)
8. [Parameter File Examples](#parameter-file-examples)
9. [Best Practices](#best-practices)
10. [Troubleshooting](#troubleshooting)

## Introduction

LISFLOOD-FP can simulate rainfall-driven flooding through multiple methods of increasing complexity:

- Uniform rainfall across the entire domain
- Time-varying uniform rainfall
- Spatially distributed constant rainfall
- Spatially distributed and time-varying rainfall (dynamic rainfall)
- Masked rainfall to specific areas
- Rainfall with routing (for large domains/low resolutions)

The appropriate method depends on the size and resolution of your domain, available rainfall data, and the flooding process you're modeling.

## Rainfall Input Methods

### 1. Uniform Constant Rainfall

Simplest method. Applies a single rainfall rate across the entire domain.

**Parameter file setting:**
```
rain  0.00001  # Rainfall rate in m/s (36mm/hr)
```

### 2. Time-varying Uniform Rainfall

Applies uniform rainfall that changes over time according to a time series.

**Parameter file setting:**
```
rain  rainfall_timeseries.txt  # Path to time series file
```

### 3. Spatially Distributed Constant Rainfall

Applies different rainfall rates across the domain, but constant in time.

**Parameter file setting:**
```
rain  rainfall_map.asc  # Path to rainfall grid
```

### 4. Spatially Distributed and Time-varying Rainfall

Most complex approach, providing full spatial and temporal variability.

**Parameter file setting:**
```
dynamicrainfile  rainfall_data.nc  # Path to NetCDF file
```

### 5. Masked Rainfall

Restricts rainfall to specific areas using a binary mask.

**Parameter file setting:**
```
rain          rainfall_timeseries.txt  # Rainfall time series
rainfallmask  urban_mask.asc           # Binary mask (0/1)
```

## File Formats

### Time Series (.txt format)

```
# Time(s) RainfallRate(m/s)
0.0       0.0
3600.0    5.556e-6       # 20mm/hr
7200.0    8.333e-6       # 30mm/hr
10800.0   2.778e-6       # 10mm/hr
14400.0   0.0
```

Important notes:
- First column is time in seconds from simulation start
- Second column is rainfall rate in meters per second
- Comments can be added with #
- No header row (except comments)
- Linear interpolation is used between time points

### Rainfall Grid (.asc format)

```
ncols        200
nrows        200
xllcorner    80000.0
yllcorner    40000.0
cellsize     5.0
NODATA_value -9999
0.0000 0.0000 0.0000 0.0000 ...
0.0000 1.1e-6 1.1e-6 1.1e-6 ...
...
```

Important notes:
- Standard ESRI ASCII grid format
- Values must be in m/s
- Grid dimensions and cell size must match the DEM
- NODATA cells will not receive rainfall

### NetCDF Dynamic Rainfall (.nc format)

The NetCDF file must contain:
- A 3D variable named "rainfall" or "rain" with dimensions (time, y, x)
- Units in m/s
- Time dimension with units attribute (e.g., "seconds since 2023-01-01 00:00:00")
- Coordinate variables for x, y, and time
- Grid dimensions that match or can be resampled to the DEM

Example NetCDF structure:
```
dimensions:
  time = UNLIMITED ;
  y = 200 ;
  x = 200 ;
variables:
  double time(time) ;
    time:units = "seconds since 2023-01-01 00:00:00" ;
  double y(y) ;
    y:units = "m" ;
  double x(x) ;
    x:units = "m" ;
  float rainfall(time, y, x) ;
    rainfall:units = "m s-1" ;
    rainfall:_FillValue = -9999.f ;
```

### Rainfall Mask (.asc format)

```
ncols        200
nrows        200
xllcorner    80000.0
yllcorner    40000.0
cellsize     5.0
NODATA_value -9999
0 0 0 0 0 0 0 ...
0 1 1 1 1 1 0 ...
0 1 1 1 1 1 0 ...
...
```

Important notes:
- Binary grid (0 = no rainfall, 1 = rainfall)
- Must match DEM dimensions exactly
- Cannot be used with dynamicrainfile option

## Unit Conversions

LISFLOOD-FP uses m/s for rainfall rates internally. Convert from common units:

| From | To m/s | Conversion Factor |
|------|--------|-------------------|
| mm/hr | m/s | Multiply by 2.778 × 10⁻⁷ |
| mm/day | m/s | Multiply by 1.157 × 10⁻⁸ |
| inch/hr | m/s | Multiply by 7.056 × 10⁻⁶ |
| inch/day | m/s | Multiply by 2.940 × 10⁻⁷ |

Common rainfall intensities in m/s:
- Light rain (1 mm/hr): 2.778 × 10⁻⁷ m/s
- Moderate rain (5 mm/hr): 1.389 × 10⁻⁶ m/s
- Heavy rain (20 mm/hr): 5.556 × 10⁻⁶ m/s
- Extreme rain (100 mm/hr): 2.778 × 10⁻⁵ m/s

## Rainfall Masks

Rainfall masks allow you to limit rainfall to specific areas of your domain. This is useful for:
- Urban drainage modeling (apply rainfall only to urban areas)
- Modeling specific catchments within a larger domain
- Testing scenarios with localized rainfall

The mask is a binary grid (0/1) where:
- 1: Rainfall is applied at the rate specified
- 0: No rainfall is applied

## Dynamic Rainfall

Dynamic rainfall combines both temporal and spatial variability, offering the most realistic representation of actual rainfall events. This is particularly important for:
- Flash flooding events
- Convective storms with high spatial variability
- Large domains with rainfall moving across the area

### Using Multiple ASC Files (Workaround)

LISFLOOD-FP does not natively support reading a folder of ASC files for time-varying spatial rainfall. However, you can work around this limitation by:

1. **Converting to NetCDF**: Convert your collection of ASC files into a single NetCDF file. This is the recommended approach.

2. **Pre-processing script**: Create a script that modifies your parameter file between time steps, changing the rainfall input file and restarting the simulation.

Example Python script to convert multiple ASC files to NetCDF:
```python
import numpy as np
import xarray as xr
import glob
import rasterio
import os
import re

# Find all rainfall ASC files
rainfall_files = sorted(glob.glob("rainfall_*.asc"))

# Extract timestamps from filenames (assuming format rainfall_HHMMSS.asc)
timestamps = []
for filename in rainfall_files:
    match = re.search(r'rainfall_(\d+)\.asc', filename)
    if match:
        time_str = match.group(1)
        # Convert HHMMSS to seconds - adjust this parsing as needed
        h = int(time_str[0:2])
        m = int(time_str[2:4])
        s = int(time_str[4:6])
        timestamps.append(h * 3600 + m * 60 + s)

# Read the first file to get dimensions
with rasterio.open(rainfall_files[0]) as src:
    height = src.height
    width = src.width
    transform = src.transform
    nodata = src.nodata
    rainfall_data = np.zeros((len(rainfall_files), height, width))
    
    # Create coordinate arrays
    y_coords = np.arange(src.bounds.bottom, src.bounds.top, src.res[1])
    if len(y_coords) > height:
        y_coords = y_coords[:height]
    y_coords = np.flip(y_coords)  # ASC files start from top
    
    x_coords = np.arange(src.bounds.left, src.bounds.right, src.res[0])
    if len(x_coords) > width:
        x_coords = x_coords[:width]

# Read all rainfall files
for i, file in enumerate(rainfall_files):
    with rasterio.open(file) as src:
        rain = src.read(1)
        rainfall_data[i, :, :] = rain
        
        # Convert from mm/hr to m/s if needed (adjust factor as needed)
        # rainfall_data[i, :, :] = rainfall_data[i, :, :] * 2.778e-7

# Create xarray dataset
ds = xr.Dataset(
    data_vars=dict(
        rainfall=(["time", "y", "x"], rainfall_data)
    ),
    coords=dict(
        time=timestamps,
        y=y_coords,
        x=x_coords
    ),
)

# Add variable attributes
ds.rainfall.attrs["units"] = "m s-1"
ds.rainfall.attrs["long_name"] = "Rainfall rate"
ds.time.attrs["units"] = "seconds since simulation start"

# Save to NetCDF
ds.to_netcdf("dynamic_rainfall.nc")
print(f"Converted {len(rainfall_files)} ASC files to NetCDF")
```

### Creating NetCDF Rainfall Files

You can create NetCDF rainfall files from various sources:

1. Using Python with xarray:
```python
import xarray as xr
import numpy as np

# Create time and space dimensions
times = np.arange(0, 86400, 3600)  # 24 hours, hourly
y = np.arange(200)
x = np.arange(200)

# Create rainfall data (example with moving storm)
rainfall = np.zeros((len(times), len(y), len(x)))
for t in range(len(times)):
    # Create a moving rain pattern
    center_x = int(x.max() * t / len(times))
    center_y = int(y.max() / 2)
    for i in range(len(y)):
        for j in range(len(x)):
            dist = np.sqrt((i - center_y)**2 + (j - center_x)**2)
            if dist < 20:  # 20-cell radius storm
                rainfall[t, i, j] = 5.556e-6 * np.exp(-dist/10)  # 20mm/hr at center

# Create dataset
ds = xr.Dataset(
    data_vars=dict(
        rainfall=(["time", "y", "x"], rainfall)
    ),
    coords=dict(
        time=times,
        y=y,
        x=x
    ),
    attrs=dict(description="Example rainfall data for LISFLOOD-FP"),
)

# Add variable attributes
ds.rainfall.attrs["units"] = "m s-1"
ds.rainfall.attrs["long_name"] = "Rainfall rate"
ds.time.attrs["units"] = "seconds since 2023-01-01 00:00:00"

# Save to NetCDF
ds.to_netcdf("rainfall_data.nc")
```

2. Converting from weather radar data:
   - Weather radar data is often available in NetCDF or HDF5 format
   - Convert radar reflectivity to rainfall rates using Z-R relationships
   - Reproject and resample to match your DEM

3. From climate model outputs:
   - Downscale GCM/RCM outputs to appropriate resolution
   - Convert to NetCDF format with the correct dimensions and units

## Rainfall with Routing

For large domains or coarse resolutions, rainfall shouldn't be assumed to immediately contribute to surface water. The routing option allows rainfall to flow following topography.

**Parameter file settings:**
```
rain       rainfall_data.txt  # Rainfall input
routing    ON                 # Enable routing
routespeed 0.1                # Flow speed (m/s)
```

Additional routing options:
```
dist_routing  ON              # Slope-dependent flow speeds
RouteSfThresh 0.1             # Friction slope threshold
```

How routing works:
1. A flow direction map is created based on DEM slopes
2. Rainfall is added to water depths
3. Water moves between cells based on surface slopes and specified routing speed
4. Water only enters the full hydraulic model when depths exceed the specified threshold

## Parameter File Examples

### Example 1: Simple Uniform Rainfall

```
# Basic simulation with uniform rainfall
DEMfile       terrain.asc
resroot       ./results
saveint       3600
sim_time      86400
fpfric        0.035
rain          0.00001         # 36 mm/hr uniform rainfall
```

### Example 2: Time-varying Rainfall with Mask

```
# Urban flood simulation with time-varying rainfall
DEMfile       city_dem.asc
resroot       ./results/urban
saveint       300
sim_time      43200
fpfric        0.05
solver        FV1             # Finite volume solver
rain          storm_hyetograph.txt  # Time-varying rainfall
rainfallmask  urban_area.asc  # Only apply to urban areas
```

### Example 3: Dynamic Rainfall with Routing

```
# Large catchment with radar rainfall
DEMfile         catchment_dem.asc
resroot         ./results/catchment
saveint         1800
sim_time        259200         # 3 days
fpfric          0.035
solver          ACC            # Acceleration solver
dynamicrainfile radar_rainfall.nc  # NetCDF radar data
routing         ON             # Enable routing
routespeed      0.1            # 0.1 m/s routing speed
dist_routing    ON             # Slope-dependent routing
```

### Example 4: Combined Rainfall and Boundary Inflows

```
# Coastal flooding with rainfall and tidal boundary
DEMfile       coastal_dem.asc
bcifile       tidal_boundary.bci
resroot       ./results/coastal
saveint       1800
sim_time      172800          # 2 days
fpfric        0.03
rain          coastal_storm.txt  # Rainfall time series
```

## Best Practices

1. **Match spatial resolution**:
   - Ensure rainfall grid resolution matches or can be resampled to your DEM
   - For dynamic rainfall, consider the appropriate spatial resolution for your application

2. **Temporal resolution**:
   - For flash flooding in small catchments: 5-15 minute intervals
   - For larger river basins: 1-hour intervals may suffice
   - Should be less than or equal to your output interval (saveint)

3. **Rainfall rates**:
   - Verify units carefully (m/s in LISFLOOD-FP)
   - Check that rates are reasonable for your region and event type
   - Consider applying a scaling factor for calibration

4. **When to use routing**:
   - Grid cell size > 50m
   - Large catchments where travel time is significant
   - When simulating both hillslope and channel processes

5. **Rainfall masks**:
   - Useful for urban drainage modeling
   - Can improve efficiency by limiting computation to areas of interest
   - Should be based on actual catchment boundaries

6. **Testing**:
   - Start with simple uniform rainfall before using complex inputs
   - For dynamic rainfall, visualize the NetCDF data before simulation
   - Check mass balance to ensure rainfall volumes are correctly applied

## Troubleshooting

| Issue | Possible Causes | Solution |
|-------|-----------------|----------|
| Rainfall not appearing in simulation | Wrong units | Convert to m/s |
| | File path issue | Check paths are correct |
| | NetCDF format incorrect | Verify variable names and dimensions |
| Mass balance errors | Boundary losses | Check if water is leaving the domain |
| | Infiltration losses | Check infiltration parameters |
| | Numerical errors | Reduce timestep (smaller initial_tstep) |
| Rainfall input too heavy/light | Unit conversion error | Verify conversion from mm/hr to m/s |
| | Scaling issue | Apply appropriate scaling factor |
| Model instability with rainfall | Timestep too large | Reduce initial_tstep |
| | Rainfall too intense | Consider routing or smoother rainfall input |
| Dynamic rainfall not matching time | Time units issue in NetCDF | Check "units" attribute for time variable |
| | Time offset issue | Ensure time starts at 0 for simulation start |

For assistance with rainfall data preparation, contact the LISFLOOD-FP development team or refer to the user manual and example files provided with the software.