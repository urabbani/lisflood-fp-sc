# LISFLOOD-FP 8.2 User Manual

## Table of Contents

1. [Introduction](#introduction)
2. [Model Overview](#model-overview)
3. [Getting Started](#getting-started)
4. [Input Files](#input-files)
5. [Running Simulations](#running-simulations)
6. [Output Files](#output-files)
7. [Solver Options](#solver-options)
8. [Advanced Features](#advanced-features)
9. [Model Parameters](#model-parameters)
10. [Troubleshooting](#troubleshooting)
11. [References](#references)

## Introduction

LISFLOOD-FP is a two-dimensional hydrodynamic model designed for simulating floodplain inundation over complex topography. It was originally developed at the University of Bristol and has been continuously enhanced to include advanced solvers and features.

Version 8.2 includes significant improvements such as:
- Enhanced CUDA-based solvers for GPU acceleration
- Dynamic grid adaptation capabilities
- Improved dam operations
- Comprehensive hazard assessment
- Advanced weir functionality

## Model Overview

LISFLOOD-FP solves the shallow water equations (SWE) using a variety of numerical schemes:

| Solver | Description | Best Use Case |
|--------|-------------|---------------|
| ACC | Acceleration solver (simplified momentum) | Large domains with gradual flow changes |
| FV1 | First-order finite volume | General purpose with good stability |
| FV2 | Second-order finite volume | Higher accuracy needs with smooth flows |
| DG2 | Discontinuous Galerkin (2nd order) | Supercritical flows and shock capturing |
| ACC_NUGRID | Non-uniform grid acceleration | Domains requiring variable resolution |

The model operates on a raster grid and can simulate:
- Fluvial (river) flooding
- Pluvial (rainfall) flooding
- Coastal flooding
- Dam breaks and levee failures
- Urban flooding scenarios

## Getting Started

### System Requirements

- 64-bit operating system (Linux, Windows)
- Modern multi-core CPU
- 8+ GB RAM (16+ GB recommended for large domains)
- NVIDIA GPU with CUDA support (optional, for acceleration)
- Storage space for input/output data

### Basic Workflow

1. Prepare your terrain data (DEM) and other inputs
2. Create a parameter file (.par)
3. Run LISFLOOD-FP with your parameter file
4. Analyze the results using the output files

## Input Files

LISFLOOD-FP requires various input files depending on your simulation requirements:

### Essential Files

| File Type | Description | Format |
|-----------|-------------|--------|
| DEM file | Digital Elevation Model | ASCII raster or binary |
| Parameter file | Simulation configuration | Text file (.par) |

### Optional Files

| File Type | Description | Format |
|-----------|-------------|--------|
| Boundary condition | Inflows/outflows | Text file (.bci) |
| Initial conditions | Starting water depths | ASCII raster |
| River network | Channel definitions | Text file (.river) |
| Rainfall | Precipitation data | ASCII/binary/NetCDF |
| Manning's n | Roughness values | ASCII raster |
| Weir data | Structure definitions | Text file |
| Dam data | Dam operations | Text file |

### Parameter File Structure

The parameter file (.par) controls all aspects of your simulation. It consists of key-value pairs:

```
# Comments start with #
DEMfile         path/to/dem.asc    # Digital Elevation Model
resroot         results/sim1        # Output location prefix
manningn        0.035              # Global Manning's n value
simulate_time   86400              # Simulation time in seconds
massint         3600               # Mass balance interval
saveint         3600               # Results save interval
```

A full example parameter file with explanations is available in the [PARAMETER_FILE.md](PARAMETER_FILE.md) document.

## Running Simulations

### Basic Command Line Usage

```bash
./lisflood -v parameter_file.par
```

Where:
- `-v` enables verbose output
- `parameter_file.par` is your parameter file

### Command Line Options

| Option | Description |
|--------|-------------|
| `-v` | Verbose mode |
| `-cuda` | Enable CUDA acceleration |
| `-noclip` | Disable optimization of wet/dry regions |
| `-overwrite` | Allow overwriting of existing results |
| `-version` | Display version information |

### Example

```bash
./lisflood -v -cuda simulation.par
```

## Output Files

LISFLOOD-FP produces various output files based on your configuration:

### Standard Outputs

| Extension | Description | Format |
|-----------|-------------|--------|
| .max | Maximum water depths | ASCII/binary |
| .wd | Water depths at save intervals | ASCII/binary |
| .elev | Water surface elevations | ASCII/binary |
| .mass | Mass balance information | Text file |
| .log | Simulation log | Text file |

### Optional Outputs

| Extension | Description | Format |
|-----------|-------------|--------|
| .Vx/.Vy | Velocity components | ASCII/binary |
| .speed | Flow speed | ASCII/binary |
| .hazard | Flood hazard rating | ASCII/binary |
| .stage | Water stage at gauge points | Text file |
| .discharge | Flow discharge at monitoring points | Text file |
| .checkpoint | Simulation state for restart | Binary |

## Solver Options

LISFLOOD-FP 8.2 includes several numerical solvers that can be selected in the parameter file:

### Acceleration Solver (ACC)
```
solver         ACC
```
- Fast and efficient for large domains
- Uses the local inertial approximation
- Good balance of performance and accuracy

### Finite Volume Solvers (FV1/FV2)
```
solver         FV1   # or FV2
```
- FV1: First-order scheme with good stability
- FV2: Second-order scheme with higher accuracy
- Solve the full shallow water equations

### Discontinuous Galerkin (DG2)
```
solver         DG2
```
- Higher-order accuracy
- Better handling of steep gradients and shocks
- More computationally intensive

### Non-uniform Grid (ACC_NUGRID)
```
solver         ACC_NUGRID
```
- Variable resolution based on topography and flow features
- Computational efficiency for large domains
- Based on the acceleration solver

## Advanced Features

### GPU Acceleration
Enable CUDA acceleration for significantly faster simulations:
```
cuda           1       # Enable CUDA
```

### Sub-grid Channels
Model river channels below the DEM resolution:
```
SGCwidth       subgrid_width.asc   # Channel width file
SGCbank        subgrid_bank.asc    # Bank elevation file
SGCbed         subgrid_bed.asc     # Bed elevation file
```

### Dynamic Rainfall
Apply spatially and temporally varying rainfall:
```
rainfall       rainfall.nc         # NetCDF rainfall file
rainfallmask   mask.asc            # Optional mask file
```

### Weirs and Structures
Model hydraulic structures and levee breaches:
```
weirfile       weir_data.txt       # Weir definitions
```

### Dam Operations
Simulate dam operations and reservoir routing:
```
damfile        dam_params.txt      # Dam parameters
dammask        dam_mask.asc        # Dam locations
```

## Model Parameters

### Physical Parameters

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `manningn` | Global Manning's roughness | 0.035 | s/m^(1/3) |
| `gravity` | Gravitational acceleration | 9.8 | m/s² |
| `depth_thresh` | Minimum depth threshold | 0.001 | m |
| `max_Froude` | Maximum Froude number | 10000 | - |

### Numerical Parameters

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `cfl` | Courant-Friedrichs-Lewy number | 0.7 | - |
| `theta` | Time-stepping parameter | 1.0 | - |
| `eps` | Adaptive grid threshold | 0.1 | - |
| `max_step` | Maximum time step | 10.0 | s |

### Output Control

| Parameter | Description | Default | Units |
|-----------|-------------|---------|-------|
| `saveint` | Output interval | 3600 | s |
| `massint` | Mass balance interval | 3600 | s |
| `overwrite` | Overwrite existing files | 0 (off) | - |
| `ascii_out` | ASCII output format | 1 (on) | - |
| `binary_out` | Binary output format | 0 (off) | - |

## Troubleshooting

### Common Issues

| Problem | Possible Cause | Solution |
|---------|----------------|----------|
| Simulation crashes | Memory limitations | Reduce domain size or use 64-bit version |
| Unstable results | Time step too large | Reduce CFL number |
| Slow performance | Large grid size | Enable CUDA or use coarser resolution |
| Water appears stuck | Poor DEM quality | Check for sinks and connectivity |
| Mass balance errors | Boundary conditions | Check BCI file and inflow definitions |

### Diagnostics

Use the verbose (-v) option to get detailed information during the simulation:
```bash
./lisflood -v simulation.par
```

Check the log file and mass balance file for warnings and errors.

## References

1. Bates, P.D., Horritt, M.S., Fewtrell, T.J., 2010. A simple inertial formulation of the shallow water equations for efficient two-dimensional flood inundation modelling. Journal of Hydrology 387, 33-45.
2. Neal, J., Schumann, G., Bates, P.D., 2012. A subgrid channel model for simulating river hydraulics and floodplain inundation over large and data sparse areas. Water Resources Research 48, W11506.
3. Shaw, J., Kesserwani, G., Neal, J., Bates, P.D., Sharifian, M.K., 2021. LISFLOOD-FP 8.0: the new discontinuous Galerkin shallow-water solver for multi-core CPUs and GPUs. Geoscientific Model Development 14, 3577–3602.

---

For additional information and updates, please visit the [LISFLOOD-FP repository](https://github.com/your-organization/lisflood-fp) or contact the development team.