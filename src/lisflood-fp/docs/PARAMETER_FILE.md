# LISFLOOD-FP Parameter File Guide

## Overview

The parameter file is the primary method for configuring LISFLOOD-FP simulations. It uses a simple key-value format with each parameter on a new line. Comments can be added using the `#` character.

This document provides a comprehensive reference for all available parameters in LISFLOOD-FP 8.2.

## Example Parameter File

Below is a fully annotated example parameter file with explanations for each parameter:

```
# LISFLOOD-FP Parameter File Example
# ==============================
# This is a comprehensive parameter file covering all common options

# --- Input Files ---
DEMfile         terrain.asc     # Digital Elevation Model
startfile       initial.asc     # Initial water depths (optional)
bcifile         boundary.bci    # Boundary condition file
weirfile        weirs.txt       # Weir structure definition file
riverfile       channel.river   # River channel definition file

# --- Output Control ---
resroot         ./results/sim1  # Output directory and prefix
dirroot         ./results       # Output directory
resname         sim1            # Output file prefix
saveint         3600            # Output interval (seconds)
massint         3600            # Mass balance calculation interval (seconds)
fpfric          0.035           # Floodplain Manning's n value
#elevoff                        # Turn off .elevation file output
#depthoff                       # Turn off .depth file output
#voutput                        # Include velocity output (.Vx, .Vy)
#startfile2d                    # Output initial 2D flow field (DG2/FV2)
#startelev                      # Output start elevation
#startq2d                       # Output start flows (DG2)
#adaptoff                       # Turn off adaptive time stepping
#binaryraster                   # Use binary raster format
netcdf_out                      # Output in NetCDF format

# --- Time Control ---
sim_time        86400           # Simulation duration (seconds)
initial_tstep   10.0            # Initial time step size (seconds)
solver_max_step 10.0            # Maximum solver time step (seconds)

# --- Solver Options ---
solver          ACC             # Solver type: ACC, FV1, FV2, DG2, ACC_NUGRID
cuda            1               # Enable CUDA acceleration (1=on, 0=off)
acceleration    1               # Use acceleration formulation (1=on, 0=off)
froude          1               # Enable Froude number limitation (1=on, 0=off)
hk_max          10              # Maximum depth for kinematic assumption (m)
krivodonovathreshold 10.0       # Krivodonova threshold (DG2)
epsilon         0.1             # AMR threshold (DG2)
courant         0.7             # CFL number for adaptive time stepping

# --- Sub-Grid Channel Model ---
SGCwidth        channel_width.asc   # Sub-grid channel width grid
SGCbank         bank_elev.asc       # Sub-grid channel bank elevation
SGCbed          bed_elev.asc        # Sub-grid channel bed elevation
SGCn            0.035               # Sub-grid channel Manning's n
SGCr            0.15                # Sub-grid channel multiplier
SGCp            0.78                # Sub-grid channel exponent
SGCm            1.0                 # Sub-grid channel meander coefficient
SGCchan         1                   # Sub-grid channel type (1=rectangular)

# --- Infiltration and Rainfall ---
rain            rainfall.asc        # Rainfall input file
rainfallmask    mask.asc            # Rainfall mask
dynamicrainfile rainfall_dyn.nc     # Dynamic rainfall NetCDF file
evaporation     evap.nc             # Evaporation input file
infinity        0.0001              # Infiltration rate (m/s)
dist_inf        inf_grid.asc        # Distributed infiltration grid

# --- Dam Operations ---
damfile         dams.txt            # Dam parameters file
dammask         dam_locs.asc        # Dam location mask

# --- Special Modes ---
diffusive                           # Enable diffusive mode
routing                             # Enable flow routing
steady                              # Enable steady-state checking
checkpoint      3600                # Checkpoint interval (seconds)
checkpointdir   ./checkpoint        # Checkpoint directory

# --- Misc Parameters ---
logfile         simulation.log      # Log file path
stagefile       gauges.stage        # Stage monitoring locations
gaugefile       gauges.txt          # Gauge locations file
profiles        1                   # Enable profiling (1=on, 0=off)
debug           0                   # Debug level (0-3)
verbose                             # Enable verbose output
```

## Parameter Descriptions

### Input Files

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `DEMfile` | Digital Elevation Model | Yes | - |
| `startfile` | Initial water depths | No | All dry |
| `bcifile` | Boundary conditions | No | No boundaries |
| `weirfile` | Weir definitions | No | No weirs |
| `riverfile` | River channel definitions | No | No rivers |
| `manningsn` | Manning's n grid | No | Global value |
| `SGCn` | Sub-grid channel Manning's n | No | Global value |
| `bdyfile` | Legacy boundary file | No | No boundaries |
| `porfile` | Porosity file | No | No porosity |

### Output Control

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `resroot` | Output path and prefix | Yes | - |
| `dirroot` | Output directory | No | Current dir |
| `resname` | Output file prefix | No | "res" |
| `saveint` | Output interval (seconds) | No | 3600 |
| `massint` | Mass balance interval (seconds) | No | 3600 |
| `elevoff` | Disable elevation output | No | Output on |
| `depthoff` | Disable depth output | No | Output on |
| `voutput` | Enable velocity output | No | Output off |
| `hazard` | Enable hazard output | No | Output off |
| `ascii_out` | Enable ASCII output | No | On |
| `binary_out` | Enable binary output | No | Off |
| `netcdf_out` | Enable NetCDF output | No | Off |

### Time Control

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `sim_time` | Simulation time (seconds) | Yes | 3600 |
| `initial_tstep` | Initial time step (seconds) | No | 10.0 |
| `solver_max_step` | Maximum solver time step | No | Input DEM cell size |
| `adaptoff` | Disable adaptive time stepping | No | Adaptive on |
| `checkpoint` | Checkpoint interval (seconds) | No | No checkpoints |
| `loadcheck` | Load checkpoint file | No | Start from scratch |

### Solver Options

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `solver` | Solver type (ACC, FV1, FV2, DG2, etc.) | No | ACC |
| `cuda` | Enable CUDA acceleration | No | Off |
| `acceleration` | Use acceleration solver | No | On |
| `froude` | Enable Froude number limitation | No | On |
| `courant` | CFL number | No | 0.7 |
| `theta` | Time-stepping parameter | No | 1.0 |
| `eps` | AMR threshold | No | 0.1 |
| `krivodonovathreshold` | DG2 slope limiter threshold | No | 10.0 |
| `diffusive` | Enable diffusive mode | No | Off |

### Sub-Grid Channel Model

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `SGCwidth` | Channel width grid | Yes (for SGC) | - |
| `SGCbank` | Bank elevation grid | Yes (for SGC) | - |
| `SGCbed` | Bed elevation grid | Yes (for SGC) | - |
| `SGCn` | Channel Manning's n | No | 0.035 |
| `SGCr` | Channel multiplier | No | 0.12 |
| `SGCp` | Channel exponent | No | 0.78 |
| `SGCm` | Meander coefficient | No | 1.0 |
| `SGCchan` | Channel type (1=rect, 2=trap, 3=para) | No | 1 |
| `SGClevee` | Enable levees | No | Off |

### Physical Parameters

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `fpfric` | Floodplain Manning's n | No | 0.035 |
| `gravity` | Gravitational acceleration | No | 9.8 |
| `depth_thresh` | Minimum depth threshold | No | 0.001 |
| `max_Froude` | Maximum Froude number | No | 10000 |

### Rainfall and Infiltration

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `rain` | Rainfall input file | No | No rain |
| `rainfallmask` | Rainfall mask | No | Uniform |
| `dynamicrainfile` | Dynamic rainfall NetCDF | No | No rain |
| `evaporation` | Evaporation input | No | No evaporation |
| `infinity` | Infiltration rate (m/s) | No | 0.0 |
| `dist_inf` | Distributed infiltration | No | Uniform |
| `routing` | Enable flow routing | No | Off |
| `routespeed` | Routing flow speed | No | 0.1 |

### Dam Operations

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `damfile` | Dam parameters file | No | No dams |
| `dammask` | Dam location mask | No | No dams |

### Monitoring

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `stagefile` | Stage monitoring locations | No | No monitoring |
| `gaugefile` | Gauge locations file | No | No gauges |
| `profiles` | Enable profiling | No | Off |
| `debug` | Debug level (0-3) | No | 0 |
| `logfile` | Log file path | No | No log |
| `verbose` | Enable verbose output | No | Off |

## Parameter File Format

- Each parameter should be on a separate line
- Format is `parameter_name parameter_value`
- Lines beginning with `#` are treated as comments
- Parameter names are case-insensitive
- Boolean flags (without values) enable the feature
- File paths can be relative or absolute

## Important Notes

1. **Required Parameters**:
   - At minimum, you must specify `DEMfile` and `resroot`
   - `sim_time` should be specified or defaults to 3600 seconds

2. **Parameter Order**:
   - Parameters can be specified in any order
   - Later parameters override earlier ones if duplicated

3. **File Paths**:
   - Use forward slashes (/) even on Windows
   - Consider using absolute paths for reliability
   - Maximum path length is 255 characters

4. **Numeric Values**:
   - Use decimal points for floating-point values
   - Scientific notation (e.g., 1.0e-6) is supported

5. **Special Characters**:
   - Avoid spaces in file paths
   - If needed, enclose paths in quotes: `"path with spaces"`

## Command Line Overrides

Most parameters can be overridden via command line using the parameter name with a hyphen prefix:

```bash
./lisflood simulation.par -saveint 1800 -fpfric 0.04
```

This overrides the `saveint` and `fpfric` values in the parameter file.

## Best Practices

1. **Organize by Section**:
   - Group related parameters together
   - Use comments to indicate sections

2. **Document Your Choices**:
   - Comment on non-standard parameter values
   - Note the rationale for specific choices

3. **Version Control**:
   - Keep parameter files in version control
   - Document changes between runs

4. **Validate Before Running**:
   - Check file paths and existence
   - Verify parameter values are in expected ranges