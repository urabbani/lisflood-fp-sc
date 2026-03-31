# LISFLOOD-FP 8.2 Quick Reference Guide

This quick reference guide provides examples for common LISFLOOD-FP tasks and use cases. For detailed information, refer to the [User Manual](USER_MANUAL.md).

## Table of Contents

1. [Basic Commands](#basic-commands)
2. [Common Parameters](#common-parameters)
3. [Example Command Lines](#example-command-lines)
4. [Solver Selection Guide](#solver-selection-guide)
5. [Output Files Reference](#output-files-reference)
6. [Common Workflows](#common-workflows)
7. [Tips and Best Practices](#tips-and-best-practices)

## Basic Commands

### Running a Simulation

```bash
# Basic run
./lisflood simulation.par

# Run with verbose output
./lisflood -v simulation.par

# Run with CUDA acceleration
./lisflood -cuda simulation.par

# Override parameters from command line
./lisflood simulation.par -fpfric 0.045 -sim_time 3600
```

### Utility Commands

```bash
# Display version information
./lisflood -version

# Run without optimizing wet/dry area calculation
./lisflood -noclip simulation.par

# Allow overwriting of existing results
./lisflood -overwrite simulation.par

# Load a checkpoint file
./lisflood simulation.par -loadcheck checkpoint.chk
```

## Common Parameters

### Essential Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `DEMfile` | Digital Elevation Model | `DEMfile    terrain.asc` |
| `resroot` | Results directory and prefix | `resroot    ./results/sim1` |
| `sim_time` | Simulation duration (seconds) | `sim_time    86400` |
| `saveint` | Output save interval (seconds) | `saveint     3600` |
| `fpfric` | Manning's n value | `fpfric      0.035` |

### Solver Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `solver` | Solver type | `solver      ACC` |
| `acceleration` | Use acceleration solver | `acceleration  1` |
| `cuda` | Enable CUDA | `cuda        1` |
| `courant` | CFL number | `courant     0.7` |

### Output Control

| Parameter | Description | Example |
|-----------|-------------|---------|
| `voutput` | Enable velocity output | `voutput` |
| `hazard` | Enable hazard output | `hazard` |
| `binary_out` | Binary output format | `binary_out` |
| `netcdf_out` | NetCDF output format | `netcdf_out` |

## Example Command Lines

### Basic Floodplain Simulation

```bash
./lisflood -v -cuda examples/simple_floodplain.par
```

### Dam Breach Analysis

```bash
./lisflood -v examples/dam_breach.par -sim_time 21600 -saveint 300
```

### Urban Flood with Rainfall

```bash
./lisflood -v -cuda examples/urban_flood.par -fpfric 0.045
```

### Levee Breach Simulation

```bash
./lisflood -v examples/weir_example.par -solver FV1
```

### Checkpointing Example

```bash
# Run initial simulation with checkpointing
./lisflood -v simulation.par -checkpoint 3600

# Continue from checkpoint
./lisflood -v simulation.par -loadcheck checkpoint_000001.chk
```

## Solver Selection Guide

| Solver | Best For | Features | Performance |
|--------|----------|----------|-------------|
| `ACC` | Large domains | Local inertial approximation | Fastest |
| `FV1` | General purpose | Full shallow water equations | Fast |
| `FV2` | Steep gradients | 2nd-order accuracy | Medium |
| `DG2` | Complex flows | Higher-order accuracy | Slower |
| `ACC_NUGRID` | Variable resolution | Adaptive mesh refinement | Fast |

### Recommendation Guidelines

1. **For large domains (>10^6 cells)**: 
   - Use `ACC` with CUDA
   - Example: `solver ACC` + `cuda 1`

2. **For accurate flow velocities**:
   - Use `FV1` or `FV2`
   - Example: `solver FV1` + `acceleration 0`

3. **For dam breaks or supercritical flows**:
   - Use `FV2` or `DG2`
   - Example: `solver DG2` + `cuda 1`

4. **For complex urban environments**:
   - Use `ACC_NUGRID` for efficiency
   - Example: `solver ACC_NUGRID` + `epsilon 0.1`

## Output Files Reference

| Extension | Description | Format | Parameter |
|-----------|-------------|--------|-----------|
| `.max` | Maximum water depths | ASCII/binary | Always produced |
| `.wd` | Water depths at intervals | ASCII/binary | Toggle with `depthoff` |
| `.elev` | Water surface elevations | ASCII/binary | Toggle with `elevoff` |
| `.Vx/.Vy` | Velocity components | ASCII/binary | Enable with `voutput` |
| `.hazard` | Hazard rating | ASCII/binary | Enable with `hazard` |
| `.mass` | Mass balance info | Text | Always produced |
| `.stage` | Water levels at points | Text | Enable with `stagefile` |
| `.check` | Checkpoint files | Binary | Enable with `checkpoint` |

## Common Workflows

### Flood Inundation Mapping

1. Prepare DEM and boundary conditions
2. Run simulation with ACC solver:
   ```bash
   ./lisflood -v -cuda inundation.par
   ```
3. Process maximum depths:
   ```bash
   # Convert to GeoTIFF (using GDAL)
   gdal_translate -a_srs EPSG:27700 results/inundation.max results/inundation_max.tif
   ```

### Dam Breach Analysis

1. Set up DEM with dam and initial water level
2. Define breach using weir file
3. Run with FV1 or FV2 solver:
   ```bash
   ./lisflood -v dam_breach.par -solver FV2
   ```
4. Analyze arrival times and maximum depths

### Levee Scenario Testing

1. Create base DEM with levees
2. Define potential breach points as weirs
3. Run multiple scenarios with different weir parameters
4. Compare results using:
   ```bash
   python postprocess/diff.py results/scenario1.max results/scenario2.max
   ```

## Tips and Best Practices

### Performance Optimization

1. **Use GPU acceleration** where available:
   ```bash
   ./lisflood -cuda simulation.par
   ```

2. **Optimize domain size**:
   - Clip DEM to area of interest
   - Use appropriate cell size for required detail

3. **Balance accuracy and speed**:
   - ACC solver for rapid assessments
   - FV1/FV2/DG2 for detailed analysis

### Stability Considerations

1. **Adapt CFL number** for challenging topography:
   ```bash
   ./lisflood simulation.par -courant 0.5
   ```

2. **Use checkpointing** for long simulations:
   ```bash
   ./lisflood simulation.par -checkpoint 3600
   ```

3. **Monitor mass balance** for simulation health:
   - Check `.mass` file for large errors
   - Reduce time step if mass errors are significant

### Data Management

1. **Use binary output** for large simulations:
   ```
   binary_out
   ```

2. **Organize results** with descriptive prefixes:
   ```
   resroot ./results/scenario1_Q100
   ```

3. **Save key parameters** in filenames:
   ```
   resname sim_n035_Q100
   ```

---

This quick reference guide covers the most common LISFLOOD-FP operations. For comprehensive information, refer to the [User Manual](USER_MANUAL.md).