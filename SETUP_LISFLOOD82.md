# Setup Guide for LISFLOOD-FP 8.2

This guide will help you set up the Self-Calibrating LISFLOOD-FP framework with your specific version: **urabbani/lisflood-fp_8.2_update**

## 🚀 Quick Setup

### 1. Clone Your LISFLOOD-FP Repository

```bash
# Clone your LISFLOOD-FP 8.2 version
cd /home/umair/.openclaw/workspace
git clone https://github.com/urabbani/lisflood-fp_8.2_update.git
cd lisflood-fp_8.2_update

# Build LISFLOOD-FP 8.2
mkdir build && cd build
cmake .. -D_CONFIG=../config/local
make -j  # Build with all cores
```

**Note:** After building, the executable will be at `build/lisflood`

### 2. Copy Flood-AutoCalib to Workspace

```bash
# Navigate to workspace
cd /home/umair/.openclaw/workspace

# If you just created the framework, it's already here
# Otherwise, clone/download it
```

### 3. Create Configuration

```bash
# Copy LISFLOOD-FP 8.2 specific config
cp config/project_lisflood82.yaml.example config/project.yaml

# Edit with your paths
nano config/project.yaml
```

### 4. Edit Configuration File

Update these key paths in `config/project.yaml`:

```yaml
model:
  type: "lisflood82"  # Use LISFLOOD-FP 8.2 adapter
  executable: "/home/umair/.openclaw/workspace/lisflood-fp_8.2_update/build/lisflood"
  param_template: "data/lisflood/template.par"
  output_dir: "data/outputs/lisflood82/"
  
  # LISFLOOD-FP 8.2 options
  use_cuda: true  # Enable GPU acceleration
  solver: "ACC"  # Use ACC solver (fastest)
```

### 5. Prepare Input Data

#### Required Files

1. **Digital Elevation Model (DEM)**
   ```bash
   # Format: ASCII raster (.asc)
   # Location: data/terrain/dem.asc
   ```

2. **Rainfall Data**
   ```bash
   # Options:
   # 1. Static ASCII: rainfall.asc
   # 2. Dynamic NetCDF: rainfall.nc (recommended)
   # Location: data/rainfall/rainfall_202606.nc
   ```

3. **Boundary Conditions**
   ```bash
   # Format: BCI file (LISFLOOD-FP format)
   # Location: data/boundary/upstream.bci
   
   # BCI file format (example):
   # TIME    FLOW
   # 0        1000
   # 3600     1250
   # 7200     1500
   # ...
   ```

4. **LISFLOOD-FP Parameter Template**
   ```bash
   # Create a basic .par template
   # Location: data/lisflood/template.par
   
   # Example template.par:
   DEMfile ../terrain/dem.asc
   resroot ./results/sim1
   nsim_time 604800
   fpfric 0.035
   infinity 0.0001
   courant 0.7
   theta 1.0
   depth_thresh 0.001
   saveint 3600
   massint 3600
   ```

#### Validation Data

5. **Gauge Observations**
   ```bash
   # Format: CSV
   # Location: data/observations/gauges/gauge_001.csv
   
   # CSV format:
   timestamp,discharge
   2026-06-01T00:00:00Z,1250.5
   2026-06-01T01:00:00Z,1265.2
   ...
   ```

6. **Satellite Inundation**
   ```bash
   # Format: GeoTIFF (binary mask or depth grid)
   # Location: data/observations/satellite/flood_mask.tif
   ```

## 📁 Project Structure (Complete)

```
flood-autocalib/
├── config/
│   ├── project.yaml                        # Your main configuration
│   └── project_lisflood82.yaml.example    # LISFLOOD-FP 8.2 template
├── data/
│   ├── rainfall/
│   │   └── rainfall_202606.nc           # Dynamic rainfall (NetCDF)
│   ├── terrain/
│   │   └── dem.asc                      # Digital Elevation Model
│   ├── boundary/
│   │   └── upstream.bci                 # Boundary conditions
│   ├── lisflood/
│   │   ├── template.par                 # LISFLOOD-FP parameter template
│   │   └── gauges.txt                   # Gauge monitoring points
│   ├── observations/
│   │   ├── gauges/
│   │   │   └── gauge_001.csv         # Discharge observations
│   │   └── satellite/
│   │       └── flood_mask.tif          # Flood extent mask
│   └── outputs/
│       └── lisflood82/                # Model outputs
├── models/
│   ├── lisflood82/
│   │   ├── adapter.py                  # LISFLOOD-FP 8.2 adapter
│   │   └── __init__.py
│   └── ...
└── ...
```

## 🔧 LISFLOOD-FP 8.2 Features

### Supported Solvers

| Solver | Description | Best Use Case | Speed |
|---------|-------------|----------------|-------|
| **ACC** | Acceleration solver (simplified momentum) | Large domains | ⚡⚡⚡ Fastest |
| **FV1** | First-order finite volume | General purpose | ⚡⚡ Fast |
| **FV2** | Second-order finite volume | Higher accuracy | ⚡ Medium |
| **DG2** | Discontinuous Galerkin (2nd order) | Supercritical flows | 🐢 Slowest |
| **ACC_NUGRID** | Non-uniform grid acceleration | Variable resolution | ⚡⚡ Fast |

### GPU Acceleration

Enable CUDA in configuration:

```yaml
model:
  use_cuda: true  # Requires NVIDIA GPU with CUDA support
```

**Speedup:** 5-20x faster than CPU (depends on grid size)

### Calibratable Parameters

The framework can automatically calibrate these LISFLOOD-FP 8.2 parameters:

| Parameter | Description | Range | Unit |
|-----------|-------------|--------|-------|
| `fpfric` | Floodplain Manning's n | 0.01-0.1 | s/m^(1/3) |
| `SGCn` | Sub-grid channel Manning's n | 0.02-0.06 | s/m^(1/3) |
| `infinity` | Infiltration rate | 0.0-0.001 | m/s |
| `courant` | CFL number | 0.3-0.9 | - |
| `theta` | Time-stepping parameter | 0.5-1.0 | - |
| `depth_thresh` | Minimum depth threshold | 0.0005-0.01 | m |
| `gravity` | Gravitational acceleration | 9.7-9.9 | m/s² |
| `max_Froude` | Maximum Froude number | 1.0-10000 | - |
| `SGCr` | Sub-grid channel multiplier | 0.1-0.2 | - |
| `SGCp` | Sub-grid channel exponent | 0.7-0.85 | - |

## 🚀 Running Calibration

### Basic Calibration (CPU)

```bash
python main.py --config config/project.yaml
```

### GPU-Accelerated Calibration

```bash
# Make sure CUDA is enabled in config.yaml
python main.py --config config/project.yaml
```

### With Plots

```bash
python main.py --config config/project.yaml --save-plots
```

### Override Iterations

```bash
python main.py --config config/project.yaml --max-iterations 200 --save-plots
```

## 📊 Output Files

After calibration completes:

```
data/outputs/lisflood82/
├── simulation.par              # Last simulation parameters
├── calibrated.par             # Best calibrated parameters
├── calibration_results.json   # Full calibration history
├── simulation.max            # Maximum water depths (ASCII)
├── simulation.wd             # Water depths at intervals
├── simulation.elev           # Water surface elevations
├── simulation.discharge       # Discharge at gauge points
└── simulation.mass           # Mass balance information

visualization/
├── calibration_progress.png  # Metrics over iterations
├── pareto_front.png         # Multi-objective trade-offs
├── hydrograph_comparison.png # Observed vs simulated discharge
└── inundation_comparison.png # Observed vs simulated flooding
```

## 🎯 Calibration Results

### Metrics Interpretation

- **NSE (Nash-Sutcliffe Efficiency)**: 
  - ≥ 0.85: Excellent
  - 0.65-0.85: Good
  - < 0.65: Poor

- **KGE (Kling-Gupta Efficiency)**:
  - ≥ 0.80: Excellent
  - 0.60-0.80: Good
  - < 0.60: Poor

- **IoU (Intersection over Union)**:
  - ≥ 0.70: Excellent spatial accuracy
  - 0.50-0.70: Moderate accuracy
  - < 0.50: Poor spatial accuracy

### Choosing from Pareto Front

The calibration produces a Pareto front of non-dominated solutions. Choose based on your priorities:

1. **Peak flow accuracy** → Choose solution with highest NSE/KGE
2. **Flood extent accuracy** → Choose solution with highest IoU
3. **Balanced** → Choose middle of Pareto front

## 🔍 Troubleshooting

### LISFLOOD-FP Build Errors

```bash
# Check CMake version
cmake --version  # Should be 3.10+

# Check CUDA
nvcc --version  # Should show CUDA version

# Build with verbose output
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j VERBOSE=1
```

### CUDA Not Found

```bash
# Check GPU
nvidia-smi

# Check CUDA installation
nvcc --version

# Disable CUDA in config.yaml if not available
model:
  use_cuda: false
```

### Simulation Crashes

**Symptom:** LISFLOOD-FP exits with error

**Possible causes:**
1. DEM has NaN values or incorrect CRS
2. Boundary condition file has wrong format
3. Time step too large (reduce `courant` parameter)
4. Memory limits (use 64-bit version)

**Solutions:**
1. Validate DEM: `gdalinfo data/terrain/dem.asc`
2. Check BCI file format
3. Reduce CFL number: `courant: 0.5`
4. Increase swap or use GPU

### No Output Files

**Check:**
1. Output directory permissions
2. `resroot` parameter in .par file
3. LISFLOOD-FP return code

**Debug:**
```bash
# Run with verbose mode
cd data/outputs/lisflood82/
/home/umair/.openclaw/workspace/lisflood-fp_8.2_update/build/lisflood -v simulation.par
```

### Poor Calibration Results

**Symptoms:** Low NSE/KGE/IoU values

**Solutions:**
1. Increase `max_iterations` (e.g., 200)
2. Adjust `objective_weights` in config
3. Check that calibration data matches simulation period
4. Verify boundary conditions are realistic
5. Try different solver (FV1, FV2)
6. Check DEM quality (sinks, artifacts)

## 🚀 Performance Tips

### Use GPU Acceleration

```yaml
model:
  use_cuda: true  # 5-20x faster
```

### Choose Fast Solver

```yaml
model:
  solver: "ACC"  # Fastest for large domains
```

### Reduce Output Frequency

```yaml
lisflood82:
  output:
    saveint: 7200  # Save every 2 hours (instead of 1 hour)
```

### Disable Unnecessary Outputs

```yaml
lisflood82:
  output:
    v_output: false  # No velocity output
    elev_off: true  # No elevation output
```

## 📚 Additional Resources

- [LISFLOOD-FP User Manual](https://github.com/urabbani/lisflood-fp_8.2_update/blob/main/docs/USER_MANUAL.md)
- [Parameter File Guide](https://github.com/urabbani/lisflood-fp_8.2_update/blob/main/docs/PARAMETER_FILE.md)
- [Installation Guide](https://github.com/urabbani/lisflood-fp_8.2_update/blob/main/docs/INSTALL.md)

## 🤝 Getting Help

1. Check logs: `logs/calibration82.log`
2. Verify all paths in `config/project.yaml`
3. Test LISFLOOD-FP manually:
   ```bash
   cd data/outputs/lisflood82/
   /path/to/lisflood -v simulation.par
   ```
4. Check LISFLOOD-FP GitHub issues for known problems

---

Ready to calibrate! 🌊
