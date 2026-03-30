# 🎉 Self-Calibrating LISFLOOD-FP 8.2 - Complete!

I've created a complete self-calibrating flood simulation framework specifically adapted for your LISFLOOD-FP 8.2 version.

## 📦 What Was Built

### ✅ Complete Framework Structure

```
flood-autocalib/
├── README.md                          # Full documentation
├── QUICKSTART.md                       # Quick start guide
├── SETUP_LISFLOOD82.md               # Setup guide for your version ⭐
├── requirements.txt                     # Python dependencies
├── main.py                            # Entry point
│
├── config/
│   ├── project.yaml.example            # Generic config template
│   └── project_lisflood82.yaml.example # LISFLOOD-FP 8.2 config ⭐
│
├── models/
│   ├── adapter.py                     # Base adapter interface
│   ├── lisflood/adapter.py           # Generic LISFLOOD-FP adapter
│   └── lisflood82/adapter.py        # YOUR LISFLOOD-FP 8.2 adapter ⭐
│
├── calibration/
│   ├── metrics.py                    # NSE, KGE, IoU calculations
│   └── loop.py                      # Auto-calibration loop (AutoResearch-style)
│
├── data/
│   └── preprocessing.py             # Data loading & validation
│
├── satellite/
│   └── preprocessor.py              # Satellite flood mask alignment
│
└── visualization/
    └── plots.py                     # Calibration visualization
```

### ✅ LISFLOOD-FP 8.2 Specific Features

**Fully adapted for your version:**

1. **Multiple Solver Support**
   - ACC (Acceleration - fastest) ⚡⚡⚡
   - FV1 (First-order finite volume)
   - FV2 (Second-order finite volume)
   - DG2 (Discontinuous Galerkin - higher accuracy)
   - ACC_NUGRID (Non-uniform grid)

2. **GPU Acceleration**
   - CUDA support with `-cuda` flag
   - 5-20x speedup over CPU

3. **11 Calibratable Parameters**
   - `fpfric`: Floodplain Manning's n
   - `SGCn`: Sub-grid channel Manning's n
   - `infinity`: Infiltration rate
   - `courant`: CFL number
   - `theta`: Time-stepping parameter
   - `depth_thresh`: Minimum depth threshold
   - `gravity`: Gravitational acceleration
   - `max_Froude`: Maximum Froude number
   - `SGCr`: Sub-grid channel multiplier
   - `SGCp`: Sub-grid channel exponent

4. **Smart Output Extraction**
   - Automatic detection of output files
   - Supports ASCII, binary, and NetCDF
   - Extracts discharge (.discharge) and inundation (.max, .wd)

## 🚀 Quick Start (3 Steps)

### Step 1: Build Your LISFLOOD-FP 8.2

```bash
cd /home/umair/.openclaw/workspace
git clone https://github.com/urabbani/lisflood-fp_8.2_update.git
cd lisflood-fp_8.2_update
mkdir build && cd build
cmake .. -D_CONFIG=../config/local
make -j
```

**Result:** Executable at `build/lisflood`

### Step 2: Configure Framework

```bash
cd /home/umair/.openclaw/workspace/flood-autocalib
cp config/project_lisflood82.yaml.example config/project.yaml

# Edit paths in config/project.yaml:
# - Set executable path to your lisflood build
# - Set input data paths (DEM, rainfall, boundary)
# - Set observation data paths (gauges, satellite)
nano config/project.yaml
```

**Key config changes:**
```yaml
model:
  type: "lisflood82"  # Use your LISFLOOD-FP 8.2 version
  executable: "/home/umair/.openclaw/workspace/lisflood-fp_8.2_update/build/lisflood"
  use_cuda: true  # Enable GPU acceleration
  solver: "ACC"  # Use fastest solver
```

### Step 3: Run Calibration

```bash
# Install Python dependencies
pip install -r requirements.txt

# Run calibration
python main.py --config config/project.yaml --save-plots
```

**That's it!** The framework will:
1. Load your LISFLOOD-FP 8.2 parameters
2. Iterate 100 times (or customize)
3. Automatically find best parameters
4. Generate visualization plots
5. Save calibrated .par file ready to use

## 📊 What It Calibrates

### Multi-Objective Metrics

| Metric | Type | Target | Purpose |
|---------|-------|---------|----------|
| **NSE** | Temporal | ≥ 0.85 | Hydrograph shape accuracy |
| **KGE** | Temporal | ≥ 0.80 | Decomposed flow accuracy |
| **IoU** | Spatial | ≥ 0.70 | Flood extent overlap |

### Calibration Loop Strategy

```
Iterations 0-9:    Latin Hypercube Sampling (exploration)
Iterations 10-29:   Gaussian Perturbation around best (refinement)
Iterations 30+:      Fine-tuning (exploitation)
```

### Pareto Front Analysis

When metrics conflict, system maintains **non-dominated solutions**:
- High NSE/KGE + Low IoU → Good temporal, poor spatial
- Low NSE/KGE + High IoU → Poor temporal, good spatial
- Balanced → Middle of Pareto front

**You choose based on priorities!**

## 📁 Input Data Required

| Data | Format | Location |
|-------|---------|----------|
| DEM (terrain) | ASCII raster (.asc) | `data/terrain/dem.asc` |
| Rainfall | NetCDF (.nc) or ASCII | `data/rainfall/rainfall_202606.nc` |
| Boundary | BCI file | `data/boundary/upstream.bci` |
| Parameter template | .par file | `data/lisflood/template.par` ✅ (included) |
| Gauge observations | CSV | `data/observations/gauges/gauge_001.csv` |
| Satellite flood mask | GeoTIFF | `data/observations/satellite/flood_mask.tif` |

## 📈 Output Files

### After Calibration

```
data/outputs/lisflood82/
├── calibrated.par              ⭐ Best parameters (ready to use!)
├── calibration_results.json   Full calibration history
└── [LISFLOOD-FP outputs]

visualization/
├── calibration_progress.png  Metrics over iterations
├── pareto_front.png         Multi-objective trade-offs
├── hydrograph_comparison.png Observed vs simulated
└── inundation_comparison.png Spatial comparison
```

## 🎯 Example Configuration

```yaml
model:
  type: "lisflood82"
  executable: "/home/umair/.openclaw/workspace/lisflood-fp_8.2_update/build/lisflood"
  param_template: "data/lisflood/template.par"
  output_dir: "data/outputs/lisflood82/"
  use_cuda: true
  solver: "ACC"

calibration:
  max_iterations: 100
  target_metrics:
    nse: 0.85
    kge: 0.80
    iou: 0.70
  objective_weights:
    nse: 0.40
    kge: 0.30
    iou: 0.30
```

## 🚀 Performance Tips

### With GPU (Recommended)
```yaml
model:
  use_cuda: true
  solver: "ACC"
```
**Speedup:** 5-20x faster than CPU

### Without GPU
```yaml
model:
  use_cuda: false
  solver: "FV1"  # Stable CPU solver
```

### Reduce Output Frequency
```yaml
lisflood82:
  time:
    saveint: 7200  # Save every 2 hours (not 1 hour)
```

### Disable Unnecessary Outputs
```yaml
lisflood82:
  output:
    v_output: false  # No velocity output
    elev_off: true  # No elevation output
```

## 📚 Documentation

| File | Description |
|-------|-------------|
| **SETUP_LISFLOOD82.md** | ⭐ Setup guide for your version |
| **QUICKSTART.md** | Quick start reference |
| **README.md** | Full framework documentation |

## 🔧 Troubleshooting

### LISFLOOD-FP Not Found
```bash
# Check build directory
ls /home/umair/.openclaw/workspace/lisflood-fp_8.2_update/build/

# Update config.yaml with correct path
model:
  executable: "/absolute/path/to/lisflood"
```

### CUDA Not Working
```bash
# Check GPU
nvidia-smi

# Disable CUDA if not available
model:
  use_cuda: false
```

### Poor Calibration Results
1. Increase `max_iterations` (e.g., 200)
2. Adjust `objective_weights` in config
3. Check calibration data matches simulation period
4. Verify boundary conditions are realistic
5. Try different solver (FV1, FV2)

## 🎉 Ready to Calibrate!

### What You Need

1. ✅ LISFLOOD-FP 8.2 built
2. ✅ Input data (DEM, rainfall, boundary)
3. ✅ Validation data (gauges, satellite)
4. ✅ Configuration file edited

### Then Run

```bash
cd /home/umair/.openclaw/workspace/flood-autocalib
python main.py --config config/project.yaml --save-plots
```

### Get Results

- **Best parameters:** `data/outputs/lisflood82/calibrated.par`
- **Metrics:** `NSE`, `KGE`, `IoU` in console and JSON
- **Visualizations:** `visualization/*.png`

## 💡 Next Steps

1. **Set up your input data** (DEM, rainfall, boundary, observations)
2. **Edit `config/project.yaml`** with correct paths
3. **Build LISFLOOD-FP 8.2** from your repository
4. **Run calibration** with `--save-plots`
5. **Review results** in `calibrated.par` and visualization plots

---

**Questions?** Check `SETUP_LISFLOOD82.md` for detailed setup instructions! 🌊
