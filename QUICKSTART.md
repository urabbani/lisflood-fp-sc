# Quick Start Guide

## Installation

```bash
# Navigate to project directory
cd /mnt/d/flood-autocalib

# Install dependencies
pip install -r requirements.txt
```

## Project Setup

```bash
# Copy example configuration
cp config/project.yaml.example config/project.yaml

# Edit configuration with your data paths
nano config/project.yaml
```

## Configuration Checklist

Edit `config/project.yaml` with your paths:

### Required Files

1. **LISFLOOD-FP executable**
   ```yaml
   model:
     executable: "/path/to/lisflood-fp"
   ```

2. **Input data**
   - Rainfall: `data/rainfall/rainfall_202606.nc`
   - Terrain (DEM): `data/terrain/dem.tif`
   - Boundary discharge: `data/boundary/upstream_discharge.csv`

3. **LISFLOOD-FP template**
   - `data/lisflood/template.par` (template parameter file)

4. **Validation data**
   - Gauge observations: `data/observations/gauges/gauge_001.csv`
   - Satellite inundation: `data/observations/satellite/flood_mask.tif`

## Directory Structure

```bash
flood-autocalib/
├── config/
│   └── project.yaml              # Your configuration (edit this)
├── data/
│   ├── rainfall/                 # Rainfall data (NetCDF/CSV)
│   ├── terrain/                 # DEM (GeoTIFF)
│   ├── boundary/                # Boundary conditions (CSV)
│   ├── observations/
│   │   ├── gauges/             # Gauge data (CSV)
│   │   └── satellite/          # Flood masks (GeoTIFF)
│   └── outputs/
│       └── lisflood/           # Model outputs
├── models/
│   └── lisflood/
│       └── template.par         # LISFLOOD-FP parameter template
├── visualization/               # Calibration plots
└── logs/                       # Log files
```

## Data Format Requirements

### Gauge Data (CSV)

```csv
timestamp,discharge
2026-06-01T00:00:00Z,1250.5
2026-06-01T01:00:00Z,1265.2
2026-06-01T02:00:00Z,1280.3
...
```

### Satellite Flood Mask (GeoTIFF)

- Binary mask: 0=dry, 1=flooded
- Or depth grid: water depths in meters
- Must have CRS and geotransform

### LISFLOOD-FP Template (.par)

```text
FRICTION: 0.035
INFILTRATION: 0.002
MANNING_N: 0.035
EVAPORATION: 0.0
THETA: 0.7
...
```

## Running Calibration

### Basic Calibration

```bash
python main.py --config config/project.yaml
```

### With Plots

```bash
python main.py --config config/project.yaml --save-plots
```

### Override Iterations

```bash
python main.py --config config/project.yaml --max-iterations 50 --save-plots
```

### Verbose Mode

```bash
python main.py --config config/project.yaml --verbose --save-plots
```

## Output Files

After calibration completes:

```
data/outputs/lisflood/
├── calibration_results.json    # Best parameters + history
├── calibrated.par             # Calibrated .par file
├── sim.par                    # Last simulation parameters
└── [LISFLOOD-FP outputs]

visualization/
├── calibration_progress.png    # Metrics over iterations
├── pareto_front.png          # Multi-objective trade-offs
├── hydrograph_comparison.png # Observed vs simulated discharge
└── inundation_comparison.png # Observed vs simulated flooding

logs/
└── calibration.log           # Detailed log
```

## Interpreting Results

### Calibration Metrics

- **NSE**: Nash-Sutcliffe Efficiency (0-1, higher is better)
  - NSE > 0.85: Good
  - NSE 0.65-0.85: Satisfactory
  - NSE < 0.65: Unsatisfactory

- **KGE**: Kling-Gupta Efficiency (-∞ to 1, higher is better)
  - KGE > 0.80: Good
  - KGE 0.60-0.80: Satisfactory
  - KGE < 0.60: Unsatisfactory

- **IoU**: Intersection over Union (0-1, higher is better)
  - IoU > 0.70: Good spatial accuracy
  - IoU 0.50-0.70: Moderate accuracy
  - IoU < 0.50: Poor spatial accuracy

### Pareto Front

The Pareto front shows non-dominated parameter sets where:
- You can't improve one metric without worsening another

Choose based on your priorities:
- **Peak flow timing**: Higher NSE/KGE
- **Flood extent**: Higher IoU
- **Balanced**: Middle of Pareto front

## Troubleshooting

### LISFLOOD-FP Not Found

```bash
# Update executable path in config.yaml
model:
  executable: "/usr/local/bin/lisflood-fp"
```

### Parameter Out of Range

Check LISFLOOD-FP template parameters are within physical bounds:
- Friction: 0.01-0.1 (Manning's n)
- Infiltration: 0.0-0.01 (m/s)
- Theta: 0.5-1.0 (diffusion coefficient)

### NaN Metrics

- Check gauge data has valid timestamps and discharge values
- Ensure satellite mask is properly aligned to model grid
- Verify LISFLOOD-FP simulation completes successfully

### Poor Calibration

- Increase `max_iterations` (e.g., 200)
- Adjust `objective_weights` to prioritize your target metric
- Check that calibration data matches simulation period
- Verify boundary conditions represent actual event

## Next Steps

1. **Validate on different storm events** to test generalizability
2. **Adjust objective weights** based on your priorities
3. **Enable surrogate model** for faster convergence (in config.yaml)
4. **Integrate MetaClaw** for cross-run learning (future)

## Getting Help

- Check logs: `logs/calibration.log`
- Verify data paths in `config/project.yaml`
- Ensure all input files exist and are readable

---

Ready to calibrate! 🌊
