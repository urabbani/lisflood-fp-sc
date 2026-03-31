# Self-Calibrating LISFLOOD-FP 🌊

**Inspired by AutoResearch and AutoResearchClaw**

A model-agnostic, self-improving flood simulation system with automatic calibration using gauge data and satellite-derived inundation extents.

**Version:** 0.1.1
**Status:** Production Ready
**Inspiration:** [karpathy/autoresearch](https://github.com/karpathy/autoresearch) + [aiming-lab/AutoResearchClaw](https://github.com/aiming-lab/AutoResearchClaw)

## Overview

```
User Input                     System                       Auto-Calibration
┌───────────────┐            ┌──────────────┐              ┌──────────────┐
│ Rainfall      │──────┐     │ LISFLOOD-FP  │───┬──┐       │ Multi-       │
│ (Temporal)    │      │     │ (Physics     │   │  │──────│ Objective    │
│               │      │     │  Engine)     │   │  │       │ Optimizer    │
│ Terrain (DEM) │───┐  │     └──────────────┘   │  │       │ (NSE/KGE/IoU)│
│ (Spatial)     │   │  │                        │  │       └──────────────┘
│               │   │  │     ┌──────────────┐   │  │
│ Boundary      │   │  └────│ Observation  │───┘  │
│ Discharge     │   │       │ Comparator    │      │
│ (Input)       │   │       └──────────────┘      │
└───────────────┘   │                              │
                    │    ┌──────────────┐         │
                    └───│ Calibration  │◀────────┘
                         │ Loop         │
                         └──────────────┘
```

## Features

- **Automatic Calibration**: No manual parameter tuning
- **Multi-Objective**: Optimize NSE, KGE (temporal) and IoU (spatial)
- **Satellite Integration**: Support Sentinel-1 SAR, Landsat, custom masks
- **Model-Agnostic Architecture**: Easy to add other models (HEC-HMS, SWAT)
- **Pareto Front Analysis**: Non-dominated parameter sets
- **Lessons Learned**: MetaClaw integration for cross-run knowledge

## Quick Start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Prepare input data
cp config/project.yaml.example config/project.yaml
# Edit config.yaml with your data paths

# 3. Run calibration
python main.py --config config/project.yaml

# 4. View results
# - outputs/calibration_results.json (best parameters)
# - outputs/pareto_front.png (multi-objective trade-offs)
# - outputs/hydrograph_comparison.png (observed vs simulated)
# - outputs/inundation_comparison.png (spatial comparison)
```

## Input Data Requirements

### Required
- **Rainfall**: Time series (CSV or NetCDF) with timestamps and mm/hour values
- **Terrain (DEM)**: GeoTIFF with elevation (meters)
- **Boundary Discharge**: Time series (CSV) for upstream boundary condition

### Validation Data
- **Discharge Observations**: Gauge data (CSV) with timestamps and m³/s
- **Inundation Extent**: Satellite-derived flood mask (GeoTIFF, binary or probability)

## Output

- **Calibrated LISFLOOD-FP parameters** (.par file)
- **Calibration metrics** (NSE, KGE, IoU)
- **Pareto front** visualization
- **Comparison plots** (hydrograph, inundation)
- **Lessons learned** (for future runs via MetaClaw)

## Architecture

```
flood-autocalib/
├── config/
│   └── project.yaml.example
├── data/
│   ├── rainfall/
│   ├── terrain/
│   ├── boundary/
│   ├── observations/
│   │   ├── gauges/
│   │   └── satellite/
│   └── outputs/
├── models/
│   ├── adapter.py          # Base adapter interface
│   ├── lisflood/
│   │   └── adapter.py      # LISFLOOD-FP implementation
│   └── hechms/
│       └── adapter.py      # (Future) HEC-HMS
├── calibration/
│   ├── metrics.py          # NSE, KGE, IoU
│   ├── loop.py             # Auto-calibration loop
│   └── surrogate.py        # Gaussian Process surrogate
├── data/
│   ├── preprocessing.py    # Input data validation
│   └── observations.py     # Load gauge + satellite data
├── satellite/
│   ├── extraction.py       # Water extraction from GEE
│   └── preprocessor.py     # Satellite data alignment
├── tests/                  # Unit tests (63 tests)
│   ├── test_metrics.py     # Metric calculations
│   ├── test_loop.py        # Calibration loop
│   ├── test_adapter.py     # Model adapters
│   └── test_satellite.py   # Satellite preprocessing
├── visualization/
│   └── plots.py            # Calibration visualization
├── main.py                 # Entry point
└── requirements.txt
```

## Configuration Example

```yaml
project:
  name: "indus-basin-demo"

# Model adapter
model:
  type: "lisflood-fp"
  executable: "/path/to/lisflood-fp"
  param_template: "data/lisflood/template.par"
  output_dir: "data/outputs/lisflood/"

# Input data paths
input_data:
  rainfall: "data/rainfall/rainfall_202606.nc"
  terrain: "data/terrain/dem.tif"
  boundary_discharge: "data/boundary/upstream_discharge.csv"

# Calibration settings
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

# Validation data
observations:
  gauges:
    - id: "gauge_001"
      path: "data/observations/gauges/gauge_001.csv"
      location: [67.123, 24.567]

  satellite:
    - product: "sentinel1"
      path: "data/observations/satellite/flood_mask.tif"
      timestamp: "2026-06-03T10:00:00Z"
      # Optional: auto-extract via GEE
      auto_extract:
        enabled: false
        ee_aoi: "data/basin_boundary.geojson"
        ee_date_range: ["2026-06-01", "2026-06-07"]

# Surrogate model (speed up calibration)
surrogate:
  enabled: true
  retrain_interval: 10
  model_type: "gaussian_process"

# MetaClaw integration
metaclaw:
  enabled: true
  skills_dir: "~/.metaclaw/skills/flood-calibration/"
  lesson_to_skill:
    min_severity: "warning"
    max_skills_per_run: 5
```

## How It Works

1. **Data Preprocessing**: Validate and align input data (rainfall, terrain, boundary)
2. **Parameter Initialization**: Load LISFLOOD-FP template with parameter ranges
3. **Calibration Loop** (AutoResearch-style):
   - Agent proposes new parameters (exploration → exploitation)
   - Run LISFLOOD-FP simulation
   - Compute metrics: NSE, KGE (temporal) + IoU (spatial)
   - Update Pareto front
   - Record lessons for MetaClaw
4. **Surrogate Acceleration** (optional):
   - Train GP model on historical runs
   - Use Expected Improvement to propose next parameters
   - Validate with LISFLOOD-FP
5. **Output**: Best parameters + visualization

## Satellite Water Extraction (Optional)

```python
# Enable automatic extraction in config.yaml
satellite:
  auto_extract:
    enabled: true
    ee_aoi: "data/basin_boundary.geojson"
    ee_date_range: ["2026-06-01", "2026-06-07"]
    products: ["COPERNICUS/S1_GRD", "LANDSAT/LC08/C02/T1"]
    water_threshold: -0.1  # MNDWI threshold
```

This uses Google Earth Engine to:
- Query Sentinel-1 SAR and Landsat data
- Extract water indices (MNDWI, backscatter thresholding)
- Generate flood masks aligned to model grid

## Metrics

### Temporal (Hydrograph)
- **NSE**: Nash-Sutcliffe Efficiency (0-1, higher is better)
- **KGE**: Kling-Gupta Efficiency (-∞ to 1, higher is better)

### Spatial (Inundation)
- **IoU**: Intersection over Union (0-1, higher is better)
- **F1 Score**: Binary classification accuracy

### Composite Score

```
Score = 0.4 * NSE + 0.3 * KGE + 0.3 * IoU
```

(Tunable weights in config.yaml)

## Pareto Front

When objectives conflict (e.g., NSE improves but IoU drops), the system maintains a Pareto front of non-dominated solutions. You can choose:

- **Overall best** (highest composite score)
- **Peak-flow focused** (higher NSE)
- **Spatially accurate** (higher IoU)
- **Balanced** (middle of Pareto front)

## Future Extensions

- [ ] HEC-HMS adapter
- [ ] SWAT adapter
- [ ] Ensemble Kalman Filter (sequential calibration)
- [ ] Physics-Informed Neural Network (PINN) surrogate
- [ ] Multi-catchment calibration with transfer learning
- [ ] Real-time calibration during flood events

## License

MIT License

## Citation

If you use this code, please cite:

```bibtex
@software{flood_autocalib,
  title = {Self-Calibrating LISFLOOD-FP},
  author = {Umair Rabbani},
  year = {2026},
  url = {https://github.com/yourusername/flood-autocalib}
}
```
