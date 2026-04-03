---
layout: default
title: Self-Calibrating LISFLOOD-FP
description: A model-agnostic, self-improving flood simulation system with automatic calibration
---

# Self-Calibrating LISFLOOD-FP 🌊

**A model-agnostic, self-improving flood simulation system with automatic calibration using gauge data and satellite-derived inundation extents.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/urabbani/self-calibrating-lisflood-fp/actions)

## Overview

![System Overview Diagram](https://via.placeholder.com/800x400/4A90E2/FFFFFF?text=Self-Calibrating+LISFLOOD-FP+System+Overview)

Self-Calibrating LISFLOOD-FP automates the calibration of flood models by automatically adjusting model parameters to optimize against both temporal (hydrograph) and spatial (inundation) observations, eliminating manual parameter tuning.

**Inspired by:** [karpathy/autoresearch](https://github.com/karpathy/autoresearch) + [aiming-lab/AutoResearchClaw](https://github.com/aiming-lab/AutoResearchClaw)

## Key Features

### 🔧 Automatic Calibration
No manual parameter tuning required - the system automatically finds optimal parameters.

### 🎯 Multi-Objective Optimization
Simultaneously optimizes:
- **NSE** (Nash-Sutcliffe Efficiency) - temporal accuracy
- **KGE** (Kling-Gupta Efficiency) - temporal accuracy  
- **IoU** (Intersection over Union) - spatial accuracy

### 🛰️ Satellite Integration
Supports Sentinel-1 SAR, Landsat, and custom satellite-derived flood masks.

### 🔄 Model-Agnostic Architecture
Designed to easily add other hydrological models (HEC-HMS, SWAT, etc.).

### ⚖️ Pareto Front Analysis
Maintains non-dominated parameter sets when objectives conflict, showing trade-offs.

### 🧠 MetaClaw Integration
Learns from previous calibration runs to improve efficiency in future similar basins.

### 🚀 Optional Surrogate Modeling
Uses Gaussian Process models to dramatically speed up calibration.

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

- Calibrated LISFLOOD-FP parameters (.par file)
- Calibration metrics (NSE, KGE, IoU)
- Pareto front visualization
- Comparison plots (hydrograph, inundation)
- Lessons learned (for future runs via MetaClaw)

## Architecture

```
flood-autocalib/
├── src/
│   └── lisflood-fp/        # LISFLOOD-FP 8.2 source (integrated)
├── scripts/
│   └── build_lisflood.py   # Compile integrated source
├── config/
│   └── project.yaml.example
├── models/
│   ├── adapter.py          # Base adapter interface
│   ├── lisflood/
│   │   └── adapter.py      # LISFLOOD-FP implementation
│   └── lisflood82/
│       └── adapter.py       # LISFLOOD-FP 8.2 (GPU, ACC/FV1/FV2/DG2)
├── calibration/
│   ├── metrics.py          # NSE, KGE, IoU
│   └── loop.py             # Auto-calibration loop
├── data/
│   └── preprocessing.py    # Input data validation
├── satellite/
│   └── preprocessor.py     # Satellite data alignment
├── tests/                  # Unit tests (63 tests)
├── visualization/
│   └── plots.py            # Calibration visualization
├── main.py                 # Entry point
└── requirements.txt
```

## Building LISFLOOD-FP

The LISFLOOD-FP 8.2 source is integrated directly. The framework auto-compiles it when `executable: "auto"` is set in config.

```bash
# Manual build (optional — auto-builds if needed)
python scripts/build_lisflood.py           # CPU-only
python scripts/build_lisflood.py --cuda    # With GPU acceleration
python scripts/build_lisflood.py --force   # Rebuild
```

**Build prerequisites:**
- CMake >= 3.13
- C++14 compiler (GCC 5+, Clang 3.4+, or MSVC 2017+)
- Optional: CUDA Toolkit >= 9.0 (for GPU acceleration)
- Optional: NetCDF libraries (for NetCDF I/O)

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

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use this code, please cite:

```bibtex
@software{flood_autocalib,
  title = {Self-Calibrating LISFLOOD-FP},
  author = {Umair Rabbani},
  year = {2026},
  url = {https://github.com/urabbani/self-calibrating-lisflood-fp}
}
```

## Contact

**Umair Rabbani**  
📧 [urabbani@example.com](mailto:urabbani@example.com)  
🐙 [github.com/urabbani](https://github.com/urabbani)  
💼 [linkedin.com/in/urabbani](https://linkedin.com/in/urabbani)

---

*Ready to revolutionize your flood modeling workflow? Try Self-Calibrating LISFLOOD-FP today!*