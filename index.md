---
layout: default
title: Self-Calibrating LISFLOOD-FP
description: Automatic calibration of flood models using multi-objective optimization with gauge data and satellite observations
---

<link rel="stylesheet" href="{{ '/assets/css/style.scss' | relative_url }}">

<div class="hero-buttons">
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp" class="btn-primary">
    <svg width="18" height="18" viewBox="0 0 16 16" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
    View on GitHub
  </a>
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/blob/main/QUICKSTART.md" class="btn-secondary">
    Quick Start Guide
  </a>
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/archive/refs/heads/main.zip" class="btn-secondary">
    Download ZIP
  </a>
</div>

## Why Self-Calibrating?

Flood model calibration is traditionally **manual, time-consuming, and subjective**. Hydrologists spend days to weeks tuning Manning's roughness coefficients, infiltration rates, and other parameters by trial and error. Self-Calibrating LISFLOOD-FP eliminates this bottleneck by automating the entire process with multi-objective optimization.

**What it does:**

- Automatically finds optimal Manning's n, infiltration, CFL, and 8 other parameters
- Simultaneously calibrates against **discharge gauges** and **satellite flood extents**
- Maintains a **Pareto front** so you can choose between peak-flow accuracy, spatial accuracy, or balanced solutions
- Learns from previous runs to accelerate future calibrations

<div class="callout">
  <p><strong>Built for LISFLOOD-FP 8.2</strong> — supports all solvers (ACC, FV1, FV2, DG2, ACC_NUGRID), GPU acceleration via CUDA, and 11 calibratable hydraulic parameters. The architecture is model-agnostic, so adapters for HEC-HMS, SWAT, and others can be added.</p>
</div>

## Key Features

<div class="feature-grid">
  <div class="feature-card">
    <div class="feature-icon">🎯</div>
    <h4>Multi-Objective Optimization</h4>
    <p>Simultaneously optimizes NSE, KGE (temporal/hydrograph) and IoU (spatial/inundation) with configurable weights.</p>
  </div>
  <div class="feature-card">
    <div class="feature-icon">🛰️</div>
    <h4>Satellite Integration</h4>
    <p>Direct comparison against Sentinel-1 SAR, Landsat, and custom satellite-derived flood masks via MNDWI and backscatter thresholding.</p>
  </div>
  <div class="feature-card">
    <div class="feature-icon">⚡</div>
    <h4>GPU Acceleration</h4>
    <p>LISFLOOD-FP 8.2 with CUDA support runs simulations orders of magnitude faster, enabling 100+ calibration iterations.</p>
  </div>
  <div class="feature-card">
    <div class="feature-icon">⚖️</div>
    <h4>Pareto Front Analysis</h4>
    <p>Maintains non-dominated parameter sets when objectives conflict. Choose peak-flow focused, spatially accurate, or balanced solutions.</p>
  </div>
  <div class="feature-card">
    <div class="feature-icon">🧠</div>
    <h4>Lessons Learned</h4>
    <p>Records parameter correlations and trade-offs from each run. Future calibrations in similar basins converge faster.</p>
  </div>
  <div class="feature-card">
    <div class="feature-icon">🚀</div>
    <h4>Surrogate Modeling</h4>
    <p>Optional Gaussian Process surrogate trained on historical runs uses Expected Improvement to propose better parameters with fewer simulations.</p>
  </div>
</div>

## How It Works

<div class="workflow-steps">
  <div class="step">
    <div class="step-content">
      <h4>Prepare Input Data</h4>
      <p>Provide rainfall time series, DEM terrain, upstream boundary discharge, and validation data (gauge observations + satellite flood masks).</p>
    </div>
  </div>
  <div class="step">
    <div class="step-content">
      <h4>Configure Parameters</h4>
      <p>Define parameter ranges for 11 calibratable variables (Manning's n, infiltration, CFL, sub-grid channel parameters, etc.) in a YAML config file.</p>
    </div>
  </div>
  <div class="step">
    <div class="step-content">
      <h4>Run Auto-Calibration Loop</h4>
      <p>The system iterates through exploration (Latin Hypercube Sampling), refinement (Gaussian Perturbation), and exploitation (fine-tuning) phases, running LISFLOOD-FP at each step.</p>
    </div>
  </div>
  <div class="step">
    <div class="step-content">
      <h4>Evaluate Against Observations</h4>
      <p>Each simulation is scored against gauge data (NSE, KGE) and satellite flood extents (IoU, F1). The Pareto front is updated with non-dominated solutions.</p>
    </div>
  </div>
  <div class="step">
    <div class="step-content">
      <h4>Extract Calibrated Parameters</h4>
      <p>Select from the Pareto front: overall best, peak-flow focused, spatially accurate, or balanced. Output includes calibrated <code>.par</code> file, metrics, and visualization plots.</p>
    </div>
  </div>
</div>

## Quick Start

```bash
# Clone the repository
git clone https://github.com/urabbani/self-calibrating-lisflood-fp.git
cd self-calibrating-lisflood-fp

# Install dependencies
pip install -r requirements.txt

# LISFLOOD-FP 8.2 auto-compiles on first run, or build manually:
python scripts/build_lisflood.py --cuda

# Configure your project
cp config/project.yaml.example config/project.yaml
# Edit config.yaml with your data paths and calibration targets

# Run calibration
python main.py --config config/project.yaml

# With visualization output
python main.py --config config/project.yaml --save-plots
```

<div class="callout">
  <p><strong>Build prerequisites:</strong> CMake >= 3.13, C++14 compiler (GCC 5+, Clang 3.4+, MSVC 2017+). Optional: CUDA Toolkit >= 9.0 for GPU acceleration, NetCDF libraries for NetCDF I/O.</p>
</div>

## Calibration Metrics

<div class="metric-badges">
  <div class="metric-badge">
    <strong>NSE (Nash-Sutcliffe)</strong>
    <span>Temporal accuracy of simulated hydrograph vs. gauge observations (0-1)</span>
  </div>
  <div class="metric-badge">
    <strong>KGE (Kling-Gupta)</strong>
    <span>Decomposes bias, correlation, and variability of discharge time series (-inf to 1)</span>
  </div>
  <div class="metric-badge">
    <strong>IoU (Intersection over Union)</strong>
    <span>Spatial overlap between simulated and satellite-derived inundation extent (0-1)</span>
  </div>
  <div class="metric-badge">
    <strong>F1 Score</strong>
    <span>Binary classification accuracy of flooded vs. dry pixels (0-1)</span>
  </div>
</div>

**Composite Score:** `0.4 * NSE + 0.3 * KGE + 0.3 * IoU` (weights configurable in YAML)

## Calibratable Parameters

The adapter exposes 11 LISFLOOD-FP parameters for automatic optimization:

| Parameter | Description | Range |
|-----------|-------------|-------|
| `fpfric` | Floodplain Manning's roughness coefficient | 0.01 - 0.10 |
| `SGCn` | Sub-grid channel Manning's roughness | 0.02 - 0.06 |
| `infinity` | Infiltration rate (m/s) | 0.0 - 0.001 |
| `courant` | CFL number for adaptive time-stepping | 0.3 - 0.9 |
| `theta` | Time-stepping weight parameter | 0.5 - 1.0 |
| `depth_thresh` | Minimum depth threshold (m) | 0.0005 - 0.01 |
| `gravity` | Gravitational acceleration (m/s2) | 9.7 - 9.9 |
| `max_Froude` | Maximum Froude number | 1.0 - 10000.0 |
| `SGCr` | Sub-grid channel width multiplier | 0.1 - 0.2 |
| `SGCp` | Sub-grid channel exponent | 0.7 - 0.85 |
| `SGCd` | Sub-grid channel depth | 0.5 - 5.0 |

## Calibration Strategy

The AutoResearch-style loop progresses through three phases:

1. **Exploration** (iterations 0-9): Latin Hypercube Sampling across the full parameter space to build an initial understanding of the response surface.

2. **Refinement** (iterations 10-29): Gaussian Perturbation around the best-known solutions, narrowing the search to promising regions.

3. **Exploitation** (iterations 30+): Fine-tuning within the best region, making small adjustments to converge on optimal parameters.

Convergence is detected when target metrics are met or the composite score plateaus for a configurable number of iterations.

## Configuration

All settings are defined in a single YAML file:

```yaml
project:
  name: "my-basin-study"

model:
  type: "lisflood-fp"
  executable: "auto"          # auto-compiles from integrated source
  param_template: "data/lisflood/template.par"

input_data:
  rainfall: "data/rainfall/event_2026.nc"
  terrain: "data/terrain/dem.tif"
  boundary_discharge: "data/boundary/upstream.csv"

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

observations:
  gauges:
    - id: "gauge_001"
      path: "data/observations/gauges/gauge_001.csv"
      location: [67.123, 24.567]
  satellite:
    - product: "sentinel1"
      path: "data/observations/satellite/flood_mask.tif"
      timestamp: "2026-06-03T10:00:00Z"
```

See the full [configuration example](https://github.com/urabbani/self-calibrating-lisflood-fp/blob/main/config/project.yaml.example) for all options including surrogate modeling and MetaClaw integration.

## Input Data Requirements

**Required inputs:**
- **Rainfall** — Time series (CSV or NetCDF) with timestamps and mm/hour values
- **Terrain (DEM)** — GeoTIFF with elevation in meters
- **Boundary Discharge** — Time series (CSV) for upstream boundary condition

**Validation data:**
- **Discharge Observations** — Gauge data (CSV) with timestamps and m3/s
- **Inundation Extent** — Satellite-derived flood mask (GeoTIFF, binary or probability)

## Output

- **Calibrated parameters** — `.par` file ready for LISFLOOD-FP production runs
- **Calibration metrics** — NSE, KGE, IoU, F1, and composite score history
- **Pareto front** — Multi-objective trade-off visualization
- **Hydrograph comparison** — Observed vs. simulated discharge over time
- **Inundation comparison** — Simulated vs. satellite flood extent maps
- **Lessons learned** — Parameter correlations recorded for future calibration runs

## Architecture

The framework follows a modular adapter pattern:

- **`models/adapter.py`** — Base adapter interface (plug in any hydraulic model)
- **`models/lisflood82/adapter.py`** — LISFLOOD-FP 8.2 adapter (ACC, FV1, FV2, DG2 solvers, GPU)
- **`calibration/loop.py`** — AutoResearch-style calibration engine
- **`calibration/metrics.py`** — NSE, KGE, IoU, F1, composite scoring
- **`satellite/preprocessor.py`** — Flood mask alignment and water extraction
- **`visualization/plots.py`** — Diagnostic and output plotting

**Test suite:** 63 unit tests covering metrics, calibration loop, adapter validation, and satellite preprocessing.

## Roadmap

- HEC-HMS adapter
- SWAT adapter
- Ensemble Kalman Filter (sequential calibration)
- Physics-Informed Neural Network (PINN) surrogate
- Multi-catchment calibration with transfer learning
- Real-time calibration during flood events

## Citation

If you use this software in your research, please cite:

<div class="citation-box">
<pre><code>@software{flood_autocalib,
  title = {Self-Calibrating {LISFLOOD-FP}},
  author = {Umair Rabbani},
  year = {2026},
  url = {https://github.com/urabbani/self-calibrating-lisflood-fp}
}</code></pre>
</div>

## Links

<div class="link-grid">
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg>
    Source Code
  </a>
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/issues">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M8 1.5a6.5 6.5 0 100 13 6.5 6.5 0 000-13zM0 8a8 8 0 1116 0A8 8 0 010 8zm9 3a1 1 0 11-2 0 1 1 0 012 0zm-.25-6.25a.75.75 0 00-1.5 0v3.5a.75.75 0 001.5 0v-3.5z"/></svg>
    Report Issues
  </a>
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/blob/main/QUICKSTART.md">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M0 1.75A.75.75 0 01.75 1h4.253c1.227 0 2.317.59 3 1.501A3.744 3.744 0 0111.006 1h4.245a.75.75 0 01.75.75v10.5a.75.75 0 01-.75.75h-4.507a2.25 2.25 0 00-1.591.659l-.622.621a.75.75 0 01-1.06 0l-.622-.621A2.25 2.25 0 005.258 13H.75a.75.75 0 01-.75-.75V1.75z"/></svg>
    Quick Start Guide
  </a>
  <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/blob/main/SETUP_LISFLOOD82.md">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M11.28 3.22a.75.75 0 00-1.06 1.06L12.94 7H7.75a.75.75 0 000 1.5h5.19l-2.72 2.72a.75.75 0 101.06 1.06l4-4a.75.75 0 000-1.06l-4-4z"/></svg>
    LISFLOOD-FP Setup
  </a>
  <a href="https://github.com/urabbani">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M10.561 8.073a6.005 6.005 0 013.432 5.142.75.75 0 11-1.498.07 4.5 4.5 0 00-8.99 0 .75.75 0 11-1.498-.07 6.004 6.004 0 013.431-5.142 3.999 3.999 0 115.123 0zM10.5 5a2.5 2.5 0 10-5 0 2.5 2.5 0 005 0z"/></svg>
    Umair Rabbani
  </a>
  <a href="https://www.linkedin.com/in/drumairrabbani/">
    <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor"><path d="M0 1.146C0 .513.526 0 1.175 0h13.65C15.474 0 16 .513 16 1.146v13.708c0 .633-.526 1.146-1.175 1.146H1.175C.526 16 0 15.487 0 14.854V1.146zM4.943 12.57V6.3H2.508v6.27h2.435zM3.726 5.302c.846 0 1.372-.561 1.372-1.261-.016-.716-.526-1.261-1.355-1.261-.83 0-1.371.545-1.371 1.261 0 .7.526 1.261 1.34 1.261h.014zm5.04 7.268h2.435V9.04c0-.215.016-.43.079-.583.173-.43.567-.876 1.228-.876.867 0 1.213.661 1.213 1.631v3.358h2.435V8.913c0-2.254-1.204-3.3-2.81-3.3-1.316 0-1.893.737-2.211 1.237h.016v-1.06H8.766c.032.682 0 6.27 0 6.27h-.001z"/></svg>
    LinkedIn
  </a>
</div>

---

<p align="center"><strong>Umair Rabbani</strong> &mdash; <a href="mailto:umairrs@gmail.com">umairrs@gmail.com</a></p>

<p align="center">Licensed under the <a href="https://github.com/urabbani/self-calibrating-lisflood-fp/blob/main/LICENSE">MIT License</a></p>
