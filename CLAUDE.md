# CLAUDE.md - Instructions for AI Assistants

This file provides context and instructions for Claude, Codex, and other AI assistants working on this project.

---

## Project Overview

**Self-Calibrating LISFLOOD-FP** is a model-agnostic flood simulation calibration system that automatically optimizes hydraulic model parameters using:

1. **Gauge data** - Temporal discharge observations (NSE, KGE metrics)
2. **Satellite data** - Spatial inundation extents (IoU metric)
3. **AutoResearch-style loop** - Iterative parameter search (100+ iterations)
4. **Pareto front analysis** - Multi-objective trade-offs
5. **Lessons learned system** - Records parameter correlations for future runs

---

## Quick Context

- **Main model:** LISFLOOD-FP 8.2 (GPU-accelerated)
- **Integrated source:** src/lisflood-fp/ (LISFLOOD-FP 8.2 with auto-compilation)
- **Inspiration:** AutoResearch (karpathy) + AutoResearchClaw (aiming-lab)
- **Calibration metrics:** NSE, KGE (temporal) + IoU (spatial)
- **Calibratable parameters:** 11 total (friction, infiltration, CFL, etc.)

---

## Project Structure

```
flood-autocalib/
├── src/                       # Integrated LISFLOOD-FP 8.2 source
│   └── lisflood-fp/
├── config/                    # Configuration files
│   ├── project.yaml.example
│   └── project_lisflood82.yaml.example  # LISFLOOD-FP 8.2 specific
├── models/                    # Model adapters
│   ├── adapter.py            # Base interface
│   ├── lisflood/adapter.py  # Generic LISFLOOD-FP
│   └── lisflood82/adapter.py # LISFLOOD-FP 8.2 (with GPU)
├── calibration/               # Calibration logic
│   ├── metrics.py            # NSE, KGE, IoU calculations
│   └── loop.py               # Auto-calibration loop
├── data/                      # Data processing
│   └── preprocessing.py     # Data loading & validation
├── satellite/                 # Satellite data handling
│   └── preprocessor.py     # Flood mask alignment
├── tests/                     # Unit tests (pytest)
│   ├── test_metrics.py       # Metric calculations
│   ├── test_loop.py          # Calibration loop components
│   ├── test_adapter.py       # Model adapter & parameter validation
│   └── test_satellite.py     # Satellite preprocessing
├── visualization/             # Plotting
│   └── plots.py             # Calibration visualizations
├── main.py                    # Entry point
└── requirements.txt           # Python dependencies
```

---

## Key Components to Understand

### 1. Model Adapter Pattern

**Base class:** `ModelAdapter` in `models/adapter.py`

All model adapters must implement:
- `get_parameters()` - Return parameters with ranges
- `set_parameters(params)` - Update model parameters
- `run_simulation(storm_event)` - Run and return outputs
- `get_outputs()` - Extract discharge + inundation data

**Implemented adapters:**
- `LISFLOOD82Adapter` - For urabbani/lisflood-fp_8.2_update
- `LISFLOODAdapter` - Generic LISFLOOD-FP

### 2. Calibration Loop

**File:** `calibration/loop.py`

**Core algorithm:**
```
For each iteration:
1. Propose new parameters (exploration → refinement → exploitation)
2. Run LISFLOOD-FP simulation
3. Compute metrics: NSE, KGE, IoU
4. Update Pareto front (non-dominated solutions)
5. Record lessons (parameter correlations, trade-offs)
6. Check convergence (tolerance + target metrics)
```

**Strategies:**
- Iterations 0-9: Latin Hypercube Sampling (exploration)
- Iterations 10-29: Gaussian Perturbation (refinement)
- Iterations 30+: Fine-tuning (exploitation)

### 3. Multi-Objective Metrics

**File:** `calibration/metrics.py`

**Temporal (hydrograph):**
- **NSE** (Nash-Sutcliffe Efficiency): 0-1, higher is better
- **KGE** (Kling-Gupta Efficiency): -∞ to 1, higher is better

**Spatial (inundation):**
- **IoU** (Intersection over Union): 0-1, higher is better
- **F1 Score**: 0-1, higher is better

**Composite score:** 0.4×NSE + 0.3×KGE + 0.3×IoU (configurable)

### 4. LISFLOOD-FP 8.2 Adapter

**File:** `models/lisflood82/adapter.py`

**Key features:**
- Supports 5 solvers: ACC, FV1, FV2, DG2, ACC_NUGRID
- GPU acceleration via CUDA (`-cuda` flag)
- 11 calibratable parameters
- Automatic output extraction (.max, .wd, .discharge files)

**Calibratable parameters:**
| Parameter | Description | Range |
|-----------|-------------|--------|
| fpfric | Floodplain Manning's n | 0.01-0.1 |
| SGCn | Sub-grid channel Manning's n | 0.02-0.06 |
| infinity | Infiltration rate (m/s) | 0.0-0.001 |
| courant | CFL number | 0.3-0.9 |
| theta | Time-stepping parameter | 0.5-1.0 |
| depth_thresh | Minimum depth threshold (m) | 0.0005-0.01 |
| gravity | Gravitational acceleration (m/s²) | 9.7-9.9 |
| max_Froude | Maximum Froude number | 1.0-10000.0 |
| SGCr | Sub-grid channel multiplier | 0.1-0.2 |
| SGCp | Sub-grid channel exponent | 0.7-0.85 |

---

## Usage for AI Assistants

### When Adding New Model Adapters

1. Create new directory under `models/` (e.g., `models/hechms/`)
2. Implement `ModelAdapter` base class
3. Define parameter ranges in `__init__`
4. Implement `run_simulation()` to return `SimulationResult`
5. Add model type to `main.py` initialization logic
6. Create config example in `config/`

### When Modifying Calibration Loop

- Keep lessons learned system intact (`_record_lessons()`)
- Maintain Pareto front logic (`update_pareto_front()`)
- Preserve convergence detection (`check_convergence()`)
- Document new strategies in comments

### When Adding New Metrics

1. Add metric function to `calibration/metrics.py`
2. Add to `composite_score()` normalization
3. Update config validation in `validate_metrics()`
4. Document in `CalibrationMetrics` docstring

### When Fixing Bugs

- Check logs: `logs/calibration.log`
- Verify LISFLOOD-FP executable path
- Validate input data paths in config
- Test LISFLOOD-FP manually with verbose mode: `./lisflood -v simulation.par`

---

## Common Tasks

### Task: Run Calibration

```bash
# Basic
python main.py --config config/project.yaml

# With plots
python main.py --config config/project.yaml --save-plots

# Override iterations
python main.py --config config/project.yaml --max-iterations 200 --save-plots
```

### Task: Add New Model Adapter

1. Create `models/newmodel/adapter.py`
2. Implement `ModelAdapter` interface
3. Update `main.py` to recognize model type
4. Add config example: `config/project_newmodel.yaml.example`
5. Update documentation

### Task: Debug Calibration Issues

1. Check config paths: all files exist?
2. Run LISFLOOD-FP manually: `/path/to/lisflood -v simulation.par`
3. Review logs: `logs/calibration.log`
4. Check gauge data format (CSV with timestamp, discharge)
5. Verify satellite mask (GeoTIFF with CRS)

---

## Dependencies

**Core:**
- numpy, scipy - Numerical computations
- pandas - Data handling
- yaml - Configuration

**Geospatial:**
- rasterio - Raster I/O
- geopandas - Vector data
- shapely - Geometry operations

**Visualization:**
- matplotlib, seaborn - Plotting

**Machine Learning (future):**
- scikit-learn, scikit-optimize - Surrogate modeling

---

## Configuration System

**Primary file:** `config/project.yaml`

**Key sections:**
- `model`: LISFLOOD-FP executable, solver options
- `input_data`: DEM, rainfall, boundary paths
- `observations`: Gauge + satellite data paths
- `calibration`: Iterations, targets, weights
- `lisflood82`: Solver settings, output control
- `visualization`: Plotting options
- `surrogate`: Surrogate model settings (optional)
- `metaclaw`: Cross-run learning (optional)

---

## Documentation Files

- **README.md** - Full framework documentation
- **SUMMARY.md** - Quick overview and features
- **QUICKSTART.md** - Quick start guide
- **SETUP_LISFLOOD82.md** - Setup guide for LISFLOOD-FP 8.2
- **CLAUDE.md** - This file (AI assistant instructions)

---

## Inspiration and Design Philosophy

### AutoResearch (karpathy/autoresearch)
- Simple loop: modify → run → evaluate → keep/discard
- Fixed time budget (5 minutes)
- Single-file scope for reviewability
- ~100 experiments overnight

### AutoResearchClaw (aiming-lab/AutoResearchClaw)
- 23-stage research pipeline
- PROCEED/REFINE/PIVOT decisions
- MetaClaw integration (lessons → skills)
- Multi-agent debate
- Anti-fabrication system

### This Framework
- AutoResearch-style loop for calibration (simpler than 23 stages)
- Multi-objective optimization (NSE, KGE, IoU)
- Pareto front analysis
- Lessons learned (ready for MetaClaw integration)
- Model-agnostic (plug-and-play adapters)

---

## Testing

### Automated Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test modules
python -m pytest tests/test_metrics.py -v
python -m pytest tests/test_loop.py -v
```

Test coverage includes:
- **test_metrics.py** (32 tests): NSE, KGE, IoU, F1, timing_error, composite_score, validate_metrics
- **test_loop.py** (10 tests): LHS sampling, local search, Pareto front, convergence, parameter application
- **test_adapter.py** (9 tests): SimulationResult defaults, parameter validation, parameter bounds
- **test_satellite.py** (13 tests): MNDWI (including div-by-zero), SAR, depth thresholding, mask smoothing, validation

### Manual Test with Sample Data

1. Use provided example files in `data/`
2. Run: `python main.py --config config/project.yaml --save-plots`
3. Check outputs in `data/outputs/lisflood82/`
4. Review plots in `visualization/`

### Expected Output

- `calibrated.par` - Best parameter set
- `calibration_results.json` - Full history
- `simulation.max`, `.wd`, `.discharge` - Model outputs
- Plots: progress, pareto, hydrograph, inundation

---

## Contributing Guidelines

- Follow PEP 8 for Python code
- Add docstrings to all functions
- Update this file when adding new features
- Test with example data before committing
- Update documentation (README, guides)
- Use descriptive commit messages

---

## Version History

**0.1.1** (2026-03-30)
- Fix: Calibration loop now correctly applies proposed parameters to model before simulation
- Fix: Convergence detection requires targets met or extended score plateau (was premature)
- Fix: Pareto front update no longer mutates list during iteration
- Fix: Proper Latin Hypercube Sampling with stratified permutations
- Fix: Local RNG via `np.random.default_rng()` (no global seed pollution)
- Fix: `timing_error()` uses consistent NaN-filtered indexing
- Fix: MNDWI division-by-zero guard for zero denominator
- Fix: Bare `except` replaced with `except Exception`
- Fix: Safe `config.get("metaclaw", {})` access prevents KeyError
- Fix: ASCII raster header parsing is case-insensitive
- Fix: `re.escape()` on parameter names in regex
- Fix: Unknown parameters now correctly rejected by `validate_parameters()`
- Fix: Cross-package relative imports converted to absolute imports
- Fix: `validate_mask()` edge pixel calculation shape mismatch
- Add: Python logging module with config-based setup
- Add: 63 unit tests across 4 test files
- Add: `LISFLOOD82Adapter` to package exports
- Remove: Unused `scipy.stats` import, duplicate `numpy` import
- Move: `shutil` imports to module level in adapter files

**0.1.0** (2026-03-30)
- Initial release
- LISFLOOD-FP 8.2 adapter with GPU support
- Multi-objective calibration (NSE, KGE, IoU)
- AutoResearch-style loop
- Pareto front analysis
- Visualization system
- Complete documentation

---

## Contact and Support

- **Main repository:** https://github.com/urabbani/self-calibrating-lisflood-fp
- **LISFLOOD-FP 8.2:** https://github.com/urabbani/lisflood-fp_8.2_update
- **Issues:** Open issue on GitHub repository

---

**Note:** This project is actively developed. Feel free to experiment, improve, and extend!

