# LISFLOOD-FP Implementation Improvements

This document outlines the improvements made to the LISFLOOD-FP codebase to address scientific and computational issues identified in the analysis.

## 1. Scientific Improvements

### 1.1 Advanced Friction Model for Shallow Flows

The friction model has been enhanced to better handle very shallow flows where traditional Manning's formula can become unstable or physically unrealistic. The new model:

- Implements laminar flow physics for very shallow depths (< 1 cm)
- Uses Darcy-Weisbach equation with Reynolds number-dependent friction factors for laminar flows
- Creates a smooth transition between laminar and turbulent flow regimes (1-5 cm)
- Maintains the standard Manning's formula for deeper flows (> 5 cm)
- Improves numerical stability for shallow flow conditions

Implementation files:
- `/src/swe/dg2/friction.cpp`

### 1.2 Green-Ampt Infiltration Scheme

A physically-based Green-Ampt infiltration model has been implemented to replace the simple constant-rate infiltration model. The new scheme:

- Tracks soil moisture content dynamically
- Models infiltration based on hydraulic conductivity, soil properties, and moisture deficit
- Accounts for cumulative infiltration effects
- Supports spatially variable soil properties
- Provides more realistic simulation of soil-water interactions
- Maintains continuity with the groundwater system

Implementation files:
- `/src/infevap.cpp`
- `/src/lisflood.h` (structure definitions)

## 2. Computational Improvements

### 2.1 CUDA Memory Management Enhancement

A modern, safe CUDA memory management system has been implemented using C++ RAII principles and smart pointers. This system:

- Replaces raw pointer allocation/deallocation with smart pointers
- Adds type safety for CUDA memory operations
- Provides automatic resource cleanup to prevent memory leaks
- Improves error detection through systematic error checking
- Simplifies memory management code

Implementation files:
- `/src/cuda/cuda_safe_memory.h`
- `/src/cuda/cuda_memory_example.cpp`

## 3. Infrastructure Updates

### 3.1 Docker Containerization

Docker support has been added to simplify deployment and enhance reproducibility:

- Multi-stage Docker build for smaller image size
- Support for both CPU and GPU execution
- Docker Compose configuration for easy deployment
- Runtime configuration via environment variables
- Volume mounts for input/output data
- Comprehensive documentation of container usage

Implementation files:
- `/Dockerfile`
- `/docker-compose.yml`
- `/docker-entrypoint.sh`
- `/run-docker.sh`
- `/.dockerignore`

## 4. Usage Examples

### 4.1 Advanced Friction Model

The improved friction model should be automatically applied in all simulations. It is especially beneficial for:

- Urban flood modeling where shallow flows are common
- Flood recession phases where depths decrease below typical modeling thresholds
- Rainfall-runoff simulations over complex terrain

### 4.2 Green-Ampt Infiltration

To use the Green-Ampt infiltration model, provide the following parameters in your simulation file:

```
# Required parameters
soildepth       1.0             # Soil depth in meters
soilporosity    soil_poros.asc  # Soil porosity grid (can be uniform)
suction_head    0.1             # Wetting front suction head (m)

# Optional parameters
init_moisture   0.2             # Initial soil moisture content
dist_inf        hyd_cond.asc    # Hydraulic conductivity grid
```

### 4.3 Docker Container

To run a simulation using the Docker container:

```bash
# Build the Docker image
docker build -t lisflood-fp .

# Run a simulation with GPU acceleration
./run-docker.sh -i /path/to/input -o /path/to/output -g -p simulation.par

# Run without GPU
./run-docker.sh -i /path/to/input -o /path/to/output -p simulation.par
```

## 5. Future Work

The following enhancements are planned for future development:

1. Multi-GPU support for very large domains
2. Python bindings for easier integration into data science workflows
3. Automated regression testing framework
4. Sediment transport model
5. Advanced evapotranspiration model
6. Groundwater coupling