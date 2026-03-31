# LISFLOOD-FP 8.2 Technical Documentation

## Table of Contents

1. [Model Overview](#model-overview)
2. [Numerical Solvers](#numerical-solvers)
3. [GPU Implementation](#gpu-implementation)
4. [Dynamic Grid Adaptation](#dynamic-grid-adaptation)
5. [Sub-Grid Channel Model](#sub-grid-channel-model)
6. [Hydraulic Structures](#hydraulic-structures)
7. [Dam Operations](#dam-operations)
8. [Hazard Assessment](#hazard-assessment)
9. [Input/Output Formats](#inputoutput-formats)
10. [Performance Considerations](#performance-considerations)
11. [Code Structure](#code-structure)
12. [Development History](#development-history)

## Model Overview

LISFLOOD-FP is a two-dimensional hydraulic model for simulating floodplain inundation. It solves the shallow water equations (SWE) or simplified versions of them, depending on the selected solver. The model operates on a regular grid and implements multiple computational methods to balance accuracy and computational efficiency.

### Core Features

- Multiple numerical solvers with varying levels of physical representation
- GPU acceleration using CUDA for high-performance computing
- Sub-grid channel representation for efficient river modeling
- Dynamic grid adaptation for multi-resolution simulations
- Complex boundary condition handling
- Comprehensive hydraulic structure modeling
- Integrated hazard assessment

### Physical Basis

LISFLOOD-FP solves various forms of the shallow water equations:

1. **Full Shallow Water Equations** (FV1, FV2, DG2 solvers):
   
   Momentum equation:
   
   $$\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u} = -g \nabla h + \mathbf{S_f}$$
   
   Continuity equation:
   
   $$\frac{\partial h}{\partial t} + \nabla \cdot (h\mathbf{u}) = 0$$

2. **Local Inertial Approximation** (ACC solver):
   
   $$\frac{\partial \mathbf{q}}{\partial t} + g h \nabla h = -g h \mathbf{S_f}$$
   
   $$\frac{\partial h}{\partial t} + \nabla \cdot \mathbf{q} = 0$$

3. **Diffusive Wave Approximation** (Diffusive mode):
   
   $$\mathbf{q} = -\frac{h^{5/3}}{\mathbf{n}} \nabla h \cdot |\nabla h|^{-1/2}$$
   
   $$\frac{\partial h}{\partial t} + \nabla \cdot \mathbf{q} = 0$$

Where:
- $h$ is water depth
- $\mathbf{u}$ is velocity vector
- $\mathbf{q}$ is discharge per unit width
- $g$ is gravitational acceleration
- $\mathbf{S_f}$ is friction slope
- $\mathbf{n}$ is Manning's roughness coefficient

## Numerical Solvers

LISFLOOD-FP 8.2 implements multiple numerical solvers, each with different characteristics:

### Acceleration Solver (ACC)

- **Mathematical Basis**: Local inertial approximation of the shallow water equations
- **Discretization**: Finite difference on a staggered grid
- **Time Integration**: Semi-implicit scheme
- **Stability**: Conditionally stable with CFL condition
- **Key Feature**: Efficient computation with good balance of accuracy and speed
- **Implementation**: Both CPU and GPU versions available
- **Best For**: Large-scale simulations where full momentum is not critical

The ACC solver uses a simplified momentum equation that neglects the convective acceleration term but retains the local acceleration and pressure terms. The semi-implicit time-stepping scheme provides good stability while allowing larger time steps than explicit schemes.

#### Key Equations:

Flow calculation:

$$q_{i+1/2,j}^{n+1} = \frac{q_{i+1/2,j}^{n} - g h_{flow} \Delta t \frac{\partial z}{\partial x}}{1 + g \Delta t n^2 |q| / h_{flow}^{7/3}}$$

Where:
- $q$ is the flow per unit width
- $h_{flow}$ is the maximum depth between adjacent cells
- $n$ is Manning's roughness coefficient
- $\Delta t$ is the time step

### Finite Volume Solvers (FV1/FV2)

- **Mathematical Basis**: Full shallow water equations
- **Discretization**: Finite volume method
- **Flux Calculation**: HLL approximate Riemann solver
- **Time Integration**: 
  - FV1: Euler explicit
  - FV2: Second-order Runge-Kutta
- **Slope Reconstruction**: 
  - FV1: Piecewise constant
  - FV2: Piecewise linear with slope limiters
- **Implementation**: CUDA-based GPU implementation
- **Best For**: Scenarios requiring full dynamic wave modeling

The FV solvers use the Godunov-type finite volume method to solve the full shallow water equations. FV2 extends FV1 with higher-order accuracy through piecewise linear reconstruction and slope limiting.

#### Key Features:

- Well-balanced for preserving steady states
- Conservative treatment of mass and momentum
- Robust handling of wetting and drying
- Shock-capturing capabilities

### Discontinuous Galerkin Solver (DG2)

- **Mathematical Basis**: Full shallow water equations
- **Discretization**: Discontinuous Galerkin (modal basis functions)
- **Flux Calculation**: HLL approximate Riemann solver
- **Time Integration**: Second-order Runge-Kutta
- **Slope Limiting**: Slope-based limiter with Krivodonova shock detector
- **Implementation**: CUDA-based GPU implementation
- **Best For**: Problems requiring high accuracy and sharp shock capturing

The DG2 solver represents the solution within each cell using a local polynomial expansion (up to second order). This provides higher-order accuracy while maintaining element-wise conservation.

#### Key Features:

- Higher-order accuracy on smooth solutions
- Excellent conservation properties
- Better representation of complex flow structures
- Efficient implementation on GPUs

### Non-Uniform Grid Solver (ACC_NUGRID)

- **Mathematical Basis**: Local inertial approximation with multi-resolution grid
- **Grid Representation**: Quadtree data structure using Morton codes
- **Refinement Criteria**: Wavelet-based error indicators
- **Implementation**: CUDA-based GPU implementation
- **Best For**: Domains with varying resolution requirements

The ACC_NUGRID solver extends the acceleration solver to support adaptive mesh refinement. It dynamically refines the grid based on solution features, providing high resolution where needed while maintaining computational efficiency.

## GPU Implementation

LISFLOOD-FP 8.2 features comprehensive GPU implementation using NVIDIA CUDA, enabling significant performance gains over CPU-only execution.

### Implementation Strategy

- **Multiple Kernels**: Different computational tasks are executed as separate CUDA kernels
- **Memory Management**: Careful management of device memory to minimize transfers
- **Grid-Block Structure**: Optimized for CUDA execution model
- **Solver-Specific Optimizations**: Each solver has tailored GPU implementations

### Key Aspects

1. **Memory Layout**: 
   - Structured grid data aligned for coalesced memory access
   - Padded arrays to avoid bank conflicts
   - Unified memory for simplified data management in newer GPUs

2. **Computational Kernels**:
   - Flux calculation kernels
   - Time integration kernels
   - Boundary condition kernels
   - Maximum field calculation kernels
   - Wetting/drying treatment kernels

3. **Optimization Techniques**:
   - Shared memory usage for frequently accessed data
   - Minimized thread divergence
   - Overlapping computation and memory transfers
   - Persistent thread blocks for iterative operations

4. **Load Balancing**:
   - Dynamic load balancing for non-uniform problems
   - Work stealing for adaptive grids
   - Efficient handling of sparse computation domains

### Performance Considerations

- **Memory Bandwidth**: LISFLOOD-FP performance on GPUs is primarily memory-bandwidth limited
- **Compute Intensity**: Higher-order methods (DG2) have better compute-to-memory ratio
- **GPU Architecture**: Performance scales with memory bandwidth and compute capability
- **Precision**: Double precision calculations are supported but may reduce performance

## Dynamic Grid Adaptation

The dynamic grid adaptation feature in LISFLOOD-FP 8.2 allows for multi-resolution simulations, focusing computational resources where they are most needed.

### Key Components

1. **Grid Representation**:
   - Quadtree-based grid structure
   - Morton code (Z-order curve) indexing for efficient neighbor finding
   - Maximum refinement level limit configurable by user

2. **Refinement Criteria**:
   - Wavelet coefficient magnitude
   - Gradient-based error indicators
   - Feature detection (wet/dry boundaries, hydraulic jumps)
   - Predefined regions (boundaries, structures)

3. **Grid Operations**:
   - Refinement: Subdivision of cells based on criteria
   - Coarsening: Merging of cells when high resolution no longer needed
   - Load balancing: Redistribution of computational load for parallel efficiency

4. **Multi-Resolution Flux Calculation**:
   - Conservative flux transfer between different resolution levels
   - Special handling of resolution transitions

### Implementation Details

- **Data Structure**: Specialized classes for managing non-uniform grid topology
- **Memory Management**: Dynamic allocation and compaction of grid data
- **Algorithm Adaptations**: Modified numerical schemes for handling resolution transitions
- **Performance Optimization**: Specialized kernels for efficient GPU execution

## Sub-Grid Channel Model

The Sub-Grid Channel (SGC) model in LISFLOOD-FP provides an efficient way to represent river channels that are smaller than the DEM resolution.

### Model Formulation

The SGC model parameterizes sub-grid scale flow using:

1. **Channel Geometry**:
   - Channel width
   - Bank elevation
   - Bed elevation

2. **Flow Parameterization**:
   - Channel cross-section shape (rectangular, trapezoidal, parabolic)
   - Manning's equation for friction
   - Cross-section area / width relationship

3. **Connection to Floodplain**:
   - Bankfull elevation as threshold for floodplain interaction
   - Conservative mass transfer between channel and floodplain

### Key Equations

Channel flow calculation:

$$Q = \frac{A^{5/3} R^{2/3}}{n} \sqrt{S}$$

Where:
- $A$ is the cross-sectional area
- $R$ is the hydraulic radius
- $n$ is Manning's roughness coefficient
- $S$ is the water surface slope

The relationship between channel depth $h_c$ and channel area $A$ depends on the channel type:

1. **Rectangular** ($t=1$):
   $A = w \cdot h_c$

2. **Trapezoidal** ($t=2$):
   $A = (w + s \cdot h_c) \cdot h_c$

3. **Parabolic** ($t=3$):
   $A = \frac{2}{3} \cdot w \cdot h_c$

Where:
- $w$ is the channel width
- $s$ is the side slope for trapezoidal channels

## Hydraulic Structures

LISFLOOD-FP 8.2 includes comprehensive modeling of hydraulic structures such as weirs, bridges, and culverts.

### Weir Modeling

Weirs are modeled using standard weir equations for both free flow and drowned (submerged) flow conditions:

1. **Free Flow**:
   $$Q = C_d \cdot w \cdot h_u^{1.5}$$

2. **Drowned Flow**:
   $$Q = C_d \cdot w \cdot h_u \cdot \frac{\sqrt{h_u - h_d}}{\sqrt{m}}$$

Where:
- $Q$ is the discharge
- $C_d$ is the discharge coefficient
- $w$ is the weir width
- $h_u$ is the upstream head above the weir crest
- $h_d$ is the downstream head above the weir crest
- $m$ is the modular limit (typically 0.7)

### Bridge Modeling

Bridges are modeled as a combination of weir flow (over the bridge deck) and orifice flow (through the bridge opening):

1. **Orifice Flow**:
   $$Q = C_d \cdot A \cdot \sqrt{2g \cdot \Delta h}$$

2. **Weir Flow** (when overtopped):
   Same as weir equations above

Where:
- $A$ is the bridge opening area
- $\Delta h$ is the head difference
- $g$ is gravitational acceleration

### Implementation Details

- **Structure Representation**: Point-based locations with attributes
- **Flow Direction Control**: Optional one-way flow constraint
- **Integration with Grid Flow**: Special handling at structure locations
- **Subgrid Implementation**: Compatible with subgrid channel model

## Dam Operations

LISFLOOD-FP 8.2 includes capabilities for modeling dam operations and reservoir routing.

### Dam Model Components

1. **Physical Representation**:
   - Dam location
   - Reservoir volume-elevation relationship
   - Spillway characteristics
   - Operating rules

2. **Operating Logic**:
   - Multiple operation schemes available
   - Rule-based releases
   - Target level management
   - Flood control operations

3. **Flow Calculations**:
   - Controlled releases based on operating rules
   - Spillway flow using weir equations
   - Mass conservation for reservoir volume

### Operation Types

The model supports multiple operation schemes:

1. **Constant Release**: Fixed discharge if available
2. **Storage-Based Release**: Release based on current storage
3. **Linear Storage Release**: Linear relationship with storage
4. **WaterGAP Model**: Based on Döll et al. (2003)
5. **PCR-GLOBWB Model**: Based on Wada et al. (2014)
6. **WBMplus Model**: Based on Wisser et al. (2010)
7. **Hanasaki Model**: Based on Hanasaki et al. (2005)
8. **Custom Rule Curve**: User-defined monthly release targets
9. **Forecast-Based**: Operations based on inflow forecasts

### Implementation Details

- **Volume Accounting**: Precise tracking of inflows, outflows, and storage
- **Integration with Grid**: Special treatment of dam cells
- **Cascade Systems**: Support for multiple dams in a river system
- **Output Metrics**: Comprehensive reporting of dam operations

## Hazard Assessment

LISFLOOD-FP 8.2 includes integrated hazard assessment capabilities for flood risk analysis.

### Hazard Metrics

1. **Depth-Velocity Product**:
   $$H = d \cdot (v + 0.5)$$
   Where:
   - $H$ is the hazard rating
   - $d$ is the water depth
   - $v$ is the flow velocity

2. **Debris Factor Inclusion**:
   $$H = d \cdot (v + 0.5) + DF$$
   Where $DF$ is the debris factor:
   - $DF = 0$ for $d < 0.25m$
   - $DF = 0.5$ for $0.25m \leq d < 0.75m$
   - $DF = 1.0$ for $d \geq 0.75m$

3. **Hazard Classification**:
   - Class 1 (Low): $H < 0.75$
   - Class 2 (Moderate): $0.75 \leq H < 1.25$
   - Class 3 (Significant): $1.25 \leq H < 2.0$
   - Class 4 (Extreme): $H \geq 2.0$

### Implementation Details

- **Velocity Calculation**: Derived from flow and depth fields
- **Maximum Hazard Tracking**: Records maximum hazard during simulation
- **Output Options**: Hazard maps in standard output formats
- **Alternative Formulations**: Multiple hazard calculation methods available

## Input/Output Formats

LISFLOOD-FP 8.2 supports multiple file formats for input and output data.

### Input Formats

1. **Raster Data**:
   - ASCII Grid (.asc)
   - Binary Grid
   - NetCDF (.nc)
   - GeoTIFF (.tif) via GDAL

2. **Time Series**:
   - Text files (.csv, .txt)
   - NetCDF (.nc)

3. **Vector Data**:
   - Text-based formats for structures, gauges, etc.
   - Simple coordinate lists

4. **Configuration**:
   - Parameter file (.par)
   - Command line options

### Output Formats

1. **Raster Results**:
   - ASCII Grid (.asc)
   - Binary Grid
   - NetCDF (.nc)

2. **Time Series**:
   - Text files (.stage, .discharge)
   - NetCDF (.nc)

3. **Checkpoint Files**:
   - Binary format for model state

### Special Format Features

- **Headers**: ASCII files contain projection information
- **Compression**: Options for compressed outputs
- **Precision Control**: Configurable precision for outputs
- **Time-Varying Outputs**: Support for time-indexed series

## Performance Considerations

LISFLOOD-FP 8.2 is designed for high-performance computing, with several factors affecting simulation speed:

### Key Performance Factors

1. **Domain Size**: Number of computational cells
2. **Solver Choice**: 
   - ACC is fastest but least physically complete
   - FV1 balances speed and physics
   - FV2/DG2 provide higher accuracy at greater computational cost
3. **Hardware Selection**:
   - GPU acceleration provides 10-100x speedup
   - CPU performance varies with core count and clock speed
4. **Simulation Duration**: Linear scaling with time steps
5. **Optimization Settings**:
   - Wet/Dry area tracking
   - CFL number selection
   - Output frequency

### Optimization Techniques

1. **Domain Decomposition**:
   - Grid partitioning for parallel execution
   - Load balancing for non-uniform workloads

2. **Memory Optimization**:
   - Efficient data structures
   - Minimized data transfers
   - Cache-friendly access patterns
   - Safe memory management practices:
     - Null pointer checking before deletion
     - RAII pattern for automatic resource cleanup
     - Consistent allocation/deallocation patterns
     - Smart pointer usage for newer code

3. **Computational Optimizations**:
   - Vectorization (SIMD)
   - Thread parallelism (OpenMP)
   - GPU acceleration (CUDA)

4. **Algorithm Selection**:
   - Adaptive time stepping
   - Solution-adaptive spatial discretization

### Benchmarks

Typical performance benchmarks on modern hardware:

| Domain Size | Solver | Hardware | Performance |
|-------------|--------|----------|-------------|
| 1000x1000 | ACC | CPU (8 cores) | ~15 time steps/second |
| 1000x1000 | ACC | GPU (RTX 2080) | ~300 time steps/second |
| 1000x1000 | FV1 | GPU (RTX 2080) | ~80 time steps/second |
| 1000x1000 | DG2 | GPU (RTX 2080) | ~20 time steps/second |

Actual performance will vary based on specific simulation characteristics and hardware configurations.

## Code Structure

LISFLOOD-FP 8.2 is organized into several major components:

### Directory Structure

```
lisflood-fp_umair/
├── src/                    # Source code
│   ├── cuda/               # CUDA GPU solvers
│   │   ├── acc/            # Acceleration solver
│   │   ├── acc_nugrid/     # Non-uniform grid solver
│   │   ├── adaptive/       # Adaptive solver
│   │   ├── common/         # Common CUDA utilities
│   │   ├── dg2/            # Discontinuous Galerkin solver
│   │   ├── fv1/            # First-order finite volume
│   │   └── fv2/            # Second-order finite volume
│   ├── lisflood2/          # Core lisflood functionality
│   ├── rain/               # Rainfall handling
│   └── swe/                # Shallow water equation solvers
├── cmake/                  # CMake modules
├── config/                 # Configuration files
├── docs/                   # Documentation
├── examples/               # Example simulations
└── tests/                  # Test cases
```

### Key Components

1. **Core Library**:
   - Parameter handling
   - Input/output management
   - Memory allocation
   - Utility functions

### Coding Standards

The LISFLOOD-FP codebase follows these coding conventions:
1. **Naming Conventions**:
   - camelCase for variables and functions
   - PascalCase for classes and types
   - ALL_CAPS for constants and macros
2. **Memory Management**:
   - All pointers checked for nullptr before deletion
   - Smart pointers preferred for new code
   - RAII pattern used where appropriate
3. **Error Handling**:
   - Comprehensive error checks and meaningful messages
   - Graceful handling of exceptional conditions
4. **Platform Independence**:
   - Abstraction of platform-specific code
   - Conditional compilation using feature detection
   - Cross-platform filesystem and system operations

2. **Solvers**:
   - CPU-based implementations
   - GPU (CUDA) implementations
   - Numerical algorithms

3. **Special Features**:
   - Sub-grid channel model
   - Hydraulic structures
   - Dam operations
   - Hazard assessment

4. **Build System**:
   - CMake configuration
   - Platform-specific settings
   - Dependency management

## Development History

LISFLOOD-FP has evolved significantly since its first release:

### Major Versions

1. **LISFLOOD-FP 1.0** (2000):
   - Original diffusive wave model
   - CPU-only implementation

2. **LISFLOOD-FP 3.0** (2010):
   - Introduction of the acceleration solver
   - Improved numerical stability

3. **LISFLOOD-FP 5.0** (2013):
   - Sub-grid channel model
   - OpenMP parallelization

4. **LISFLOOD-FP 7.0** (2017):
   - Initial CUDA GPU implementation
   - Finite volume solvers

5. **LISFLOOD-FP 8.0** (2021):
   - Discontinuous Galerkin solver
   - Enhanced GPU implementation
   - Multi-core CPU implementation

6. **LISFLOOD-FP 8.2** (2023):
   - Dynamic grid adaptation
   - Enhanced dam operations
   - Comprehensive hazard assessment
   - Advanced weir functionality
   - CUDA solver consistency improvements
   - Cross-platform enhancements

### Future Directions

Ongoing and planned developments include:

1. **Performance Enhancements**:
   - Multi-GPU support
   - Enhanced vectorization

2. **Physical Process Extensions**:
   - Sediment transport
   - Water quality modeling
   - Urban drainage integration

3. **User Interface Improvements**:
   - GUI for model setup
   - Real-time visualization

4. **Integration Capabilities**:
   - Coupling with hydrological models
   - Integration with GIS frameworks
   - Web service interfaces

---

This technical documentation provides a comprehensive overview of the LISFLOOD-FP 8.2 model. For specific implementation details, refer to the source code and inline documentation.