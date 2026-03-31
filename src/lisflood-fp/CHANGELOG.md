# Changelog

All notable changes to LISFLOOD-FP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [8.2.0] - 2023-08-15

### Added
- Dynamic grid adaptation capabilities
- Advanced dam operation models
- Comprehensive hazard assessment
- Improved weir functionality
- NetCDF input/output support
- DG2 solver enhancements

### Fixed
- Memory management issues in several components
- Numerical stability in wetting/drying front
- Platform-specific code for better cross-platform support
- CUDA solver consistency improvements

### Changed
- Improved documentation
- Enhanced build system
- Streamlined configuration options
- Better error handling and messaging

## [8.1.0] - 2022-10-20

### Added
- CUDA FV2 solver
- Multi-core CPU optimizations
- Extended hydraulic structure support
- Enhanced checkpoint capabilities
- New diagnostics and statistics

### Fixed
- Multiple numerical stability issues
- Memory leaks in long simulations
- Boundary condition handling
- GPU resource management

## [8.0.0] - 2021-06-30

### Added
- Discontinuous Galerkin (DG2) solver
- CUDA GPU acceleration
- Sub-grid channel model improvements
- Finite volume solvers (FV1)
- Enhanced visualization options

### Changed
- Complete code restructuring
- New build system using CMake
- Improved memory management
- Enhanced documentation

### Removed
- Legacy ATS solver
- Deprecated visualization tools