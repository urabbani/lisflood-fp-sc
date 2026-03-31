# LISFLOOD-FP 8.2

LISFLOOD-FP is a two-dimensional hydraulic model designed for simulating floodplain inundation over complex topography. This version includes significant enhancements including GPU acceleration, dynamic grid adaptation, and improved hydraulic structure modeling.

## Features

- Multiple numerical solvers with varying levels of physical representation:
  - Acceleration solver (ACC) for simplified momentum modeling
  - Finite volume solvers (FV1/FV2) for full shallow water equations
  - Discontinuous Galerkin solver (DG2) for higher-order accuracy
  - Non-uniform grid solver (ACC_NUGRID) for adaptive resolution

- GPU acceleration using NVIDIA CUDA for high-performance computing

- Sub-grid channel model for efficient river network representation

- Comprehensive hydraulic structure modeling (weirs, bridges, culverts)

- Dam operations and reservoir routing capabilities

- Dynamic grid adaptation for multi-resolution simulations

- Integrated hazard assessment for flood risk analysis

- Support for multiple input/output formats (ASCII, binary, NetCDF)

## Installation

See [docs/INSTALL.md](docs/INSTALL.md) for detailed installation instructions.

Quick start:

```bash
# Linux/macOS
mkdir build && cd build
cmake .. -D_CONFIG=../config/local
make -j

# Windows with Visual Studio
mkdir build && cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

## Documentation

- [Installation Guide](docs/INSTALL.md)
- [User Manual](docs/USER_MANUAL.md)
- [Parameter File Guide](docs/PARAMETER_FILE.md)
- [Technical Details](docs/TECHNICAL_DETAILS.md)

## Examples

The `examples` directory contains several example parameter files and input data for common simulation scenarios:

- Basic floodplain simulation
- Weir/levee configuration
- Dam breach modeling
- Urban flooding with rainfall

## Usage

Basic command line usage:

```bash
./lisflood -v parameter_file.par
```

Where:
- `-v` enables verbose output
- `parameter_file.par` is your parameter file

Common options:
- `-cuda` - Enable CUDA acceleration
- `-noclip` - Disable optimization of wet/dry regions
- `-overwrite` - Allow overwriting of existing results

## Docker Container

To run a simulation using the Docker container:

```bash
# Build the Docker image
docker build -t lisflood-fp .

# Run a simulation with GPU acceleration
./run-docker.sh -i /path/to/input -o /path/to/output -g -p simulation.par

# Run without GPU
./run-docker.sh -i /path/to/input -o /path/to/output -p simulation.par

## License

LISFLOOD-FP is distributed under the terms of the GNU GPL license. See the LICENSE file for details.

## Acknowledgements

LISFLOOD-FP was originally developed at the University of Bristol. This distribution includes significant enhancements including GPU acceleration, dynamic grid adaptation, and improved hydraulic structure modeling.

## Citation

If you use LISFLOOD-FP in your research, please cite:

Shaw, J., Kesserwani, G., Neal, J., Bates, P.D., Sharifian, M.K. (2021). LISFLOOD-FP 8.0: the new discontinuous Galerkin shallow-water solver for multi-core CPUs and GPUs. Geoscientific Model Development 14, 3577–3602.

## Contributing

We welcome contributions from the community! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on how to contribute to this project.

## Contact

For questions, support, or to report issues, please open an issue on the GitHub repository.