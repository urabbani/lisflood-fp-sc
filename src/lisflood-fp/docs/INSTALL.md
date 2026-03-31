# LISFLOOD-FP 8.2 Installation Guide

LISFLOOD-FP is a two-dimensional hydrodynamic flood inundation model designed for simulating floodplain inundation. This installation guide provides comprehensive instructions for building and installing LISFLOOD-FP on both Linux and Windows systems.

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Linux Installation](#linux-installation)
3. [Windows Installation](#windows-installation)
4. [Configuration Options](#configuration-options)
5. [Troubleshooting](#troubleshooting)

## Prerequisites

### Common Dependencies
- **CMake** (version 3.13 or higher)
- **C++ Compiler** with C++14 support (GCC 5+, Clang 3.4+, or MSVC 2017+)
- **Git** (optional, for version control)

### Optional Dependencies
- **CUDA Toolkit** (version 9.0 or higher) - For GPU acceleration
- **NetCDF library** - For NetCDF file support
- **OpenMP** - For multi-core CPU utilization
- **NUMA library** (Linux only) - For non-uniform memory access

## Linux Installation

### Getting the Source Code
Clone the repository or extract the source archive to your preferred location:

```bash
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair
```

### Using CMake (Recommended)

1. Create a build directory:
```bash
mkdir build
cd build
```

2. Configure with CMake:
```bash
cmake .. -D_CONFIG=../config/local
```

3. Build:
```bash
make -j$(nproc)
```

4. Install (optional):
```bash
sudo make install
```

### Environment Setup
You may want to add the LISFLOOD-FP executable directory to your PATH:

```bash
export PATH=$PATH:/path/to/lisflood-fp_umair/build
```

Add this line to your `~/.bashrc` or `~/.profile` file to make it permanent.

## Windows Installation

### Getting the Source Code
1. Download and extract the source code archive, or clone using Git:
```cmd
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair
```

### Using CMake with Visual Studio

1. Install Visual Studio with C++ development tools
2. Install CMake (3.13 or later)
3. For CUDA support, install NVIDIA CUDA Toolkit
4. For NetCDF support, the required libraries are included in the repository

#### Command Line Build
```cmd
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

#### Using Visual Studio GUI
1. Open Visual Studio
2. Select "Open a local folder" and navigate to the LISFLOOD-FP directory
3. Visual Studio should detect the CMakeLists.txt file and configure the project
4. Build using the Visual Studio interface

### Environment Setup
To add the LISFLOOD-FP executable to your PATH:
1. Go to System Properties > Advanced > Environment Variables
2. Edit the PATH variable and add the build directory path
3. Click OK to save changes

## Configuration Options

The default configuration options are in `config.default.cmake`. You can customize these by:

1. Creating your own configuration file (e.g., in the `config` directory)
2. Specifying it with the `-D_CONFIG=path/to/config` option when running CMake

### Common Configuration Options

| Option | Description | Default |
|--------|-------------|---------|
| `_NETCDF` | Enable NetCDF support | 1 (ON) |
| `_NUMERIC_MODE` | Numeric precision mode | 1 (double) |
| `_ONLY_RECT` | Only use rectangular grids | 1 (ON) |
| `_PROFILE_MODE` | Enable profiling | 0 (OFF) |
| `_DISABLE_WET_DRY` | Disable wet/dry optimization | 0 (OFF) |
| `_CALCULATE_Q_MODE` | Calculate discharge | 1 (ON) |
| `_SGM_BY_BLOCKS` | Use SGC by block mode | 0 (OFF) |
| `_BALANCE_TYPE` | Load balancing type | 0 (default) |

## Troubleshooting

### CUDA Issues
- Ensure you have compatible CUDA drivers installed
- Check CUDA_ARCHITECTURES setting in CMakeLists.txt matches your GPU
- Use `nvidia-smi` to verify your GPU is detected and working properly

### NetCDF Issues
- On Linux, ensure NetCDF libraries are properly installed
- On Windows, the required NetCDF libraries are included in the windep directory
- Check for environment variables like `NETCDF_DIR` or `NETCDF_HOME`

### Common Error: "NUMA library not found" (Linux)
- Install libnuma-dev package:
```bash
sudo apt-get install libnuma-dev  # For Debian/Ubuntu
sudo yum install numactl-devel    # For CentOS/RHEL/Fedora
```

### Compilation Errors
- Ensure your compiler supports C++14
- Update your build tools if necessary
- Check for any missing dependencies

### Performance Issues
- Ensure you're using the appropriate solver for your system
- Enable CUDA if you have a compatible GPU
- Check OpenMP settings for multi-core performance

## Getting Help

If you encounter issues not covered in this guide, please submit an issue on our GitHub repository or contact the development team.