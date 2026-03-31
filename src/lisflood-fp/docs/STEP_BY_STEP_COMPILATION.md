# Step-by-Step Compilation Guide for LISFLOOD-FP 8.2

This guide provides detailed, step-by-step instructions for compiling LISFLOOD-FP 8.2 on various operating systems, including Linux, Windows, macOS, and Windows Subsystem for Linux (WSL).

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Linux (Ubuntu/Debian)](#linux-ubuntudebian)
3. [Linux (CentOS/RHEL/Fedora)](#linux-centosrhelfedora)
4. [Windows with Visual Studio](#windows-with-visual-studio)
5. [Windows with MinGW](#windows-with-mingw)
6. [macOS](#macos)
7. [Windows Subsystem for Linux (WSL)](#windows-subsystem-for-linux-wsl)
8. [Common Issues and Troubleshooting](#common-issues-and-troubleshooting)

## Prerequisites

Regardless of your platform, you'll need the following:

- **C++ compiler** with C++14 support
- **CMake** (version 3.13 or higher)
- **Git** (optional, for version control)

Optional but recommended components:

- **CUDA Toolkit** (for GPU acceleration)
- **NetCDF library** (for NetCDF file support)
- **OpenMP** (for multi-core CPU utilization)

## Linux (Ubuntu/Debian)

### 1. Install Required Packages

```bash
# Update package lists
sudo apt update

# Install essential build tools
sudo apt install build-essential cmake git

# Install NUMA library
sudo apt install libnuma-dev

# Install NetCDF (optional)
sudo apt install libnetcdf-dev

# Install CUDA (optional, for GPU acceleration)
# First add the NVIDIA repository
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"
sudo apt update
sudo apt install cuda
```

### 2. Clone or Extract LISFLOOD-FP

If using Git:
```bash
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair
```

Or extract from the provided archive:
```bash
tar -xzf lisflood-fp_umair.tar.gz
cd lisflood-fp_umair
```

### 3. Configure with CMake

```bash
# Create a build directory
mkdir -p build
cd build

# Configure with CMake
cmake .. -D_CONFIG=../config/local
```

If you want to use a specific compiler:
```bash
cmake .. -D_CONFIG=../config/local -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_CXX_COMPILER=g++-10
```

### 4. Build

```bash
# Build using all available cores
make -j$(nproc)
```

### 5. Verify the Build

```bash
# Check if the executable was created
ls -l lisflood

# Try running with version flag
./lisflood -version
```

## Linux (CentOS/RHEL/Fedora)

### 1. Install Required Packages

```bash
# For CentOS/RHEL
sudo yum install gcc gcc-c++ make cmake git
sudo yum install numactl-devel

# For Fedora
sudo dnf install gcc gcc-c++ make cmake git
sudo dnf install numactl-devel

# Install NetCDF (optional)
sudo yum install netcdf-devel   # CentOS/RHEL
# or
sudo dnf install netcdf-devel   # Fedora

# Install CUDA (optional, for GPU acceleration)
# Go to https://developer.nvidia.com/cuda-downloads and follow instructions for your distribution
```

### 2-5. Follow the same steps as in Ubuntu/Debian section above

## Windows with Visual Studio

### 1. Install Required Software

1. **Install Visual Studio**:
   - Download and install [Visual Studio](https://visualstudio.microsoft.com/downloads/)
   - During installation, select the "Desktop development with C++" workload
   - Ensure C++14 support is enabled

2. **Install CMake**:
   - Download and install [CMake](https://cmake.org/download/)
   - Ensure to add CMake to your system PATH during installation

3. **Install CUDA (optional, for GPU acceleration)**:
   - Download and install [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads)
   - Follow NVIDIA's installation guide

4. **Install Git (optional)**:
   - Download and install [Git for Windows](https://git-scm.com/download/win)

### 2. Clone or Extract LISFLOOD-FP

Using Git Bash or Command Prompt:
```bash
# If using Git
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair

# Or extract the ZIP archive using Windows Explorer
```

### 3. Configure with CMake

Using Command Prompt or PowerShell:
```bash
# Create a build directory
mkdir build
cd build

# Configure with CMake (change the generator to match your Visual Studio version)
cmake .. -G "Visual Studio 17 2022" -A x64 -D_CONFIG=../config/local
```

For older versions of Visual Studio, use the appropriate generator:
- Visual Studio 2019: `"Visual Studio 16 2019"`
- Visual Studio 2017: `"Visual Studio 15 2017"`

### 4. Build

**Option 1**: Build from command line:
```bash
cmake --build . --config Release
```

**Option 2**: Build using Visual Studio:
1. Open the generated solution file (`.sln`) in Visual Studio
2. Select "Release" configuration from the dropdown
3. Right-click on the "lisflood" project and select "Build"

### 5. Verify the Build

The executable will be in the `Release` directory:
```bash
# Navigate to the Release directory
cd Release

# Run with version flag
lisflood.exe -version
```

## Windows with MinGW

### 1. Install Required Software

1. **Install MinGW-w64**:
   - Download [MinGW-w64](https://www.mingw-w64.org/downloads/) or use [MSYS2](https://www.msys2.org/)
   - Ensure you install a version that supports C++14
   - Add the MinGW `bin` directory to your system PATH

2. **Install CMake**:
   - Download and install [CMake](https://cmake.org/download/)
   - Add CMake to your system PATH during installation

3. **Install Git (optional)**:
   - Download and install [Git for Windows](https://git-scm.com/download/win)

### 2. Clone or Extract LISFLOOD-FP

Same as in the Visual Studio section.

### 3. Configure with CMake

```bash
# Create a build directory
mkdir build
cd build

# Configure with CMake
cmake .. -G "MinGW Makefiles" -D_CONFIG=../config/local
```

### 4. Build

```bash
# Build
mingw32-make -j4  # Adjust the number based on your CPU cores
```

### 5. Verify the Build

```bash
# Check if the executable was created
dir lisflood.exe

# Try running with version flag
lisflood.exe -version
```

## macOS

### 1. Install Required Packages

Using [Homebrew](https://brew.sh/):

```bash
# Install Homebrew if you don't have it
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install required packages
brew install cmake
brew install gcc@11  # or any other GCC version
brew install netcdf  # optional
brew install libomp  # for OpenMP support

# CUDA is not directly supported on newer Macs with Apple Silicon
# For Intel Macs, you can install CUDA following NVIDIA's instructions
```

### 2. Clone or Extract LISFLOOD-FP

```bash
# Using Git
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair

# Or extract from the provided archive
tar -xzf lisflood-fp_umair.tar.gz
cd lisflood-fp_umair
```

### 3. Configure with CMake

```bash
# Create a build directory
mkdir -p build
cd build

# Configure with CMake (using GCC instead of default Clang)
cmake .. -D_CONFIG=../config/local -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11
```

### 4. Build

```bash
# Build
make -j$(sysctl -n hw.ncpu)
```

### 5. Verify the Build

```bash
# Check if the executable was created
ls -l lisflood

# Try running with version flag
./lisflood -version
```

## Windows Subsystem for Linux (WSL)

WSL allows you to run a Linux environment directly on Windows. This is often the easiest way to build Linux software on a Windows machine.

### 1. Install WSL

Open PowerShell as Administrator and run:

```powershell
# Install WSL with Ubuntu
wsl --install -d Ubuntu
```

This will download and install Ubuntu on WSL. After installation completes, a Ubuntu terminal will open.

### 2. Set Up the Development Environment in WSL

In the Ubuntu terminal:

```bash
# Update package lists
sudo apt update

# Install development tools
sudo apt install build-essential cmake git

# Install NUMA library
sudo apt install libnuma-dev

# Install NetCDF (optional)
sudo apt install libnetcdf-dev

# Install CUDA (optional, for GPU acceleration)
# Note: NVIDIA drivers must be installed in Windows for this to work
# Go to https://developer.nvidia.com/cuda-downloads for specific WSL installation
```

### 3. Navigate to Your Windows Drives

WSL mounts Windows drives under `/mnt/`. For example, `C:` is mounted as `/mnt/c/`.

```bash
# Navigate to your Windows user directory
cd /mnt/c/Users/YourUsername/

# Create or navigate to your project directory
mkdir -p Projects
cd Projects
```

### 4. Clone or Extract LISFLOOD-FP

```bash
# Clone the repository (if using Git)
git clone https://github.com/your-username/lisflood-fp_umair.git
cd lisflood-fp_umair

# Or extract from archive (if you have downloaded it)
tar -xzf lisflood-fp_umair.tar.gz
cd lisflood-fp_umair
```

### 5. Configure and Build

```bash
# Create a build directory
mkdir -p build
cd build

# Configure with CMake
cmake .. -D_CONFIG=../config/local

# Build
make -j$(nproc)
```

### 6. Verify the Build

```bash
# Check if the executable was created
ls -l lisflood

# Try running with version flag
./lisflood -version
```

### 7. Running LISFLOOD-FP from Windows

You can run the WSL-built LISFLOOD-FP executable directly from Windows using:

```powershell
wsl ~/Projects/lisflood-fp_umair/build/lisflood -version
```

Or for running a simulation:

```powershell
wsl ~/Projects/lisflood-fp_umair/build/lisflood ~/Projects/lisflood-fp_umair/examples/simple_floodplain.par
```

## Common Issues and Troubleshooting

### CUDA-Related Issues

1. **"CUDA not found" error**:
   - Check if CUDA is properly installed: `nvcc --version`
   - Ensure the CUDA path is in your system PATH
   - Try specifying the CUDA path explicitly: `-DCUDA_TOOLKIT_ROOT_DIR=/path/to/cuda`

2. **CUDA architecture mismatch**:
   - Edit `CMakeLists.txt` and set `CMAKE_CUDA_ARCHITECTURES` to match your GPU
   - You can find your CUDA architecture using: `nvidia-smi --query-gpu=compute_cap --format=csv`

### Build Errors

1. **Compiler version issues**:
   - Ensure your compiler supports C++14
   - For GCC, version 5.0 or higher is required
   - For Visual Studio, 2017 or newer is required

2. **Missing libraries**:
   - Install the required development libraries
   - Check for specific error messages indicating missing headers or libraries

3. **CMake configuration errors**:
   - Ensure you're using CMake 3.13 or higher: `cmake --version`
   - Delete the build directory and start fresh if you've made changes to CMake files

### NetCDF Issues

1. **NetCDF library not found**:
   - Install the NetCDF development package
   - Set the NetCDF path explicitly: `-DNETCDF_DIR=/path/to/netcdf`
   - Disable NetCDF support if not needed: set `_NETCDF=0` in your config file

### NUMA Issues (Linux)

1. **"NUMA library not found" error**:
   - Install the NUMA development package:
     - Ubuntu/Debian: `sudo apt install libnuma-dev`
     - CentOS/RHEL: `sudo yum install numactl-devel`

### WSL-Specific Issues

1. **File permission problems**:
   - WSL and Windows have different file permission systems
   - Try running the build process entirely within the WSL filesystem
   - For data files, ensure they have appropriate read permissions

2. **GPU access from WSL**:
   - Only WSL2 supports GPU passthrough
   - You need Windows 10 build 20150 or higher
   - Install NVIDIA WSL-compatible drivers on Windows
   - Install CUDA toolkit for WSL inside the WSL environment

---

If you encounter any issues not covered here, check the error messages for specific details and search for solutions online or contact the development team for support.