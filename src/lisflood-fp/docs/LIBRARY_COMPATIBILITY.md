# LISFLOOD-FP Library Compatibility Guide

This document provides information about compatibility with prerequisite libraries, recommended versions, and how to handle library updates.

## Table of Contents

1. [Library Dependencies](#library-dependencies)
2. [Version Compatibility Matrix](#version-compatibility-matrix)
3. [Handling Library Updates](#handling-library-updates)
4. [Troubleshooting Library Issues](#troubleshooting-library-issues)
5. [Contributing Library Updates](#contributing-library-updates)

## Library Dependencies

LISFLOOD-FP relies on several external libraries for various functions. This section outlines the main dependencies and their roles.

### Core Dependencies

| Library | Purpose | Required | Notes |
|---------|---------|----------|-------|
| CMake | Build system | Yes | Version 3.13+ required |
| C++ Standard Library | Core functionality | Yes | C++14 support required |
| NUMA | Non-uniform memory access | Linux only | Required on Linux systems |
| OpenMP | Multi-threading | Optional | Enables CPU parallelization |

### Optional Dependencies

| Library | Purpose | Required | Notes |
|---------|---------|----------|-------|
| CUDA | GPU acceleration | Optional | Significant performance benefits |
| NetCDF | File I/O | Optional | For NetCDF file format support |
| GDAL | Raster processing | Optional | For GeoTIFF and other format support |

## Version Compatibility Matrix

This matrix shows tested and compatible versions of key dependencies across different platforms.

### CUDA Compatibility

| LISFLOOD-FP Version | Minimum CUDA | Recommended CUDA | Maximum Tested | Notes |
|---------------------|--------------|------------------|----------------|-------|
| 8.2 | 9.0 | 11.x | 12.1 | CUDA 9.0+ required for C++14 support |
| 8.1 | 8.0 | 10.x | 11.6 | |
| 8.0 | 8.0 | 10.x | 11.0 | |

### Compiler Compatibility

| Platform | Compiler | Minimum Version | Recommended | Notes |
|----------|----------|-----------------|-------------|-------|
| Linux | GCC | 5.0 | 9.0+ | C++14 support required |
| Linux | Clang | 3.4 | 10.0+ | C++14 support required |
| Windows | MSVC | 2017 (19.10) | 2019+ | VS2017 or newer required |
| Windows | MinGW | 7.0 | 8.0+ | C++14 support required |
| macOS | Apple Clang | 8.0 | 12.0+ | C++14 support required |
| macOS | GCC | 6.0 | 10.0+ | Via Homebrew/MacPorts |

### CMake Compatibility

| LISFLOOD-FP Version | Minimum CMake | Recommended | Notes |
|---------------------|---------------|-------------|-------|
| 8.2 | 3.13 | 3.20+ | Older versions may not support CUDA properly |
| 8.1 | 3.10 | 3.15+ | |
| 8.0 | 3.8 | 3.15+ | |

### NetCDF Compatibility

| LISFLOOD-FP Version | Minimum NetCDF | Recommended | Notes |
|---------------------|----------------|-------------|-------|
| 8.2 | 4.4.0 | 4.8.0+ | API compatibility maintained |
| 8.1 | 4.3.0 | 4.6.0+ | |
| 8.0 | 4.3.0 | 4.5.0+ | |

## Handling Library Updates

LISFLOOD-FP includes mechanisms to handle library updates gracefully. This section explains how the code adapts to different library versions.

### Version Detection

The build system automatically detects the versions of installed libraries and adjusts compilation settings accordingly:

```cmake
# In CMakeLists.txt
if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 9.0)
  message(WARNING "CUDA 9.0 or newer is recommended for full C++14 support")
endif()

# OpenMP version-specific settings
if(OpenMP_CXX_VERSION VERSION_GREATER_EQUAL 4.5)
  add_compile_definitions(OPENMP_4_5_FEATURES)
endif()
```

### API Compatibility Layers

For libraries with changing APIs, we implement compatibility layers:

1. **NetCDF API Changes**: The `netcdf_io.cpp` file includes conditional compilation for different NetCDF API versions.

2. **CUDA API Evolution**: The CUDA implementation uses wrapper functions to handle API differences between versions:

```cpp
// In cuda_util.h
template<typename T>
void setDeviceMemory(T* dest, T value, size_t count) {
#if CUDA_VERSION >= 11000
  // Use newer API in CUDA 11+
  cudaMemset(dest, value, count * sizeof(T));
#else
  // Use older method for CUDA < 11
  cudaMemcpy(dest, &value, sizeof(T), cudaMemcpyHostToDevice);
  for(size_t i = 1; i < count; i++) {
    cudaMemcpy(dest + i, dest, sizeof(T), cudaMemcpyDeviceToDevice);
  }
#endif
}
```

3. **Compiler Differences**: The code handles compiler-specific features using preprocessor macros:

```cpp
#ifdef _MSC_VER
  // Microsoft Visual C++ specific code
#elif defined(__GNUC__)
  // GCC specific code
#elif defined(__clang__)
  // Clang specific code
#endif
```

### Dynamic Feature Detection

For optional features that depend on library capabilities, LISFLOOD-FP uses runtime checks:

```cpp
bool checkCUDACapabilities() {
  int deviceCount = 0;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) return false;
  
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, 0);
  
  // Check for minimum compute capability
  if (deviceProp.major < 3) return false;
  
  // Check for specific features
  return true;
}
```

## Troubleshooting Library Issues

This section provides guidance for diagnosing and resolving issues related to library updates.

### Common Symptoms and Solutions

1. **Build Failures After Library Update**

   Symptom: Compilation errors mentioning missing functions or incompatible types
   
   Solutions:
   - Clear CMake cache: `rm -rf build/* && cmake ..`
   - Check error messages for specific API changes
   - Verify library versions match compatibility matrix
   - Consider downgrading to a known compatible version

2. **Runtime Errors with Library Loading**

   Symptom: "Library not found" or "Symbol not found" errors
   
   Solutions:
   - Check library paths: `ldd lisflood` (Linux) or `dumpbin /dependents lisflood.exe` (Windows)
   - Ensure runtime libraries match build-time libraries
   - Set LD_LIBRARY_PATH (Linux) or PATH (Windows) appropriately

3. **Performance Degradation After Update**

   Symptom: Slower execution after updating libraries
   
   Solutions:
   - Check for debug vs. release builds of libraries
   - Verify optimization settings in CMake
   - Compare with previous versions using benchmarks
   - Check for library-specific performance issues in release notes

### Diagnosing Version Mismatches

Use these commands to verify library versions:

```bash
# CUDA version
nvcc --version

# GCC version
gcc --version

# NetCDF version
nc-config --version  # Linux/macOS
ncdump --version     # Windows

# Check dynamic libraries (Linux)
ldd lisflood | grep cuda
ldd lisflood | grep netcdf

# Check dynamic libraries (macOS)
otool -L lisflood | grep cuda
otool -L lisflood | grep netcdf

# Check dynamic libraries (Windows)
dumpbin /dependents lisflood.exe
```

## Contributing Library Updates

When updating library compatibility, please follow these guidelines:

1. **Testing Protocol**:
   - Test with the minimum supported version
   - Test with the latest stable version
   - Document any issues or workarounds

2. **Update Documentation**:
   - Update the compatibility matrix in this document
   - Document any API changes that required code changes
   - Note any performance implications

3. **Code Changes**:
   - Use conditional compilation for version-specific code
   - Add clear comments explaining version differences
   - Minimize changes to core algorithms
   - Prefer compatibility layers over extensive code changes

4. **Pull Request Process**:
   - Include information about the library version tested
   - Include before/after performance metrics if relevant
   - Reference any relevant library documentation

---

By following these guidelines, we can maintain LISFLOOD-FP's compatibility with evolving library ecosystems while ensuring stability and performance.