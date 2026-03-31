#include "cuda_safe_memory.h"
#include <iostream>

// Example implementation showing how to use the safe memory management classes
void cudaMemoryExample() {
    const size_t dataSize = 1024;
    
    try {
        // Create device memory using the RAII wrapper
        CudaMemoryManager<float> deviceData(dataSize, CudaMemoryManager<float>::DEVICE);
        
        // Initialize device memory to zero
        deviceData.memset(0);
        
        // Create host memory for input and output
        auto hostInput = makeCudaHostPtr<float>(dataSize);
        auto hostOutput = makeCudaHostPtr<float>(dataSize);
        
        // Initialize host input data
        for (size_t i = 0; i < dataSize; ++i) {
            hostInput.get()[i] = static_cast<float>(i);
        }
        
        // Copy data from host to device
        deviceData.copyFromHost(hostInput.get());
        
        // Execute kernel (would be called here)
        
        // Copy results back to host
        deviceData.copyToHost(hostOutput.get());
        
        // Managed memory example (easier to use)
        CudaMemoryManager<float> managedData(dataSize, CudaMemoryManager<float>::MANAGED);
        
        // With managed memory, we can access it directly from CPU code
        float* managedPtr = managedData.get();
        for (size_t i = 0; i < dataSize; ++i) {
            managedPtr[i] = static_cast<float>(i * 2);
        }
        
        // No explicit copy needed for managed memory - synchronize before kernel launch
        cudaDeviceSynchronize();
        
        // Execute kernel on managed memory (would be called here)
        
        // Synchronize after kernel to ensure CPU can access the data
        cudaDeviceSynchronize();
        
        // Access results directly
        float sum = 0.0f;
        for (size_t i = 0; i < dataSize; ++i) {
            sum += managedPtr[i];
        }
        
        std::cout << "Sum of managed memory array: " << sum << std::endl;
        
    } catch (const std::runtime_error& e) {
        std::cerr << "CUDA memory error: " << e.what() << std::endl;
    }
}

// Example of using the low-level smart pointers directly
void cudaSmartPointerExample() {
    try {
        // Allocate device memory using smart pointer
        auto deviceArray = makeCudaDevicePtr<int>(100);
        
        // Allocate host memory using smart pointer
        auto hostArray = makeCudaHostPtr<int>(100);
        
        // Initialize host data
        for (int i = 0; i < 100; ++i) {
            hostArray.get()[i] = i;
        }
        
        // Copy to device
        cudaSafeMemcpy(deviceArray.get(), hostArray.get(), 100, cudaMemcpyHostToDevice);
        
        // Memory is automatically freed when smart pointers go out of scope
    } catch (const std::runtime_error& e) {
        std::cerr << "CUDA error: " << e.what() << std::endl;
    }
}

// Demonstrates how to use these classes to replace legacy raw pointer methods
void migrateFromLegacyCode() {
    int* oldDevicePtr = nullptr;
    int* oldHostPtr = nullptr;
    
    // Legacy code would do:
    // cudaMalloc(&oldDevicePtr, 100 * sizeof(int));
    // cudaMallocHost(&oldHostPtr, 100 * sizeof(int));
    // ... use memory ...
    // cudaFree(oldDevicePtr);
    // cudaFreeHost(oldHostPtr);
    
    // New code with better memory safety:
    auto devicePtr = makeCudaDevicePtr<int>(100);
    auto hostPtr = makeCudaHostPtr<int>(100);
    
    // Use the raw pointers for compatibility with existing CUDA kernel calls
    int* rawDevicePtr = devicePtr.get();
    int* rawHostPtr = hostPtr.get();
    
    // Memory is freed automatically when smart pointers go out of scope
}