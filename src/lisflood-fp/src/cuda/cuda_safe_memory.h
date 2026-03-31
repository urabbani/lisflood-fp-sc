#pragma once

#include <memory>
#include <utility>
#include <stdexcept>
#include <string>
#include <cuda_runtime.h>

// CUDA Error checking helper
inline void checkCudaError(cudaError_t result, const char* file, int line) {
    if (result != cudaSuccess) {
        throw std::runtime_error(
            std::string("CUDA Error: ") + 
            cudaGetErrorString(result) + 
            " at " + file + ":" + std::to_string(line)
        );
    }
}

#define CUDA_CHECK(call) checkCudaError(call, __FILE__, __LINE__)

// Deleter functors for smart pointers
struct CudaDeleter {
    void operator()(void* ptr) {
        if (ptr) {
            cudaFree(ptr);
        }
    }
};

struct CudaHostDeleter {
    void operator()(void* ptr) {
        if (ptr) {
            cudaFreeHost(ptr);
        }
    }
};

struct CudaManagedDeleter {
    void operator()(void* ptr) {
        if (ptr) {
            cudaFree(ptr); // Same as device memory
        }
    }
};

// Smart pointer type aliases for different CUDA memory types
template <typename T>
using CudaDevicePtr = std::unique_ptr<T, CudaDeleter>;

template <typename T>
using CudaHostPtr = std::unique_ptr<T, CudaHostDeleter>;

template <typename T>
using CudaManagedPtr = std::unique_ptr<T, CudaManagedDeleter>;

// Factory functions for creating CUDA memory
template <typename T>
CudaDevicePtr<T> makeCudaDevicePtr(size_t count = 1) {
    T* ptr = nullptr;
    CUDA_CHECK(cudaMalloc(&ptr, count * sizeof(T)));
    return CudaDevicePtr<T>(ptr);
}

template <typename T>
CudaHostPtr<T> makeCudaHostPtr(size_t count = 1) {
    T* ptr = nullptr;
    CUDA_CHECK(cudaMallocHost(&ptr, count * sizeof(T)));
    return CudaHostPtr<T>(ptr);
}

template <typename T>
CudaManagedPtr<T> makeCudaManagedPtr(size_t count = 1) {
    T* ptr = nullptr;
    CUDA_CHECK(cudaMallocManaged(&ptr, count * sizeof(T)));
    return CudaManagedPtr<T>(ptr);
}

// Memory operations with error checking
template <typename T>
void cudaSafeMemset(T* ptr, int value, size_t count) {
    CUDA_CHECK(cudaMemset(ptr, value, count * sizeof(T)));
}

template <typename T>
void cudaSafeMemcpy(T* dst, const T* src, size_t count, cudaMemcpyKind kind) {
    CUDA_CHECK(cudaMemcpy(dst, src, count * sizeof(T), kind));
}

// RAII class for CUDA memory
template <typename T>
class CudaMemoryManager {
public:
    enum MemoryType {
        DEVICE,
        HOST,
        MANAGED
    };

    CudaMemoryManager(size_t count, MemoryType type = DEVICE) 
        : m_count(count), m_type(type) {
        allocate();
    }

    ~CudaMemoryManager() = default; // Smart pointers handle cleanup

    // Get raw pointer
    T* get() const {
        switch (m_type) {
            case DEVICE: return m_devicePtr.get();
            case HOST: return m_hostPtr.get();
            case MANAGED: return m_managedPtr.get();
            default: return nullptr;
        }
    }

    // Reset with new memory of same size and type
    void reset() {
        allocate();
    }

    // Set memory to a value
    void memset(int value) {
        cudaSafeMemset(get(), value, m_count);
    }

    // Copy memory from host
    void copyFromHost(const T* src) {
        switch (m_type) {
            case DEVICE:
                cudaSafeMemcpy(m_devicePtr.get(), src, m_count, cudaMemcpyHostToDevice);
                break;
            case HOST:
                cudaSafeMemcpy(m_hostPtr.get(), src, m_count, cudaMemcpyHostToHost);
                break;
            case MANAGED:
                cudaSafeMemcpy(m_managedPtr.get(), src, m_count, cudaMemcpyHostToDevice);
                break;
        }
    }

    // Copy memory to host
    void copyToHost(T* dst) const {
        switch (m_type) {
            case DEVICE:
                cudaSafeMemcpy(dst, m_devicePtr.get(), m_count, cudaMemcpyDeviceToHost);
                break;
            case HOST:
                cudaSafeMemcpy(dst, m_hostPtr.get(), m_count, cudaMemcpyHostToHost);
                break;
            case MANAGED:
                cudaSafeMemcpy(dst, m_managedPtr.get(), m_count, cudaMemcpyDeviceToHost);
                break;
        }
    }

    // Get the count
    size_t size() const { return m_count; }

    // Get memory type
    MemoryType type() const { return m_type; }

private:
    void allocate() {
        switch (m_type) {
            case DEVICE:
                m_devicePtr = makeCudaDevicePtr<T>(m_count);
                break;
            case HOST:
                m_hostPtr = makeCudaHostPtr<T>(m_count);
                break;
            case MANAGED:
                m_managedPtr = makeCudaManagedPtr<T>(m_count);
                break;
        }
    }

    size_t m_count;
    MemoryType m_type;
    CudaDevicePtr<T> m_devicePtr;
    CudaHostPtr<T> m_hostPtr;
    CudaManagedPtr<T> m_managedPtr;
};