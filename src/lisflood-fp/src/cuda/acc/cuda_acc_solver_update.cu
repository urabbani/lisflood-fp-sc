#include "cuda_acc_solver.cuh"
#include "cuda_boundary.cuh"
#include "cuda_hll.cuh"
#include "cuda_solver.cuh"
#include <algorithm>

// Kernel to zero velocity gradients in shallow water to improve stability
__global__ void zero_thin_depth_slopes_kernel
(
    Flow U,
    NUMERIC_TYPE depth_thresh
)
{
    int global_i = blockIdx.x * blockDim.x + threadIdx.x;
    int global_j = blockIdx.y * blockDim.y + threadIdx.y;
    
    for (int j = global_j; j < cuda::geometry.ysz + 2; j += blockDim.y * gridDim.y) {
        for (int i = global_i; i < cuda::pitch; i += blockDim.x * gridDim.x) {
            const int idx = j * cuda::pitch + i;
            
            // Check if water depth is below threshold for velocity calculation
            if (U.H[idx] <= depth_thresh * C(1.1)) {
                // Zero the velocity to prevent instability
                U.Qx[idx] = C(0.0);
                U.Qy[idx] = C(0.0);
                
                // Zero velocity fields if present
                if (cuda::solver_params.voutput == ON) {
                    U.Vx[idx] = C(0.0);
                    U.Vy[idx] = C(0.0);
                }
            }
        }
    }
}