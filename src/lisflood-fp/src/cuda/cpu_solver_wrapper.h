/*
 * cpu_solver_wrapper.h
 *
 * CUDA wrapper for CPU-based solvers to ensure all solvers are CUDA-based
 */

#ifndef CPU_SOLVER_WRAPPER_H
#define CPU_SOLVER_WRAPPER_H

#ifdef __CUDACC__

#include "../lisflood.h"
#include "cuda_solver.cuh"
#include "../swe/fv1.h"
#include "../swe/dg2.h"

namespace lis {
namespace cuda {
namespace cpu_wrapper {

/**
 * CUDA wrapper for CPU-based acceleration solver
 */
class AccelerationSimulation {
public:
    /**
     * Run simulation using CPU-based acceleration solver wrapped in CUDA interface
     * 
     * @param params Simulation parameters
     * @param states Simulation states
     * @param solver Solver parameters
     * @param verbose Verbose mode flag
     * @return Maximum water depth
     */
    template<typename F>
    NUMERIC_TYPE run(Fnames& params, States& states, Pars& pars, Solver<F>& solver, int verbose);

private:
    /**
     * Initialize CUDA memory for CPU solver data
     */
    void initialize_cuda_memory();

    /**
     * Copy data from CPU to GPU
     */
    void copy_to_device();

    /**
     * Copy results from GPU back to CPU
     */
    void copy_from_device();

    /**
     * Clean up CUDA memory
     */
    void cleanup();
};

/**
 * CUDA wrapper for CPU-based FV1 solver
 */
class FV1Simulation {
public:
    /**
     * Run simulation using CPU-based FV1 solver wrapped in CUDA interface
     * 
     * @param params Simulation parameters
     * @param states Simulation states
     * @param solver Solver parameters
     * @param verbose Verbose mode flag
     * @return Maximum water depth
     */
    template<typename F>
    NUMERIC_TYPE run(Fnames& params, States& states, Pars& pars, Solver<F>& solver, int verbose);

private:
    /**
     * Initialize CUDA memory for CPU solver data
     */
    void initialize_cuda_memory();

    /**
     * Copy data from CPU to GPU
     */
    void copy_to_device();

    /**
     * Copy results from GPU back to CPU
     */
    void copy_from_device();

    /**
     * Clean up CUDA memory
     */
    void cleanup();
};

/**
 * CUDA wrapper for CPU-based DG2 solver
 */
class DG2Simulation {
public:
    /**
     * Run simulation using CPU-based DG2 solver wrapped in CUDA interface
     * 
     * @param params Simulation parameters
     * @param states Simulation states
     * @param solver Solver parameters
     * @param verbose Verbose mode flag
     * @return Maximum water depth
     */
    template<typename F>
    NUMERIC_TYPE run(Fnames& params, States& states, Pars& pars, Solver<F>& solver, int verbose);

private:
    /**
     * Initialize CUDA memory for CPU solver data
     */
    void initialize_cuda_memory();

    /**
     * Copy data from CPU to GPU
     */
    void copy_to_device();

    /**
     * Copy results from GPU back to CPU
     */
    void copy_from_device();

    /**
     * Clean up CUDA memory
     */
    void cleanup();
};

} // namespace cpu_wrapper
} // namespace cuda
} // namespace lis

#endif // CPU_SOLVER_WRAPPER_H
#endif // __CUDACC__
