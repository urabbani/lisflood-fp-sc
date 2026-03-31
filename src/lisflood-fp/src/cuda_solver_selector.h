/*
 * cuda_solver_selector.h
 *
 * CUDA solver selection to ensure all solvers are CUDA-based
 */

#ifndef CUDA_SOLVER_SELECTOR_H
#define CUDA_SOLVER_SELECTOR_H

#include "lisflood.h"

#ifdef CUDA
#include "cuda/acc/cuda_acc_simulate.cuh"
#include "cuda/fv1/cuda_fv1_simulate.cuh"
#include "cuda/dg2/cuda_dg2_simulate.cuh"
#include "cuda/fv2/cuda_fv2_simulate.cuh"
#include "cuda/acc_nugrid/cuda_acc_nugrid_simulate.cuh"
#include "cuda/adaptive/cuda_adaptive_simulate.cuh"
#include "cuda/cpu_solver_wrapper.h"
#endif

namespace lis {
namespace solver {

/**
 * Selects the appropriate solver based on simulation states and runs the simulation
 * 
 * @param params Simulation parameters
 * @param states Simulation states
 * @param pars Physical parameters
 * @param solver Solver parameters
 * @param verbose Verbose mode flag
 * @param argc Command line argument count
 * @param argv Command line arguments
 * @return Simulation result code (0 for success)
 */
int select_and_run_solver(Fnames& params, States& states, Pars& pars, Solver& solver, int verbose, int argc, char** argv);

} // namespace solver
} // namespace lis

#endif // CUDA_SOLVER_SELECTOR_H