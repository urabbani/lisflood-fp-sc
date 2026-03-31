/*
 * cuda_solver_selector.cpp
 *
 * CUDA solver selection implementation to ensure all solvers are CUDA-based
 */

#include "cuda_solver_selector.h"
#include <stdio.h>

namespace lis {
namespace solver {

int select_and_run_solver(Fnames& params, States& states, Pars& pars, Solver& solver, int verbose, int argc, char** argv) {
#ifdef CUDA
    if (states.cuda == ON) {
        // CUDA is enabled, use CUDA-based solvers
        if (states.fv1 == ON) {
            lis::cuda::fv1::Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.fv2 == ON) {
            lis::cuda::fv2::Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.dg2 == ON) {
            lis::cuda::dg2::Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.acceleration == ON) {
            lis::cuda::acc::Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.acc_nugrid == ON) {
            lis::cuda::acc_nugrid::Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.mwdg2 == ON || states.hwfv1 == ON) {
            lis::cuda::adaptive::Simulation simulation;
            simulation.run(argc, argv);
        }
        else {
            fprintf(stderr, "CUDA mode is enabled, but no CUDA-compatible solver was selected.\n");
            return -1;
        }
    }
    else {
        // CUDA is available but not enabled, wrap CPU solvers with CUDA interface for consistent behavior
        if (states.fv1 == ON) {
            lis::cuda::cpu_wrapper::FV1Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.dg2 == ON) {
            lis::cuda::cpu_wrapper::DG2Simulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else if (states.acceleration == ON) {
            lis::cuda::cpu_wrapper::AccelerationSimulation simulation;
            simulation.run(params, states, pars, solver, verbose);
        }
        else {
            fprintf(stderr, "No CUDA-compatible solver was selected.\n");
            return -1;
        }
    }
    return 0;
#else
    // CUDA is not available, use CPU-based solvers
    fprintf(stderr, "CUDA support is not available. Recompile with CUDA=ON for GPU acceleration.\n");
    return -1;
#endif
}

} // namespace solver
} // namespace lis