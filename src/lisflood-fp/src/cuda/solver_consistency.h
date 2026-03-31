/*
 * solver_consistency.h
 *
 * Utilities for ensuring consistent behavior across different solver implementations
 */

#ifndef SOLVER_CONSISTENCY_H
#define SOLVER_CONSISTENCY_H

#include "cuda_solver.cuh"
#include "acc/cuda_acc_solver.cuh"
#include "fv1/cuda_fv1_solver.cuh"
#include "fv2/cuda_fv2_solver.cuh"
#include "dg2/cuda_dg2_solver.cuh"
#include "acc_nugrid/cuda_acc_nugrid_simulate.cuh"

namespace lis {
namespace cuda {
namespace solver_consistency {

/**
 * Checks if all solvers have consistent handling of boundary conditions
 * 
 * @return True if boundary conditions are consistently handled, false otherwise
 */
bool check_boundary_consistency();

/**
 * Checks if all solvers have consistent implementations of essential numerical operations
 * 
 * @return True if numerical operations are consistently implemented, false otherwise
 */
bool check_numerical_implementation_consistency();

/**
 * Checks if all solvers handle edge cases consistently (e.g., thin water)
 * 
 * @return True if edge cases are consistently handled, false otherwise
 */
bool check_edge_case_consistency();

/**
 * Checks for consistent mass conservation across solvers
 * 
 * @return True if mass conservation is consistently implemented, false otherwise
 */
bool check_mass_conservation_consistency();

/**
 * Utility class for validating solver implementations
 */
class SolverValidator {
public:
    /**
     * Validates that a solver calculates flux terms correctly
     * 
     * @param solver The solver to validate
     * @return True if flux calculations are valid, false otherwise
     */
    template <typename SolverType>
    static bool validate_flux_calculations(SolverType& solver);
    
    /**
     * Validates that a solver applies boundary conditions correctly
     * 
     * @param solver The solver to validate
     * @return True if boundary handling is valid, false otherwise
     */
    template <typename SolverType>
    static bool validate_boundary_handling(SolverType& solver);
    
    /**
     * Validates that a solver handles thin water depths correctly
     * 
     * @param solver The solver to validate
     * @return True if thin water handling is valid, false otherwise
     */
    template <typename SolverType>
    static bool validate_thin_water_handling(SolverType& solver);
    
    /**
     * Validates that a solver conserves mass correctly
     * 
     * @param solver The solver to validate
     * @return True if mass conservation is valid, false otherwise
     */
    template <typename SolverType>
    static bool validate_mass_conservation(SolverType& solver);
};

} // namespace solver_consistency
} // namespace cuda
} // namespace lis

#endif // SOLVER_CONSISTENCY_H