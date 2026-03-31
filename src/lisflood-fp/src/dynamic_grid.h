/*
 * dynamic_grid.h
 *
 * Header file for dynamic grid functionality in LISFLOOD-FP
 * Provides a flexible framework for adaptive mesh refinement (AMR)
 * 
 */

#ifndef DYNAMIC_GRID_H
#define DYNAMIC_GRID_H

#include "lisflood.h"
#include "cuda/acc_nugrid/MortonCode.h"
#include "cuda/acc_nugrid/AssembledSolution.h"
#include "cuda/acc_nugrid/NonUniformNeighbours.h"
#include "cuda/acc_nugrid/Details.h"
#include "cuda/acc_nugrid/NonUniformInterfaces.h"

namespace dynamic_grid {

/**
 * Grid refinement criteria
 */
enum RefinementCriteria {
    DEPTH_GRADIENT = 0,      // Refine based on depth gradient
    VELOCITY_GRADIENT = 1,   // Refine based on velocity gradient
    COMBINED_GRADIENT = 2,   // Refine based on combined gradients
    WAVELET_COEFFICIENT = 3, // Refine based on wavelet coefficients
    FIXED_MASK = 4,          // Refine based on fixed regions
    TIME_VARYING_MASK = 5    // Refine based on time-varying regions
};

/**
 * Configuration parameters for dynamic grid
 */
struct DynamicGridConfig {
    int max_refinement_level;             // Maximum refinement level
    RefinementCriteria refinement_criteria;  // Refinement criteria to use
    NUMERIC_TYPE refinement_threshold;    // Threshold for refinement
    NUMERIC_TYPE coarsening_threshold;    // Threshold for coarsening
    bool enable_dynamic_load_balancing;   // Enable dynamic load balancing
    bool refine_boundaries;               // Refine grid around boundaries
    bool refine_point_sources;            // Refine grid around point sources
    bool refine_gauge_points;             // Refine grid around gauge points
    int buffer_zone_width;                // Width of buffer zone around refined areas
    int update_frequency;                 // Grid update frequency (in timesteps)
};

/**
 * Initialize dynamic grid
 * @param params Simulation parameters
 * @param config Dynamic grid configuration
 * @param topo_grid Topography grid
 * @param h_grid Water depth grid
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 * @return Assembled solution for dynamic grid
 */
lis::cuda::acc_nugrid::AssembledSolution initializeDynamicGrid(
    Pars* params,
    DynamicGridConfig* config,
    NUMERIC_TYPE* topo_grid,
    NUMERIC_TYPE* h_grid,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Update dynamic grid based on refinement criteria
 * @param current_time Current simulation time
 * @param config Dynamic grid configuration
 * @param assembled_solution Assembled solution
 * @param h_grid Water depth grid
 * @param v_grid Velocity grid
 * @param timestep_count Current timestep count
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 * @return True if grid was updated, false otherwise
 */
bool updateDynamicGrid(
    NUMERIC_TYPE current_time,
    DynamicGridConfig* config,
    lis::cuda::acc_nugrid::AssembledSolution* assembled_solution,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* v_grid,
    int timestep_count,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Calculate refinement criterion for a cell
 * @param config Dynamic grid configuration
 * @param h_grid Water depth grid
 * @param v_grid Velocity grid
 * @param topo_grid Topography grid
 * @param i X-index
 * @param j Y-index
 * @param grid_cols_padded Number of padded columns in grid
 * @return Refinement criterion value
 */
NUMERIC_TYPE calculateRefinementCriterion(
    DynamicGridConfig* config,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* v_grid,
    NUMERIC_TYPE* topo_grid,
    int i,
    int j,
    int grid_cols_padded
);

/**
 * Mark cells for refinement
 * @param config Dynamic grid configuration
 * @param h_grid Water depth grid
 * @param v_grid Velocity grid
 * @param topo_grid Topography grid
 * @param refinement_flags Refinement flags grid
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 */
void markCellsForRefinement(
    DynamicGridConfig* config,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* v_grid,
    NUMERIC_TYPE* topo_grid,
    int* refinement_flags,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Apply buffer zone around refined cells
 * @param config Dynamic grid configuration
 * @param refinement_flags Refinement flags grid
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 */
void applyBufferZone(
    DynamicGridConfig* config,
    int* refinement_flags,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Rebuild grid hierarchy based on refinement flags
 * @param config Dynamic grid configuration
 * @param assembled_solution Assembled solution
 * @param refinement_flags Refinement flags grid
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 */
void rebuildGridHierarchy(
    DynamicGridConfig* config,
    lis::cuda::acc_nugrid::AssembledSolution* assembled_solution,
    int* refinement_flags,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Map data from old grid to new grid
 * @param assembled_solution New assembled solution
 * @param old_assembled_solution Old assembled solution
 * @param h_grid Water depth grid
 * @param v_grid Velocity grid
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 */
void mapDataToNewGrid(
    lis::cuda::acc_nugrid::AssembledSolution* assembled_solution,
    lis::cuda::acc_nugrid::AssembledSolution* old_assembled_solution,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* v_grid,
    int grid_cols,
    int grid_rows,
    int grid_cols_padded
);

/**
 * Update grid neighbors and interfaces
 * @param assembled_solution Assembled solution
 * @param neighbors Non-uniform neighbors
 * @param interfaces Non-uniform interfaces
 */
void updateNeighborsAndInterfaces(
    lis::cuda::acc_nugrid::AssembledSolution* assembled_solution,
    lis::cuda::acc_nugrid::NonUniformNeighbours* neighbors,
    lis::cuda::acc_nugrid::NonUniformInterfaces* interfaces
);

} // namespace dynamic_grid

#endif // DYNAMIC_GRID_H
