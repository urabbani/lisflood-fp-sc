/*
 * hazard_calc.h
 *
 * Hazard calculation functions for LISFLOOD-FP
 * This module implements standardized hazard calculation metrics
 *
 */

#ifndef HAZARD_CALC_H
#define HAZARD_CALC_H

#include "lisflood.h"

namespace hazard {

/**
 * Calculate hazard rating based on velocity and depth
 * Uses the HR = d * (v + 0.5) equation from Defra/EA
 * @param depth Water depth in meters
 * @param velocity Velocity in m/s
 * @return Hazard rating value
 */
NUMERIC_TYPE calculateDefraDvHazard(NUMERIC_TYPE depth, NUMERIC_TYPE velocity);

/**
 * Calculate hazard rating based on velocity and depth
 * Based on the method in Foudi et al. (2015)
 * @param depth Water depth in meters
 * @param velocity Velocity in m/s
 * @return Hazard class (1-4) with 4 being most severe
 */
int calculateFoudiHazardClass(NUMERIC_TYPE depth, NUMERIC_TYPE velocity);

/**
 * Calculate debris factor based on water depth
 * @param depth Water depth in meters
 * @return Debris factor
 */
NUMERIC_TYPE calculateDebrisFactor(NUMERIC_TYPE depth);

/**
 * Calculate full hazard rating with debris factor
 * Based on the Australian Rainfall and Runoff Guideline (2019)
 * HR = d * (v + 0.5) + DF
 * @param depth Water depth in meters
 * @param velocity Velocity in m/s
 * @return Hazard rating with debris factor
 */
NUMERIC_TYPE calculateARRHazard(NUMERIC_TYPE depth, NUMERIC_TYPE velocity);

/**
 * Update maximum hazard in grid cells
 * @param grid_cols Number of columns in grid
 * @param grid_rows Number of rows in grid
 * @param grid_cols_padded Number of padded columns in grid
 * @param depth_thresh Depth threshold for hazard calculation
 * @param h_grid Water depth grid
 * @param Vx_grid X-velocity grid
 * @param Vy_grid Y-velocity grid
 * @param maxHazard_grid Maximum hazard grid
 */
void updateMaxHazard(
    int grid_cols, int grid_rows, int grid_cols_padded,
    NUMERIC_TYPE depth_thresh,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* Vx_grid, NUMERIC_TYPE* Vy_grid,
    NUMERIC_TYPE* maxHazard_grid);

} // namespace hazard

#endif // HAZARD_CALC_H