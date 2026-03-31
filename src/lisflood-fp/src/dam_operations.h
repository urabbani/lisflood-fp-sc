/*
 * dam_operations.h
 *
 * Header file for dam operations module in LISFLOOD-FP
 * This module provides flexible dam operation capabilities
 * 
 */

#ifndef DAM_OPERATIONS_H
#define DAM_OPERATIONS_H

#include "lisflood.h"

namespace dam_ops {

// Dam operation types enum for easier identification
enum DamOperationType {
    CONSTANT_RELEASE = 1,         // Basic constant release
    STORAGE_BASED_RELEASE = 2,    // Release based on storage
    LINEAR_STORAGE_RELEASE = 3,   // Linear relationship with storage
    WATERGAP_MODEL = 4,           // Based on Doll et al. (2003)
    PCR_GLOBWB_MODEL = 5,         // Based on Wada et al. (2014)
    WBMPLUS_MODEL = 6,            // Based on Wisser et al. (2010)
    HANASAKI_MODEL = 7,           // Based on Hanasaki et al. (2005)
    CUSTOM_RULE_CURVE = 8,        // User-defined rule curve
    FORECAST_BASED = 9,           // Forecast-based operation
    CASCADE_DAM = 10              // Part of a cascade of dams
};

/**
 * Operation rule curve for dams
 * Defines monthly release targets
 */
struct RuleCurve {
    NUMERIC_TYPE monthly_targets[12];  // Target releases for each month
};

/**
 * Calculate dam operational release
 * This is the original implementation from sgm_fast.cpp
 * @param delta_time Time step (s)
 * @param curr_time Current simulation time (s)
 * @param h_grid Water depth grid
 * @param Qx_grid X-direction discharge grid
 * @param Qy_grid Y-direction discharge grid
 * @param Damptr Pointer to dam data structure
 * @param g Gravitational acceleration
 * @param dam_id Dam index
 */
void calculateDamRelease(
    const NUMERIC_TYPE delta_time, 
    const NUMERIC_TYPE curr_time, 
    const NUMERIC_TYPE* h_grid, 
    NUMERIC_TYPE* Qx_grid, 
    NUMERIC_TYPE* Qy_grid, 
    DamData* Damptr, 
    const NUMERIC_TYPE g, 
    const int dam_id
);

/**
 * Initialize rule curves for dams
 * @param Damptr Pointer to dam data structure
 * @param rule_curves Array of rule curves
 * @param dam_id Dam index
 */
void initializeRuleCurve(
    DamData* Damptr,
    RuleCurve* rule_curves,
    int dam_id
);

/**
 * Calculate dam operational release using rule curves
 * @param delta_time Time step (s)
 * @param curr_time Current simulation time (s)
 * @param Damptr Pointer to dam data structure
 * @param rule_curves Array of rule curves
 * @param dam_id Dam index
 * @return Operational release (m³/s)
 */
NUMERIC_TYPE calculateRuleCurveRelease(
    const NUMERIC_TYPE delta_time,
    const NUMERIC_TYPE curr_time,
    DamData* Damptr,
    RuleCurve* rule_curves,
    int dam_id
);

/**
 * Calculate spillway flow
 * @param Damptr Pointer to dam data structure
 * @param g Gravitational acceleration
 * @param dam_id Dam index
 * @return Spillway flow (m³/s)
 */
NUMERIC_TYPE calculateSpillwayFlow(
    DamData* Damptr,
    const NUMERIC_TYPE g,
    int dam_id
);

/**
 * Update dam state (volume, height) after operations
 * @param delta_time Time step (s)
 * @param Damptr Pointer to dam data structure
 * @param dam_id Dam index
 */
void updateDamState(
    const NUMERIC_TYPE delta_time,
    DamData* Damptr,
    int dam_id
);

/**
 * Process all dams in the simulation
 * @param delta_time Time step (s)
 * @param curr_time Current simulation time (s)
 * @param grid_cols Number of columns in grid
 * @param grid_cols_padded Number of padded columns in grid
 * @param h_grid Water depth grid
 * @param volume_grid Volume grid
 * @param Qx_grid X-direction discharge grid
 * @param Qy_grid Y-direction discharge grid
 * @param Damptr Pointer to dam data structure
 * @param g Gravitational acceleration
 */
void processDams(
    const NUMERIC_TYPE delta_time,
    const NUMERIC_TYPE curr_time,
    const int grid_cols_padded,
    const NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* volume_grid,
    NUMERIC_TYPE* Qx_grid,
    NUMERIC_TYPE* Qy_grid,
    DamData* Damptr,
    const NUMERIC_TYPE g
);

} // namespace dam_ops

#endif // DAM_OPERATIONS_H