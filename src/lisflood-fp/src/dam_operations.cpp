/*
 * dam_operations.cpp
 *
 * Implementation of dam operations module for LISFLOOD-FP
 * Based on original implementation in sgm_fast.cpp
 */

#include "dam_operations.h"
#include <math.h>
#include <time.h>

namespace dam_ops {

void calculateDamRelease(
    const NUMERIC_TYPE delta_time, 
    const NUMERIC_TYPE curr_time, 
    const NUMERIC_TYPE* h_grid, 
    NUMERIC_TYPE* Qx_grid, 
    NUMERIC_TYPE* Qy_grid, 
    DamData* Damptr, 
    const NUMERIC_TYPE g, 
    const int dam_id)
{
    // Get operation code
    int op_code = Damptr->DamOperationCode[dam_id];
    
    // Initialize with default value
    Damptr->DamOperationQ[dam_id] = C(0.0);
    
    // Calculate release based on operation code
    switch (op_code) {
        case CONSTANT_RELEASE: 
            // Constant release (original code = 1)
            Damptr->DamOperationQ[dam_id] = Damptr->DamMeanQ[dam_id];
            
            // Check available volume
            if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
            }
            break;
            
        case STORAGE_BASED_RELEASE:
            // Storage-based release (original code = 2)
            {
                NUMERIC_TYPE upper_limit = getmin((Damptr->Volmax[dam_id] / delta_time), Damptr->DamMeanQ[dam_id]);
                NUMERIC_TYPE lower_limit = getmax((Damptr->Volmax[dam_id] - Damptr->DamVol[dam_id]), 0) / delta_time;
                Damptr->DamOperationQ[dam_id] = getmax(lower_limit, upper_limit);
                
                // Check available volume
                if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                    Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
                }
            }
            break;
            
        case LINEAR_STORAGE_RELEASE:
            // Linear storage-based release (original code = 3)
            {
                NUMERIC_TYPE tmp = Damptr->DamMeanQ[dam_id] * ((Damptr->DamVol[dam_id] / Damptr->Volmax[dam_id]) + C(0.5));
                NUMERIC_TYPE upper_limit = getmin((Damptr->Volmax[dam_id] / delta_time), tmp);
                NUMERIC_TYPE lower_limit = getmax((Damptr->DamVol[dam_id] - Damptr->Volmax[dam_id]), 0) / delta_time;
                Damptr->DamOperationQ[dam_id] = getmax(lower_limit, upper_limit);
                
                // Check available volume
                if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                    Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
                }
            }
            break;
            
        case WATERGAP_MODEL:
            // WaterGAP model (original code = 4)
            // Q= kS[S/Smax]^1.5, where k=0.01/d
            Damptr->DamOperationQ[dam_id] = (C(0.01) * Damptr->DamVol[dam_id] * 
                pow((Damptr->DamVol[dam_id] / Damptr->Volmax[dam_id]), C(1.5))) / delta_time;
                
            // Check available volume
            if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
            }
            break;
            
        case PCR_GLOBWB_MODEL:
            // PCR-GLOBWB model (original code = 5)
            // Q=(S-Smin)/(Smax-Smin)*Qmean; Smin=10% of Smax
            Damptr->DamOperationQ[dam_id] = (((Damptr->DamVol[dam_id] - (Damptr->Volmax[dam_id] * C(0.2))) / 
                (Damptr->Volmax[dam_id] * C(0.7)))) * Damptr->DamMeanQ[dam_id];
                
            // Check available volume
            if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
            }
            break;
            
        case WBMPLUS_MODEL:
            // WBMplus model (original code = 6)
            {
                NUMERIC_TYPE Kappa = C(0.16);
                NUMERIC_TYPE Lambda = C(0.6);
                
                Damptr->DamOperationQ[dam_id] = Kappa * (Damptr->DamVin[dam_id] / delta_time);
                
                if ((Damptr->DamVin[dam_id] / delta_time) < Damptr->DamMeanQ[dam_id]) {
                    Damptr->DamOperationQ[dam_id] = Lambda * (Damptr->DamVin[dam_id] / delta_time) + 
                        (Damptr->DamMeanQ[dam_id] - (Damptr->DamVin[dam_id] / delta_time));
                }
                
                // Check available volume
                if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                    Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
                }
            }
            break;
            
        case HANASAKI_MODEL:
            // Hanasaki et al. (2005) model (original code = 7)
            {
                NUMERIC_TYPE Alpha = C(0.85);
                NUMERIC_TYPE Storage = Damptr->Volmax[dam_id] / Damptr->DamMeanQ[dam_id];
                
                // Update yearly calculation
                if (curr_time <= Damptr->DamYear[dam_id]) {
                    Damptr->OP7_Kappa[dam_id] = Damptr->DamVol[dam_id] / (Alpha * Damptr->Volmax[dam_id]);
                    Damptr->DamYear[dam_id] += (C(365.0) * C(86400.0)); // One year in seconds
                }
                
                // Daily release
                NUMERIC_TYPE DayRelease = Damptr->DamMeanQ[dam_id];
                
                // Calculate release
                Damptr->DamOperationQ[dam_id] = Damptr->OP7_Kappa[dam_id] * DayRelease;
                
                if (Storage < C(0.5)) {
                    Damptr->DamOperationQ[dam_id] = pow((Storage / 0.5), 2) * Damptr->OP7_Kappa[dam_id] * 
                        DayRelease + (1 - pow((Storage / 0.5), 2)) * (Damptr->DamVin[dam_id] / delta_time);
                }
                
                // Check available volume
                if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                    Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
                }
            }
            break;
            
        case CUSTOM_RULE_CURVE:
            // Custom rule curve - implemented in calculateRuleCurveRelease
            // Implementation would be added here when rule curves are available
            break;
            
        case FORECAST_BASED:
            // Forecast-based operation - placeholder for future implementation
            // Would consider forecast inflows for optimized operation
            break;
            
        case CASCADE_DAM:
            // Cascade dam operation - placeholder for future implementation
            // Would consider upstream/downstream dam operations
            break;
            
        default:
            // Default to constant release if unknown operation code
            Damptr->DamOperationQ[dam_id] = Damptr->DamMeanQ[dam_id];
            
            // Check available volume
            if (Damptr->DamOperationQ[dam_id] * delta_time > Damptr->DamVol[dam_id]) {
                Damptr->DamOperationQ[dam_id] = Damptr->DamVol[dam_id] / delta_time;
            }
            break;
    }
}

void initializeRuleCurve(
    DamData* Damptr,
    RuleCurve* rule_curves,
    int dam_id)
{
    // Default initialization sets all months to mean flow
    for (int month = 0; month < 12; month++) {
        rule_curves[dam_id].monthly_targets[month] = Damptr->DamMeanQ[dam_id];
    }
    
    // Optionally could read from file in future implementation
}

NUMERIC_TYPE calculateRuleCurveRelease(
    const NUMERIC_TYPE delta_time,
    const NUMERIC_TYPE curr_time,
    DamData* Damptr,
    RuleCurve* rule_curves,
    int dam_id)
{
    // Get current month (0-11) from simulation time
    time_t sim_time = (time_t)curr_time;
    struct tm* timeinfo = gmtime(&sim_time);
    int month = timeinfo->tm_mon;
    
    // Get target release for current month
    NUMERIC_TYPE target_release = rule_curves[dam_id].monthly_targets[month];
    
    // Adjust release based on storage
    NUMERIC_TYPE storage_ratio = Damptr->DamVol[dam_id] / Damptr->Volmax[dam_id];
    
    // Adjustment factor based on storage
    NUMERIC_TYPE adjustment_factor;
    
    if (storage_ratio > 0.9) {
        adjustment_factor = 1.2; // Release more when storage is high
    } else if (storage_ratio < 0.3) {
        adjustment_factor = 0.8; // Release less when storage is low
    } else {
        adjustment_factor = 1.0; // Normal release
    }
    
    NUMERIC_TYPE release = target_release * adjustment_factor;
    
    // Check available volume
    if (release * delta_time > Damptr->DamVol[dam_id]) {
        release = Damptr->DamVol[dam_id] / delta_time;
    }
    
    return release;
}

NUMERIC_TYPE calculateSpillwayFlow(
    DamData* Damptr,
    const NUMERIC_TYPE g,
    int dam_id)
{
    // Calculate spillway flow if water level exceeds spillway crest
    if (Damptr->InitialHeight[dam_id] <= Damptr->SpillHeight[dam_id]) {
        return C(0.0); // No spillway flow
    } else {
        // Standard weir equation for spillway flow
        return Damptr->Spill_Cd[dam_id] * Damptr->SpillWidth[dam_id] * 
            pow(g, C(0.5)) * pow((Damptr->InitialHeight[dam_id] - Damptr->SpillHeight[dam_id]), C(1.5));
    }
}

void updateDamState(
    const NUMERIC_TYPE delta_time,
    DamData* Damptr,
    int dam_id)
{
    // Update dam volume
    Damptr->DamLoss = Damptr->DamLoss + Damptr->DamVol[dam_id];
    Damptr->DamVol[dam_id] += Damptr->DamVin[dam_id];
    Damptr->DamVol[dam_id] -= (Damptr->DamTotalQ[dam_id] * delta_time);
    Damptr->DamLoss = Damptr->DamLoss - Damptr->DamVol[dam_id];
    
    // Update dam height
    Damptr->InitialHeight[dam_id] = (Damptr->DamVol[dam_id] / Damptr->Volmax[dam_id]) * 
        Damptr->DamHeight[dam_id] + (Damptr->SpillHeight[dam_id] - Damptr->DamHeight[dam_id]);
}

void processDams(
    const NUMERIC_TYPE delta_time,
    const NUMERIC_TYPE curr_time,
    const int grid_cols_padded,
    const NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* volume_grid,
    NUMERIC_TYPE* Qx_grid,
    NUMERIC_TYPE* Qy_grid,
    DamData* Damptr,
    const NUMERIC_TYPE g)
{
    // Process all dams
    for (int dam_id = 0; dam_id < Damptr->NumDams; dam_id++) {
        // Initialize dam flows
        Damptr->DamOperationQ[dam_id] = C(0.0);
        Damptr->DamTotalQ[dam_id] = C(0.0);
        
        // Calculate operational release
        calculateDamRelease(delta_time, curr_time, h_grid, Qx_grid, Qy_grid, Damptr, g, dam_id);
        Damptr->DamTotalQ[dam_id] += Damptr->DamOperationQ[dam_id];
        
        // Calculate spillway flow
        Damptr->SpillQ[dam_id] = calculateSpillwayFlow(Damptr, g, dam_id);
        Damptr->DamTotalQ[dam_id] += Damptr->SpillQ[dam_id];
        
        // Update dam state (volume and height)
        updateDamState(delta_time, Damptr, dam_id);
        
        // Output dam release to model grid
        int output_cell = Damptr->OutputCellX[dam_id] + Damptr->OutputCellY[dam_id] * grid_cols_padded;
        volume_grid[output_cell] += (Damptr->DamTotalQ[dam_id] * delta_time);
    }
}

} // namespace dam_ops