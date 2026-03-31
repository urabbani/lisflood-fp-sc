/*
 * hazard_calc.cpp
 *
 * Implementation of hazard calculation functions for LISFLOOD-FP
 */

#include "hazard_calc.h"
#include <math.h>

namespace hazard {

NUMERIC_TYPE calculateDefraDvHazard(NUMERIC_TYPE depth, NUMERIC_TYPE velocity) {
    if (depth < C(0.001)) return C(0.0);
    
    // HR = d * (v + 0.5) equation from Defra/EA
    NUMERIC_TYPE hazard = depth * (velocity + C(0.5));
    
    return hazard;
}

int calculateFoudiHazardClass(NUMERIC_TYPE depth, NUMERIC_TYPE velocity) {
    if (depth < C(0.001)) return 0;
    
    // Foudi et al. (2015) hazard classification
    if (depth < C(0.4) && velocity < C(1.0)) return 1;
    if (depth < C(0.8) && velocity < C(1.0)) return 2;
    if (depth < C(1.2) && velocity < C(1.0)) return 3;
    if (depth >= C(1.2) || velocity >= C(1.0)) return 4;
    
    return 0; // Default case
}

NUMERIC_TYPE calculateDebrisFactor(NUMERIC_TYPE depth) {
    // Based on Australian Rainfall and Runoff Guidelines
    if (depth < C(0.5)) return C(0.0);
    if (depth < C(1.2)) return C(0.5);
    return C(1.0);
}

NUMERIC_TYPE calculateARRHazard(NUMERIC_TYPE depth, NUMERIC_TYPE velocity) {
    if (depth < C(0.001)) return C(0.0);
    
    // HR = d * (v + 0.5) + DF
    NUMERIC_TYPE debris_factor = calculateDebrisFactor(depth);
    NUMERIC_TYPE hazard = depth * (velocity + C(0.5)) + debris_factor;
    
    return hazard;
}

void updateMaxHazard(
    int grid_cols, int grid_rows, int grid_cols_padded,
    NUMERIC_TYPE depth_thresh,
    NUMERIC_TYPE* h_grid,
    NUMERIC_TYPE* Vx_grid, NUMERIC_TYPE* Vy_grid,
    NUMERIC_TYPE* maxHazard_grid) {
    
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < grid_rows; j++) {
        for (int i = 0; i < grid_cols; i++) {
            int index = i + j * grid_cols_padded;
            
            // Only calculate for wet cells
            if (h_grid[index] > depth_thresh) {
                // Get velocity components
                int Vx_index1 = i + j * grid_cols_padded;
                int Vx_index2 = (i + 1) + j * grid_cols_padded;
                int Vy_index1 = i + j * grid_cols_padded;
                int Vy_index2 = i + (j + 1) * grid_cols_padded;
                
                // Calculate average velocity for cell
                NUMERIC_TYPE vx = C(0.5) * (Vx_grid[Vx_index1] + Vx_grid[Vx_index2]);
                NUMERIC_TYPE vy = C(0.5) * (Vy_grid[Vy_index1] + Vy_grid[Vy_index2]);
                
                // Calculate velocity magnitude
                NUMERIC_TYPE velocity = SQRT(vx * vx + vy * vy);
                
                // Calculate hazard rating using Defra method
                NUMERIC_TYPE hazard = calculateARRHazard(h_grid[index], velocity);
                
                // Update max hazard if current hazard is higher
                if (hazard > maxHazard_grid[index]) {
                    maxHazard_grid[index] = hazard;
                }
            }
        }
    }
}

} // namespace hazard