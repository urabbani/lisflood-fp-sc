#include "find_interfaces.cuh"
#include <algorithm>

__global__ void lis::cuda::acc_nugrid::find_interfaces
(
	NonUniformNeighbours d_non_uniform_nghbrs,
	AssembledSolution d_assem_sol,
	NonUniformInterfaces d_non_uniform_itrfaces
)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= d_non_uniform_nghbrs.length) return;

	int cumu_interface_count = d_non_uniform_nghbrs.cumu_itrface_counts[idx];

	int nghbr_dir = d_non_uniform_nghbrs.dirs[idx]; // MKS

	NUMERIC_TYPE dx_sum = C(0.0);
	NUMERIC_TYPE dx     = C(0.0);

	int count = 0;

	// FINDING OWNER INTERFACES //

	int owner_elem_idx = d_non_uniform_nghbrs.owner_elem_idcs[idx];

	//
	//int level = d_assem_sol.levels[d_non_uniform_nghbrs.owner_elem_idcs[idx]];
	//MortonCode code = d_assem_sol.act_idcs[d_non_uniform_nghbrs.owner_elem_idcs[idx]] - get_lvl_idx(level);
	//int x = compact(code);
	//int y = compact(code >> 1);
	//

	for (int i = 0; i < d_assem_sol.nghbr_counts[owner_elem_idx]; i++)
	{
		int owner_itrfaces_idx = d_assem_sol.cumu_nghbr_counts[owner_elem_idx] + i;

		int interface_dir = d_non_uniform_nghbrs.dirs[owner_itrfaces_idx];

		bool vertical =
			(
				(
					nghbr_dir == I_EAST
					||
					nghbr_dir == I_WEST
					)
				&&
				(
					interface_dir == I_NORTH
					||
					interface_dir == I_SOUTH
					)
				);

		bool horizontal =
			(
				(
					nghbr_dir == I_NORTH
					||
					nghbr_dir == I_SOUTH
					)
				&&
				(
					interface_dir == I_EAST
					||
					interface_dir == I_WEST
					)
				);

		// ------------------------ //

		index_1D nghbr_act_idx = d_non_uniform_nghbrs.nghbr_act_idcs[idx];

		// Handle ghost interfaces properly
		if (nghbr_act_idx < 0)
		{
			// Mark this as a ghost cell boundary
			d_non_uniform_nghbrs.nghbr_elem_idcs[idx] = -1;
			
			// For boundary conditions, handle based on boundary type
			int boundary_type = -nghbr_act_idx - 1; // Convert negative index to boundary type
			
			if (boundary_type >= 0 && boundary_type < 4) // Valid boundary types
			{
				// Setup for different boundary types (0=closed, 1=open, 2=specified Q, 3=specified H)
				switch (boundary_type)
				{
					case 0: // Closed boundary
						// Zero flow across boundary
						d_non_uniform_itrfaces.dx[cumu_interface_count] = C(0.0);
						d_non_uniform_itrfaces.q_vol[cumu_interface_count] = C(0.0);
						break;
					
					case 1: // Open boundary (free flow)
						// Setup for normal flow calculation
						d_non_uniform_itrfaces.dx[cumu_interface_count] = d_non_uniform_nghbrs.dx[owner_itrfaces_idx];
						d_non_uniform_itrfaces.q_vol[cumu_interface_count] = C(0.0);
						break;
					
					case 2: // Specified Q
					case 3: // Specified H
						// Setup for specified boundary conditions
						// These will be handled in flow calculation
						d_non_uniform_itrfaces.dx[cumu_interface_count] = d_non_uniform_nghbrs.dx[owner_itrfaces_idx];
						d_non_uniform_itrfaces.q_vol[cumu_interface_count] = C(0.0);
						break;
					
					default:
						// Default to closed boundary
						d_non_uniform_itrfaces.dx[cumu_interface_count] = C(0.0);
						d_non_uniform_itrfaces.q_vol[cumu_interface_count] = C(0.0);
				}
				
				// Update interface information
				d_non_uniform_itrfaces.load_idcs[cumu_interface_count] = owner_itrfaces_idx;
				cumu_interface_count++;
				count++;
			}
			
			return;
		}
		// ----------------------------------- //


		if (vertical || horizontal)
		{
			dx = d_non_uniform_nghbrs.dx[owner_itrfaces_idx];

			dx_sum += dx;

			d_non_uniform_itrfaces.dx[cumu_interface_count] = dx;
			d_non_uniform_itrfaces.load_idcs[cumu_interface_count] = owner_itrfaces_idx;
			cumu_interface_count++;
			count++;
		}
	}

	// FINDING NEIGHBOUR INTERFACES //
	int nghbr_elem_idx = d_non_uniform_nghbrs.nghbr_elem_idcs[idx];

	for (int i = 0; i < d_assem_sol.nghbr_counts[nghbr_elem_idx]; i++)
	{
		int nghbr_itrfaces_idx = d_assem_sol.cumu_nghbr_counts[nghbr_elem_idx] + i;

		int interface_dir = d_non_uniform_nghbrs.dirs[nghbr_itrfaces_idx];

		bool vertical =
			(
				(
					nghbr_dir == I_EAST
					||
					nghbr_dir == I_WEST
					)
				&&
				(
					interface_dir == I_NORTH
					||
					interface_dir == I_SOUTH
					)
				);

		bool horizontal =
			(
				(
					nghbr_dir == I_NORTH
					||
					nghbr_dir == I_SOUTH
					)
				&&
				(
					interface_dir == I_EAST
					||
					interface_dir == I_WEST
					)
				);

		if (vertical || horizontal)
		{
			dx = d_non_uniform_nghbrs.dx[nghbr_itrfaces_idx];

			dx_sum += dx;

			d_non_uniform_itrfaces.dx[cumu_interface_count] = dx;
			d_non_uniform_itrfaces.load_idcs[cumu_interface_count] = nghbr_itrfaces_idx;
			cumu_interface_count++;
			count++;
		}
	}

	// ---------------------------- //

	d_non_uniform_nghbrs.dx_sum[idx] = dx_sum;

	int nghbr = idx;

	int interface_count = d_non_uniform_nghbrs.itrface_counts[idx];

	idx = d_non_uniform_nghbrs.cumu_itrface_counts[idx];

	for (int i = 0; i < interface_count; i++)
	{
		d_non_uniform_itrfaces.q_vol[idx + i] = C(0.0);
		d_non_uniform_itrfaces.nghbrs[idx + i] = nghbr;
		d_non_uniform_itrfaces.nghbr_counts[idx + i] = count;
	}
}