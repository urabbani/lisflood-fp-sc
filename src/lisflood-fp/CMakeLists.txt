# LISFLOOD-FP CMake Build Configuration
# Provides build instructions for the LISFLOOD-FP 8.2 model

# Version information
set(LISFLOOD_VERSION_MAJOR 8)
set(LISFLOOD_VERSION_MINOR 2)
set(LISFLOOD_VERSION_PATCH 0)
set(LISFLOOD_VERSION "${LISFLOOD_VERSION_MAJOR}.${LISFLOOD_VERSION_MINOR}.${LISFLOOD_VERSION_PATCH}")

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
include(CheckLanguage)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80) # Default to compute capability 8.0 if not specified
endif()

set(CUDAToolkit_ROOT "/usr/local/cuda")
project(lisflood-fp)

# _CONFIG will be overridden by -D_CONFIG=../config/local during docker build
set(_CONFIG config.default.cmake CACHE FILEPATH "CMake config file")
include(${_CONFIG}) # Include config/local here, which sets CMAKE_CXX_FLAGS

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release) # Default to Release build if not specified
endif()


enable_language(CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# --- Base library (collects common sources and dependencies) ---
add_library(lisflood-base INTERFACE)
target_sources(lisflood-base INTERFACE
	src/boundary.cpp             #<-- Added src/
	src/ch_flow.cpp              #<-- Added src/
	src/chkpnt.cpp               #<-- Added src/
	src/fp_acc.cpp               #<-- Added src/
	src/fp_flow.cpp              #<-- Added src/
	src/fp_trent.cpp             #<-- Added src/
	src/infevap.cpp              #<-- Added src/
 	src/cuda/solver_consistency.cpp   #<-- Added src/
 	src/dynamic_grid.cpp         #<-- Added src/
 	src/dam_operations.cpp       #<-- Added src/
 	src/hazard_calc.cpp          #<-- Added src/
	src/input.cpp                #<-- Added src/
	src/iterateq.cpp             #<-- Added src/
	src/output.cpp               #<-- Added src/
	src/pars.cpp                 #<-- Added src/
	src/por_flow.cpp             #<-- Added src/
	src/sgc.cpp                  #<-- Added src/
	src/util.cpp                 #<-- Added src/
	src/utility.cpp              #<-- Added src/
	src/weir_flow.cpp            #<-- Added src/
	src/lisflood2/DataTypes.cpp  #<-- Added src/
	src/lisflood2/file_tool.cpp  #<-- Added src/
	src/lisflood2/lis2_output.cpp #<-- Added src/
	src/lisflood2/lisflood_processing.cpp #<-- Added src/
	src/lisflood2/sgm_fast.cpp   #<-- Added src/
	src/lisflood2/windows_ext.cpp #<-- Added src/
	src/rain/rain.cpp            #<-- Added src/
	src/swe/boundary.cpp         #<-- Added src/
	src/swe/dg2.cpp              #<-- Added src/
	src/swe/fields.cpp           #<-- Added src/
	src/swe/flux.cpp             #<-- Added src/
	src/swe/fv1.cpp              #<-- Added src/
	src/swe/hll.cpp              #<-- Added src/
	src/swe/input.cpp            #<-- Added src/
	src/swe/output.cpp           #<-- Added src/
	src/swe/stats.cpp            #<-- Added src/
	src/swe/fv1/modifiedvars.cpp #<-- Added src/
	src/swe/dg2/dg2_output.cpp   #<-- Added src/
	src/swe/dg2/fields.cpp       #<-- Added src/
	src/swe/dg2/friction.cpp     #<-- Added src/
	src/swe/dg2/modifiedvars.cpp #<-- Added src/
	src/swe/dg2/slope_limiter.cpp #<-- Added src/
)

# --- Main Executable ---
add_executable(lisflood)
target_sources(lisflood PRIVATE src/lisflood.cpp)
target_link_libraries(lisflood PRIVATE lisflood-base)

# Add linker options specifically for the lisflood executable
# These were previously set globally via LDFLAGS in config/local
target_link_options(lisflood PRIVATE -O3 -static)


# --- CUDA Support (Conditional) ---
set(CMAKE_CUDA_COMPILER_WORKS TRUE)
enable_language(CUDA)
if (CMAKE_CUDA_COMPILER)
	find_package(CUDA REQUIRED)
	include_directories(${CUDA_INCLUDE_DIRS})
	message(STATUS "CUDA Compiler Found: ${CMAKE_CUDA_COMPILER}")
	message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
	set_property(TARGET lisflood PROPERTY CUDA_SEPARABLE_COMPILATION ON) # Note: Repeated line removed below
	# set_property(TARGET lisflood PROPERTY CUDA_SEPARABLE_COMPILATION ON) # Removed duplicate

	add_compile_definitions(CUDA USE_CUDA_FOR_ALL_SOLVERS)
		
	# Add extra flags to make CUDA compilation more permissive
	target_compile_options(lisflood PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--Wno-deprecated-gpu-targets -x cu>)

	target_sources(lisflood PRIVATE
		src/cuda/ghostraster.cpp
		src/cuda/io.cpp
		src/cuda/sample.cpp
		src/cuda/stats.cpp
		src/cuda/cpu_solver_wrapper.cu
		src/cuda_solver_selector.cpp
	)

	set_property(SOURCE
		src/cuda/cuda_boundary.cu
		src/cuda/cuda_dem.cu
		src/cuda/cuda_flow.cu
		src/cuda/cuda_geometry.cu
		src/cuda/cuda_hll.cu
		src/cuda/cuda_max_field.cu
		src/cuda/cuda_rain.cu
		src/cuda/cuda_sample.cu
		src/cuda/cuda_simulate.cu
		src/cuda/cuda_snapshot.cu
		src/cuda/cuda_solver.cu
		src/cuda/cuda_stats.cu
		src/cuda/cuda_unifiedallocator.cu
		src/cuda/cuda_util.cu
		src/cuda/fv1/cuda_fv1_flow.cu
		src/cuda/fv1/cuda_fv1_simulate.cu
		src/cuda/fv1/cuda_fv1_snapshot.cu
		src/cuda/fv1/cuda_fv1_solver.cu
		src/cuda/fv2/cuda_fv2_dem.cu
		src/cuda/fv2/cuda_fv2_flow.cu
		src/cuda/fv2/cuda_fv2_simulate.cu
		src/cuda/fv2/cuda_fv2_slope_limit.cu
		src/cuda/fv2/cuda_fv2_snapshot.cu
		src/cuda/fv2/cuda_fv2_solver.cu
		src/cuda/dg2/cuda_dg2_dem.cu
		src/cuda/dg2/cuda_dg2_flow.cu
		src/cuda/dg2/cuda_dg2_simulate.cu
		src/cuda/dg2/cuda_dg2_slope_limit.cu
		src/cuda/dg2/cuda_dg2_snapshot.cu
		src/cuda/dg2/cuda_dg2_solver.cu
		src/cuda/acc/cuda_acc_flow.cu
		src/cuda/acc/cuda_acc_simulate.cu
		src/cuda/acc/cuda_acc_snapshot.cu
		src/cuda/acc/cuda_acc_solver.cu
		src/cuda/acc/cuda_acc_max.cu
		src/cuda/acc_nugrid/cuda_acc_nugrid_simulate.cu
		src/cuda/acc_nugrid/get_num_blocks.cpp
		src/cuda/acc_nugrid/cuda_utils.cu
		src/cuda/acc_nugrid/read_gauge_points.cu
		src/cuda/acc_nugrid/read_point_srcs.cu
		src/cuda/acc_nugrid/read_bdy_conds.cu
		src/cuda/acc_nugrid/gen_bdy_morton_codes.cu
		src/cuda/acc_nugrid/preflag_details.cu
		src/cuda/acc_nugrid/get_modal_values.cu
		src/cuda/acc_nugrid/read_raster_file.cpp
		src/cuda/acc_nugrid/generate_all_morton_codes.cu
		src/cuda/acc_nugrid/sort_finest_scale_coefficients_z_order.cu
		src/cuda/acc_nugrid/copy_finest_coefficients.cu
		src/cuda/acc_nugrid/get_max_scale_coefficients.cu
		src/cuda/acc_nugrid/preflag_topo.cu
		src/cuda/acc_nugrid/encode_and_thresh_topo.cu
		src/cuda/acc_nugrid/get_reg_tree.cu
		src/cuda/acc_nugrid/reg_new.cu
		src/cuda/acc_nugrid/project_assem_sol.cu
		src/cuda/acc_nugrid/reinsert_assem_sol.cu
		src/cuda/acc_nugrid/traverse_tree_of_sig_details.cu
		src/cuda/acc_nugrid/reverse_z_order_assembled_solution.cu
		src/cuda/acc_nugrid/write_reals_to_file.cu
		src/cuda/acc_nugrid/write_raster_maps.cu
		src/cuda/acc_nugrid/write_soln_vtk.cu
		src/cuda/acc_nugrid/reverse_z_order_act_idcs.cu
		src/cuda/acc_nugrid/get_compaction_flags.cu
		src/cuda/acc_nugrid/compaction.cu
		src/cuda/acc_nugrid/count_neighbours.cu
		src/cuda/acc_nugrid/find_nonuniform_neighbours.cu
		src/cuda/acc_nugrid/count_interfaces_per_neighbours.cu
		src/cuda/acc_nugrid/sort_nghbr_assem_sol_row_major.cu
		src/cuda/acc_nugrid/init_neighbours.cu
		src/cuda/acc_nugrid/init_interfaces.cu
		src/cuda/acc_nugrid/find_interfaces.cu
		src/cuda/acc_nugrid/load_interface_q_vol.cu
		src/cuda/acc_nugrid/load_neighbour_h.cu
		src/cuda/acc_nugrid/compute_q.cu
		src/cuda/acc_nugrid/update_h.cu
		src/cuda/acc_nugrid/calculate_dt.cu
		src/cuda/acc_nugrid/get_dt_CFL.cu
		src/cuda/acc_nugrid/init_h.cu
		src/cuda/acc_nugrid/init_q.cu
		src/cuda/acc_nugrid/read_all_bdy_conds.cu
		src/cuda/acc_nugrid/read_num_stage_points.cpp
		src/cuda/acc_nugrid/read_num_bdy_cells.cpp
		src/cuda/acc_nugrid/read_num_point_srcs.cpp
		src/cuda/acc_nugrid/write_gauge_point_data.cu
		src/cuda/acc_nugrid/unifiedallocator.cu
		src/cuda/acc_nugrid/rain.cu
		src/cuda/acc_nugrid/write_mass_data.cu
		src/cuda/acc_nugrid/rain_uniform.cu
		src/cuda/acc_nugrid/update_max.cu
		src/cuda/acc_nugrid/write_max_maps.cu
		src/cuda/acc_nugrid/init_max.cu
		src/cuda/acc_nugrid/write_chkpnt.cu
		src/cuda/acc_nugrid/read_chkpnt.cu
		PROPERTY LANGUAGE CUDA)


	add_subdirectory(src/cuda/adaptive)

	# Enable line info generation for CUDA code (useful for debugging/profiling)
	target_compile_options(lisflood PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--generate-line-info>)

	target_include_directories(lisflood-base INTERFACE
		/usr/local/cuda/include
		${CMAKE_CURRENT_SOURCE_DIR}/src/cuda
		${CMAKE_CURRENT_SOURCE_DIR}/src/cuda/common
	)
	target_link_libraries(lisflood-base INTERFACE cuda_runtime)

else()
	message(STATUS "CUDA Compiler Not Found. Building CPU-only version.")
endif()


# --- Helper Executables ---
add_executable(generateDG2DEM src/preprocess/generateDG2DEM.cpp)
target_link_libraries(generateDG2DEM PRIVATE lisflood-base)
# Note: Link options (-O3, -static) are NOT added to helpers unless required.

add_executable(generateDG2start src/preprocess/generateDG2start.cpp)
target_link_libraries(generateDG2start PRIVATE lisflood-base)
# Note: Link options (-O3, -static) are NOT added to helpers unless required.

add_executable(DG2downscale src/postprocess/DG2downscale.cpp)
target_link_libraries(DG2downscale PRIVATE lisflood-base)
# Note: Link options (-O3, -static) are NOT added to helpers unless required.


# --- Dependencies applied to lisflood-base (inherited by linked targets) ---
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
	message(STATUS "OpenMP Found. Enabling OpenMP support.")
	target_link_libraries(lisflood-base INTERFACE OpenMP::OpenMP_CXX) # Links -lgomp or equivalent
else()
	message(WARNING "OpenMP Not Found. Building without OpenMP parallelization.")
endif()

if (UNIX AND NOT APPLE) # NUMA is typically Linux-specific
	find_package(NUMA REQUIRED)
	message(STATUS "NUMA Found. Linking against ${NUMA_LIBRARIES}")
	target_include_directories(lisflood-base INTERFACE ${NUMA_INCLUDE_DIRS})
	target_link_libraries(lisflood-base INTERFACE ${NUMA_LIBRARIES}) # Links -lnuma
endif()

# Check if NetCDF support is requested via the _NETCDF definition
# This definition comes from the CMAKE_CXX_FLAGS set in config/local
# We need to check the flag string itself, not a CMake variable _NETCDF
string(FIND "${CMAKE_CXX_FLAGS}" "-D_NETCDF=1" _NETCDF_ENABLED)
if (NOT _NETCDF_ENABLED EQUAL -1)
	message(STATUS "NetCDF support requested (_NETCDF=1 found in CXX_FLAGS).")
	if (WIN32)
		# Support both x86 and x64 architectures on Windows
		message(STATUS "Configuring NetCDF for Windows.")
		if(CMAKE_SIZEOF_VOID_P EQUAL 8)
			# 64-bit architecture
			message(STATUS "Using 64-bit NetCDF libraries.")
			target_include_directories(lisflood-base INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/windep/netCDF4-64/include)
			target_link_directories(lisflood-base INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/windep/netCDF4-64/lib)
		else()
			# 32-bit architecture
			message(STATUS "Using 32-bit NetCDF libraries.")
			target_include_directories(lisflood-base INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/windep/netCDF4-32/include)
			target_link_directories(lisflood-base INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/windep/netCDF4-32/lib)
		endif()
		target_link_libraries(lisflood-base INTERFACE netcdf)
		# Copy DLLs - Consider if this is the best approach; might be better to ensure they are in PATH
		# file(GLOB DLLS "${CMAKE_CURRENT_SOURCE_DIR}/DLL's/*.dll") # Use CMAKE_CURRENT_SOURCE_DIR
		# if(DLLS)
		# 	message(STATUS "Copying NetCDF DLLs to build directory.")
		# 	file(COPY ${DLLS} DESTINATION ${CMAKE_BINARY_DIR})
		# else()
		# 	message(WARNING "NetCDF DLLs not found in ${CMAKE_CURRENT_SOURCE_DIR}/DLL's/")
		# endif()
	else() # Non-Windows (Linux/macOS)
		message(STATUS "Searching for NetCDF using find_package.")
		find_package(NetCDF REQUIRED)
		message(STATUS "NetCDF Found: Includes=${NetCDF_INCLUDE_DIRS}, Libraries=${NetCDF_LIBRARIES}")
		target_include_directories(lisflood-base INTERFACE ${NetCDF_INCLUDE_DIRS})
		target_link_libraries(lisflood-base INTERFACE NetCDF::NetCDF) # Links -lnetcdf
	endif()
	# Add the compile definition if NetCDF is enabled (already done via CXX_FLAGS, but belts and suspenders)
	target_compile_definitions(lisflood-base INTERFACE _NETCDF=1)
else()
	message(STATUS "NetCDF support not requested (_NETCDF=1 not found in CXX_FLAGS).")
	# Ensure the definition is explicitly 0 if not enabled
	target_compile_definitions(lisflood-base INTERFACE _NETCDF=0)
endif()

# Debug definition (already handled by CMAKE_BUILD_TYPE implicitly in many cases, but can be explicit)
if (CMAKE_BUILD_TYPE MATCHES Debug)
	target_compile_definitions(lisflood-base INTERFACE _DEBUG)
endif()

# --- Installation (Optional - uncomment and adjust if needed) ---
# install(TARGETS lisflood generateDG2DEM generateDG2start DG2downscale
#	 RUNTIME DESTINATION bin
#	 LIBRARY DESTINATION lib
#	 ARCHIVE DESTINATION lib/archive)
# install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/config.default.txt DESTINATION share/lisflood-fp RENAME lisflood_config.txt)
