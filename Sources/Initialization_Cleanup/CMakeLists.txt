# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_funct3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_ns.f
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_nunv.f
	${CMAKE_CURRENT_SOURCE_DIR}/close_all_files.f
	${CMAKE_CURRENT_SOURCE_DIR}/fixaray.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_funct3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_ns.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_nunv.f
	${CMAKE_CURRENT_SOURCE_DIR}/free_persistent_mem.f
	${CMAKE_CURRENT_SOURCE_DIR}/guess_axis.f
	${CMAKE_CURRENT_SOURCE_DIR}/heading.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_geometry.f90
	${CMAKE_CURRENT_SOURCE_DIR}/initialize_radial.f
	${CMAKE_CURRENT_SOURCE_DIR}/initialize_vmec_arrays.f
	${CMAKE_CURRENT_SOURCE_DIR}/load_xc_from_wout.f
	${CMAKE_CURRENT_SOURCE_DIR}/magnetic_fluxes.f
	${CMAKE_CURRENT_SOURCE_DIR}/open_output_files.f
	${CMAKE_CURRENT_SOURCE_DIR}/profil1d.f
	${CMAKE_CURRENT_SOURCE_DIR}/profil3d.f
	#${CMAKE_CURRENT_SOURCE_DIR}/profile_functions.f
	${CMAKE_CURRENT_SOURCE_DIR}/reset_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/vsetup.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
