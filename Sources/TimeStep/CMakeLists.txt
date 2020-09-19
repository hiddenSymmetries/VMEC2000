# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/eqsolve.f
	${CMAKE_CURRENT_SOURCE_DIR}/evolve.f
	${CMAKE_CURRENT_SOURCE_DIR}/gmres_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/parallel_include_module.f90
	${CMAKE_CURRENT_SOURCE_DIR}/parallel_vmec_module.f90
	${CMAKE_CURRENT_SOURCE_DIR}/interp.f
	${CMAKE_CURRENT_SOURCE_DIR}/restart.f
	${CMAKE_CURRENT_SOURCE_DIR}/runvmec.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_history.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
