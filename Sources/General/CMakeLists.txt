# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/add_fluxes.f90
	${CMAKE_CURRENT_SOURCE_DIR}/reinit.f90
	${CMAKE_CURRENT_SOURCE_DIR}/alias.f
	${CMAKE_CURRENT_SOURCE_DIR}/angle_constraints.f90
	${CMAKE_CURRENT_SOURCE_DIR}/aspectratio.f
	${CMAKE_CURRENT_SOURCE_DIR}/bcovar.f
	${CMAKE_CURRENT_SOURCE_DIR}/blocktridiagonalsolver_bst.f90
	${CMAKE_CURRENT_SOURCE_DIR}/blocktridiagonalsolver.f90
	${CMAKE_CURRENT_SOURCE_DIR}/convert.f
	${CMAKE_CURRENT_SOURCE_DIR}/csplinx.f
	${CMAKE_CURRENT_SOURCE_DIR}/directaccess.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fbal.f
	${CMAKE_CURRENT_SOURCE_DIR}/forces.f
	${CMAKE_CURRENT_SOURCE_DIR}/funct3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/getfsq.f
	${CMAKE_CURRENT_SOURCE_DIR}/jacobian.f
	${CMAKE_CURRENT_SOURCE_DIR}/lamcal.f90
	${CMAKE_CURRENT_SOURCE_DIR}/precondn.f
	${CMAKE_CURRENT_SOURCE_DIR}/realspace.f
	${CMAKE_CURRENT_SOURCE_DIR}/residue.f90
	${CMAKE_CURRENT_SOURCE_DIR}/scalfor.f
	${CMAKE_CURRENT_SOURCE_DIR}/spectrum.f
	${CMAKE_CURRENT_SOURCE_DIR}/symforce.f
	${CMAKE_CURRENT_SOURCE_DIR}/symrzl.f
	${CMAKE_CURRENT_SOURCE_DIR}/tomnsp_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/totzsp_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/vforces.f
	#${CMAKE_CURRENT_SOURCE_DIR}/vmec_c_interface.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_dim.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_main.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_params.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_persistent.f
	${CMAKE_CURRENT_SOURCE_DIR}/xstuff.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/fbal.f APPEND PROPERTY OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/vforce.f")
