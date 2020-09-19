# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/add_tension.f
	${CMAKE_CURRENT_SOURCE_DIR}/getspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/gety2.f
	${CMAKE_CURRENT_SOURCE_DIR}/initspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/jacprod.f
	${CMAKE_CURRENT_SOURCE_DIR}/set_dual.f
	${CMAKE_CURRENT_SOURCE_DIR}/setspline.f
	${CMAKE_CURRENT_SOURCE_DIR}/setup_int.f
	${CMAKE_CURRENT_SOURCE_DIR}/sort_data.f
	${CMAKE_CURRENT_SOURCE_DIR}/splinint.f
	${CMAKE_CURRENT_SOURCE_DIR}/splints.f
	${CMAKE_CURRENT_SOURCE_DIR}/splintx.f
	${CMAKE_CURRENT_SOURCE_DIR}/vspline.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
