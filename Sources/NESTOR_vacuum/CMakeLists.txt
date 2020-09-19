# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/analysum.f
	${CMAKE_CURRENT_SOURCE_DIR}/analysum2.f
	${CMAKE_CURRENT_SOURCE_DIR}/analyt.f
	${CMAKE_CURRENT_SOURCE_DIR}/becoil.f
	${CMAKE_CURRENT_SOURCE_DIR}/belicu.f
	${CMAKE_CURRENT_SOURCE_DIR}/bextern.f
	${CMAKE_CURRENT_SOURCE_DIR}/fouri.f
	${CMAKE_CURRENT_SOURCE_DIR}/fourp.f
	${CMAKE_CURRENT_SOURCE_DIR}/greenf.f
	${CMAKE_CURRENT_SOURCE_DIR}/precal.f
	${CMAKE_CURRENT_SOURCE_DIR}/scalpot.f
	${CMAKE_CURRENT_SOURCE_DIR}/surface.f
	${CMAKE_CURRENT_SOURCE_DIR}/tolicu.f
	${CMAKE_CURRENT_SOURCE_DIR}/vac_persistent.f
	${CMAKE_CURRENT_SOURCE_DIR}/vacmod.f
	${CMAKE_CURRENT_SOURCE_DIR}/vacmod0.f
	${CMAKE_CURRENT_SOURCE_DIR}/vacuum.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
