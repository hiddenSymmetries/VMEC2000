# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/precon2d.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
