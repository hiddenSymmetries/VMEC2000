	# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bss.f
	${CMAKE_CURRENT_SOURCE_DIR}/dump_output.f
	${CMAKE_CURRENT_SOURCE_DIR}/elongation.f
	${CMAKE_CURRENT_SOURCE_DIR}/eqfor.f
	${CMAKE_CURRENT_SOURCE_DIR}/fileout.f
	${CMAKE_CURRENT_SOURCE_DIR}/freeb_data.f
	${CMAKE_CURRENT_SOURCE_DIR}/fsym_fft.f
	${CMAKE_CURRENT_SOURCE_DIR}/fsym_invfft.f
	${CMAKE_CURRENT_SOURCE_DIR}/getbrho.f
	${CMAKE_CURRENT_SOURCE_DIR}/getcurmid.f
	${CMAKE_CURRENT_SOURCE_DIR}/jxbforce.f
	${CMAKE_CURRENT_SOURCE_DIR}/mercier.f
	#${CMAKE_CURRENT_SOURCE_DIR}/perf2likwid.f90
	${CMAKE_CURRENT_SOURCE_DIR}/printout.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_indata.f
	${CMAKE_CURRENT_SOURCE_DIR}/readin.f
	${CMAKE_CURRENT_SOURCE_DIR}/TestWout.f
	${CMAKE_CURRENT_SOURCE_DIR}/timer_sub.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_io.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmercier.f
	${CMAKE_CURRENT_SOURCE_DIR}/write_dcon.f
	${CMAKE_CURRENT_SOURCE_DIR}/wrout.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
