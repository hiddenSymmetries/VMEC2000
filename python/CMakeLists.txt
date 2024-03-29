#cmake_minimum_required(VERSION 3.17)

#project(cylinder VERSION 1.0 LANGUAGES C CXX Fortran)
#list(APPEND CMAKE_MODULE_PATH "/home/mbkumar/anaconda3/envs/testv/lib/python3.8/site-packages/skbuild/resources/cmake")
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#list(APPEND CMAKE_MODULE_PATH /u/bmedasan/.conda/envs/stlopt_v/lib/python3.8/site-packages/skbuild/resources/cmake)
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
find_package(F90Wrap REQUIRED)
#if(NOT APPLE)
#	set(CMAKE_INSTALL_RPATH $ORIGIN)
#endif()

# Fortran preprocessing compiler 
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
	set(FPP_COMPILER fpp)
    set(FPP_COMPILE_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
	set(FPP_COMPILER gfortran)
    set(FPP_COMPILE_FLAGS "-E -cpp")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
	set(FPP_COMPILER ftn)
    set(FPP_COMPILE_FLAGS "-e Z")
endif()

get_directory_property(COMP_DEFS COMPILE_DEFINITIONS)
message(STATUS "Compile definitions for preprocessor are ${COMP_DEFS}")
string(REPLACE ";" " " COMP_DEF_STR "${COMPILE_DEFINITIONS}")
string(REPLACE ";" " " SCALAPACK_LIB_STR "${SCALAPACK_LIB}")

	
function(preprocess_fortran outvar)
    message(STATUS "preprocess_fortran arguments: ${outvar}, followed by ${ARGN}")
    set(srcs)
    foreach(f ${ARGN})
        # is it a Fortran file?
        if(f MATCHES "\\.[Ff](9[05])?")
          	message(STATUS "Got fortran file: ${f}")
          	# construct output filename
          	if(NOT IS_ABSOLUTE "${f}")
            	get_filename_component(f "${f}" ABSOLUTE)
       		endif()
        	file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
        	get_filename_component(e "${r}" EXT)
        	get_filename_component(n "${r}" NAME_WE)
        	get_filename_component(p "${r}" PATH)
        	set(of "${CMAKE_CURRENT_BINARY_DIR}/${n}_fpp${e}")
        	message(STATUS "Output name: ${of}")
        	# preprocess the thing
        	if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
				add_custom_command(OUTPUT "${of}"
					#COMMAND ${FPP_EXE} ${incflags} "${defflags}" "${f}" "${of}"
					#COMMAND ${CMAKE_Fortran_COMPILER} -E -cpp "${f}" -o "${of}"
					COMMAND ${FPP_COMPILER} ${FPP_COMPILE_FLAGS} ${COMP_DEF_STR} "${f}" "${of}"
					IMPLICIT_DEPENDS Fortran "${f}"
					COMMENT "Preprocessing ${f}"
					VERBATIM
				)
        	else()
				add_custom_command(OUTPUT "${of}"
					COMMAND ${FPP_COMPILER} ${FPP_COMPILE_FLAGS} ${COMP_DEF_STR} "${f}" -o "${of}"
					IMPLICIT_DEPENDS Fortran "${f}"
					COMMENT "Preprocessing ${f}"
					VERBATIM
				)
        	endif()
        	list(APPEND srcs "${of}")
      		#else()
       		#  list(APPEND srcs "${f}")
     	endif()
    endforeach()
    # return the (preprocessed) sources
    set(${outvar} "${srcs}" PARENT_SCOPE)
endfunction()

set(fortran_src_files ${CMAKE_CURRENT_SOURCE_DIR}/fortran_src/runvmec.f
                      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_src/reinit.f90
                      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_src/cleanup.f90
                      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_src/vmec_input.f
                      ${CMAKE_CURRENT_SOURCE_DIR}/fortran_src/read_wout_mod.f90
)
message(STATUS "fortran_src_files is ${fortran_src_files}")
preprocess_fortran(fpp_files ${fortran_src_files})
message(STATUS "fpp_files is ${fpp_files}")

#add_library(cyl_lib SHARED ${fpp_files})

#function(f90wrap_outputs outvar)
#    message(STATUS "process f90 files with f90wrap")
#    set(srcs)
#    foreach(f ${ARGN})
#        # is it a Fortran file?
#        if(f MATCHES "\\.[Ff](9[05])?")
#            message(STATUS "Got fortran file: ${f}")
#            # construct output filename
#            if(NOT IS_ABSOLUTE "${f}")
#                get_filename_component(f "${f}" ABSOLUTE)
#            endif()
#            file(RELATIVE_PATH r "${CMAKE_CURRENT_BINARY_DIR}" "${f}")
#            get_filename_component(e "${r}" EXT)
#            get_filename_component(n "${r}" NAME_WE)
#            get_filename_component(p "${r}" PATH)
#            set(of "${CMAKE_CURRENT_BINARY_DIR}/f90wrap_${n}${e}")
#            message(STATUS "Output name: ${of}")
       
#            list(APPEND srcs "${of}")
#        #else()
#        #  list(APPEND srcs "${f}")
#        endif()
#    endforeach()
#    # return the (preprocesysed) sources
#    set(${outvar} "${srcs}" PARENT_SCOPE)
#endfunction()

#f90wrap_outputs(f90wrap_output_files ${fpp_files})
#f90wrap_outputs(f90wrap_output_files ${fortran_src_files})

# ----------------------------------------------------------------------------
# NOTE: There is no way to identify the f90wrap---.f90 files ahead of running f90wrap
# NOTE: The files produced have no one->one relation with the source files.
# NOTE: So giving the names of f90wrap_---.f90 files manually 
#-----------------------------------------------------------------------------
set(f90wrap_output_files ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90
			             ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_read_wout_mod.f90  
						 ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_vmec_input.f90 
)

set(kind_map_file ${CMAKE_CURRENT_SOURCE_DIR}/kind_map)
set(python_mod_name vmec)
set(python_mod_file ${CMAKE_CURRENT_BINARY_DIR}/${python_mod_name}.py)

add_custom_command(OUTPUT ${python_mod_file} ${f90wrap_output_files}
    COMMAND "${F90Wrap_EXECUTABLE}" -m "${python_mod_name}" ${fortran_src_files} -k "${kind_map_file}"
    #IMPLICIT_DEPENDS Fortran ${fpp_files}
    DEPENDS ${fortran_src_files} ${kind_map_file}
    COMMENT "Executing F90Wrap for" ${fortran_src_files}
    VERBATIM
)
#execute_process(
#	COMMAND ${F90Wrap_EXECUTABLE} -m ${python_mod_name} ${fortran_src_files} -k ${kind_map_file}
#	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#)
#FILE(GLOB f90wrap_output_files ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_*.f*)

add_custom_target("${python_mod_name}_pymod"
  DEPENDS ${python_mod_file} ${f90wrap_output_files} 
)

#file(GLOB f90wrap_output_files CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/f90wrap_*.f90")
set(f2py_module_name "_${python_mod_name}")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
message(STATUS "Python exten suffix expansion: ${PYTHON_EXTENSION_MODULE_SUFFIX}")
message(STATUS "f90_wrap_output_files: "  ${f90wrap_output_files})
message(STATUS "f2py_module_name: ${f2py_module_name}")
message(STATUS "generated_module_name: ${generated_module_file}")

#set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)
#set(CMAKE_SKIP_BUILD_PATH FALSE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
include_directories("${NumPy_INCLUDE_DIRS}" "${F2PY_INCLUDE_DIRS}" "${CMAKE_CURRENT_BINARY_DIR}")
add_custom_target(${f2py_module_name} ALL
  #DEPENDS ${generated_module_file} ${python_mod_name} cyl_lib ${f90wrap_output_files}
  DEPENDS ${generated_module_file} vmec ${f90wrap_output_files}
)

add_custom_command(
    OUTPUT ${generated_module_file}
    COMMAND ${F2PY_EXECUTABLE} 
        -m ${f2py_module_name} 
        --build-dir ${CMAKE_CURRENT_BINARY_DIR} 
        --f90exec=${CMAKE_Fortran_COMPILER}
        --f77exec=${CMAKE_Fortran_COMPILER}
        -c 
        #${SCALAPACK_LIB} ${NETCDF_F} ${NETCDF_C}
        ${f90wrap_output_files} 
        -I${CMAKE_BINARY_DIR}/build/modules/vmec
		--verbose
        ${CMAKE_BINARY_DIR}/build/lib/libvmec.a  
        ${VMEC_LINK_LIB}
    #IMPLICIT_DEPENDS Fortran ${f90wrap_output_files}
    DEPENDS vmec ${f90wrap_output_files}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    #VERBATIM
    COMMAND_EXPAND_LISTS
)

python_extension_module(${generated_module_file})

install(FILES ${python_mod_file} ${generated_module_file} #${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
        DESTINATION python/vmec
)
install(TARGETS xvmec vmec)
#        LIBRARY DESTINATION ${CMAKE_INSTALL_DIR}/.
#        RUNTIME DESTINATION bin
#)


