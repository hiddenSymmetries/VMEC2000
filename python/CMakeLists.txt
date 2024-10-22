set(Python_FIND_STRATEGY LOCATION)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_program(F2PY_EXECUTABLE NAMES f2py${PYTHON_VERSION_MAJOR} f2py)
find_package(F90Wrap REQUIRED)
find_package(PythonExtensions REQUIRED)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Python headers are at ${PYTHON_INCLUDE_DIRS}")
message(STATUS "NumPy headers are at ${NumPy_INCLUDE_DIRS}")
message(STATUS "F2PY headers are at ${F2PY_INCLUDE_DIRS}")

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

add_custom_target("${python_mod_name}_pymod"
    DEPENDS ${python_mod_file} ${f90wrap_output_files} 
)

set(f2py_module_name "_${python_mod_name}")
set(generated_module_file ${CMAKE_BINARY_DIR}/build/lib/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
message(STATUS "Python exten suffix expansion: ${PYTHON_EXTENSION_MODULE_SUFFIX}")
message(STATUS "f90_wrap_output_files: "  ${f90wrap_output_files})
message(STATUS "f2py_module_name: ${f2py_module_name}")
message(STATUS "generated_module_name: ${generated_module_file}")

add_custom_target(${f2py_module_name_ct} ALL
  DEPENDS ${generated_module_file} vmec ${f90wrap_output_files}
)

set(f2py_module_c "${f2py_module_name}module.c")
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # "${CMAKE_CURRENT_BINARY_DIR}/${f2py_f90_wrapper}"
  COMMAND ${F2PY_EXECUTABLE}
          -m "${f2py_module_name}"
          --build-dir ${CMAKE_CURRENT_BINARY_DIR}
          --lower # Important
          "${f90wrap_output_files}"
  DEPENDS ${fortran_src_files} ${f2py_module_name_ct} # Fortran source
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND_EXPAND_LISTS
)
add_library(${f2py_module_name} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIRS}/fortranobject.c"
            ${f90wrap_output_files}
            # ${f2py_f90_wrapper}
)
target_include_directories(${f2py_module_name} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${NumPy_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS}
                           ${CMAKE_BINARY_DIR}/build/modules/vmec)
target_link_libraries(${f2py_module_name} vmec)
set_target_properties(${f2py_module_name} PROPERTIES SUFFIX "${PYTHON_EXTENSION_MODULE_SUFFIX}")
set_target_properties(${f2py_module_name} PROPERTIES PREFIX "")

# Linker fixes
if (UNIX)
  if (APPLE)
    set_target_properties(${f2py_module_name} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${f2py_module_name} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

install(FILES ${python_mod_file} # ${generated_module_file} #${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
       DESTINATION python/vmec
)
install(TARGETS ${f2py_module_name} DESTINATION python/vmec)

