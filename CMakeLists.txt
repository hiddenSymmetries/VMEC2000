# cmake build system was created by Mark Cianciosa. Current this only supports Mac OS X and Linux.
#
# All options can be define by using -DOPTION_NAME="Definition of option"
#
# Options:
# EXTRA_FLAGS         -- Compiler flags added to all builds.
# EXTRA_RELEASE_FLAGS -- Compiler flags added to release builds.
# EXTRA_DEBUG_FLAGS   -- Compiler flags added to debug builds.
#
# NETCDF_LIB_PATH     -- Path to the netcdf libraries.
# NETCDF_INC_PATH     -- Path to the netcdf include files.
# NETCDF_F_LIB_NAME   -- Name of Netcdf_f library if nonstandard
# NETCDF_C_LIB_NAME   -- Name of netcdf_c library if nonstandard
#
# SCALAPACK_LIB_NAME  -- Name of Sscalapack library if nonstandard
#
# Some builds of SCALAPACK require linking to external BLACS LIBS
# These options allow you to define nonstandard BLACS LIB NAMES if needed
#
# BLACS_LIB_NAME      -- Name of BLACS library
# BLACS_CINIT_NAME    -- Name of BLACS Cinit library
# BLACS_F77INIT_NAME  -- Name of BLACS F77init library
#
# BUILD_GTOVMI        -- Optionally builds GTOVMI

cmake_minimum_required (VERSION 2.8.7)

message(STATUS "CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

# Start Project and note the language used. The BLAS and LAPACK libraries for ifort require a working C and C++ compiler.
project (VMEC Fortran C CXX)

enable_testing ()


# Activate or Deactivate code profiling. The profiler maybe activated by setting -DUSE_PROFILER=ON on the cmake command line.
OPTION (USE_PROFILER "Activate code profiling." OFF)
if (USE_PROFILER)
        add_definitions (-DPROFILE_ON)
endif ()

# Activate or Deactivate fast copies. The fast copies maybe activated by setting -DUSE_FAST_COPY=ON on the cmake command line.
OPTION (USE_FAST_COPY "Activate code fast copies." OFF)
if (USE_FAST_COPY)
	add_definitions (-DFAST_COPY)
endif ()

# Add the ability to auto generate API documentation with Doxygen
#find_package (Doxygen)
#if (DOXYGEN_FOUND)
#	configure_file (${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/build/Doxyfile @ONLY)
#	add_custom_target (doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/build/Doxyfile 
#	                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build
#	                   COMMENT "Generating API documentation with Doxygen" VERBATIM)
#endif ()

# Configure MPI
# If using the compiler wrapper, there is no need to find the MPI libraries.
get_filename_component (compiler ${CMAKE_Fortran_COMPILER} NAME)
if (${compiler} STREQUAL mpiifort OR ${compiler} STREQUAL mpifort)
	add_definitions (-DMPI_OPT)
else ()
	find_package (MPI REQUIRED)
endif ()

if (MPI_Fortran_FOUND)
	add_definitions (-DMPI_OPT)
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
	include_directories (${MPI_Fortran_INCLUDE_PATH})

#  Some of the tests need to be oversubscribbed to run. There is a problem where
#  depending on the current platform or implementation of MPI the oversubscribe
#  flag changes. MPIEXEC may not write to stdout correctly. Add the ability to
#  set this variable manually.
	if (NOT DEFINED ${MPI_OVERSUBSCRIBE_FLAG}) 
		execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE MPI_VERSION)
		if (${MPI_VERSION} MATCHES "Open MPI" OR
		    ${MPI_VERSION} MATCHES "OpenRTE"  OR
		    ${MPI_VERSION} MATCHES "slurm"    OR
		    ${MPI_VERSION} MATCHES "aprun")
			set (MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
	 	endif ()
	endif ()
endif ()

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

# Add extra compiler flags to various builds.
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_DEBUG_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS}")

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor. The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/build/modules)

# Set Apple specific options. This adds a predefine for DARWIN and adds lapack from the Accelerate framework.
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	add_definitions (-DDARWIN)
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
endif ()

# Set Linux specific options. This adds a predefine for LINUX.
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	add_definitions (-DLINUX)
	find_package (BLAS REQUIRED)
	find_package (LAPACK REQUIRED)
endif ()

# Use netcdf if libraries are found. The values of NETCDF_LIB_PATH and NETCDF_INC_PATH should be set from the command line or Cmake GUI. Convert path from native format to cmake internal format.
if (NETCDF_INC_PATH)
	file (TO_NATIVE_PATH ${NETCDF_INC_PATH} NETCDF_INC_PATH)
	find_file (NETCDF_I netcdf.inc ${NETCDF_INC_PATH} NO_DEFAULT_PATH)
endif ()
message(STATUS "NETCDF_I value is ${NETCDF_I}")

if (NETCDF_LIB_PATH)
	file (TO_NATIVE_PATH ${NETCDF_LIB_PATH} NETCDF_LIB_PATH)
        find_library (NETCDF_F NAMES netcdff ${NETCDF_F_LIB_NAME} PATHS ${NETCDF_LIB_PATH} NO_DEFAULT_PATH)
        find_library (NETCDF_C NAMES netcdf ${NETCDF_C_LIB_NAME} PATHS ${NETCDF_LIB_PATH} NO_DEFAULT_PATH)
endif ()
message(STATUS "NETCDF_F value is ${NETCDF_F}")
message(STATUS "NETCDF_C value is ${NETCDF_C}")
if (NETCDF_C_LIB_PATH)
        find_library (NETCDF_C NAMES netcdf ${NETCDF_C_LIB_NAME} PATHS ${NETCDF_C_LIB_PATH} NO_DEFAULT_PATH)
endif ()

if (NETCDF_I AND (NETCDF_F OR NETCDF_C))
	message (STATUS "Using netcdf")
	add_definitions (-DNETCDF)
	include_directories (${NETCDF_INC_PATH})

	if (NETCDF_F)
		list (APPEND NETCDF_LIBRARIES ${NETCDF_F})
	endif ()
	if (NETCDF_C)
		list (APPEND NETCDF_LIBRARIES ${NETCDF_C})
	endif ()
else ()
	message (FATAL_ERROR "Failed to find the required netcdf libraries. These must be found by setting the NETCDF_INC_PATH and NETCDF_LIB_PATH variables.")
endif ()

SET(SCALAPACK_PREFIX "/usr/local" CACHE PATH "scalapack is installed under this directory")
SET(BLACS_PREFIX "/usr/local" CACHE PATH "blacs is installed under this directory")
#ENDIF ()
#FIND_PATH(yaml_cpp_INCLUDE_DIRS "yaml-cpp/yaml.h" HINTS "${YAML_CPP_PREFIX}/include")
#    IF ( NOT yaml_cpp_INCLUDE_DIRS )
#        MESSAGE(FATAL_ERROR "Could not find the yaml-cpp headers")
#    ENDIF ()
#FIND_LIBRARY(yaml_cpp_LIBRARIES "yaml-cpp" HINTS "${YAML_CPP_PREFIX}/lib" "${YAML_CPP_PREFIX}/lib" ENV LD_LIBRARY_PATH)
#    IF ( NOT yaml_cpp_LIBRARIES )
#        MESSAGE(FATAL_ERROR "Could not find the yaml-cpp library")
#    ENDIF ()

find_library (SCALAPACK_LIB NAMES scalapack ${SCALAPACK_LIB_NAME} HINTS "${SCALAPACK_PREFIX}/lib" ENV LD_LIBRARY_PATH)
message(STATUS "scalapack_lib is ${SCALAPACK_LIB}")
if (SCALAPACK_LIB)
# Some builds of SCALAPACK require linking to external BLACS Libraries
	find_library (BLACS_LIB NAMES blacs ${BLACS_LIB_NAME} HINTS "${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACS_LIB)
		list (APPEND SCALAPACK_LIB ${BLACS_LIB})
	endif ()
	find_library (BLACSCinit_LIB NAMES blacsCinit ${BLACS_CINIT_NAME} HINTS "${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACSCinit_LIB)
		list (APPEND SCALAPACK_LIB ${BLACSCinit_LIB})
	endif ()
	find_library (BLACSF77init_LIB NAMES blacsF77init ${BLACS_F77INIT_NAME} HINTS "${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACSF77init_LIB)
		list (APPEND SCALAPACK_LIB ${BLACSF77init_LIB})
	endif ()
else ()
	if (USE_PARVMEC)
		message (FATAL_ERROR "PARVMEC requires a SCALAPACK library. A suitable library was not found.")
	endif ()
endif ()
message(STATUS "updated scalapack_lib is ${SCALAPACK_LIB}")

add_definitions (-DMRC)

set (vmec_sources "")
# List all the subdirectories these are in library dependancy order.
add_subdirectory (Sources)

add_library (vmec STATIC ${vmec_sources})
#add_dependencies (vmec stell)
set_target_properties (vmec PROPERTIES COMPILE_DEFINITIONS "SKS")
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (vmec PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
else ()
	set_target_properties (vmec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
endif ()
	

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/Sources/TimeStep/vmec.f)
set_target_properties (xvmec PROPERTIES COMPILE_DEFINITIONS "SKS")
#add_dependencies (xvmec vmec)
target_link_libraries (xvmec vmec)

if (CMAKE_VERSION VERSION_GREATER "2.8.8")
	set_target_properties (xvmec PROPERTIES INCLUDE_DIRECTORIES "${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
	#target_include_directories (xvmec PRIVATE "${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
else ()
	set_target_properties (xvmec PROPERTIES COMPILE_FLAGS "-I${CMAKE_Fortran_MODULE_DIRECTORY}/vmec")
endif ()

target_link_libraries (xvmec ${NETCDF_LIBRARIES})

if (SCALAPACK_LIB)
	target_link_libraries (xvmec ${SCALAPACK_LIB})
endif ()
if (BLAS_FOUND)
	target_link_libraries (xvmec ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (xvmec ${LAPACK_LIBRARIES})
endif ()
if (MPI_Fortran_FOUND)
        target_link_libraries (xvmec ${MPI_Fortran_LIBRARIES})
endif ()

