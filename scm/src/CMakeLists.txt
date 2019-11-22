CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(scm)
set(PROJECT "CCPP-SCM")

####################################################################
# Begin CCPP prebuild step                                         #
####################################################################
if(STATIC)
  # Start with empty list of suites
  message (STATUS "Generating list of suites to compile for static CCPP")
  set(SUITES "")
  # Get list of all suite definition files (with full path)
  file(GLOB SUITE_DEFINITION_FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/../../ccpp/suites/suite_*.xml"
  )
  # Extract file name and suite name and append to SUITES
  foreach(suite_definition_filepath IN LISTS SUITE_DEFINITION_FILES)
    get_filename_component(suite_definition_filename ${suite_definition_filepath} NAME)
    string(REGEX REPLACE "^suite_(.+)\\.xml$" "\\1" suite_name ${suite_definition_filename})
    set(SUITES ${SUITES}${suite_name},)
    message (STATUS "  adding suite ${suite_name}")
  endforeach(suite_definition_filepath IN LISTS SUITE_DEFINITION_FILES)
  # Abort if no suites found
  if ("${SUITES}" STREQUAL "")
    message(FATAL_ERROR "No suites found to compile for static CCPP")
  endif("${SUITES}" STREQUAL "")
  # Remove trailing comma from list of suites
  string(REGEX REPLACE "(.+),$" "\\1" SUITES ${SUITES})
  # Run CCPP prebuild.py
  message (STATUS "Running ccpp_prebuild.py for static CCPP")
  execute_process(
    COMMAND ccpp/framework/scripts/ccpp_prebuild.py --config=ccpp/config/ccpp_prebuild_config.py --static --suites=${SUITES} --builddir=${CMAKE_CURRENT_BINARY_DIR}/..
    OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
    ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
    RESULT_VARIABLE return_code
  )
else(STATIC)
  # Run CCPP prebuild.py
  message (STATUS "Running ccpp_prebuild.py for dynamic CCPP")
  execute_process(
    COMMAND ccpp/framework/scripts/ccpp_prebuild.py --config=ccpp/config/ccpp_prebuild_config.py --builddir=${CMAKE_CURRENT_BINARY_DIR}/..
    OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
    ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
    RESULT_VARIABLE return_code
  )
endif(STATIC)
# Check return code from CCPP prebuild.py
if(return_code EQUAL 0)
  message (STATUS "CCPP prebuild step completed successfully")
else(return_code EQUAL 0)
  message (FATAL_ERROR "CCPP prebuild step failed, check ccpp_prebuild.out/ccpp_prebuild.err")
endif(return_code EQUAL 0)
####################################################################
# End CCPP prebuild step                                           #
####################################################################

ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

IF(DEFINED ENV{NETCDF})
  MESSAGE(STATUS "The NETCDF environment variable is $ENV{NETCDF}")
  SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})
ELSE(DEFINED ENV{NETCDF})
  MESSAGE(FATAL_ERROR "The NETCDF environement variable must be set to point to your NetCDF installation before building. Stopping...")
ENDIF(DEFINED ENV{NETCDF})

IF(DEFINED ENV{NCEPLIBS_DIR})
  MESSAGE(STATUS "The NCEPLIBS_DIR environment variable is $ENV{NCEPLIBS_DIR}")
  set(NCEPLIBS_DIR $ENV{NCEPLIBS_DIR})
ELSE(DEFINED ENV{NCEPLIBS_DIR})
  MESSAGE(FATAL_ERROR "The NCEPLIBS_DIR environment variable must be set to point to your NCEPlibs installation before building. Stopping...")
ENDIF(DEFINED ENV{NCEPLIBS_DIR})

set(BACIO_LIB4   ${NCEPLIBS_DIR}/lib/libbacio_4.a)
set(SP_LIBd      ${NCEPLIBS_DIR}/lib/libsp_v2.0.2_d.a)
set(W3NCO_LIBd   ${NCEPLIBS_DIR}/lib/libw3nco_d.a)

if(STATIC)
  ADD_DEFINITIONS(-DSTATIC)
endif(STATIC)

SET(CCPP_SRC ${CMAKE_SOURCE_DIR}/../../ccpp/framework)
SET(GFSPHYSICS_SRC ${CMAKE_SOURCE_DIR}/../../ccpp/physics)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/physics)

#------------------------------------------------------------------------------
# Add required preprocessor flags for build type
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  ADD_DEFINITIONS(-DDEBUG)
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

#------------------------------------------------------------------------------
# Add -DCCPP preprocessor flag (needed to preprocess GFS_typedefs.F90 from FV3)
ADD_DEFINITIONS(-DCCPP)

#------------------------------------------------------------------------------
# Add host-model specific preprocessor flag (needed for some physics schemes)
ADD_DEFINITIONS(-DSCM)

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
#ADD_DEFINITIONS(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
ADD_DEFINITIONS(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM -DGFSL64)  # add GFSL64 to use correct GFS ak and bks
#ADD_DEFINITIONS(-DINTERNAL_FILE_NML)

#------------------------------------------------------------------------------
# Detect OS and set preprocessor flags
if(APPLE)
  ADD_DEFINITIONS(-DMACOSX)
elseif(UNIX)
  ADD_DEFINITIONS(-DLINUX)
else (APPLE)
  message (FATAL_ERROR "Unsupported platform, only Linux and MacOSX are supported at this time.")
endif()

#------------------------------------------------------------------------------
# Set compile options
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  # Basic settings
  #if CFLAGS, CXX_FLAGS, or FFLAGS are defined via environment variables, use them
  if (DEFINED ENV{CFLAGS})
    set (CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -O0 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{CFLAGS})
    set (CMAKE_C_FLAGS       "-O0 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{CFLAGS})
  if (DEFINED ENV{CXXFLAGS})
    set (CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} -O0 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{CXXFLAGS})
    set (CMAKE_CXX_FLAGS     "-O0 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{CXXFLAGS})
  if (DEFINED ENV{FFLAGS})
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{FFLAGS})
    set (CMAKE_Fortran_FLAGS "-O0 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{FFLAGS})
  # Compiler-dependent settings
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set (CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -ggdb")
    set (CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} -ggdb")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all -fbacktrace -ffpe-trap=zero -finit-real=nan -ggdb -ffree-line-length-none -cpp -fdefault-real-8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set (CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -traceback")
    set (CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} -traceback")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpe0 -warn -traceback -debug extended -fpp -r8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -traceback")
    set (CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} -traceback")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mnoipa -traceback -Mpreprocess -r8")
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  #if CFLAGS, CXX_FLAGS, or FFLAGS are defined via environment variables, use them
  if (DEFINED ENV{CFLAGS})
    set (CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} -O2 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{CFLAGS})
    set (CMAKE_C_FLAGS       "-O2 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{CFLAGS})
  if (DEFINED ENV{CXXFLAGS})
    set (CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} -O2 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{CXXFLAGS})
    set (CMAKE_CXX_FLAGS     "-O2 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{CXXFLAGS})
  if (DEFINED ENV{FFLAGS})
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -g -fPIC" CACHE STRING "" FORCE)
  else (DEFINED ENV{FFLAGS})
    set (CMAKE_Fortran_FLAGS "-O2 -g -fPIC" CACHE STRING "" FORCE)
  endif (DEFINED ENV{FFLAGS})
  # Compiler-dependent settings
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -cpp -fdefault-real-8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ftz -debug minimal -fpp -r8")
  elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mipa=noconst -Mpreprocess -r8")
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

#------------------------------------------------------------------------------
# Set netCDF flags for preprocessor, compiler and linker (if defined)
set(NETCDF_DIR $ENV{NETCDF})
if(NETCDF_DIR)
  set (NETCDF_INC "-I${NETCDF_DIR}/include")
  set (NETCDF_LIB "-L${NETCDF_DIR}/lib -lnetcdff -lnetcdf")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_INC} ${NETCDF_LIB}")
  ADD_DEFINITIONS(-DNETCDF)
  message (STATUS "Enable netCDF support")
else(NETCDF_DIR)
  message (STATUS "Disable netCDF support")
endif(NETCDF_DIR)

#------------------------------------------------------------------------------
# Set SIONlib flags for preprocessor, compiler and linker (if defined)
set(SIONLIB $ENV{SIONLIB})
if (SIONLIB)
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --cflags --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_FLAGS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --cflags --f90' failed")
  endif (NOT exit_code EQUAL "0")
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --libs --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_LIBS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --libs --f90' failed")
  endif (NOT exit_code EQUAL "0")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SIONlib_Fortran_FLAGS} ${SIONlib_Fortran_LIBS}")
  ADD_DEFINITIONS(-DSION)
  message (STATUS "Enable SIONlib support")
else (SIONLIB)
  message (STATUS "Disable SIONlib support")
endif (SIONLIB)

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
# CMAKE_Fortran_COMPILER_ID = {"Intel", "PGI", "GNU", "Clang", "MSVC", ...}
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../ccpp/framework/cmake")
#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
else(OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif()

#------------------------------------------------------------------------------
# Configure sources
ADD_SUBDIRECTORY(${CCPP_SRC} ${CMAKE_BINARY_DIR}/ccpp/framework)
ADD_SUBDIRECTORY(${GFSPHYSICS_SRC} ${CMAKE_BINARY_DIR}/ccpp/physics)
ADD_DEPENDENCIES(ccppphys ccpp)

SET(scm_source_files gmtb_scm.F90
            gmtb_scm_input.f90
            gmtb_scm_utils.f90
            gmtb_scm_vgrid.F90
            gmtb_scm_setup.f90
            gmtb_scm_forcing.f90
            gmtb_scm_time_integration.f90
            gmtb_scm_output.f90
            gmtb_scm_kinds.f90
            gmtb_scm_physical_constants.f90
            gmtb_scm_type_defs.f90
)

if(STATIC)
  ADD_EXECUTABLE(gmtb_scm ${scm_source_files} ccpp_static_api.F90)
  TARGET_LINK_LIBRARIES(gmtb_scm ccppphys ccpp ${BACIO_LIB4} ${SP_LIBd} ${W3NCO_LIBd})
else(STATIC)
  ADD_EXECUTABLE(gmtb_scm ${scm_source_files})
  TARGET_LINK_LIBRARIES(gmtb_scm ccppphys ccpp)
endif(STATIC)

set_target_properties(gmtb_scm PROPERTIES
                               COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                               LINK_FLAGS "${CMAKE_Fortran_FLAGS}")

add_custom_command(
        TARGET gmtb_scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_SOURCE_DIR}/run_gmtb_scm.py
                ${CMAKE_CURRENT_BINARY_DIR}/run_gmtb_scm.py)

add_custom_command(
        TARGET gmtb_scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_SOURCE_DIR}/multi_run_gmtb_scm.py
                ${CMAKE_CURRENT_BINARY_DIR}/multi_run_gmtb_scm.py)
