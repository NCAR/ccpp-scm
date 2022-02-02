cmake_minimum_required(VERSION 3.14)

project(scm
        VERSION 5.0.0
        LANGUAGES C CXX Fortran)
set(PROJECT "CCPP-SCM")

####################################################################
# Begin CCPP prebuild step                                         #
####################################################################
# Start with empty list of suites
message (STATUS "Generating list of suites to compile for CCPP")
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
  message(FATAL_ERROR "No suites found to compile for CCPP")
endif("${SUITES}" STREQUAL "")
# Remove trailing comma from list of suites
string(REGEX REPLACE "(.+),$" "\\1" SUITES ${SUITES})
# Run CCPP prebuild.py
message (STATUS "Running ccpp_prebuild.py for CCPP")
# Make the directories where the ccpp_prebuild.py script wants to write caps and make snippets
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics/physics)
execute_process(
  COMMAND ccpp/framework/scripts/ccpp_prebuild.py --config=ccpp/config/ccpp_prebuild_config.py --suites=${SUITES} --builddir=${CMAKE_CURRENT_BINARY_DIR}
  OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
  ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
  RESULT_VARIABLE return_code
)
  
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

# First, need to figure out if using the new NCEPLIBS-external/NCEPLIBS packages
# or the legacy/contrib netCDF and NCEPLIBS installations
IF(DEFINED ENV{bacio_ROOT})
  message(STATUS "Configuring SCM build for new NCEPLIBS-external/NCEPLIBS software packages")
  set(NEW_NCEPLIBS ON)
ELSE()
  message(STATUS "Configuring SCM build for legacy/contrib netCDF/NCEPLIBS software packages")
  set(NEW_NCEPLIBS OFF)
ENDIF()

IF(NEW_NCEPLIBS)
  MESSAGE(STATUS "Using FindNetCDF.cmake module for $ENV{CMAKE_Platform}")
  FIND_PACKAGE(NetCDF REQUIRED COMPONENTS C Fortran)
else()
  IF(DEFINED ENV{NETCDF})
    MESSAGE(STATUS "The NETCDF environment variable is $ENV{NETCDF}")
    SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})
  ELSE(DEFINED ENV{NETCDF})
    MESSAGE(FATAL_ERROR "The NETCDF environement variable must be set to point to your NetCDF installation before building. Stopping...")
  ENDIF(DEFINED ENV{NETCDF})
ENDIF()

IF(NEW_NCEPLIBS)
  find_package(bacio REQUIRED)
  find_package(sp    REQUIRED)
  find_package(w3nco REQUIRED)
ELSE()
  # Find bacio
  IF(DEFINED ENV{BACIO_LIB4})
    MESSAGE(STATUS "The BACIO_LIB4 environment variable is $ENV{BACIO_LIB4}")
    SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{BACIO_LIB4})
    SET(BACIO_LIB4 $ENV{BACIO_LIB4})
  ELSE()
    MESSAGE(FATAL_ERROR "The BACIO_LIB4 environment variable must be set to point to your BACIO installation (part of NCEPLIBS) before building. Stopping...")
  ENDIF()
  # Find sp
  IF(DEFINED ENV{SP_LIBd})
    MESSAGE(STATUS "The SP_LIBd environment variable is $ENV{SP_LIBd}")
    SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{SP_LIBd})
    SET(SP_LIBd $ENV{SP_LIBd})
  ELSE(DEFINED ENV{SP_LIBd})
    MESSAGE(FATAL_ERROR "The SP_LIBd environment variable must be set to point to your SP installation (part of NCEPLIBS) before building. Stopping...")
  ENDIF(DEFINED ENV{SP_LIBd})
  # Find w3nco
  IF(DEFINED ENV{W3NCO_LIBd})
    MESSAGE(STATUS "The W3NCO_LIBd environment variable is $ENV{W3NCO_LIBd}")
    SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{W3NCO_LIBd})
    SET(W3NCO_LIBd $ENV{W3NCO_LIBd})
  ELSE(DEFINED ENV{W3NCO_LIBd})
    MESSAGE(FATAL_ERROR "The W3NCO_LIBd environment variable must be set to point to your W3NCO installation (part of NCEPLIBS) before building. Stopping...")
  ENDIF(DEFINED ENV{W3NCO_LIBd})
ENDIF()

SET(CCPP_FRAMEWORK_SRC ${CMAKE_SOURCE_DIR}/../../ccpp/framework)
SET(CCPP_PHYSICS_SRC ${CMAKE_SOURCE_DIR}/../../ccpp/physics)

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

# Build Options
set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in SCM)")
set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
set(AVX             OFF CACHE BOOL "Enable AVX-I instruction set")
set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")

#------------------------------------------------------------------------------
# Set compile options
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check")
  
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  
  if(NOT 32BIT)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  endif()
  
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
  else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
  endif()
  
  set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE     "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_BITFORBIT     "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align -sox -fp-model source")
  
  if(NOT 32BIT)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  endif()
  
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv -init=snan,arrays")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ftrapuv")
  else()
    if(32BIT)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3 -no-prec-div -no-prec-sqrt")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -debug minimal")
    set(FAST "-fast-transcendentals")
    if(AVX2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx2")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx2")
    elseif(SIMDMULTIARCH)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axSSE4.2,CORE-AVX2")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -axSSE4.2,CORE-AVX2")
    elseif(AVX)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx-i")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx-i")
    endif()
  endif()
  
  set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE     "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_BITFORBIT     "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "This program has only been compiled with gfortran and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")


#------------------------------------------------------------------------------
# Set flag for 32bit dynamics build
if(32BIT)
  message(STATUS "Compile CCPP slow physics with 64-bit precision, fast physics with 32-bit precision")
  add_definitions(-DOVERLOAD_R4)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-real-size 64 -no-prec-div -no-prec-sqrt")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-fdefault-real-8 -fdefault-double-8")
  endif()
  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
else()
  message(STATUS "Compile CCPP physics with 64-bit precision")
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
  set(CMAKE_Fortran_FLAGS_PHYSICS  "")
  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
endif()

#------------------------------------------------------------------------------
# Set netCDF flags for preprocessor, compiler and linker (if defined)
if(NEW_NCEPLIBS)
  #If using hera.intel, we use the target_link_libraries(NetCDF::NetCDF_[Fortran,C,CXX]), but ccpp-phyics inherits ${CMAKE_Fortran_FLAGS} for compiling physics, so we still need to set these
  ADD_DEFINITIONS(-DNETCDF)
  message (STATUS "Enable netCDF support")
  message (STATUS "NetCDF_Fortran_INCLUDE_DIRS is ${NetCDF_Fortran_INCLUDE_DIRS}")
  message (STATUS "NetCDF_Fortran_LIBRARIES is ${NetCDF_Fortran_LIBRARIES}")
  set (NETCDF_INC "-I${NetCDF_Fortran_INCLUDE_DIRS}")
  set (NETCDF_LIB ${NetCDF_Fortran_LIBRARIES})
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_INC} ${NETCDF_LIB}")
else()
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
endif()

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
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../ccpp/framework/cmake")
#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  # OpenMP broken for clang compiler
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang*")
    message(FATAL_ERROR "OpenMP not supported for the LLVM Clang compiler")
  endif()
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
ADD_SUBDIRECTORY(${CCPP_FRAMEWORK_SRC} ${CMAKE_BINARY_DIR}/ccpp/framework)
ADD_SUBDIRECTORY(${CCPP_PHYSICS_SRC} ${CMAKE_BINARY_DIR}/ccpp/physics)
ADD_DEPENDENCIES(ccpp_physics ccpp_framework)

SET(scm_source_files scm.F90
            scm_input.F90
            scm_utils.F90
            scm_vgrid.F90
            scm_setup.F90
            scm_forcing.F90
            scm_time_integration.F90
            scm_output.F90
)


ADD_EXECUTABLE(scm ${scm_source_files} ccpp_static_api.F90)
if(NEW_NCEPLIBS)
  #the FindNetCDF.cmake module suggests to use NetCDF::NetCDF_[Fortran,C,CXX] rather than set compiler flags; you actually need both Fortan and C components to compile successfully
  TARGET_LINK_LIBRARIES(scm ccpp_physics ccpp_framework NetCDF::NetCDF_Fortran NetCDF::NetCDF_C bacio::bacio_4 sp::sp_d w3nco::w3nco_d)
else()
  TARGET_LINK_LIBRARIES(scm ccpp_physics ccpp_framework ${BACIO_LIB4} ${SP_LIBd} ${W3NCO_LIBd})
endif()

set_target_properties(scm PROPERTIES
                               COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                               LINK_FLAGS "${CMAKE_Fortran_FLAGS}")

add_custom_command(
        TARGET scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_SOURCE_DIR}/run_scm.py
                ${CMAKE_CURRENT_BINARY_DIR}/run_scm.py)

add_custom_command(
        TARGET scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_SOURCE_DIR}/multi_run_scm.py
                ${CMAKE_CURRENT_BINARY_DIR}/multi_run_scm.py)
