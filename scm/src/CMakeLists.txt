cmake_minimum_required(VERSION 3.14)

project(scm
        VERSION 5.0.0
        LANGUAGES C Fortran)
set(PROJECT "CCPP-SCM")

####################################################################
# Begin CCPP prebuild step                                         #
####################################################################
if(DEFINED CCPP_SUITES)
  if (${CCPP_SUITES} MATCHES "ALL")
    message("Calling CCPP code generator (ccpp_prebuild.py) for all available suites ...")
  else()
    set(_ccpp_suites_arg "--suites=${CCPP_SUITES}")
    message("Calling CCPP code generator (ccpp_prebuild.py) for suites ${_ccpp_suites_arg} ...")
  endif()
  unset(CCPP_SUITES CACHE)
else()
  execute_process(
    COMMAND scm/src/suite_info.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
    OUTPUT_VARIABLE suite_string
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(_ccpp_suites_arg "--suites=${suite_string}")
    message("Calling CCPP code generator (ccpp_prebuild.py) for SUPPORTED suites ${_ccpp_suites_arg} ...")
endif()
# Run CCPP prebuild.py
message (STATUS "Running ccpp_prebuild.py for CCPP")
# Make the directories where the ccpp_prebuild.py script wants to write caps and make snippets
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ccpp/physics/physics)
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  execute_process(
    COMMAND ccpp/framework/scripts/ccpp_prebuild.py --config=ccpp/config/ccpp_prebuild_config.py ${_ccpp_suites_arg} --builddir=${CMAKE_CURRENT_BINARY_DIR} --debug --verbose
    OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
    ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
    RESULT_VARIABLE return_code
    )
else()
  execute_process(
    COMMAND ccpp/framework/scripts/ccpp_prebuild.py --config=ccpp/config/ccpp_prebuild_config.py ${_ccpp_suites_arg} --builddir=${CMAKE_CURRENT_BINARY_DIR}
    OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
    ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../..
    RESULT_VARIABLE return_code
    )
endif()
  
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

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../../CMakeModules/Modules)

find_package(NetCDF REQUIRED COMPONENTS C Fortran)
find_package(bacio REQUIRED)
find_package(sp    REQUIRED)
find_package(w3emc REQUIRED)

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
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/physics)

#------------------------------------------------------------------------------
# Add required preprocessor flags for build type
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  ADD_DEFINITIONS(-DDEBUG)
endif ()

#------------------------------------------------------------------------------
# Add -DCCPP preprocessor flag (needed to preprocess GFS_typedefs.F90 from FV3)
ADD_DEFINITIONS(-DCCPP)

#------------------------------------------------------------------------------
# Add host-model specific preprocessor flag (needed for some physics schemes)
ADD_DEFINITIONS(-DSCM)

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
ADD_DEFINITIONS(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
#ADD_DEFINITIONS(-DINTERNAL_FILE_NML)

#------------------------------------------------------------------------------
# Detect OS and set preprocessor flags
if(APPLE)
  ADD_DEFINITIONS(-DMACOSX)
elseif(UNIX)
  ADD_DEFINITIONS(-DLINUX)
else ()
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
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
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
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
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
# Set netCDF flags for preprocessor
ADD_DEFINITIONS(-DNETCDF)

# Set OpenMP flags for C/C++/Fortran
if (OPENMP)
  # Clang compiler does not support OpenMP out of the box
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang*")
    message(FATAL_ERROR "OpenMP not supported for Apple/LLVM Clang compilers")
  endif()
  find_package(OpenMP REQUIRED)
endif()

# How about using proper compile targets etc?
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

if(OPENMP)
  TARGET_LINK_LIBRARIES(scm PUBLIC OpenMP::OpenMP_Fortran)
endif()
TARGET_LINK_LIBRARIES(scm NetCDF::NetCDF_Fortran)
TARGET_LINK_LIBRARIES(scm bacio::bacio_4)
TARGET_LINK_LIBRARIES(scm sp::sp_d)
TARGET_LINK_LIBRARIES(scm w3emc::w3emc_d)
TARGET_LINK_LIBRARIES(scm ccpp_framework)
TARGET_LINK_LIBRARIES(scm ccpp_physics)

SET_TARGET_PROPERTIES(scm PROPERTIES
                               COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                               LINK_FLAGS "${CMAKE_Fortran_FLAGS}")

ADD_CUSTOM_COMMAND(
        TARGET scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_SOURCE_DIR}/run_scm.py
                ${CMAKE_CURRENT_BINARY_DIR}/run_scm.py)

