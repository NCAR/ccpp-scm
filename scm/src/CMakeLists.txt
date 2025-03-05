
cmake_minimum_required(VERSION 3.14)

project(scm
        VERSION 5.0.0
        LANGUAGES C Fortran)
set(PROJECT "CCPP-SCM")

#------------------------------------------------------------------------------
#
# Set where the CCPP Framework lives
#
#------------------------------------------------------------------------------
get_filename_component(TEST_ROOT "${CMAKE_SOURCE_DIR}" DIRECTORY)
get_filename_component(CCPP_ROOT "${TEST_ROOT}" DIRECTORY)

#------------------------------------------------------------------------------
#
# Create list of SCHEME_FILES and HOST_FILES
# Paths should be relative to CMAKE_SOURCE_DIR (this file's directory)
#
#------------------------------------------------------------------------------
LIST(APPEND SCHEME_FILES "${CCPP_ROOT}/ccpp/config/scheme_files.txt")
LIST(APPEND HOST_FILES   "${CCPP_ROOT}/ccpp/config/host_files.txt")

#------------------------------------------------------------------------------
#
# Which Suites are we building with?
#
#------------------------------------------------------------------------------ 
if(DEFINED CCPP_SUITES)
  set(suite_string ${CCPP_SUITES})
  message(STATUS "Calling CCPP code generator (ccpp_capgen.py) for suites ${suite_string} ...")
  unset(CCPP_SUITES CACHE)
else()
  execute_process(
    COMMAND scm/src/suite_info.py
    WORKING_DIRECTORY ${CCPP_ROOT}
    OUTPUT_VARIABLE suite_string
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
    set(_ccpp_suites_arg "--suites=${suite_string}")
    message(STATUS "Calling CCPP code generator (ccpp_capgen.py) for suites listed in scm/src/suite_info.py: ${_ccpp_suites_arg}.")
endif()

#------------------------------------------------------------------------------
#
# End of project-specific input
#
#------------------------------------------------------------------------------

# By default, no verbose output
SET(VERBOSITY 0 CACHE STRING "Verbosity level of output (default: 0)")
# By default, generated caps go in ccpp subdir
SET(CCPP_CAP_FILES "${CMAKE_BINARY_DIR}/ccpp" CACHE
  STRING "Location of CCPP-generated cap files")

SET(CCPP_FRAMEWORK ${CCPP_ROOT}/ccpp/framework/scripts)

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CCPP_FRAMEWORK}/cmake")
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
#
# Create metadata and source file lists from ${SCHEME_FILES}
#
#------------------------------------------------------------------------------
FOREACH(FILE ${SCHEME_FILES})
  FILE(STRINGS ${FILE} FILENAMES)
  LIST(APPEND SCHEME_FILENAMES ${FILENAMES})
  FOREACH(FILE ${SCHEME_FILENAMES})
    MATH(EXPR NSCHEME_FILES "${NSCHEME_FILES}+1")
  ENDFOREACH(FILE)
ENDFOREACH(FILE)
#
MATH(EXPR FCOUNT 0)
FOREACH(FILE ${SCHEME_FILENAMES})
  # Set character seperator
  MATH(EXPR FCOUNT "${FCOUNT}+1")
  set(SEPERATOR "")
  string(CONCAT SCHEME_FILENAME ${CCPP_ROOT} "/" ${FILE} ${SEPERATOR} )
  string(REGEX REPLACE "\n$" "" SCHEME_FILENAME "${SCHEME_FILENAME}")
  list(APPEND SCHEME_SOURCE ${SCHEME_FILENAME})
  # Meta file list
  string(REGEX REPLACE "\\.[^.]*$" ".meta" METADATA_FILE ${FILE})
  string(CONCAT METADATA_FILENAME ${CCPP_ROOT} "/" ${METADATA_FILE} ${SEPERATOR} )
  list(APPEND SCHEME_METADATA ${METADATA_FILENAME})
ENDFOREACH(FILE)
string(REPLACE ";" "," SCHEME_METADATA "${SCHEME_METADATA}")

#------------------------------------------------------------------------------
#
# Create metadata and source file lists	from ${HOST_FILES}
#
#------------------------------------------------------------------------------
FOREACH(FILE ${HOST_FILES})
  FILE(STRINGS ${FILE} FILENAMES)
  LIST(APPEND HOST_FILENAMES "${FILENAMES}")
  FOREACH(FILE ${HOST_FILENAMES})
    MATH(EXPR NHOST_FILES "${NHOST_FILES}+1")
  ENDFOREACH(FILE)
ENDFOREACH(FILE)
#
MATH(EXPR FCOUNT 0)
FOREACH(FILE ${HOST_FILENAMES})
  # Set character seperator
  MATH(EXPR FCOUNT "${FCOUNT}+1")
  set(SEPERATOR "")
  string(CONCAT HOST_FILENAME ${CCPP_ROOT} "/" ${FILE} ${SEPERATOR} )
  list(APPEND HOST_SOURCE ${HOST_FILENAME})
  #
  string(REGEX REPLACE "\\.[^.]*$" ".meta" METADATA_FILE ${FILE})
  string(CONCAT METADATA_FILENAME ${CCPP_ROOT} "/" ${METADATA_FILE} ${SEPERATOR} )
  list(APPEND HOST_METADATA ${METADATA_FILENAME})
ENDFOREACH(FILE)
string(REPLACE ";" "," HOST_METADATA "${HOST_METADATA}")

#------------------------------------------------------------------------------
#
# Create SDF list from ${suite_string), which is provided via the command line,
# or comes from scm/src/suite_info.py
#
#------------------------------------------------------------------------------
string(REPLACE "," ";" suite_string "${suite_string}")
set(suite_list ${suite_string})
FOREACH(suite ${suite_list})
  list(APPEND SUITES "${CCPP_ROOT}/ccpp/suites/suite_${suite}.xml")
ENDFOREACH(suite)
string(REPLACE ";" "," SUITES "${SUITES}")

#------------------------------------------------------------------------------
#
# Run ccpp_capgen
#
#------------------------------------------------------------------------------

set(CAPGEN_CMD "${CCPP_FRAMEWORK}/ccpp_capgen.py")
list(APPEND CAPGEN_CMD "--host-files")
list(APPEND CAPGEN_CMD ${HOST_METADATA})
list(APPEND CAPGEN_CMD "--scheme-files")
list(APPEND CAPGEN_CMD ${SCHEME_METADATA})
list(APPEND CAPGEN_CMD "--suites")
list(APPEND CAPGEN_CMD "${SUITES}")
list(APPEND CAPGEN_CMD "--host-name")
list(APPEND CAPGEN_CMD "SCM")
list(APPEND CAPGEN_CMD "--output-root")
list(APPEND CAPGEN_CMD "${CCPP_CAP_FILES}")
list(APPEND CAPGEN_CMD "--kind-type")
list(APPEND CAPGEN_CMD "kind_phys=REAL64")
list(APPEND CAPGEN_CMD "--kind-type")
list(APPEND CAPGEN_CMD "kind_dyn=REAL64")
list(APPEND CAPGEN_CMD "--kind-type")
list(APPEND CAPGEN_CMD "kind_sngl_prec=REAL32")
list(APPEND CAPGEN_CMD "--kind-type")
list(APPEND CAPGEN_CMD "kind_dbl_prec=REAL64")
list(APPEND CAPGEN_CMD "--ccpp-cfgfile")
list(APPEND CAPGEN_CMD "${CCPP_ROOT}/ccpp/config/ccpp_capgen_config.py")
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  list(APPEND CAPGEN_CMD "--debug")
endif()
while (VERBOSITY GREATER 0)
  list(APPEND CAPGEN_CMD "--verbose")
  MATH(EXPR VERBOSITY "${VERBOSITY} - 1")
endwhile ()
string(REPLACE ";" " " CAPGEN_STRING "${CAPGEN_CMD}")

MESSAGE(STATUS "Running: ${CAPGEN_STRING}")
EXECUTE_PROCESS(COMMAND ${CAPGEN_CMD} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CAPGEN_OUT ERROR_VARIABLE CAPGEN_OUT RESULT_VARIABLE RES)
MESSAGE(STATUS "${CAPGEN_OUT}")
if (RES EQUAL 0)
  MESSAGE(STATUS "CCPP cap generation completed")
else(RES EQUAL 0)
  MESSAGE(FATAL_ERROR "CCPP cap generation FAILED: result = ${RES}")
endif(RES EQUAL 0)

MESSAGE(STATUS "CCPP capgen step complete")

####################################################################
# End CCPP capgen step                                             #
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

INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/physics)
include_directories(${MPI_Fortran_INCLUDE_PATH})

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

ADD_DEFINITIONS(-DMPI)


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

include(${CCPP_ROOT}/scm/bin/ccpp/physics/CCPP_CAPS.cmake)
include(${CCPP_ROOT}/scm/src/CCPP_STATIC_API.cmake)
ADD_EXECUTABLE(scm ${scm_source_files} ${API})

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

