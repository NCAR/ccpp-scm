cmake_minimum_required(VERSION 3.20)

project(scm
        VERSION 8.0.0
        LANGUAGES C Fortran)
set(PROJECT "CCPP-SCM")

enable_testing()
include(cmake/ccpp_capgen.cmake)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules/Modules)

find_package(NetCDF REQUIRED COMPONENTS C Fortran)
find_package(bacio REQUIRED)
# Use ip@5 or later if available, fall back to sp
find_package(ip 5)
if(NOT ip_FOUND)
  find_package(sp  REQUIRED)
endif()
find_package(w3emc REQUIRED)

find_package(MPI REQUIRED Fortran)
if(NOT MPI_Fortran_HAVE_F08_MODULE)
  message(FATAL_ERROR "MPI implementation does not support the Fortran 2008 mpi_f08 interface")
endif()

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
#set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
#set(AVX             OFF CACHE BOOL "Enable AVX-I instruction set")
#set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
#set(RRTMGP_32BIT    OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in rte-rrtmgp radiation)")

##------------------------------------------------------------------------------
## Set compile options
#if (32BIT)
#  add_definitions(-DSINGLE_PREC)
#  add_definitions(-DRTE_USE_SP)
#endif()

#------------------------------------------------------------------------------
# GNU
#------------------------------------------------------------------------------
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message(STATUS "Compile using GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -ggdb -fbacktrace -cpp -fcray-pointer -ffree-line-length-none -fno-range-check")
  #set(CMAKE_Fortran_FLAGS_OPENMP_OFF "-fno-openmp")

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

  #set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
  #set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  #set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
  #set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
##------------------------------------------------------------------------------
## Intel oneAPI
##------------------------------------------------------------------------------
# elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "IntelLLVM")
#  message(STATUS "Compile using Intel OneAPI")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -align array64byte -qno-opt-dynamic-align")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align -sox -fp-model source")
#  set(CMAKE_Fortran_FLAGS_OPENMP_OFF "-fno-openmp")
#
#  if(NOT 32BIT)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
#  endif()
#
#  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
##    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fstack-protector-all -fpe0 -debug -ftrapuv -init=snan,arrays")
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -debug")
#    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ftrapuv")
#  else()
#    if(32BIT)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
#    else()
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3 -no-prec-div")
#    endif()
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -debug minimal")
#    set(FAST "-fast-transcendentals")
#    if(AVX2)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx2")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx2")
#    elseif(SIMDMULTIARCH)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axSSE4.2,CORE-AVX2")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -axSSE4.2,CORE-AVX2")
#    elseif(AVX)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx-i")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx-i")
#    endif()
#  endif()
#  set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
##------------------------------------------------------------------------------
## Intel Classic
##------------------------------------------------------------------------------
# elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
#  message(STATUS "Compile using Intel Classic")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fpp -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte -qno-opt-dynamic-align")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qno-opt-dynamic-align -sox -fp-model source")
#  set(CMAKE_Fortran_FLAGS_OPENMP_OFF "-fno-openmp")
#
#  if(NOT 32BIT)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
#  endif()
#
#  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check -fstack-protector-all -fpe0 -debug -ftrapuv -init=snan,arrays")
#    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ftrapuv")
#  else()
#    if(32BIT)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
#    else()
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3 -no-prec-div")
#    endif()
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -debug minimal")
#    set(FAST "-fast-transcendentals")
#    if(AVX2)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx2")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx2")
#    elseif(SIMDMULTIARCH)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axSSE4.2,CORE-AVX2")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -axSSE4.2,CORE-AVX2")
#    elseif(AVX)
#        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx-i")
#        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx-i")
#    endif()
#  endif()
#
#  set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
#
##------------------------------------------------------------------------------
## Nvidia HPC stack
##------------------------------------------------------------------------------
#elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "NVHPC")
#  message(STATUS "Compile using Nvidia HPC Stack")
#  if(ENABLE_NVIDIA_OPENACC MATCHES "true")
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -Minfo=accel")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -Minfo=accel")
#  else()
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
#  endif()
#  set(CMAKE_Fortran_FLAGS_OPENMP_OFF "-nomp")
#
#  if(NOT 32BIT)
#      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
#  endif()
#
#  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
#    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g")
#    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g")
#  else()
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
#  endif()
#
#  set(MPI_C_COMPILER mpicc)
#  set(MPI_CXX_COMPILER mpicxx)
#  set(MPI_Fortran_COMPILER mpif90)
#
#  set(CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_C_FLAGS_BITFORBIT       "-O2 -fPIC" CACHE STRING "" FORCE)
#  set(CMAKE_Fortran_FLAGS_BITFORBIT "-O2 -fPIC" CACHE STRING "" FORCE)
else()
  message (FATAL_ERROR "This program has only been compiled with gfortran and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
endif()

##------------------------------------------------------------------------------
## Set flag for 32bit dynamics build
#
#if(32BIT)
#  message(STATUS "Compile CCPP slow physics with 64-bit precision, fast physics with 32-bit precision")
#  add_definitions(-DOVERLOAD_R4)
#  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    # set(CMAKE_Fortran_FLAGS_PHYSICS  "-real-size 64 -no-prec-div")
#  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    # set(CMAKE_Fortran_FLAGS_PHYSICS  "-fdefault-real-8 -fdefault-double-8")
#  endif()
#  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
#
#else()
#  message(STATUS "Compile CCPP physics with 64-bit precision")
#  remove_definitions(-DOVERLOAD_R8)
#  remove_definitions(-DOVERLOAD_R4)
#  set(CMAKE_Fortran_FLAGS_PHYSICS  "")
#  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
#endif()
#
##------------------------------------------------------------------------------
## Set flags for rte-rrtmgp radiation
#if(RRTMGP_32BIT)
#  message(STATUS "Compile CCPP RTE-RRTMGP with 32-bit precision")
#  add_definitions(-DRTE_USE_SP)
#  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    set(CMAKE_Fortran_FLAGS_RTERRTMGP  "-real-size 32")
#  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(CMAKE_Fortran_FLAGS_RTERRTMGP  "-fno-default-real-8 -fdefault-double-8")
#  endif()
#else()
#  message(STATUS "Compile CCPP RTE-RRTMGP with 64-bit precision")
#  remove_definitions(-DRTE_USE_SP)
#  set(CMAKE_Fortran_FLAGS_RTERRTMGP  ${CMAKE_Fortran_FLAGS_PHYSICS})
#endif()
#
##------------------------------------------------------------------------------
## Set netCDF flags for preprocessor
#ADD_DEFINITIONS(-DNETCDF)
#
## Set OpenMP flags for C/C++/Fortran
#if (OPENMP)
#  # Clang compiler does not support OpenMP out of the box
#  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang*")
#    message(FATAL_ERROR "OpenMP not supported for Apple/LLVM Clang compilers")
#  endif()
#  find_package(OpenMP REQUIRED)
#endif()

add_subdirectory(ccpp)
add_subdirectory(scm/src)
