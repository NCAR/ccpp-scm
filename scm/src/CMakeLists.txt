CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(scm)
set(PROJECT "CCPP-SCM")
ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

#MESSAGE("The NETCDF env. variable is $ENV{NETCDF}")
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})

SET(CCPP_SRC ${CMAKE_SOURCE_DIR}/../../ccpp-framework)
SET(GFSPHYSICS_SRC ${CMAKE_SOURCE_DIR}/../../ccpp-physics)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#Copy/link data to build directory
set(physics_data_dir ${CMAKE_SOURCE_DIR}/../data/GFS_physics_data)
file(GLOB physics_data_files "${physics_data_dir}/*")
file(COPY ${physics_data_files} DESTINATION ${CMAKE_BINARY_DIR}/)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    #message(STATUS "Setting build type to 'Debug' as none was specified.")
    #set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

ADD_COMPILE_OPTIONS(-O0)

ADD_SUBDIRECTORY(${CCPP_SRC} ${CMAKE_BINARY_DIR}/ccpp-framework)
ADD_SUBDIRECTORY(${GFSPHYSICS_SRC} ${CMAKE_BINARY_DIR}/ccpp-physics)

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

#this should call the custome FindNetCDF.cmake file in /src/cmake/modules and populate the NETCDF_INCLUDES and NETCDF_LIBRARIES variables with the appropriate values
SET(NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)

INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp-framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp-physics)

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
#  gfortran
#  MESSAGE("gfortran being used.")
  ADD_COMPILE_OPTIONS(-fcheck=all)
  ADD_COMPILE_OPTIONS(-fbacktrace)
  ADD_COMPILE_OPTIONS(-ffpe-trap=zero)
  ADD_COMPILE_OPTIONS(-finit-real=nan)
  ADD_COMPILE_OPTIONS(-ggdb)
  ADD_COMPILE_OPTIONS(-ffree-line-length-none)
  ADD_COMPILE_OPTIONS(-cpp)
  #set(f_flags -fcheck=all -fbacktrace -ffpe-trap=zero -finit-real=nan -ggdb -ffree-line-length-none -g -cpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
#  ifort
#  MESSAGE("ifort being used.")
  #ADD_COMPILE_OPTIONS(-check all)
  ADD_COMPILE_OPTIONS(-fpe0)
  ADD_COMPILE_OPTIONS(-warn)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-debug extended)
  ADD_COMPILE_OPTIONS(-fpp)
  #set(f_flags -g -fpe0 -warn -traceback -debug extended -fpp)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
#  pgf90
#  MESSAGE("pgf90 being used.")
  ADD_COMPILE_OPTIONS(-g)
  ADD_COMPILE_OPTIONS(-Mipa=noconst)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-Mfree)
  ADD_COMPILE_OPTIONS(-Mfptrap)
  ADD_COMPILE_OPTIONS(-Mpreprocess)
  #set(f_flags -g -Mipa=noconst -traceback -Mfree -Mfptrap -Mpreprocess)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${CMAKE_SOURCE_DIR}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#if(${CMAKE_VERSION} LESS 3.3)
#  string (REPLACE ";" " " f_flags_str "${f_flags}")
#  set(CMAKE_Fortran_FLAGS "${f_flags_str}")
#  message (STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
#else(${CMAKE_VERSION} LESS 3.3)
#  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
#endif (${CMAKE_VERSION} LESS 3.3)

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
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../ccpp-framework/cmake")
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

ADD_EXECUTABLE(gmtb_scm ${scm_source_files})
TARGET_LINK_LIBRARIES(gmtb_scm ccppphys ccpp ${NETCDF_LIBRARIES})
#set_target_properties(gmtb_scm PROPERTIES
#                               COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
#                               LINK_FLAGS "${CMAKE_Fortran_FLAGS}")

