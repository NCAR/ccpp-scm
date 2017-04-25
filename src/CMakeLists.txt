CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(scm_ipd)
ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

#if(NOT DEFINED FC)
#    set(FC $ENV{FC} CACHE STRING "Fortran Compiler")
#endif()

#CMAKE_FORCE_Fortran_COMPILER(${FC} "Fortran Compiler")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules)

# FFLAGS depend on the compiler; set this variable name so that other CMakeLists.txt use the correct flags given the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#MESSAGE("The NETCDF env. variable is $ENV{NETCDF}")
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{NETCDF})

#By default NEMS_src is set to nems_trunk_atmos, unless it is specified differently as a -D variable when cmake is invoked
#SET(NEMS_SRC NEMSLegacy_trunk_r81662 CACHE STRING "The name of the directory containing the NEMS source code to compile.")
#MESSAGE("Using ${NEMS_SRC} as the NEMS directory.")

#SET(W3LIB_SRC ${CMAKE_SOURCE_DIR}/w3nco_v2.0.6)
#SET(GOCART_PATH NEMS/src/chem/gocart/src/GMAO_Shared)
#SET(GOCART_SRC ${CMAKE_SOURCE_DIR}/${NEMS_SRC}/${GOCART_PATH})
#SET(GFS_PHYS_PATH NEMS/src/atmos/gsm/gsmphys)
#SET(GFS_PHYS_SRC ${CMAKE_SOURCE_DIR}/${NEMS_SRC}/${GFS_PHYS_PATH})
SET(CCPP_SRC ${CMAKE_SOURCE_DIR}/ccpp)

ADD_COMPILE_OPTIONS(-O0)

#ADD_SUBDIRECTORY(${W3LIB_SRC})
#ADD_SUBDIRECTORY(${GOCART_SRC})
#ADD_SUBDIRECTORY(${GFS_PHYS_SRC})
ADD_SUBDIRECTORY(${CCPP_SRC})

SET(scm_source_files gmtb_scm.f90
            gmtb_scm_input.f90
            gmtb_scm_utils.f90
            gmtb_scm_vgrid.f90
            gmtb_scm_setup.f90
            gmtb_scm_forcing.f90
            gmtb_scm_time_integration.f90
            gmtb_scm_output.f90
            gmtb_scm_kinds.f90
            gmtb_scm_physical_constants.f90
)

#this should call the custome FindNetCDF.cmake file in /src/cmake/modules and populate the NETCDF_INCLUDES and NETCDF_LIBRARIES variables with the appropriate values
SET(NETCDF_F90 "YES")
FIND_PACKAGE(NetCDF)

#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/w3nco_v2.0.6)
#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/${NEMS_SRC}/${GOCART_PATH})
#INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/${NEMS_SRC}/${GFS_PHYS_PATH})
INCLUDE_DIRECTORIES(${NETCDF_INCLUDES})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/ccpp/src)

#ADD_COMPILE_OPTIONS(-Wuninitialized)
#Not safe to use optimization yet with gmtb_scm code
#ADD_COMPILE_OPTIONS(-O0)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  # gfortran
#  MESSAGE("gfortran being used.")
  ADD_COMPILE_OPTIONS(-fcheck=all)
  ADD_COMPILE_OPTIONS(-fbacktrace)
  ADD_COMPILE_OPTIONS(-ffpe-trap=zero)
  ADD_COMPILE_OPTIONS(-finit-real=nan)
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#  # ifort (untested)
#  MESSAGE("ifort being used.")
  #ADD_COMPILE_OPTIONS(-check all)
  ADD_COMPILE_OPTIONS(-fpe0)
  ADD_COMPILE_OPTIONS(-warn)
  ADD_COMPILE_OPTIONS(-traceback)
  ADD_COMPILE_OPTIONS(-debug extended)
elseif (Fortran_COMPILER_NAME MATCHES "g77")
#  # g77
#  MESSAGE("g77 being used.")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#  message ("No optimized Fortran compiler flags are known, we just try -O2...")
#  #set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#  #set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

ADD_EXECUTABLE(gmtb_scm ${scm_source_files})
#TARGET_LINK_LIBRARIES(gmtb_scm w3 gfs_phys ${NETCDF_LIBRARIES})
TARGET_LINK_LIBRARIES(gmtb_scm ccpp ${NETCDF_LIBRARIES})
