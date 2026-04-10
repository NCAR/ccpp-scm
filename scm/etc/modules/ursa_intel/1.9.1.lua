help([[
This module loads libraries for building the CCPP Single-Column Model on
the NOAA RDHPC machine Ursa using Intel-2021.13.0
]])

whatis([===[Loads libraries needed for building the CCPP SCM on Ursa with Intel compilers ]===])

prepend_path("MODULEPATH", "/contrib/spack-stack/spack-stack-1.9.1/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

load("stack-oneapi/2024.2.1")
load("stack-intel-oneapi-mpi/2021.13")
load("stack-python/3.11.7")
load("py-f90nml")
load("py-netcdf4")
load("cmake/3.30.2")

load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.4.1")
load("sp/2.5.0")
load("w3emc/2.10.0")
load("esmf")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","ursa.intel")

execute{cmd="source /scratch3/BMC/gmtb/ccpp-scm-software/spack-stack-1.9.1/bin/activate", modeA={"load"}}
