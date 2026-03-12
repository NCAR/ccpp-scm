help([[
This module loads libraries for building the CCPP Single-Column Model on
the NOAA RDHPC machine Ursa using GNU 12.4
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Ursa with Intel compilers ]===])

prepend_path("MODULEPATH","/contrib/spack-stack/spack-stack-2.0.0/envs/ue-gcc-12.4.0/modules/Core")

load("stack-gcc/12.4.0")
load("stack-openmpi/4.1.6")
load("python/3.11.11")
load("cmake/3.31.8")

load("hdf5/1.14.5")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.6.0")
load("ip/5.4.0")
load("w3emc/2.11.0")
load("esmf/8.8.0")

load("py-f90nml")
load("py-netcdf4/1.7.2")

setenv("CMAKE_C_COMPILER","mpicc")
setenv("CMAKE_CXX_COMPILER","mpicxx")
setenv("CMAKE_Fortran_COMPILER","mpif90")
setenv("CMAKE_Platform","ursa.gnu")

execute{cmd="source /scratch3/BMC/gmtb/ccpp-scm-software/spack-stack-2.0.0-gnu/bin/activate", modeA={"load"}}
