help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using Intel oneAPI 2025.2.1
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Derecho with Intel compilers]===])


prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-2.0.0/envs/ue-oneapi-2025.2.1/modules/Core/")

load("stack-intel-oneapi-compilers/2025.2.1")
load("stack-cray-mpich/8.1.32")
load("python/3.11.11")
load("cmake/3.31.8")

load("hdf5/1.14.5")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.6.0")
load("ip/5.4.0")
load("w3emc/2.11.0")

load("py-f90nml")
load("py-netcdf4/1.7.2")

setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","derecho.intel")
