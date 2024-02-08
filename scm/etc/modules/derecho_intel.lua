help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using Intel-classic-2023.0.0
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Derecho with Intel compilers]===])

setenv("LMOD_TMOD_FIND_FIRST","yes")
load("ncarenv/23.09")
load("cmake/3.26.3")

prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/modulefiles")
prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")

load("stack-intel/2021.10.0")
load("stack-cray-mpich/8.1.25")
load("stack-python/3.10.8")
load("py-f90nml")
load("py-netcdf4/1.5.8")

load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.0")
load("bacio/2.4.1")
load("sp/2.3.3")
load("w3emc")

setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","derecho.intel")
