help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using Intel-oneAPI
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Derecho with Intel oneAPI compilers]===])

-- setenv("LMOD_TMOD_FIND_FIRST","yes")
-- load("ncarenv/24.12")

prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.3/envs/ue-oneapi-2024.2.1/install/modulefiles/Core/")


load("stack-oneapi/2024.2.1")
load("stack-cray-mpich/8.1.29")
load("stack-python/3.11.7")
load("cmake/3.27.9")

load("hdf5/1.14.3")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.4.1")
load("sp/2.5.0")
load("w3emc/2.10.0")

load("py-f90nml")
load("py-netcdf4/1.7.1.post2")

setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","derecho.intel")
