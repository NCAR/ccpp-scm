help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using GNU 12.4.0
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Derecho with GNU compilers]===])

prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.3/envs/ue-gcc-12.4.0/install/modulefiles/Core/")

load("stack-gcc/12.4.0")
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

setenv("CMAKE_C_COMPILER","mpicc")
setenv("CMAKE_CXX_COMPILER","mpicxx")
setenv("CMAKE_Fortran_COMPILER","mpif90")
setenv("CMAKE_Platform","derecho.gnu")
