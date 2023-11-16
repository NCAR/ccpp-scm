help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using Intel-classic-2023.0.0
]])

whatis([===[Loads libraries needed for building the CCPP SCM on Derecho ]===])

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.26.3"))
load(pathJoin("ncarenv", os.getenv("ncarenv_ver") or "23.06"))
load(pathJoin("craype", os.getenv("craype_ver") or "2.7.20"))

unload("conda")
prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/modulefiles")
prepend_path("MODULEPATH","/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")
load("stack-gcc/12.2.0")
load("stack-cray-mpich/8.1.25")
load("stack-python/3.10.8")
load("py-f90nml")
load("py-netcdf4/1.5.8")

load("bacio/2.4.1")
load("sp/2.3.3")
load("w3emc")

setenv("CMAKE_C_COMPILER","mpicc")
setenv("CMAKE_CXX_COMPILER","mpicxx")
setenv("CMAKE_Fortran_COMPILER","mpif90")
setenv("CMAKE_Platform","derecho.gnu")
