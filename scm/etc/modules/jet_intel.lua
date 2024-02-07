help([[
This module loads libraries for building the CCPP Single-Column Model on
the NOAA RDHPC machine Jet using Intel-2021.5.0
]])

whatis([===[Loads libraries needed for building the CCPP SCM on Jet with Intel compilers ]===])

prepend_path("MODULEPATH", "/lfs4/HFIP/hfv3gfs/spack-stack/modulefiles")

load("cmake/3.26.4")
load("miniconda/3.9.12")

prepend_path("MODULEPATH", "/lfs4/HFIP/hfv3gfs/role.epic/modulefiles")
prepend_path("MODULEPATH", "/mnt/lfs4/HFIP/hfv3gfs/role.epic/spack-stack/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")

load("stack-intel/2021.5.0")
load("stack-intel-oneapi-mpi/2021.5.1")
load("stack-python/3.10.8")
load("py-f90nml")
load("py-netcdf4/1.5.8")

load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.0")
load("bacio/2.4.1")
load("sp/2.3.3")
load("w3emc")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","jet.intel")
