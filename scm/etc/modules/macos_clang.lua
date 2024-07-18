help([[
This module loads libraries for building the CCPP Single-Column Model on
a MacOS machine with clang compilers
]])

whatis([===[Loads libraries needed for building the CCPP SCM on a MacOS machine with clang compilers]===])

local ssd=os.getenv("SPACK_STACK_DIR") or LmodError ("Environment variable SPACK_STACK_DIR is not set")
prepend_path("MODULEPATH", ssd .. "/envs/scm-test/install/modulefiles/Core")

load("stack-apple-clang")
load("stack-python/3.10.13")
load("stack-openmpi/4.1.6")


load("cmake/3.28.3")

load("py-f90nml/1.4.3")
load("py-netcdf4/1.5.8")

load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.4.1")
load("sp/2.5.0")
load("w3emc/2.10.0")

setenv("CMAKE_C_COMPILER","clang")
setenv("CMAKE_CXX_COMPILER","clang++")
setenv("CMAKE_Fortran_COMPILER","gfortran-12")
setenv("CMAKE_Platform","macos.clang")
