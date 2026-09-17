help([[
This module loads libraries for building the CCPP Single-Column Model on
the CISL machine Derecho (Cray) using Intel oneAPI 2025.3.1
]])

whatis([===[Loads spack-stack libraries needed for building the CCPP SCM on Derecho with Intel compilers]===])

local spack_root = "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-2.1.1/envs/ue-oneapi-2025.3.1"
prepend_path("MODULEPATH", pathJoin(spack_root, "modules/Core/"))
-- External compiler and Cray modules required by the stack's site/setup.sh.
prepend_path("MODULEPATH", "/opt/cray/pe/modulefiles")
prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/installs/oneapi-2025.3.1/modulefiles")
-- Expose libfabric/1.22.0 while retaining the stack's required crayenv.
append_path("MODULEPATH", "/glade/u/apps/derecho/modules/25.10/Core")


load("stack-intel-oneapi-compilers/2025.3.1")
load("stack-cray-mpich/8.1.32")
load("python/3.11.11")
load("cmake/3.31.8")

load("hdf5/1.14.5")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("bacio/2.6.0")
load("ip/5.4.0")
load("w3emc/2.13.0")

load("py-f90nml/1.4.3")
load("py-netcdf4/1.7.2")

-- emacs not in intel spack stack 2.1.1
prepend_path("PATH", "/glade/u/apps/derecho/25.10/spack/opt/spack/emacs/30.2/gcc/12.5.0/azrq/bin")

setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","derecho.intel")
