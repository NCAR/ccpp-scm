help([[
This module loads libraries for building the CCPP Single-Column Model on
the NOAA RDHPC machine Hera using Intel-2022.1.2
]])

whatis([===[Loads libraries needed for building the CCPP SCM on Hera ]===])

prepend_path("MODULEPATH","/contrib/sutils/modulefiles")
load("sutils")

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.4.1/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.1"
load(pathJoin("cmake", cmake_ver))

load_any("netcdf/4.9.2","netcdf-c/4.9.2")
load_any("netcdf/4.9.2","netcdf-fortran/4.6.0")

load("bacio/2.4.1")
load("sp/2.3.3")
load("w3emc/2.9.2")

load(pathJoin("nccmp", os.getenv("nccmp_ver") or "1.9.0.1"))
load(pathJoin("nco", os.getenv("nco_ver") or "4.9.3"))
load("ufs-pyenv")

setenv("CMAKE_C_COMPILER","mpiicc")
setenv("CMAKE_CXX_COMPILER","mpiicpc")
setenv("CMAKE_Fortran_COMPILER","mpiifort")
setenv("CMAKE_Platform","hera.intel")

prepend_path("MODULEPATH","/scratch1/NCEPDEV/nems/role.epic/miniconda3/modulefiles")
load(pathJoin("miniconda3", os.getenv("miniconda3_ver") or "4.12.0"))

if mode() == "load" then
   LmodMsgRaw([===[Please do the following to activate conda:
       > conda activate /scratch1/BMC/gmtb/SCM_anaconda/envs/pyccpp
]===])
end

