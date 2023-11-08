help([[
This module loads libraries for building the CCPP Single-Column Model on
the NOAA RDHPC machine Hera using GNU 9.2.0
]])

whatis([===[Loads libraries needed for building the UFS SRW App on Hera using GNU 9.2.0 ]===])

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.4.1/envs/unified-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/scratch1/NCEPDEV/jcsda/jedipara/spack-stack/modulefiles")

load("stack-gcc/9.2.0")
load("stack-openmpi/4.1.5")
load("stack-python/3.9.12")
load("cmake/3.23.1")

load("netcdf/4.7.2")

load("bacio/2.4.1")
load("sp/2.3.3")
load("w3emc/2.9.2")

load(pathJoin("nccmp", os.getenv("nccmp_ver") or "1.9.0.1"))
load(pathJoin("nco", os.getenv("nco_ver") or "5.0.6"))
load(pathJoin("openblas", os.getenv("openblas_ver") or "0.3.19"))
load("ufs-pyenv")

setenv("CMAKE_C_COMPILER","mpicc")
setenv("CMAKE_CXX_COMPILER","mpicxx")
setenv("CMAKE_Fortran_COMPILER","mpif90")
setenv("CMAKE_Platform","hera.gnu")

prepend_path("MODULEPATH","/scratch1/NCEPDEV/nems/role.epic/miniconda3/modulefiles")
load(pathJoin("miniconda3", os.getenv("miniconda3_ver") or "4.12.0"))

if mode() == "load" then
   LmodMsgRaw([===[Please do the following to activate conda:
       > conda activate /scratch1/BMC/gmtb/SCM_anaconda/envs/pyccpp
]===])
end

