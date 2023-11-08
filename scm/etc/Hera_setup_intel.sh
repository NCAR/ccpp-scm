#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Hera with icc/ifort"

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export SCM_ROOT=$MYDIR/../..

#load the modules in order to compile the CCPP SCM
echo "Loading spack-stack modules..."
module purge
#module use /scratch1/NCEPDEV/jcsda/jedipara/spack-stack/modulefiles
module use /scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core
module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load stack-python/3.10.8
module load netcdf-fortran/4.6.0
module load bacio/2.4.1
module load sp/2.3.3
module load w3emc/2.10.0
module load py-f90nml/1.4.3
module load cmake/3.23.1

export CMAKE_Platform=hera.intel

