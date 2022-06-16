#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Hera with icc/ifort"

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export SCM_ROOT=$MYDIR/../..

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module purge
module use /scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack
module load hpc/1.1.0
module load hpc-intel/2022.1.2
module load hpc-impi/2022.1.2
module load netcdf

echo "Setting up NCEPLIBS"
export bacio_ROOT=/scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/intel-2022.1.2/bacio/2.4.1
export sp_ROOT=/scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/intel-2022.1.2/sp/2.3.3
export w3nco_ROOT=/scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/intel-2022.1.2/w3nco/2.4.1

echo "Setting CC/CXX/FC environment variables"
export CC=icc
export CXX=icpc
export FC=ifort

echo "Loading cmake"
module load cmake/3.20.1
export CMAKE_C_COMPILER=icc
export CMAKE_CXX_COMPILER=icpc
export CMAKE_Fortran_COMPILER=ifort
export CMAKE_Platform=hera.intel

echo "Loading the SCM python environment"
. "/scratch1/BMC/gmtb/SCM_anaconda/etc/profile.d/conda.sh"
conda activate pyccpp
