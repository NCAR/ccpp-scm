#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Hera with icc/ifort"

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.7.0

echo "Setting CC/CXX/FC environment variables"
export CC=icc
export CXX=icpc
export FC=ifort

echo "Setting NCEPLIBS environment variables"
module use /scratch1/BMC/gmtb/software/NCEPLIBS-ufs-v2.0.0/intel-18.0.5.274/impi-2018.0.4/modules
module load NCEPLIBS/2.0.0

echo "Loading cmake"
module load cmake/3.16.1
export CMAKE_C_COMPILER=icc
export CMAKE_CXX_COMPILER=icpc
export CMAKE_Fortran_COMPILER=ifort
export CMAKE_Platform=hera.intel

echo "Loading the SCM python environment"
. "/scratch1/BMC/gmtb/SCM_anaconda/etc/profile.d/conda.sh"
conda activate pyccpp
