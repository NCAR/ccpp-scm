#!/bin/tcsh

echo "Setting environment variables for CCPP-SCM on Hera with icc/ifort"

setenv SCM_ROOT $PWD

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.7.0

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

echo "Setting NCEPLIBS environment variables"
##
## load modules for above compiler / MPI combination
##
module use /scratch1/BMC/gmtb/software/NCEPLIBS-ufs-v2.0.0/intel-18.0.5.274/impi-2018.0.4/modules
module load NCEPLIBS/2.0.0

echo "Loading cmake"
module load cmake/3.16.1
setenv CMAKE_C_COMPILER icc
setenv CMAKE_CXX_COMPILER icpc
setenv CMAKE_Fortran_COMPILER ifort
setenv CMAKE_Platform hera.intel

echo "Loading the SCM python environment"
source /scratch1/BMC/gmtb/SCM_anaconda/etc/profile.d/conda.csh
conda activate pyccpp
