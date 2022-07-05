#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Cheyenne with icc/ifort"

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export SCM_ROOT=$MYDIR/../..

#start with a "clean" environment; activate and deactivate ncar_pylib in order to successfully deactivate previously activated environment without errors
module load ncarenv/1.3
conda deactivate
module purge

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module load ncarenv/1.3
module load intel/2022.1
module load mpt/2.25
module load ncarcompilers/0.5.0
module load netcdf

echo "Setting up NCEPLIBS"
module use /glade/work/epicufsrt/GMTB/tools/intel/2022.1/hpc-stack-v1.2.0_6eb6/modulefiles/stack
module load hpc/1.2.0
module load hpc-intel/2022.1
module load hpc-mpt/2.25
export bacio_ROOT=/glade/work/epicufsrt/GMTB/tools/intel/2022.1/hpc-stack-v1.2.0_6eb6/intel-2022.1/bacio/2.4.1
export sp_ROOT=/glade/work/epicufsrt/GMTB/tools/intel/2022.1/hpc-stack-v1.2.0_6eb6/intel-2022.1/sp/2.3.3
export w3nco_ROOT=/glade/work/epicufsrt/GMTB/tools/intel/2022.1/hpc-stack-v1.2.0_6eb6/intel-2022.1/w3nco/2.4.1

echo "Setting CC/CXX/FC environment variables"
export CC=icc
export CXX=icpc
export FC=ifort

echo "Loading cmake"
module load cmake/3.22.0
export CMAKE_C_COMPILER=icc
export CMAKE_CXX_COMPILER=icpc
export CMAKE_Fortran_COMPILER=ifort
export CMAKE_Platform=cheyenne.intel

echo "Setting up python environment for running and plotting."
module load conda/latest

conda activate /glade/p/ral/jntp/GMTB/CCPP_SCM/conda/ccpp-scm

