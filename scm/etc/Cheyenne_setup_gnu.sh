#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Cheyenne with gcc/gfortran"

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export SCM_ROOT=$MYDIR/../..

#start with a "clean" environment; activate and deactivate ncar_pylib in order to successfully deactivate previously activated environment without errors
module load ncarenv/1.3
conda deactivate
module purge

#load the modules in order to compile the CCPP SCM
echo "Loading gnu and netcdf modules..."
module load ncarenv/1.3
module load gnu/10.1.0
module load mpt/2.22
module load ncarcompilers/0.5.0
module load netcdf

echo "Setting up NCEPLIBS"
module use /glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/modulefiles/stack
module load hpc/1.2.0
module load hpc-gnu/10.1.0
module load hpc-mpt/2.22
export bacio_ROOT=/glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/bacio/2.4.1
export sp_ROOT=/glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/sp/2.3.3
export w3emc_ROOT=/glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/w3emc/2.9.2

echo "Setting CC/CXX/FC environment variables"
export CC=gcc
export CXX=g++
export FC=gfortran

echo "Loading cmake"
module load cmake/3.22.0
export CMAKE_C_COMPILER=gcc
export CMAKE_CXX_COMPILER=g++
export CMAKE_Fortran_COMPILER=gfortran
export CMAKE_Platform=cheyenne.gnu

echo "Setting up python environment for running and plotting."
module load conda/latest

conda activate /glade/p/ral/jntp/GMTB/CCPP_SCM/conda/ccpp-scm

