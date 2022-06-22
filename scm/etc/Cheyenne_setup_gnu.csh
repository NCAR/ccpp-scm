#!/bin/tcsh

echo "Setting environment variables for CCPP-SCM on Cheyenne with gcc/gfortran"

set called=($_)

if ( "$called" != "") then  ### called by source
    set MYSCRIPT=`readlink -f -n $called[2]`
else                        ### called by direct execution of the script
    set MYSCRIPT=`readlink -f -n '$0'`
endif
set MYDIR=`dirname $MYSCRIPT`
set MYDIR=`cd $MYDIR && pwd -P`

setenv SCM_ROOT $MYDIR/../..

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
setenv bacio_ROOT /glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/bacio/2.4.1
setenv sp_ROOT /glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/sp/2.3.3
setenv w3nco_ROOT /glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/gnu-10.1.0/w3nco/2.4.1

echo "Setting CC/CXX/FC environment variables"
setenv CC gcc
setenv CXX g++
setenv FC gfortran

echo "Loading cmake"
module load cmake/3.16.4
setenv CMAKE_C_COMPILER gcc
setenv CMAKE_CXX_COMPILER g++
setenv CMAKE_Fortran_COMPILER gfortran
setenv CMAKE_Platform cheyenne.gnu

echo "Setting up python environment for running and plotting."
module load conda/latest

conda activate /glade/p/ral/jntp/GMTB/CCPP_SCM/conda/ccpp-scm

