#!/bin/tcsh

echo "Setting environment variables for CCPP-SCM on Hera with icc/ifort"

set called=($_)

if ( "$called" != "") then  ### called by source
    set MYSCRIPT=`readlink -f -n $called[2]`
else                        ### called by direct execution of the script
    set MYSCRIPT=`readlink -f -n '$0'`
endif
set MYDIR=`dirname $MYSCRIPT`
set MYDIR=`cd $MYDIR && pwd -P`

setenv SCM_ROOT $MYDIR/../..

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/2022.1.2
module load impi/2022.1.2
module use /scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/modulefiles/stack
module load hpc/1.2.0
module load hpc-intel/2022.1.2
module load hpc-impi/2022.1.2
module load netcdf

echo "Setting up NCEPLIBS"
setenv bacio_ROOT /scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/intel-2022.1.2/bacio/2.4.1
setenv sp_ROOT /scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/intel-2022.1.2/sp/2.3.3
setenv w3emc_ROOT /scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/intel-2022.1.2/w3emc/2.9.2

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

echo "Loading cmake"
module load cmake/3.20.1
setenv CMAKE_C_COMPILER icc
setenv CMAKE_CXX_COMPILER icpc
setenv CMAKE_Fortran_COMPILER ifort
setenv CMAKE_Platform hera.intel

echo "Loading the SCM python environment"
source /scratch1/BMC/gmtb/SCM_anaconda/etc/profile.d/conda.csh
conda activate pyccpp
