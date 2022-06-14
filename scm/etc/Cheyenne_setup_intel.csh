#!/bin/tcsh

echo "Setting environment variables for CCPP-SCM on Cheyenne with icc/ifort"

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
ncar_pylib
deactivate
module purge

#load the modules in order to compile the CCPP SCM
echo "Loading intel and netcdf modules..."
module load ncarenv/1.3
module load intel/19.1.1
module load mpt/2.19
module load ncarcompilers/0.5.0
module load netcdf/4.7.3

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

echo "Setting NCEPLIBS environment variables"
module use /glade/p/ral/jntp/GMTB/tools/NCEPLIBS-ufs-v2.0.0/intel-19.1.1/mpt-2.19/modules
module load  NCEPLIBS/2.0.0

echo "Loading cmake"
module load cmake/3.16.4
setenv CMAKE_C_COMPILER icc
setenv CMAKE_CXX_COMPILER icpc
setenv CMAKE_Fortran_COMPILER ifort
setenv CMAKE_Platform cheyenne.intel

echo "Setting up python environment for plotting. A NCAR Package Library for python will be cloned into /glade/work/$USER."
module load python/3.7.5
ncar_pylib
if (-d "/glade/work/$USER/ccpp_scm_python3_clone") then
    echo "ccpp_scm_python3_clone NPL exists. Loading..."
    ncar_pylib ccpp_scm_python3_clone 
else
    echo "ccpp_scm_python3_clone does not exist yet. Creating..."
    ncar_pylib -c 20200417 /glade/work/$USER/ccpp_scm_python3_clone
    ncar_pylib ccpp_scm_python3_clone
endif

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if ( $? != 0 ) then
        echo "Not found; installing f90nml"
        pip install --no-cache-dir f90nml==0.19
else
        echo "f90nml is installed"
endif
