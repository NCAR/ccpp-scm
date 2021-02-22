#!/bin/bash

echo "Setting environment variables for CCPP-SCM on Cheyenne with gcc/gfortran"

#load the modules in order to compile the CCPP SCM
echo "Loading gnu and netcdf modules..."
module purge
module load ncarenv/1.3
module load gnu/8.3.0
module load mpt/2.19
module load ncarcompilers/0.5.0
module load netcdf/4.7.3

echo "Setting CC/CXX/FC environment variables"
export CC=gcc
export CXX=g++
export FC=gfortran

echo "Setting NCEPLIBS environment variables"
module use /glade/p/ral/jntp/GMTB/tools/modulefiles/gnu-8.3.0/mpt-2.19
module load  NCEPlibs/1.1.0

echo "Loading cmake"
module load cmake/3.16.4
export CMAKE_C_COMPILER=gcc
export CMAKE_CXX_COMPILER=g++
export CMAKE_Fortran_COMPILER=gfortran
export CMAKE_Platform=cheyenne.gnu

echo "Setting up python environment for plotting. A NCAR Package Library for python will be cloned into /glade/work/$USER."
module load python/3.7.5
ncar_pylib
if [ -d "/glade/work/$USER/ccpp_scm_python3_clone" ]; then
    echo "ccpp_scm_python3_clone NPL exists. Loading..."
    ncar_pylib ccpp_scm_python3_clone
else
    echo "ccpp_scm_python3_clone does not exist yet. Creating..."
    ncar_pylib -c 20200417 /glade/work/$USER/ccpp_scm_python3_clone
    ncar_pylib ccpp_scm_python3_clone
fi

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if [ $? -ne 0 ]; then
        echo "Not found; installing f90nml"
        pip install --no-cache-dir f90nml==0.19
else
        echo "f90nml is installed"
fi
