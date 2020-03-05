#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on Cheyenne with gcc/gfortran"

#load the modules in order to compile the GMTB SCM
echo "Loading gnu and netcdf modules..."
module purge
module load ncarenv/1.3
module load gnu/8.3.0
module load mpt/2.19
module load ncarcompilers/0.5.0
module load netcdf/4.7.3

echo "Setting CC/CXX/FC environment variables"
setenv CC gcc
setenv CXX g++
setenv FC gfortran

echo "Setting NCEPLIBS environment variables"
module use -a /glade/p/ral/jntp/GMTB/tools/modulefiles/gnu-8.3.0/mpt-2.19
module load  NCEPlibs/1.0.0beta02

echo "Loading cmake"
module load cmake/3.16.4
setenv CMAKE_Platform cheyenne.gnu

echo "Setting up python environment for plotting. A NCAR Package Library for python will be cloned into /glade/work/$USER."
module load python/2.7.16
ncar_pylib
if (-d "/glade/work/$USER/gmtb_scm_python_clone") then
    echo "gmtb_scm_python_clone NPL exists. Loading..."
    ncar_pylib gmtb_scm_python_clone
else
    echo "gmtb_scm_python_clone does not exist yet. Creating..."
    ncar_pylib -c 20190627 /glade/work/$USER/gmtb_scm_python_clone
    ncar_pylib gmtb_scm_python_clone
endif

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if ( $? != 0 ) then
        echo "Not found; installing f90nml"
        cd etc/scripts/f90nml-0.19
        pip install .
        cd ../..
else
        echo "f90nml is installed"
endif
