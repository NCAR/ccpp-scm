#!/bin/bash

echo "Setting environment variables for SCM-CCPP on Cheyenne with pgcc/pgf90"

#load the modules in order to compile the GMTB SCM
echo "Loading pgi and netcdf modules..."
module purge
module load ncarenv/1.3
module load pgi/19.3
module load ncarcompilers/0.5.0
module load netcdf/4.6.3

echo "Setting CC/CXX/FC environment variables"
export CC=pgcc
export CXX=pgc++
export FC=pgf90

echo "Setting NCEPLIBS_DIR environment variable"
NCEPLIBS_DIR=/glade/p/ral/jntp/GMTB/tools/NCEPlibs/20190307/pgi-17.9/mpt-2.19
export NCEPLIBS_DIR=$NCEPLIBS_DIR

#install f90nml for the local user

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if [ $? -ne 0 ]; then
	echo "Not found; installing f90nml"
	cd etc/scripts/f90nml-0.19
	python setup.py install --user
	cd ../..
else
	echo "f90nml is installed"
fi
