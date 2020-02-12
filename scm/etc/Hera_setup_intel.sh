#!/bin/bash

echo "Setting environment variables for SCM-CCPP on Hera with icc/ifort"

#load the modules in order to compile the GMTB SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.6.1
module load pnetcdf/1.10.0

echo "Setting CC/CXX/FC environment variables"
export CC=icc
export CXX=icpc
export FC=ifort

echo "Setting NCEPLIBS_DIR environment variable"
module use -a /scratch1/BMC/gmtb/software/modulefiles/intel-18.0.5.274/impi-2018.0.4
module load NCEPlibs/9.9.9

echo "Loading the anaconda python distribution"
module load contrib
module load anaconda/anaconda2

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
