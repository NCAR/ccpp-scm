#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on Cheyenne with pgcc/pgf90"

#load the modules in order to compile the GMTB SCM
echo "Loading pgi and netcdf modules..."
module purge
module load ncarenv/1.2
module load pgi/17.9
module load ncarcompilers/0.4.1
module load netcdf/4.4.1.1

echo "Setting CC/CXX/FC environment variables"
setenv CC pgcc
setenv CXX pgc++
setenv FC pgf90

#install f90nml for the local user

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if ( $? != 0 ) then
	echo "Not found; installing f90nml"
	cd etc/scripts/f90nml-0.19
	python setup.py install --user
	cd ../..
else
	echo "f90nml is installed"
endif
