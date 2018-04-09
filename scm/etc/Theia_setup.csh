#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on Theia with icc/ifort"

#load the modules in order to compile the GMTB SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/18.1.163
module load netcdf/4.3.0

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

#prepend the anaconda installation to the path so that the anaconda version of python (with its many installed modules) is used; check if the path already contains the right path first
echo "Checking if the path to the anaconda python distribution is in PATH"
echo $PATH | grep '/contrib/ananconda/2.3.0/bin$' >&/dev/null
if ( $? != 0 ) then
	echo "anaconda path not found in PATH; prepending anaconda path to PATH environment variable"
	setenv PATH /contrib/anaconda/2.3.0/bin:$PATH
else
	echo "PATH already has the anaconda path in it"
endif

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
