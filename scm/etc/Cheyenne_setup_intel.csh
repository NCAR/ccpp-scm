#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on Cheyenne with icc/ifort"

#load the modules in order to compile the GMTB SCM
echo "Loading intel and netcdf modules..."
module purge
module load ncarenv/1.3
module load intel/19.0.2
module load ncarcompilers/0.5.0
module load mpt/2.19
module load netcdf/4.6.3
module load python/2.7.16

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

echo "Setting NCEPLIBS_DIR environment variable"
set NCEPLIBS_DIR = "/glade/p/ral/jntp/GMTB/tools/NCEPlibs/20190307/intel-19.0.2/mpt-2.19"
setenv NCEPLIBS_DIR $NCEPLIBS_DIR

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
