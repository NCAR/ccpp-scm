#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on Hera with icc/ifort"

#load the modules in order to compile the GMTB SCM
echo "Loading intel and netcdf modules..."
module purge
module load intel/18.0.5.274
module load impi/2018.0.4
module load netcdf/4.7.0

echo "Setting CC/CXX/FC environment variables"
setenv CC icc
setenv CXX icpc
setenv FC ifort

echo "Setting NCEPLIBS environment variables"
##
## load modules for above compiler / MPI combination
##
module use -a /scratch1/BMC/gmtb/software/modulefiles/intel-18.0.5.274/impi-2018.0.4
module load NCEPlibs/1.1.0

echo "Loading cmake"
module load cmake/3.16.1
setenv CMAKE_C_COMPILER icc
setenv CMAKE_CXX_COMPILER icpc
setenv CMAKE_Fortran_COMPILER ifort
setenv CMAKE_Platform hera.intel

echo "Loading the anaconda python distribution"
module use -a /contrib/anaconda/modulefiles
module load anaconda/anaconda3-4.4.0

#install f90nml for the local user

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if ( $? != 0 ) then
	echo "Not found; installing f90nml"
	pip install -e git://github.com/marshallward/f90nml.git@v0.19\#egg=f90nml --user
else
	echo "f90nml is installed"
endif

#install shapely for the local user

#check to see if shapely is installed locally
echo "Checking if shapely python module is installed"
python -c "import shapely"

if ( $? != 0 ) then
	echo "Not found; installing shapely"
	pip install --index-url http://anaconda.rdhpcs.noaa.gov/simple --trusted-host anaconda.rdhpcs.noaa.gov shapely --user
else
	echo "shapely is installed"
endif

#check to see if configobj is installed locally
echo "Checking if configobj python module is installed"
python -c "import configobj"

if ( $? != 0 ) then
	echo "Not found; installing configobj"
	pip install --index-url http://anaconda.rdhpcs.noaa.gov/simple --trusted-host anaconda.rdhpcs.noaa.gov configobj --user
else
	echo "configobj is installed"
endif

#check to see if netCDF4 is installed locally
echo "Checking if netCDF4 python module is installed"
python -c "import netCDF4"

if ( $? != 0 ) then
	echo "Not found; installing netCDF4"
	pip install --index-url http://anaconda.rdhpcs.noaa.gov/simple --trusted-host anaconda.rdhpcs.noaa.gov netCDF4 --user
else
	echo "netCDF4 is installed"
endif
