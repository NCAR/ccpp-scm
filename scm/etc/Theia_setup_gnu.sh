#!/bin/bash

echo "Setting environment variables for SCM-CCPP on Theia with gcc/gfortran"

#load the modules in order to compile the GMTB SCM
echo "Loading gnu and netcdf modules..."
module purge
module load gcc/6.2.0
# netcdf-4.5.0, compiled with gnu/6.2.0 and mvapich2-2.2, and its dependencies
export PATH="/scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2/bin:/scratch4/BMC/gmtb/parallel-netcdf-1.8.1/gnu-6.2.0/mvapich2-2.2/bin:/scratch4/BMC/gmtb/hdf5-1.8.20/gnu-6.2.0/mvapich2-2.2/bin:/scratch4/BMC/gmtb/szip-2.1.1/gnu-6.2.0/bin:/scratch4/BMC/gmtb/zlib-1.2.11/gnu-6.2.0/bin:/scratch4/BMC/gmtb/mvapich2-2.2/gnu-6.2.0/bin${PATH:+:$PATH}"
export LD_LIBRARY_PATH="/scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2/lib:/scratch4/BMC/gmtb/parallel-netcdf-1.8.1/gnu-6.2.0/mvapich2-2.2/lib:/scratch4/BMC/gmtb/hdf5-1.8.20/gnu-6.2.0/mvapich2-2.2/lib:/scratch4/BMC/gmtb/szip-2.1.1/gnu-6.2.0/lib:/scratch4/BMC/gmtb/zlib-1.2.11/gnu-6.2.0/lib:/scratch4/BMC/gmtb/mvapich2-2.2/gnu-6.2.0/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MV2_ENABLE_AFFINITY=0
export NETCDF=/scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2

echo "Setting CC/CXX/FC environment variables"
export CC=gcc
export CXX=g++
export FC=gfortran

#prepend the anaconda installation to the path so that the anaconda version of python (with its many installed modules) is used; check if the path already contains the right path first
echo "Checking if the path to the anaconda python distribution is in PATH"
echo $PATH | grep '/contrib/ananconda/2.3.0/bin$' >&/dev/null
if [ $? -ne 0 ]; then
	echo "anaconda path not found in PATH; prepending anaconda path to PATH environment variable"
	export PATH=/contrib/anaconda/2.3.0/bin:$PATH
else
	echo "PATH already has the anaconda path in it"
fi

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
