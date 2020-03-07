#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

setenv CC /usr/local/bin/gcc
setenv CXX /usr/local/bin/g++
setenv F77 /usr/local/bin/gfortran
setenv F90 /usr/local/bin/gfortran
setenv FC /usr/local/bin/gfortran

setenv NETCDF /usr/local

echo "Running NCEPLIBS installation script for SCM-CCPP"
cd ..
./contrib/build_nceplibs.sh $PWD/nceplibs
cd scm
