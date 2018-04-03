#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on UBUNTU with gcc/gfortran"

setenv CC gcc
setenv CXX g++
setenv F77 gfortran
setenv F90 gfortran
setenv FC gfortran

setenv NETCDF /usr
