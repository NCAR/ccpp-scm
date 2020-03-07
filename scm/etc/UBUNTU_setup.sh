#!/bin/bash

echo "Setting environment variables for SCM-CCPP on UBUNTU with gcc/gfortran"

export CC=gcc
export CXX=g++
export F77=gfortran
export F90=gfortran
export FC=gfortran

export NETCDF=/usr

echo "Running NCEPLIBS installation script for SCM-CCPP"
cd ..
./contrib/build_nceplibs.sh $PWD/nceplibs
cd scm
