#!/bin/bash

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

export CC=/usr/local/bin/gcc
export CXX=/usr/local/bin/g++
export F77=/usr/local/bin/gfortran
export F90=/usr/local/bin/gfortran
export FC=/usr/local/bin/gfortran

export NETCDF=/usr/local

