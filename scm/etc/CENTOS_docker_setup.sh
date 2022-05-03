#!/bin/bash

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export SCM_ROOT=$MYDIR/../..

export CC=/opt/rh/devtoolset-9/root/usr/bin/gcc
export CXX=/opt/rh/devtoolset-9/root/usr/bin/g++
export F77=/opt/rh/devtoolset-9/root/usr/bin/gfortran
export F90=/opt/rh/devtoolset-9/root/usr/bin/gfortran
export FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran

export NETCDF=/comsoftware/libs/netcdf

echo "Running NCEPLIBS installation script for SCM-CCPP"
cd ..
./contrib/build_nceplibs.sh $PWD/nceplibs
cd scm
