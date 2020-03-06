#!/bin/bash

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

export CC=/opt/rh/devtoolset-8/root/usr/bin/gcc
export CXX=/opt/rh/devtoolset-8/root/usr/bin/g++
export F77=/opt/rh/devtoolset-8/root/usr/bin/gfortran
export F90=/opt/rh/devtoolset-8/root/usr/bin/gfortran
export FC=/opt/rh/devtoolset-8/root/usr/bin/gfortran

export NETCDF=/comsoftware/libs/netcdf

echo "Running NCEPLIBS installation script for SCM-CCPP"
cd ..
./contrib/build_nceplibs.sh $PWD/nceplibs
cd scm

echo "Setting NCEPLIBS environment variables for SCM-CCPP"
export BACIO_LIB4=/comsoftware/gmtb-scm/nceplibs/lib/libbacio_v2.2.0_4.a
export SP_LIBd=/comsoftware/gmtb-scm/nceplibs/lib/libsp_v2.1.0_d.a
export W3NCO_LIBd=/comsoftware/gmtb-scm/nceplibs/lib/libw3nco_v2.1.0_d.a

