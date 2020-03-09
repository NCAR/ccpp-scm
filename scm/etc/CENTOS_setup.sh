#!/bin/bash

echo "Setting environment variables for SCM-CCPP on CENTOS with gcc/gfortran"

export CC=/usr/local/bin/gcc
export CXX=/usr/local/bin/g++
export F77=/usr/local/bin/gfortran
export F90=/usr/local/bin/gfortran
export FC=/usr/local/bin/gfortran

export NETCDF=/usr/local

if [[ -z "${BACIO_LIB4}" ]]
then
echo "The BACIO_LIB4 environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
fi

if [[ -z "${SP_LIBd}" ]]
then
echo "The SP_LIBd environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
fi

if [[ -z "${W3NCO_LIBd}" ]]
then
echo "The W3NCO_LIBd environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
fi
