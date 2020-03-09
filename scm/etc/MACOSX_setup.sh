#!/bin/bash

echo "Setting environment variables for SCM-CCPP on MACOSX with clang/gfortran"

export CC=/usr/local/bin/mpicc
export CXX=/usr/local/bin/mpicxx
export F77=/usr/local/bin/mpif77
export F90=/usr/local/bin/mpif90
export FC=/usr/local/bin/mpif90
export CPP="/usr/local/bin/mpif90 -E -x f95-cpp-input"

export LDFLAGS="-L/usr/local/opt/zlib/lib -L/usr/local/opt/llvm/lib"
export CPPFLAGS="-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
export CFLAGS="-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
export CXXFLAGS="-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
export FFLAGS="-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
export FCFLAGS="-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"

export PATH="/usr/local/opt/llvm/bin${PATH:+:$PATH}"
export LD_LIBRARY_PATH="/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib:${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"

export NETCDF=/usr/local

if [! -v BACIO_LIB4]
then
echo "The BACIO_LIB4 environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
exit 1
fi

if [! -v SP_LIBd]
then
echo "The SP_LIBd environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
exit 1
fi

if [! -v W3NCO_LIBd]
then
echo "The W3NCO_LIBd environment variable must be set before proceeding."
echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.sh script to set the appropriate environment variables."
echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
exit 1
fi
