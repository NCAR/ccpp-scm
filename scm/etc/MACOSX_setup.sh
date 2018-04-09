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
