#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on MACOSX with clang/gfortran"

setenv CC /usr/local/bin/mpicc
setenv CXX /usr/local/bin/mpicxx
setenv F77 /usr/local/bin/mpif77
setenv F90 /usr/local/bin/mpif90
setenv FC /usr/local/bin/mpif90
setenv CPP "/usr/local/bin/mpif90 -E -x f95-cpp-input"

setenv LDFLAGS "-L/usr/local/opt/zlib/lib -L/usr/local/opt/llvm/lib"
setenv CPPFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv CFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv CXXFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv FFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv FCFLAGS " -I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"

if (! $?PATH) then
  setenv PATH "/usr/local/opt/llvm/bin"
else
  setenv PATH "/usr/local/opt/llvm/bin:$PATH"
endif
if (! $?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib"
else
  setenv LD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib:$LD_LIBRARY_PATH"
endif
if (! $?DYLD_LIBRARY_PATH) then
  setenv DYLD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib"
else
  setenv DYLD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib:$DYLD_LIBRARY_PATH"
endif

setenv NETCDF /usr/local

if (! $?BACIO_LIB4) then
  echo "The BACIO_LIB4 environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
  exit 1
endif

if (! $?SP_LIBd) then
  echo "The SP_LIBd environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
  exit 1
endif

if (! $?W3NCO_LIBd) then
  echo "The W3NCO_LIBd environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
  exit 1
endif
