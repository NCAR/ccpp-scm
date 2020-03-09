#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on MACOSX with clang/gfortran"

if (! $?CC) then
  setenv CC /usr/local/bin/mpicc
endif
if (! $?CXX) then
  setenv CXX /usr/local/bin/mpicxx
endif
if (! $?F77) then
  setenv F77 /usr/local/bin/mpif77
endif
if (! $?F90) then
  setenv F90 /usr/local/bin/mpif90
endif
if (! $?FC) then
  setenv FC /usr/local/bin/mpif90
endif
setenv CPP "${FC} -E -x f95-cpp-input"

if (! $?NETCDF) then
  setenv NETCDF /usr/local
endif

if (! $?BACIO_LIB4) then
  echo "The BACIO_LIB4 environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
endif

if (! $?SP_LIBd) then
  echo "The SP_LIBd environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
endif

if (! $?W3NCO_LIBd) then
  echo "The W3NCO_LIBd environment variable must be set before proceeding."
  echo "If NCEPLIBS were installed on this machine using the umbrella build, source the generated setenv_nceplibs.csh script to set the appropriate environment variables."
  echo "If the subset of NCEPLIBS required for SCM were installed previously using the supplied build_nceplibs.sh in the contrib directory, set the environment variables as instructed by that build script."
  echo "If NCEPLIBS have not been installed, run ./contrib/build_nceplibs.sh /path/to/where/you/want/to/install/NCEPLIBS from the top level gmtb-scm directory and set the environment variables as instructed by that script."
endif
