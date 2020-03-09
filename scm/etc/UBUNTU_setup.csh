#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on UBUNTU with gcc/gfortran"

setenv CC gcc
setenv CXX g++
setenv F77 gfortran
setenv F90 gfortran
setenv FC gfortran

setenv NETCDF /usr

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
