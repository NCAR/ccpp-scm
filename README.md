# User's Guide

A more complete User's Guide can be found at https://dtcenter.org/GMTB/v3.0/scm-ccpp-guide-v3.pdf. If the instructions in this README and the more complete User's Guide differ, the linked guide should be more up-to-date and accurate.

This guide provides instructions for obtaining, compiling and running a simple
case for the GMTB single column model (SCM). The SCM code calls CCPP-compliant
physics schemes through the CCPP infrastructure code. As such, it requires the
CCPP infrastructure code and physics code, both of which are included as git
submodules within the SCM code. This package can be considered a simple example
for an atmospheric model to interact with physics through the CCPP.

## Prerequisite
There are several utility libraries as part of the NCEPlibs package that must be installed prior to building the SCM.
* bacio - Binary I/O Library
* sp - Spectral Transformation Library
* w3nco - GRIB decoder and encoder library

These libraries are prebuilt on most NOAA machines using the Intel compiler. For those needing to build the libraries themselves, GMTB recommends using the source code from GitHub at https://github.com/NCAR/NCEPlibs.git, which includes build files for various compilers and machines using OpenMP flags and which are threadsafe. Instructions for installing NCEPlibs are included on the GitHub repository webpage, but for the sake of example, execute the following for obtaining and building from source in /usr/local/NCEPlibs on a Mac:
1. `cd /usr/local/src`
2. `git clone https://github.com/NCAR/NCEPlibs.git`
3. `cd NCEPlibs`
4. `./make_ncep_libs.sh -s macosx -c gnu -d /usr/local/NCEPlibs -o 1 -m 0`

Note that the option `-m 0` can be used if MPI is not installed on the machine that is being used. The nemsio library will not be installed, however, since it requires MPI. Once NCEPlibs is built, the NCEPLIBS_DIR environment variable must be set to the location of the installation. For example, if NCEPlibs was installed in /usr/local/NCEPlibs, one would execute

`export NCEPLIBS_DIR=/usr/local/NCEPlibs`

If using Theia or Cheyenne HPC systems, this environment variable is automatically set to an appropriate installation of NCEPlibs on those machines through use of one of the setup scripts described below.

## Obtaining Code

For obtaining the last stable release, execute the following:

1. Clone the source using:
  * `git clone --recursive -b v3.0 https://github.com/NCAR/gmtb-scm`
2. Change directory into the project.
  * `cd gmtb-scm`

For working with the development branches, after executing the steps above, check out the master branches of the repository (and submodules):

1. `git checkout master`
2. `cd ccpp/physics`
3. `git checkout master`
4. `cd ../framework`
5. `git checkout master`
6. `cd ../..`

## Building and Compiling the SCM with CCPP
1. Run the CCPP prebuild script to match required physics variables with those
available from the dycore (SCM) and to generate physics caps and makefile
segments.
  * `./ccpp/framework/scripts/ccpp_prebuild.py --config=./ccpp/config/ccpp_prebuild_config.py`
  Note: add `--debug` to see the full output of the script.
2. Change directory to the top-level SCM directory.
  * `cd scm`
3. [Optional] Run the machine setup script if necessary. This script loads
compiler modules (Fortran 2003-compliant), netCDF module, etc. and sets
compiler environment variables.
  * `source etc/Theia_setup_gnu.csh` (for csh) or `. etc/Theia_setup_gnu.sh` (for bash)
  * `source etc/Theia_setup_intel.csh` (for csh) or `. etc/Theia_setup_intel.sh` (for bash)
  * `source etc/Theia_setup_pgi.csh` (for csh) or `. etc/Theia_setup_pgi.sh` (for bash)
  * `source etc/Cheyenne_setup_gnu.csh` (for csh) or `. etc/Cheyenne_setup_gnu.sh` (for bash)
  * `source etc/Cheyenne_setup_intel.csh` (for csh) or `. etc/Cheyenne_setup_intel.sh` (for bash)
  * `source etc/Cheyenne_setup_pgi.csh` (for csh) or `. etc/Cheyenne_setup_pgi.sh` (for bash)
  * `source etc/UBUNTU_setup.csh` (for csh) or `. etc/UBUNTU_setup.sh` (for bash) if following the instructions in doc/README_UBUNTU.txt
  * `source etc/CENTOS_setup.csh` (for csh) or `. etc/CENTOS_setup.sh` (for bash) if following the instructions in doc/README_CENTOS.txt
  * `source etc/MACOSX_setup.csh` (for csh) or `. etc/MACOSX_setup.sh` (for bash) if following the instructions in doc/README_MACOSX.txt
  * NOTE: The NETCDF environment variable must be set to the path of the netCDF installation that was compiled with the same compiler used in the following steps.
4. Make a build directory and change into it.
  * `mkdir bin && cd bin`
5. Invoke cmake on the source code to build.
  * `cmake ../src` (without threading/OpenMP)
  * `cmake -DOPENMP=1 ../src` (with threading/OpenMP)
  * `cmake -DCMAKE_BUILD_TYPE=Debug ../src` (debug mode)
6. Compile. Add `VERBOSE=1` to obtain more information on the build process.
  * `make`

## Running the SCM with CCPP
1. Run the SCM with a supplied case. The SCM will go through the time
 steps, applying forcing and calling the physics defined in the suite definition
 file.
  * `.run_gmtb_scm.py -c CASE_NAME [-s SUITE_NAME] [-n PHYSICS_NAMELIST_PATH] [-g]`
  * When invoking the run script, the only required argument is the name of the case to run. The case name used must match one of the case configuration files located in ../etc/case_config (without the .nml extension!). If specifying a suite other than the default, the suite name used must match the value of the suite name in one of the suite definition files located in ../../ccpp/suites, (e.g. `SCM_GFS_v15`). If specifying a namelist other than the default, the value must be an entire filename that exists in ../../ccpp/physics_namelists. The -g flag can be used to run the executable through the gdb debugger (assuming it is installed on the system).
2. A netcdf output file is generated in the location specified in the case
configuration file (is present), or in an output directory created by default in `bin` with the case name and suite name appended.

## Running the SCM with FV3GFS initial conditions
model initial conditions are needed to initialize the land surface in order to run with an interactive land model
1. Prepare model initial conditions.
  * `cd to scm/etc/scripts/`
  modify path to files you can read in extract_FV3GFS_column_ic.py, this is set up for C96
2. run extract_FV3GFS_column_ic.py, it will create fv3_model_point.nc in ../../data/processed_case_input/
  * `./extract_FV3GFS_column_ic.py`
3. cd to bin directory
  * `cd ../../bin/`
4. Run the SCM with the fv3_model_point case and C96 namelist
  * `./run_gmtb_scm.py -c fv3_model_point -n input_GFS_v15_C96.nml`
