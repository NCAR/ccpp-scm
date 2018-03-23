# User's Guide

This guide provides instructions for obtaining, compiling and running a simple
case for the GMTB single column model (SCM). The SCM code calls CCPP-compliant
physics schemes through the CCPP infrastructure code. As such, it requires the
CCPP infrastructure code and physics code, both of which are included as git
submodules within the SCM code. This package can be considered a simple example
for an atmospheric model to interact with physics through the CCPP.

## Obtaining Code
1. Download a compressed file or clone the source using
  * `git clone https://[username]@github.com/NCAR/gmtb-scm.git`
  and enter your github password when prompted.
2. Change directory into the project.
  * `cd gmtb-scm`
3. Initialize the CCPP infrastructure and physics submodules.
  * `git submodule init`
4. Update (download) the submodules.
  * `git submodule update`
  and, if asked, enter your github credentials again. If the machine is running an older
  version of git and you are denied access, you may need to configure the
  submodule URLs before repeating step 4 by executing this command:
    * `git config submodule.src/ccpp.url https://[username]@github.com/NCAR/gmtb-ccpp.git`
    * `git config submodule.gmtb-gfsphysics.url https://[username]@github.com/NCAR/gmtb-gfsphysics.git`

## Building and Compiling the SCM with CCPP
1. Run the CCPP prebuild script to match required physics variables with those
available from the dycore (SCM) and to generate physics caps and makefile
segments.
  * `cd gmtb-ccpp/scripts`
  * `./ccpp_prebuild.py`
  Note: use `./ccpp_prebuild.py --debug` to see the full output of the script.
2. Change directory to the top-level SCM directory.
  * `cd ../../gmtb-scm`
3. [Optional] Run the machine setup script if necessary. This script loads
compiler modules (Fortran 2003-compliant Intel), netCDF module, etc. and sets
compiler environment variables.
  * `source Theia_setup.csh` (for csh) or `. Theia_setup.sh` (for bash)
  * `source Cheyenne_setup.csh` (for csh) or `. Cheyenne_setup.sh` (for bash)
  * `source MACOSX_setup.csh` (for csh) or `. MACOSX_setup.sh` (for bash) if following the instructions in doc/README_MACOSX.txt
3. Make a build directory and change into it.
  * `mkdir bin && cd bin`
4. Invoke cmake on the source code to build.
  * `cmake ../src` (without threading/OpenMP)
  * `cmake -DOPENMP=1 ../src` (with threading/OpenMP)

  For extensive debugging output, add `-DCMAKE_BUILD_TYPE=Debug` to the `cmake` command.
5. Compile. Add `VERBOSE=1` to obtain more information on the build process.
  * `make`

## Running the SCM with CCPP
1. Run the SCM with the supplied case (twpice). The SCM will go through the time
 steps, applying forcing and calling the physics defined in the suite definition
 file.
  * `./gmtb_scm twpice`
2. A netcdf output file is generated in the location specified in the case
configuration file. For the twpice case, it is located in `../output_twpice/output.nc`

## Setting up the physics suite
First, a physics suite is defined using an XML file located in
src/ccpp/examples. The XML file contains the suite name, the number of suite
"parts" (suite parts exist so that an atmosphere "cap" can execute code between
calls to the physics driver), the number of subcycles for each scheme (if
physics schemes require smaller time steps than the dynamics), and the scheme
names. Scheme names found in the suite XML files must correspond to schemes
located within the physics directory.

Using a suite in the SCM framework involves specifying its name in the case
configuration file to be used. For example, for the twpice case
(case_config/twpice.nml), the variable 'physics suite' is set to the desired
suite name. NOTE: As mentioned in the 'Running' section above, since the schemes
 are in their own libraries, you must specify the path to the compiled scheme
 libraries that are being used in the suite by appending the `LD_LIBRARY_PATH`.
 Without this step, the scheme libraries will not be found at runtime.

## CCPP Implementation Notes
