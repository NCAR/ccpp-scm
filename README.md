# User's Guide

## Obtaining Code
1. Download a compressed file or clone the source using
* `git clone https://[username]@github.com/NCAR/gmtb-scm.git`
  and enter your github password when prompted.
2. Change directory into the project.
* `cd gmtb-scm`
3. Initialize the CCPP submodule.
* `git submodule init`
4. Update the CCPP submodule.
* `git submodule update`
  and enter your github credentials again. If the machine is running an older version of git and you are denied access, you may need to configure the CCPP submodule URL before repeating step 4 by executing this command:
* `git config submodule.src/ccpp.url https://[username]@github.com/NCAR/gmtb-ccpp.git`
5. Change directory into the CCPP source.
* `cd src/ccpp`
6. Checkout a branch or tag.
* `git checkout master` or
* `git checkout v0.1.0`

## Building and Compiling the SCM with CCPP
1. Change directory to the top-level directory.
* `cd ../..`
2. [Optional] Run the machine setup script if necessary. This script loads compiler modules (Fortran 2003-compliant Intel), netCDF module, etc. and sets compiler environment variables.
* `source Theia_setup` or
* `source Yellowstone_setup` or
* `source Cheyenne_setup`
3. Make a build directory and change into it.
* `mkdir bin && cd bin`
4. Invoke cmake on the source code to build.
* `cmake ../src`
5. Compile.
* `make`

## Running the SCM with CCPP
1. First append the library path to the physics suite libraries that you need within the CCPP. For example, for the dummy schemes, use
* `setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:[/path/to/build/directory/]ccpp/schemes/scm/src/scm-build`
2. Run the SCM with one of the supplied cases (twpice, arm_sgp_summer_1997, astex). The SCM will go through the time steps, applying forcing and calling the dummy physics.
* `./gmtb_scm twpice`
