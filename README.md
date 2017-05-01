# User's Guide

This code is a testing platform for the Common Community Physics Package (CCPP) and its
driver. As of v0.1.0, this code simply calls a "do-nothing" physics suite with
"do-nothing" physics schemes through the CCPP driver within the SCM framework. It
represents a proof-of-concept and a platform on which to add proper physics in the coming months.

The CCPP and its driver are separated into its [own repository](https://github.com/NCAR/gmtb-ccpp) and are pulled in to this SCM
code as a git submodule. The instructions below explain how to obtain the SCM code, bring
in the CCPP as a submodule, build the combined code into an application, and run a simple case.
Using the dummy physics, when the SCM runs, it performs all functions as if it is running
a real case (including time-stepping), but instead of calling real physics, it calls the
dummy physics, which simply print a confirmation message when called.

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

## Setting up the physics suite
First, a physics suite is defined using an XML file located in src/ccpp/examples. The XML file contains the suite name, the number of suite "parts" (suite parts exist so that an atmosphere "cap" can execute code between calls to the physics driver), the number of subcycles for each scheme (if physics schemes require smaller time steps than the dynamics), and the scheme names. Scheme names found in the suite XML files must correspond to schemes located within the ccpp/schemes directory. Schemes within this directory are compiled as their own libraries, whose names must also be specified when using a scheme.

Using a suite in the SCM framework involves specifying its name in the case configuration file to be used. For example, for the twpice case (case_config/twpice.nml), the variable 'physics suite' is set to the desired suite name. NOTE: As mentioned in the 'Running' section above, since the schemes are in their own libraries, you must specify the path to the compiled scheme libraries that are being used in the suite by appending the LD_LIBRARY_PATH (or DYLD_LIBRARY_PATH for Mac OS). Without this step, the scheme libraries will not be found at runtime. This step is likely temporary and will be handled automatically at build time in future versions.

## CCPP Implementation Notes

Interaction with the CCPP is accomplished in two source files: gmtb_scm.f90 and gmtb_scm_time_integration.f90. Within gmtb_scm.f90, before the time integration loop, the physics suite is initialized by calling 'ccpp_init' with the path to the suite XML file and the ccpp_t data type to be filled in. Immediately following, the 'ccpp_fields_add' subroutine is called to fill in the ccpp_t data type with the model data that is used (input, output, diagnostics, etc.) in the physics suite. Once the ccpp_t data type is filled, the physics suite can be called with the 'ccpp_ipd_run' routine, passing in the scheme to be called and the ccpp_t data type. For the SCM, the first time step is handled after initialization. To advance in the time loop, the 'do_time_step' subroutine is called from gmtb_scm_time_integration.f90. This subroutine handles applying SCM forcing and calling the physics. For testing purposes, the parts, subcycle loops, and schemes within the physics suite are looped through, although calls to ccpp_ipd_run can be split apart as necessary with code executed between calls. The 'atmosphere' cap has the control over how to call the schemes.
