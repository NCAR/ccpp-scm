.. _`chapter: quick`:

Quick Start Guide
=================

This chapter provides instructions for obtaining and compiling the CCPP
SCM. The SCM code calls CCPP-compliant physics schemes through the CCPP
framework code. As such, it requires the CCPP framework code and physics
code, both of which are included as submodules within the SCM git repository. This
package can be considered a simple example for an atmospheric model to
interact with physics through the CCPP.

Alternatively, if one doesn’t have access to or care to set up a machine
with the appropriate system requirements but has a working Docker
installation, it is possible to create and use a Docker container with a
pre-configured computing environment with a pre-compiled model. This is
also an avenue for running this software with a Windows PC. See section
:numref:`Section %s <docker>` for more information.

.. _obtaining_code:

Obtaining Code
--------------

The source code for the SCM, CCPP, and their required components are provided through GitHub.
The latest release branch contains the tested and supported version for
general use, while the development branch (``main``) contains the latest
developer code, but may not be as stable or consistent with existing documentation. 
Instructions for using either option are discussed here.

Release Code
~~~~~~~~~~~~

Clone the source using

.. code:: bash

   git clone --recursive -b v6.0.0 https://github.com/NCAR/ccpp-scm

By using the ``--recursive`` option, it guarantees that you are checking out the commits
of ccpp-physics and ccpp-framework that were tested with the latest
commit of the SCM main branch. If not included initially, you can always retrieve the commits of
the submodules that were intended to be used with a given commit of the
SCM by executing the following command from the SCM directory:

.. code:: bash

   git submodule update --init --recursive

The CCPP framework can be found in the ``ccpp/framework`` subdirectory at
this level. The CCPP physics parameterizations can be found in the
``ccpp/physics`` subdirectory.

.. _`development_code`:

Development Code
~~~~~~~~~~~~~~~~

If you would like to contribute as a developer to this project, please
see (in addition to the rest of this guide) the scientific and technical
documentation included with this release for both the SCM and the CCPP:

https://dtcenter.org/community-code/common-community-physics-package-ccpp/documentation

There you will find links to all of the documentation pertinent to
developers.

For working with the development branches (stability not guaranteed),
check out the ``main`` branch of the repository:

.. code:: bash

   git clone --recursive -b main https://github.com/NCAR/ccpp-scm

Recall that the ``--recursive`` option in this command clones the main ccpp-scm
repository and all subrepositories (ccpp-physics and ccpp-framework).
If not included initially, you can always retrieve the commits of
the submodules that were intended to be used with a given commit of the
SCM by executing the following command from the SCM directory:

.. code:: bash

   git submodule update --init --recursive

You can try to use the latest commits of the ccpp-physics and
ccpp-framework submodules if you wish, but this may not have been tested
(i.e. SCM development may lag ccpp-physics and/or ccpp-framework
development). To do so:

#. Navigate to the ccpp-physics directory.

   .. code:: bash

      cd ccpp-scm/ccpp/physics

#. Check out main.

   .. code:: bash

      git checkout main

#. Pull down the latest changes just to be sure.

   .. code:: bash

      git pull

#. Do the same for ccpp-framework

   .. code:: bash

      cd ../framework
      git checkout main
      git pull

#. Change back to the main directory for following the instructions in
   :numref:`Section %s <compiling>`, assuming system requirements in
   section :numref:`Section %s <systemrequirements>` are met.

   .. code:: bash

      cd ../..

.. _`systemrequirements`:

System Requirements, Libraries, and Tools
-----------------------------------------

The source code for the SCM and CCPP components is in the form of
programs written in FORTRAN 90 (with some required features from the 
FORTRAN 2008 standard), and C. In addition, the model I/O
relies on the NetCDF libraries. Beyond the standard scripts, the build
system relies on use of the Python scripting language, along with cmake,
GNU make and date.

The following software stacks have been tested with this code. Other
versions of various components will likely still work, however.

-  gfortran 12.1.0, gcc 12.1.0, cmake 3.23.2, NetCDF 4.7.4, Python
   3.9.12

-  GNU compilers 10.1.0, cmake 3.16.4, NetCDF 4.8.1, Python 3.7.12

-  GNU compilers 11.1.0, cmake 3.18.2, NetCDF 4.8.1, Python 3.8.5

-  Intel compilers 2022.0.2, cmake 3.20.1, NetCDF 4.7.4, Python 3.7.11

-  Intel compilers 2022.1.0, cmake 3.22.0, NetCDF 4.8.1, Python 3.7.12

Because these tools are typically the purview of system administrators
to install and maintain, they are considered part of the basic system
requirements. The Unified Forecast System (UFS) Short-Range Weather
Application release v1.0.0 of March 2021 provides software packages and
detailed instructions to install these prerequisites and the hpc-stack
on supported platforms (see :numref:`Section %s <setup_supported_platforms>`)

Further, there are several utility libraries as part of the hpc-stack
package that must be installed with environment variables pointing to
their locations prior to building the SCM.

-  bacio - Binary I/O Library

-  sp - Spectral Transformation Library

-  w3emc - GRIB decoder and encoder library

The following environment variables are used by the build system to
properly link these libraries: ``bacio_ROOT``, ``sp_ROOT``, and ``w3emc_ROOT`` Computational platforms on
which these libraries are prebuilt and installed in a central location
are referred to as *preconfigured* platforms. Examples of preconfigured
platforms are most NOAA high-performance computing machines (using the
Intel compiler) and the NCAR Cheyenne system (using the Intel and GNU
compilers). The machine setup scripts mentioned in
:numref:`Section %s <compiling>` load these libraries (which are identical
to those used by the UFS Short and Medium Range Weather Applications on
those machines) and set these environment variables for the user
automatically. For installing the libraries and its prerequisites on
supported platforms, existing UFS packages can be used (see
:numref:`Section %s <setup_supported_platforms>`).

Compilers
~~~~~~~~~

The CCPP and SCM have been tested on a variety of computing platforms.
Currently the CCPP system is actively supported on Linux and MacOS
computing platforms using the Intel or GNU Fortran compilers. Windows
users have a path to use this software through a Docker container that
uses Linux internally (see section `1.5 <#docker>`__). Please use
compiler versions listed in the previous section as unforeseen build
issues may occur when using older versions. Typically the best results
come from using the most recent version of a compiler. If you have
problems with compilers, please check the “Known Issues” section of the
release website
(https://dtcenter.org/community-code/common-community-physics-package-ccpp/download).

.. _`use_preconfigured_platforms`:

Using Existing Libraries on Preconfigured Platforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the SCM can be built using the so-called `"spack-stack"
libraries <https://ufs-weather-model.readthedocs.io/en/latest/Glossary.html#term-spack-stack>`__
maintained for the UFS Weather Model effort, there are many platforms
where the SCM can be built using those existing libraries. This can be
done by loading provided modules in the directory (must be done from the
top-level "ccpp-scm" directory; otherwise the command should point to
the corresponding absolute path):

.. code:: sh

   module purge
   module use scm/etc/modules
   module load [machine]_[compiler]

View the contents of the directory to see if your machine/compiler
combination is supported. As of this writing, modulefiles are provided
for Intel and GNU compilers on the NCAR machine Derecho, the NOAA
machines Hera and Jet, and the NOAA/MSU machine Orion. Loading these
modules will set up all the needed compilers, libraries, and other
programs needed for building, as well as the python libraries needed for
both building and running the SCM.

.. _`setup_supported_platforms`:

Installing Libraries on Non-preconfigured Platforms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For users on supported platforms such as generic Linux or macOS systems
that have not been preconfigured, the project is suggested for
installing prerequisite libraries. Visit
https://github.com/NOAA-EMC/hpc-stack for instructions for installing
prerequisite libraries via *hpc-stack* in their docs directory. UFS users who
already installed libraries via the *hpc-stack* package only need to set the
compiler (``CC``, ``CXX``, ``FC``), NetCDF (``NetCDF_ROOT``), and ``bacio``,
``sp`` and ``w3emc`` (``bacio_ROOT``, ``sp_ROOT``, ``w3emc_ROOT``) environment variables to point
to their installation paths in order to compile the SCM.

The SCM uses only a small part of the UFS *hpc-stack* package and has fewer
prerequisites (i.e. no ESMF or wgrib2 needed). Users who are not planning to use the
UFS can install only NetCDF/NetCDF-Fortran manually or using the
software package manager (apt, yum, brew).

The Python environment must provide the module for the SCM scripts to
function. Users can test if f90nml is installed using this command in
the shell:

::

   python -c "import f90nml"

If is installed, this command will succeed silently, otherwise an ``ImportError: No module named f90nml``
will be printed to screen. To install the ``f90nml`` (v0.19) Python module, use the
install method preferred for your Python environment (one of the
following):

-  ::

      easy_install f90nml==0.19

-  ::

      pip install f90nml==0.19

-  ::

      conda install f90nml=0.19

or perform the following steps to install it manually from source:

::

   cd /directory/with/write/priveleges
   git clone -b v0.19 https://github.com/marshallward/f90nml
   cd f90nml
   python setup.py install [--prefix=/my/install/directory or --user]

The directory ``/my/install/directory`` must exist and its subdirectory
``/my/install/directory/lib/python[version]/site-packages`` (or ``lib64``
instead of ``lib``, depending on the system) must be in the ``PYTHONPATH``
environment variable.

.. _`compiling`:

Compiling SCM with CCPP
-----------------------

The first step in compiling the CCPP and SCM is to properly setup your
user environment as described in
sections :numref:`%s <use_preconfigured_platforms>` and :numref:`Section %s <setup_supported_platforms>`. The second step is
to download the lookup tables and other large datasets (large binaries,
:math:`<`\ 1 GB) needed by the physics schemes and place them in the
correct directory: From the top-level code directory (``ccpp-scm`` by default),
execute the following scripts:

.. code:: bash

   ./contrib/get_all_static_data.sh
   ./contrib/get_thompson_tables.sh

If the download step fails, make sure that your system’s firewall does
not block access to GitHub. If it does, download the files ``comparison_data.tar.gz``,
``physics_input_data.tar.gz``, ``processed_case_input.tar.gz``, and ``raw_case_input.tar.gz``
from the GitHub release website using your browser and manually extract its
contents in the directory ``scm/data``. Similarly, do the same for 
``thompson_tables.tar.gz`` and ``MG_INCCN_data.tar.gz`` and extract
to ``scm/data/physics_input_data/``.

Following this step, the top level build system will use ``cmake`` to query system
parameters, execute the CCPP prebuild script to match the physics
variables (between what the host model – SCM – can provide and what is
needed by physics schemes in the CCPP for the chosen suites), and build
the physics caps needed to use them. Finally, ``make`` is used to compile the
components.

#. From the top-level code directory (``ccpp-scm`` by default), change directory to
   the top-level SCM directory.

   .. code:: bash

      cd scm

#. Make a build directory and change into it.

   .. code:: bash

      mkdir bin && cd bin

#. Invoke ``cmake`` on the source code to build using one of the options below.
   This step is used to identify for which suites the ccpp-framework
   will build caps and which suites can be run in the SCM without
   recompiling.

   -  Default mode

      .. code:: bash

         cmake ../src

      By default, this option uses all supported suites. The list of
      supported suites is controlled by ``scm/src/suite_info.py``.

   -  All suites mode

      .. code:: bash

         cmake -DCCPP_SUITES=ALL ../src

      All suites in ``scm/src/suite_info.py``, regardless of whether they’re supported, will be
      used. This list is typically longer for the development version of
      the code than for releases.

   -  Selected suites mode

      .. code:: bash

         cmake -DCCPP_SUITES=SCM_GFS_v16,SCM_RAP ../src

      This only compiles the listed subset of suites (which should still
      have a corresponding entry in ``scm/src/suite_info.py``)

   -  The statements above can be modified with the following options
      (put before ``../src``):

      -  Use threading with openmp (not for macOS with clang+gfortran)

         .. code:: bash

            -DOPENMP=ON

      -  Debug mode

         .. code:: bash

            -DCMAKE_BUILD_TYPE=Debug

   -  One can also save the output of this step to a log file:

      .. code:: bash

         cmake [-DCMAKE_BUILD_TYPE ...] ../src 2>&1 | tee log.cmake

   CMake automatically runs the CCPP prebuild script to match required
   physics variables with those available from the dycore (SCM) and to
   generate physics caps and makefile segments. It generates software
   caps for each physics group defined in the supplied Suite Definition
   Files (SDFs) and generates a static library that becomes part of the
   SCM executable.

   If necessary, the CCPP prebuild script can be executed manually from
   the top level directory (``ccpp-scm``). The basic syntax is

   .. code:: bash

      ./ccpp/framework/scripts/ccpp_prebuild.py --config=./ccpp/config/ccpp_prebuild_config.py --suites=SCM_GFS_v16,SCM_RAP[...] --builddir=./scm/bin [--debug]

   where the argument supplied via the ``--suites`` variable is a comma-separated
   list of suite names that exist in the directory. Note that suite
   names are the suite definition filenames minus the ``suite_`` prefix and ``.xml`` suffix.

#. Compile. Add ``VERBOSE=1`` to obtain more information on the build process.

   .. code:: bash

      make

   -  One may also use more threads for compilation and/or save the
      output of the compilation to a log file:

      .. code:: bash

         make -j4 2>&1 | tee log.make

The resulting executable may be found at ./scm (Full path of ``ccpp-scm/scm/bin/scm``).

Although ``make clean`` is not currently implemented, an out-of-source build is used,
so all that is required to clean the build directory is (from the ``bin``
directory)

.. code:: bash

   pwd #confirm that you are in the ccpp-scm/scm/bin directory before deleting files
   rm -rfd *

Note: This command can be dangerous (deletes files without confirming),
so make sure that you’re in the right directory before executing!

If you encounter errors, please capture a log file from all of the
steps, and start a thread on the support forum at:
https://dtcenter.org/forum/ccpp-user-support/ccpp-single-column-model

Run the SCM with a supplied case
--------------------------------

There are several test cases provided with this version of the SCM. For
all cases, the SCM will go through the time steps, applying forcing and
calling the physics defined in the chosen suite definition file using
physics configuration options from an associated namelist. The model is
executed through a Python run script that is pre-staged into the ``bin``
directory: ``run_scm.py``. It can be used to run one integration or several
integrations serially, depending on the command line arguments supplied.

.. _`singlerunscript`:

Run Script Usage
~~~~~~~~~~~~~~~~

Running a case requires four pieces of information: the case to run
(consisting of initial conditions, geolocation, forcing data, etc.), the
physics suite to use (through a CCPP suite definition file), a physics
namelist (that specifies configurable physics options to use), and a
tracer configuration file. As discussed in :numref:`Chapter %c <cases>`, cases are set up via their own
namelists in ``../etc/case_config``. A default physics suite is provided as a user-editable
variable in the script and default namelists and tracer configurations
are associated with each physics suite (through ``../src/suite_info.py``), so, technically, one
must only specify a case to run with the SCM when running just one
integration. For running multiple integrations at once, one need only
specify one argument (``-m``) which runs through all permutations of supported
suites from ``../src/suite_info.py`` and cases from ``../src/supported_cases.py``. The run script’s options are described
below where option abbreviations are included in brackets.

-  ``--case [-c]``

   -  **This or the ``--multirun`` option are the minimum required arguments.** The
      case should correspond to the name of a case in ``../etc/case_config`` (without the
      ``.nml`` extension).

-  ``--suite [-s]``

   -  The suite should correspond to the name of a suite in ``../ccpp/suites`` (without the
      ``.xml`` extension) that was supplied in the ``cmake`` or ``ccpp_prebuild`` step.

-  ``--namelist [-n]``

   -  The namelist should correspond to the name of a file in ``../ccpp/physics_namelists`` (WITH the
      ``.txt`` extension). If this argument is omitted, the default namelist for
      the given suite in ``../src/suite_info.py`` will be used.

-  ``--tracers [-t]``

   -  The tracers file should correspond to the name of a file in ``../etc/tracer_config`` (WITH
      the ``.txt`` extension). If this argument is omitted, the default tracer
      configuration for the given suite in ``../src/suite_info.py`` will be used.

-  ``--multirun [-m]``

   -  **This or the ``--case`` option are the minimum required arguments.** When
      used alone, this option runs through all permutations of supported
      suites from ``../src/suite_info.py`` and cases from ``../src/supported_cases.py``. When used in conjunction with the
      ``--file`` option, only the runs configured in the file will be run.

-  ``--file [-f]``

   -  This option may be used in conjunction with the ``--multirun`` argument. It
      specifies a path and filename to a python file where multiple runs
      are configured.

-  ``--gdb [-g]``

   -  Use this to run the executable through the ``gdb`` debugger (if it is
      installed on the system).

-  ``--docker [-d]``

   -  Use this argument when running in a docker container in order to
      successfully mount a volume between the host machine and the
      Docker container instance and to share the output and plots with
      the host machine.

-  ``--runtime``

   -  Use this to override the runtime provided in the case
      configuration namelist.

-  ``--runtime_mult``

   -  Use this to override the runtime provided in the case
      configuration namelist by multiplying the runtime by the given
      value. This is used, for example, in regression testing to reduce
      total runtimes.

-  ``--levels [-l]

   -  Use this to change the number of vertical levels.

-  ``--npz_type``

   -  Use this to change the type of FV3 vertical grid to produce (see
      ``src/scm_vgrid.F90`` for valid values).

-  ``--vert_coord_file``

   -  Use this to specify the path/filename of a file containing the a_k
      and b_k coefficients for the vertical grid generation code to use.

-  ``--bin_dir``

   -  Use this to specify the path to the build directory.

-  ``--run_dir``

   -  Use this to specify the path to the run directory.

-  ``--case_data_dir``

   -  Use this to specify the path to the directory containing the case
      data file (useful for using the DEPHY case repository).

-  ``--n_itt_out``

   -  Use this to specify the period of writing instantaneous output in
      timesteps (if different than the default specified in the script).

-  ``--n_itt_diag``

   -  Use this to specify the period of writing instantaneous and
      time-averaged diagnostic output in timesteps (if different than
      the default specified in the script).

-  ``--timestep [-dt]``

   -  Use this to specify the timestep to use (if different than the
      default specified in ``../src/suite_info.py``).

-  ``--verbose [-v]``

   -  Use this option to see additional debugging output from the run
      script and screen output from the executable.

When invoking the run script, the only required argument is the name of
the case to run. The case name used must match one of the case
configuration files located in ``../etc/case_config`` (*without the .nml extension!*). If
specifying a suite other than the default, the suite name used must
match the value of the suite name in one of the suite definition files
located in ``../../ccpp/suites`` (Note: not the filename of the suite definition file). As
part of the sixth CCPP release, the following suite names are valid:

#. SCM_GFS_v16

#. SCM_GFS_v17p8

#. SCM_RAP

#. SCM_HRRR

#. SCM_RRFS_v1beta

#. SCM_WoFS_v0

Note that using the Thompson microphysics scheme requires the
computation of look-up tables during its initialization phase. As of the
release, this process has been prohibitively slow with this model, so it
is HIGHLY suggested that these look-up tables are downloaded and staged
to use this scheme as described in :numref:`Section %s <compiling>`. The issue appears to be
machine/compiler-specific, so you may be able to produce the tables with
the SCM, especially when invoking ``cmake`` with the ``-DOPENMP=ON`` option.

Also note that some cases require specified surface fluxes. Special
suite definition files that correspond to the suites listed above have
been created and use the ``*_prescribed_surface`` decoration. It is not necessary to specify this
filename decoration when specifying the suite name. If the ``spec_sfc_flux`` variable in
the configuration file of the case being run is set to ``.true.``, the run script
will automatically use the special suite definition file that
corresponds to the chosen suite from the list above.

If specifying a namelist other than the default, the value must be an
entire filename that exists in ``../../ccpp/physics_namelists``. Caution should be exercised when
modifying physics namelists since some redundancy between flags to
control some physics parameterizations and scheme entries in the CCPP
suite definition files currently exists. Values of numerical parameters
are typically OK to change without fear of inconsistencies. If
specifying a tracer configuration other than the default, the value must
be an entire filename that exists in ``../../scm/etc/tracer_config``. The tracers that are used should
match what the physics suite expects, lest a runtime error will result.
Most of the tracers are dependent on the microphysics scheme used within
the suite. The tracer names that are supported as of this release are
given by the following list. Note that running without ``sphum``, ``o3mr``, and ``liq_wat`` and may
result in a runtime error in all supported suites.

#. sphum

#. o3mr

#. liq_wat

#. ice_wat

#. rainwat

#. snowwat

#. graupel

#. hailwat

#. cld_amt

#. water_nc

#. ice_nc

#. rain_nc

#. snow_nc

#. graupel_nc

#. hail_nc

#. graupel_vol

#. hail_vol

#. ccn_nc

#. sgs_tke

#. liq_aero

#. ice_aero

#. q_rimef

A NetCDF output file is generated in an output directory located named
with the case and suite within the run directory. If using a Docker
container, all output is copied to the directory in container space for
volume-mounting purposes. Any standard NetCDF file viewing or analysis
tools may be used to examine the output file (ncdump, ncview, NCL, etc).

Batch Run Script
~~~~~~~~~~~~~~~~

If using the model on HPC resources and significant amounts of processor
time is anticipated for the experiments, it will likely be necessary to
submit a job through the HPC’s batch system. An example script has been
included in the repository for running the model on Hera’s batch system
(SLURM). It is located in ``ccpp-scm/scm/etc/scm_slurm_example.py``. Edit the ``job_name``, ``account``, etc. to suit your needs and
copy to the ``bin`` directory. The case name to be run is included in the ``command``
variable. To use, invoke

.. code:: bash

   ./scm_slurm_example.py

from the ``bin`` directory.

Additional details regarding the SCM may be found in the remainder of
this guide. More information on the CCPP can be found in the CCPP
Technical Documentation available at
https://ccpp-techdoc.readthedocs.io/en/v6.0.0/.

.. _docker:

Creating and Using a Docker Container with SCM and CCPP
-------------------------------------------------------

In order to run a precompiled version of the CCPP SCM in a container,
Docker will need to be available on your machine. Please visit
https://www.docker.com to download and install the version compatible
with your system. Docker frequently releases updates to the software; it
is recommended to apply all available updates. NOTE: In order to install
Docker on your machine, you will be required to have root access
privileges. More information about getting started can be found at
https://docs.docker.com/get-started

The following tips were acquired during a recent installation of Docker
on a machine with Windows 10 Home Edition. Further help should be
obtained from your system administrator or, lacking other resources, an
internet search.

-  Windows 10 Home Edition does not support Docker Desktop due to lack
   of “Hyper-V” support, but does work with Docker Toolbox. See the
   installation guide
   (https://docs.docker.com/toolbox/toolbox_install_windows/).

-  You may need to turn on your CPU’s hardware virtualization capability
   through your system’s BIOS.

-  After a successful installation of Docker Toolbox, starting with
   Docker Quickstart may result in the following error even with
   virtualization correctly enabled: ``This computer doesn’t have VT-X/AMD-v enabled. Enabling it in the BIOS is mandatory.``
   We were able to bypass this error
   by opening a bash terminal installed with Docker Toolbox, navigating
   to the directory where it was installed, and executing the following
   command:

   .. code:: bash

      docker-machine create default --virtualbox-no-vtx-check

Building the Docker image
~~~~~~~~~~~~~~~~~~~~~~~~~

The Dockerfile builds CCPP SCM v6.0.0 from source using the GNU
compiler. A number of required codes are built and installed via the
DTC-supported common community container. For reference, the common
community container repository can be accessed here:
https://github.com/NCAR/Common-Community-Container.

The CCPP SCM has a number of system requirements and necessary libraries
and tools. Below is a list, including versions, used to create the the
GNU-based Docker image:

-  gfortran - 9.3

-  gcc - 9.3

-  cmake - 3.16.5

-  NetCDF - 4.6.2

-  HDF5 - 1.10.4

-  ZLIB - 1.2.7

-  SZIP - 2.1.1

-  Python - 3

-  NCEPLIBS subset: bacio v2.4.1_4, sp v2.3.3_d, w3emc v2.9.2_d

A Docker image containing the SCM, CCPP, and its software prerequisites
can be generated from the code in the software repository obtained by
following the instructions in :numref:`Section %s <obtaining_code>`,
and then executing the following steps:

NOTE: Windows users can execute these steps in the terminal application
that was installed as part of Docker Toolbox.

#. Navigate to the ``ccpp-scm/docker`` directory.

#. Run the ``docker build`` command to generate the Docker image, using the supplied
   Dockerfile.

   .. code:: bash

      docker build -t ccpp-scm .

   Inspect the Dockerfile if you would like to see details for how the
   image is built. The image will contain SCM prerequisite software from
   DTC, the SCM and CCPP code, and a pre-compiled executable for the SCM
   with the 6 supported suites for the SCM. A successful build will show
   two images: dtcenter/common-community-container, and ccpp-scm. To
   list images, type:

   .. code:: bash

      docker images

Using a prebuilt Docker image from Dockerhub
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A prebuilt Docker image for this release is available on Dockerhub if it
is not desired to build from source. In order to use this, execute the
following from the terminal where Docker is run:

.. code:: bash

   docker pull dtcenter/ccpp-scm:v6.0.0

To verify that it exists afterward, run

.. code:: bash

   docker images

Running the Docker image
~~~~~~~~~~~~~~~~~~~~~~~~

NOTE: Windows users can execute these steps through the Docker
Quickstart application installed with Docker Toolbox.

#. Set up a directory that will be shared between the host machine and
   the Docker container. When set up correctly, it will contain output
   generated by the SCM within the container for manipulation by the
   host machine. For Mac/Linux,

   .. code:: bash

      mkdir -p /path/to/output

   For Windows, you can try to create a directory of your choice to
   mount to the container, but it may not work or require more
   configuration, depending on your particular Docker installation. We
   have found that Docker volume mounting in Windows can be difficult to
   set up correctly. One method that worked for us was to create a new
   directory under our local user space, and specifying the volume mount
   as below. In addition, with Docker Toolbox, double check that the
   mounted directory has correct permissions. For example, open
   VirtualBox, right click on the running virtual machine, and choose
   “Settings”. In the dialog that appears, make sure that the directory
   you’re trying to share shows up in “Shared Folders" (and add it if it
   does not) and make sure that the “auto-mount" and “permanent" options
   are checked.

#. Set an environment variable to use for your SCM output directory. For
   *t/csh* shells,

   .. code:: bash

      setenv OUT_DIR /path/to/output

   For bourne/bash shells,

   .. code:: bash

      export OUT_DIR=/path/to/output

   For Windows, the format that worked for us followed this example:
   ``/c/Users/myusername/path/to/directory/to/mount``

#. To run the SCM, you can run the Docker container that was just
   created and give it the same run commands as discussed in :numref:`Section %s <singlerunscript>`
   **Be sure to remember to include the ``-d``
   include the option for all run commands**. For example,

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm ./run_scm.py -c twpice -d

   will run through the TWPICE case using the default suite and namelist
   and put the output in the shared directory. NOTE: Windows users may
   need to omit the curly braces around environment variables: use ``$OUT_DIR``
   instead of ``${OUT_DIR}``. For running through all supported cases and suites, use

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm ./run_scm.py -m -d

   The options included in the above ``run`` commands are the following:

   -  ``−−rm`` removes the container when it exits

   -  ``-it`` interactive mode with terminal access

   -  ``-v`` specifies the volume mount from host directory (outside container)
      to inside the container. Using volumes allows you to share data
      between the host machine and container. For running the SCM, the
      output is being mounted from inside the container to the on the
      host machine. Upon exiting the container, data mounted to the host
      machine will still be accessible.

   -  ``−−name`` names the container. If no name is provided, the daemon will
      autogenerate a random string name.

   NOTE: If you are using a prebuilt image from Dockerhub, substitute
   the name of the image that was pulled from Dockerhub in the commands
   above; i.e. instead of ``ccpp-scm`` above, one would have ``dtcenter/ccpp-scm:v6.0.0``.

#. To use the SCM interactively, run non-default configurations, create
   plots, or even develop code, issue the following command:

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm /bin/bash

   You will be placed within the container space and within the
   directory of the SCM with a pre-compiled executable. At this point,
   one could use the run scripts as described in previous sections
   (remembering to include the option on run scripts if output is to be
   shared with the host machine). NOTE: If developing, since the
   container is ephemeral, one should push their changes to a remote git
   repository to save them (i.e. a fork on GitHub.com).
