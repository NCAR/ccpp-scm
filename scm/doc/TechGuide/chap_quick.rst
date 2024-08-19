.. _`chapter: quick`:

Quick Start Guide
=================

This chapter provides instructions for obtaining and compiling the CCPP
SCM. We provide instructions on building the code from scratch (:numref:`Section %s <obtaining_code>`), as well as
using Docker containers for machines that have Docker software installed (:numref:`Section %s <docker>`).

.. _obtaining_code:

Obtaining Code
--------------

The source code for the SCM, CCPP, and their required components are provided through GitHub.
The latest release branch contains the tested and supported version for
general use, while the development branch (``main``) contains the latest
developer code, but may not be as stable or consistent with existing documentation.
Instructions for using either option are discussed here.

Release Code
^^^^^^^^^^^^

Clone the source using

.. code:: bash

   git clone --recursive -b v7.0.0 https://github.com/NCAR/ccpp-scm

The ``--recursive`` option is required to retrieve the ccpp-physics and ccpp-framework code,
which are stored in separate repositories and linked to the SCM repository as submodules.
If not included initially, you can always retrieve the submodules
by executing the following command from the SCM directory:

.. code:: bash

   git submodule update --init --recursive

The CCPP framework can be found in the ``ccpp/framework`` subdirectory at
this level. The CCPP physics parameterizations can be found in the
``ccpp/physics`` subdirectory.

.. _`development_code`:

Development Code
^^^^^^^^^^^^^^^^

Developers seeking to contribute code to the SCM or CCPP will need to use the most up-to-date
version of the code, which can be found on the ``main`` branch of the repository:

.. code:: bash

   git clone --recursive -b main https://github.com/NCAR/ccpp-scm

Recall that the ``--recursive`` option in this command clones the main ccpp-scm
repository and all subrepositories (ccpp-physics and ccpp-framework).
If not included initially, you can always retrieve the commits of
the submodules that were intended to be used with a given commit of the
SCM by executing the following command from the SCM directory:

.. code:: bash

   git submodule update --init --recursive

While the ``main`` branch is tested regularly for compilation and basic functionality (as described in :numref:`Section %s <testing>`),
it may not be as stable or scientifically vetted as the latest release code, and may be lacking in up-to-date documentation.

If you would like to contribute as a developer to this project, please
see (in addition to the rest of this guide) the scientific and technical
documentation included with this release for both the SCM and the CCPP:

https://dtcenter.org/community-code/common-community-physics-package-ccpp/documentation

There you will find links to all of the documentation pertinent to
developers.


While the SCM is updated with the latest commits to the CCPP submodules (ccpp-physics and ccpp-framework)
on a fairly regular basis, it may be behind by a few commits at times. You can try to use the latest commits of the ccpp-physics and
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
relies on the NetCDF libraries, as well as the NCEP libraries ``bacio``, ``sp`` and ``w3emc``.

Beyond the standard shell scripts, the build
system relies on use of the Python scripting language, along with cmake,
GNU make and date.

For the latest release, the minimum required Python version is 3.10, and CMake requires a minimum version of 3.23.
While exact minimum required versions of other prerequisites have not been established, users can reference the
list of Continuous Integration tests run on the CCPP SCM repository (see :numref:`Section %s <continuous integration>`)
for examples of known working configurations.

Spack-stack
^^^^^^^^^^^^

This is a joint effort between NOAA's Unified Forecast System (UFS) and Joint Effort for Data assimilation Integration (JEDI).
It is designed to be a comprehensive, all-in-one package containing prerequisite libraries and tools needed for all
software in the UFS ecosystem, including the CCPP SCM. As of the version 7, installing spack-stack is the main
supported method of installing the prerequisites needed for building the SCM. The latest version of the SCM is meant
to be built with spack-stack v1.6.0. Older versions may work, but are not guaranteed. Version 1.6.0 of spack-stack
contains the following set of libraries needed for building the SCM:

 - Netcdf-c (v4.9.2)

 - Netcdf-FORTRAN (v4.6.1)

 - BACIO (v2.4.1) - Binary I/O Library

 - SP (v2.5.0) - Spectral Transformation Library

 - W3EMC (2.10.0) - GRIB decoder and encoder library

Instructions for installing spack-stack can be found in the `spack-stack documentation <https://spack-stack.readthedocs.io/en/latest/>`__.
Spack-stack is already installed and maintained on many HPC platforms, including NSF NCAR's Derecho, NOAA's Hera and
Jet, and MSU's Orion.

Compilers
^^^^^^^^^
The CCPP and SCM have been tested on a variety of computing platforms.
Currently the CCPP system is actively supported on Linux and MacOS
computing platforms using the Intel or GNU Fortran compilers. Windows
users have a path to use this software through a Docker container that
uses Linux internally (see :numref:`Section %s <docker>`). Typically the best chance of successfully building and
running the SCM on a new machine comes from using the most recent version of a compiler. If you have
problems with compilers, please check the “Known Issues” section of the
release website
(https://dtcenter.org/community-code/common-community-physics-package-ccpp/download).

.. _`use_preconfigured_platforms`:

Using Existing Libraries on Preconfigured Platforms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For users on supported platforms such as generic Linux or macOS systems
that have not been preconfigured, installing ``spack-stack`` (see :numref:`Section %s <spack-stack>`)
is highly recommended, as it provides all the necessary prerequisite libraries needed for installing the SCM.

The CCPP/SCM team does not support spack-stack, so users with questions or requiring help with spack-stack installation
should reference the `spack-stack documentation <https://spack-stack.readthedocs.io/en/latest/>`__.
However, we have provided an example procedure in
`this GitHub discussion <https://github.com/NCAR/ccpp-scm/discussions/464>`__.

The main downside to spack-stack is that it contains a large number of libraries and utilities used by the whole
Unified Forecast System and related applications, only a minority of which are required for the SCM. Users may
install libraries manually if they wish, but they will need to make sure the appropriate environment variables
are set to the correct values so that the build system can find them, as described in the following paragraphs.


<<<<<<< HEAD
Setting up compilation environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For users on a pre-configured platform, the spack-stack environment can be loaded via one of the provided modules in ``scm/etc/modules/`` as described in :numref:`Section %s <use_preconfigured_platforms>`.
=======
For users on a pre-configured platform, you can load the spack-stack environment via one of the provided modules in ``scm/etc/modules/``.
For example, users on the NSF NCAR machine Derecho who wish to use Intel compilers can do the following:

::

   cd [path/to/ccpp-scm/]
   module use scm/etc/modules/
   module load derecho_intel
>>>>>>> feature/modulefile_updates

Additionally, for users who have installed spack-stack on their own MacOS or Linux machine can use the provided ``macos_clang``
or ``linux_gnu`` modules.

.. note::

  The provided modules assume ``clang``/``gfortran`` compilers on MacOS and GNU compilers for Linux.
  If you are using a different set of compilers, you may need to modify the module file.

If libraries were installed manually, users will need to set some environment variables
needed for specifying the location of the various prerequisites. Users will need to set variables for the
compilers (``CC``, ``CXX``, ``FC``), as well as the root directories for the library installs of NetCDF (``NetCDF_ROOT``),
``bacio`` (``bacio_ROOT``), ``sp`` (``sp_ROOT``), and ``w3emc`` (``w3emc_ROOT``). This is the procedure used in the
provided Dockerfile in ``ccpp-scm/docker/``, so users can reference that file for guidance on how to install this software
and set these variables.

If libraries were installed via spack-stack, users can load modules similarly to those available on pre-configured platforms.
For a user on MacOS, who has installed spack-stack with ``clang``/``gfortran`` compilers, they can set up the build environment
by setting the SPACK_STACK_DIR variable to the appropriate path, and loading the module as on pre-configured platforms described above.

::

   export SPACK_STACK_DIR=[/path/to/spack-stack]
   cd [path/to/ccpp-scm/]
   module use scm/etc/modules/
   module load macos_clang

A module file is also provided for a generic linux platform with gnu compilers. For other platforms/combinations, you may be able
to modify the provided modulefiles to work with your spack-stack install, otherwise reference the above procedure for manually installed libraries.

Python requirements
"""""""""""""""""""""

The SCM build system invokes the ``ccpp_prebuild.py`` script, and so the Python environment must be set up prior to building.
As mentioned earlier, a minimum Python version of 3.10 is required. Additionally, there are a few non-default modules required for the SCM to
function: ``f90nml`` (`documentation <https://f90nml.readthedocs.io/en/latest/index.html>`__) and
``netcdf4`` (`documentation <https://unidata.github.io/netcdf4-python/>`__). Users can test if these are installed using this command in
the shell:

::

   python -c "import f90nml; import netcdf4"

If is installed, this command will succeed silently, otherwise an ``ImportError: No module named f90nml``
will be printed to screen. To install the ``f90nml`` (v1.4.4; ) and ``netcdf4`` (v1.6.5) Python modules, use the
install method preferred for your Python environment (one of the following):

-  ::

      easy_install f90nml==1.4.4 netcdf4==1.6.5

-  ::

      pip install f90nml==1.4.4 netcdf4==1.6.5

-  ::

      conda install -c conda-forge f90nml==1.4.4 netcdf4==1.6.5


.. _`compiling`:

Compiling SCM with CCPP
-----------------------

The first step in compiling the CCPP and SCM is to properly setup your
user environment as described in :numref:`Section %s <use_preconfigured_platforms>`
and :numref:`Section %s <setup_supported_platforms>`.

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

      -  Debug mode, which compiles with lower optimization and additional compile-time checks. Only
         recommended for development and debugging, because code compiled in this mode will run slower.

         .. code:: bash

            -DCMAKE_BUILD_TYPE=Debug

      -  Single Precision, lowers the default precision of variables from double to single precision.
         A very few calculations are done in double precision where it is crucial to achieve results comparable to the default double precision.

         .. code:: bash

            -D32BIT=ON

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

.. warning::
  This command can be dangerous (deletes files without confirming),
  so make sure that you’re in the right directory before executing!

If you encounter errors, please capture a log file from all of the
steps, and start a thread on the Github Discussions support forum at:
https://github.com/NCAR/ccpp-scm/discussions

Run the SCM with a supplied case
--------------------------------

There are several test cases provided with this version of the SCM. For
all cases, the SCM will go through the time steps, applying forcing and
calling the physics defined in the chosen suite definition file using
physics configuration options from an associated namelist. The model is
executed through a Python run script that is pre-staged into the ``bin``
directory: ``run_scm.py``. It can be used to run one integration or several
integrations serially, depending on the command line arguments supplied.

Downloading input data
^^^^^^^^^^^^^^^^^^^^^^
The various SCM cases require staged input data in order to run. This includes
input data for cases and lookup tables for runtime use. This is a large dataset
(:math:`<`\ 1 GB) so it is not stored in the SCM repository, and must be downloaded
separately. To download this data place it in the correct directories,
execute the following scripts:

.. code:: bash

   ./contrib/get_all_static_data.sh
   ./contrib/get_thompson_tables.sh

If the download step fails, make sure that your system’s firewall does
not block access to GitHub. If it does, download the files ``comparison_data.tar.gz``,
``physics_input_data.tar.gz``, ``processed_case_input.tar.gz``, and ``raw_case_input.tar.gz``
from the `SCM release page <https://github.com/NCAR/ccpp-scm/releases/tag/v7.0.0>`__ using your browser and manually extract its
contents in the directory ``scm/data``. Similarly, do the same for
``thompson_tables.tar.gz`` and ``MG_INCCN_data.tar.gz`` and extract
to ``scm/data/physics_input_data/``.

New with the SCM v7 release, static data is available for running cases with GOCART climatological aerosols (where the value of ``iaer`` in the ``&gfs_physics_nml`` namelist starts with 1; see the `CCPP Scientific Documentation <https://dtcenter.ucar.edu/GMTB/v7.0.0/sci_doc/_c_c_p_psuite_nml_desp.html>`__ for more information); one example of this is with the default namelist settings for the GFS_v17_p8_ugwpv1 scheme. This dataset is very large (~12 GB), so it is recommended only to download it if you will be using it.

.. code:: bash

   ./contrib/get_aerosol_climo.sh

.. _`singlerunscript`:

Run Script Usage
^^^^^^^^^^^^^^^^

Running a case requires four pieces of information: the case to run
(consisting of initial conditions, geolocation, forcing data, etc.), the
physics suite to use (through a CCPP suite definition file), a physics
namelist (that specifies configurable physics options to use), and a
tracer configuration file. As discussed in :numref:`Chapter %c <cases>`, cases are set up via their own
namelists in ``../etc/case_config``. A default physics suite is provided as a user-editable
variable in the script and default namelists and tracer configurations
are associated with each physics suite (through ``../src/suite_info.py``), so, technically, one
must only specify a case to run with the SCM when running just one
integration. For example, to run the "BOMEX" case:

.. code:: bash

  ./run_scm.py -c bomex

For running multiple integrations at once, the run script can accept a file that contains a list of tests to run.
The file ``ccpp-scm/test/rt_test_cases.py`` contains the full list of regression test cases, so you could run that list
of tests with the following command:

.. code:: bash

 ./run_scm.py -f ../../test/rt_test_cases.py

To see the full list of available options, use the ``--help`` flag:

.. code:: bash

  ./run_scm.py --help


The run script’s full set of options are described below, where optional abbreviations are included in brackets.
If using the main branch, you should run the above command to ensure you have the most up-to-date list of options.

-  ``--case [-c]``

   -  The provided argument should correspond to the name of a case in
      ``../etc/case_config`` (without the ``.nml`` extension).

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

-  ``--file [-f]``

   -  This option may be used to specify a list of tests to run; see ../../test/rt_test_cases.py for an example.

-  ``--gdb [-g]``

   -  Use this to run the executable through the ``gdb`` debugger (if it is
      installed on the system).

-  ``--docker [-d]``

   -  Use this argument when running in a docker container in order to
      successfully mount a volume between the host machine and the
      Docker container instance, allowing the container to share the output and plots with
      the host machine.

-  ``--runtime``

   -  Use this to override the runtime provided in the case
      configuration namelist.

-  ``--runtime_mult``

   -  Use this to override the runtime provided in the case
      configuration namelist by multiplying the runtime by the given
      value. This is used, for example, in regression testing to reduce
      total runtimes (e.g., ``--runtime_mult 0.1``).

-  ``--levels [-l]``

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
<<<<<<< HEAD
part of the CCPP SCM v7.0.0 release, the following suite names are supported:
=======
part of the seventh CCPP release, the following suite names are supported:
>>>>>>> feature/modulefile_updates

#. SCM_GFS_v16

#. SCM_GFS_v16_RRTMGP

#. SCM_GFS_v17_p8_ugwpv1

#. SCM_HRRR_gf

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
been created and use the ``*_ps`` decoration. It is not necessary to specify this
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

A NetCDF output file is generated in an output directory named
with the case and suite within the run directory. If using a Docker
container, all output is copied to the directory in container space for
volume-mounting purposes. Any standard NetCDF file viewing or analysis
tools may be used to examine the output file (ncdump, ncview, NCL, etc).

Batch Run Script
^^^^^^^^^^^^^^^^

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
https://ccpp-techdoc.readthedocs.io/en/v7.0.0/.

.. _docker:

Creating and Using a Docker Container with SCM and CCPP
-------------------------------------------------------

In order to run a precompiled version of the CCPP SCM in a container,
Docker will need to be available on your machine. Please visit
https://www.docker.com to download and install the version compatible
with your system. Docker frequently releases updates to the software; it
is recommended to apply all available updates.

.. note::
  In order to install Docker on your machine, you will be required to have root access
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
^^^^^^^^^^^^^^^^^^^^^^^^^

The Dockerfile builds CCPP SCM v7.0.0 from source using the GNU
compiler.

The CCPP SCM has a number of system requirements and necessary libraries
and tools. Below is a list, including versions, used to create the the
GNU-based Docker image. These are included for reference, but recall that
the Docker container contains all of this software built-in, you do not need to install them separately!

-  gfortran - 12.2.0

-  gcc - 12.2.0

-  cmake - 3.25.1

-  NetCDF - 4.9.0

-  Python - 3.11.2

-  NCEPLIBS BACIO - v2.4.1

-  NCEPLIBS SP - v2.3.3

-  NCEPLIBS W3EMC  - v2.11.0

A Docker image containing the SCM, CCPP, and its software prerequisites
can be generated from the code in the software repository obtained by
following the instructions in :numref:`Section %s <obtaining_code>`,
and then executing the following steps:

.. note::
  Windows users can execute these steps in the terminal application
  that was installed as part of Docker Toolbox.

#. Navigate to the ``ccpp-scm/docker`` directory.

#. Run the ``docker build`` command to generate the Docker image, using the supplied
   Dockerfile.

   .. code:: bash

      docker build -t ccpp-scm .

   Inspect the Dockerfile if you would like to see details for how the
   image is built. The image will contain SCM prerequisite software from
   DTC, the SCM and CCPP code, and a pre-compiled executable for the SCM
   with the 5 supported suites for the SCM. To view

   .. code:: bash

      > docker images

      REPOSITORY           TAG       IMAGE ID       CREATED       SIZE
      ccpp-scm             latest    1b2e0a0afdf9   2 days ago    3.21GB


Using a prebuilt Docker image from Dockerhub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A prebuilt Docker image for this release is available on Dockerhub if it
is not desired to build from source. In order to use this, execute the
following from the terminal where Docker is run:

.. code:: bash

   docker pull dtcenter/ccpp-scm:v7.0.0

To verify that it exists afterward, run

.. code:: bash

   docker images

Running the Docker image
^^^^^^^^^^^^^^^^^^^^^^^^

.. note::
  Windows users can execute these steps through the Docker
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

#. Set an environment variable to use for your SCM output directory.

   .. code:: bash

      export OUT_DIR=/path/to/output

   For Windows, the format that worked for us followed this example:
   ``/c/Users/myusername/path/to/directory/to/mount``

#. To run the SCM, you can run the Docker container that was just
   created and give it the same run commands as discussed in :numref:`Section %s <singlerunscript>`
   **Be sure to remember to include the ``-d`` and ``--mpi_command "mpirun -np 1 --allow-run-as-root"`` 
   options for all run commands**. For example,

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm ./run_scm.py -c twpice --mpi_command "mpirun -np 1 --allow-run-as-root" -d

   will run through the TWPICE case using the default suite and namelist
   and put the output in the shared directory.

   .. note::
     Windows users may need to omit the curly braces around environment variables: use ``$OUT_DIR``
     instead of ``${OUT_DIR}``.

   For running through all supported cases and suites, use

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm ./run_scm.py -f ../../test/rt_test_cases.py --runtime_mult 0.1 --mpi_command "mpirun -np 1 --allow-run-as-root" -d

   The options included in the above ``run`` commands are the following:

   -  ``−−rm`` removes the container when it exits

   -  ``-it`` interactive mode with terminal access

   -  ``-v`` specifies the volume mount from host directory (outside container)
      to inside the container. Using volumes allows you to share data
      between the host machine and container. For running the SCM, the
      output is being mounted from inside the container to the
      host machine. Upon exiting the container, data mounted to the host
      machine will still be accessible.

   -  ``−−name`` names the container. If no name is provided, the daemon will
      autogenerate a random string name.

   .. note::
     If you are using a prebuilt image from Dockerhub, substitute
     the name of the image that was pulled from Dockerhub in the commands
     above; i.e. instead of ``ccpp-scm`` above, one would have ``dtcenter/ccpp-scm:v7.0.0``.

#. To use the SCM interactively, run non-default configurations, create
   plots, or even develop code, issue the following command:

   .. code:: bash

      docker run --rm -it -v ${OUT_DIR}:/home --name run-ccpp-scm ccpp-scm /bin/bash

   You will be placed within the container space and within the
   directory of the SCM with a pre-compiled executable. At this point,
   one could use the run scripts as described in previous sections
   (remembering to include the option on run scripts if output is to be
   shared with the host machine).

   .. warning::

     If developing or modifying code, since the container is ephemeral, one should push their changes to a remote git
     repository to save them (i.e. a fork on GitHub.com).
