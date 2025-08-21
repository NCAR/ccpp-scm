.. _`chapter: repository`:

Repository
==========

What is included in the repository?
-----------------------------------

The repository contains all code required to build the CCPP SCM and
scripts that can be used to obtain data to run it (e.g. downloading
large initialization tables for the Thompson microphysics schemes
discussed in :numref:`Subsection %s <singlerunscript>` and
processed case data). It is functionally separated into 3 subdirectories
representing the SCM model infrastructure (``scm`` directory), the CCPP
infrastructure (``ccpp/framework`` directory), and the CCPP physics schemes
(``ccpp/physics`` directory). The entire ``ccpp-scm`` repository resides on
Github’s NCAR space, and the ``ccpp/framework`` and ``ccpp/physics`` directories
are git submodules that point to repositories ``ccpp-framework`` and ``ccpp-physics`` on the
same space. The structure of the entire repository is represented below.
Note that the ``ccpp-physics`` repository also contains files needed for using the CCPP
with the UFS Atmosphere host model that uses the Finite-Volume
Cubed-Sphere (FV3) dynamical core.

| ``ccpp-scm/``
| ``├── CMakeModules``
| ``├── CODEOWNERS`` - List of code maintainers/developers who are automatically assigned to review Pull Requests on GitHub
| ``├── LICENSE``
| ``├── README.md``
| ``├── ccpp``
| ``│   ├── config`` - Contains the CCPP prebuild configuration file
| ``│   ├── framework`` - Contains CCPP framework submodule. See https://github.com/NCAR/ccpp-framework for contents
| ``│   ├── physics`` - Contains CCPP physics submodule. See https://github.com/NCAR/ccpp-physics for contents
| ``│   ├── physics_namelists`` - Contains physics namelist files associated with suites
| ``│   └── suites`` - Contains suite definition files
| ``├── contrib``
| ``│   ├── get_all_static_data.sh`` - Script for downloading/extracting the processed SCM case data
| ``│   ├── get_mg_inccn_data.sh`` - Script for downloading/extracting the Morrison-Gettelman data
| ``│   └── get_thompson_tables.sh`` - Script for downloading/extracting the Thompson lookup tables
| ``│   └── get_aerosol_climo.sh`` - Script for downloading/extracting the GOCART climatological aerosol data
| ``├── docker``
| ``│   └── Dockerfile`` - Contains Docker instructions for building the CCPP SCM image
| ``├── environment-suite-sim.yml`` - Python environment dependency file for the CCPP Suite Simulator
| ``├── environment-ufscasegen.yml`` - Python environment dependency file for the UFS Case Generator capability
| ``├── environment.yml`` - Python environment dependency file for the SCM
| ``├── scm``
| ``│   ├── LICENSE.txt`` - Contains licensing information
| ``│   ├── data`` - Directory where data is staged by scripts in the ``ccpp/contrib/`` directory
| ``│   │   └── vert_coord_data`` - Contains data to calculate vertical coordinates (from GSM-based GFS only)
| ``│   ├── doc``
| ``│   │   └── TechGuide`` - Contains source code and other files for this User’s/Technical Guide
| ``│   ├── etc`` - Contains case configuration, machine setup scripts, and plotting scripts. These might be machine dependent and not all scripts are supported
| ``│   │   ├── CENTOS_docker_setup.sh`` - Contains machine setup for Docker container
| ``│   │   ├── case_config`` - Contains case configuration files
| ``│   │   ├── modules`` - Contains module files for loading build environments on both pre-configured and custom platforms
| ``│   │   ├── scm_qsub_example.py`` - Example ``qsub`` (LSF) run script
| ``│   │   ├── scm_slurm_example.py`` - Example ``srun`` (SLURM) run script
| ``│   │   ├── scripts`` - Python scripts for setting up cases, plotting, and the CCPP Suite Simulator
| ``│   │   │   ├── ccpp_suite_sim`` - Python scripts for the CCPP Suite Simulator
| ``│   │   │   ├── plot_configs`` - Plot configuration files
| ``│   │   └── tracer_config`` - Tracer configuration files
| ``│   └── src`` - Source code for SCM infrastructure, Python run script, CMakeLists.txt for the SCM, example multirun setup files, suite_info.py
| ``└── test`` - Contains scripts for regression testing, Continuous Integration tests

Testing
-----------------

Regression Testing
^^^^^^^^^^^^^^^^^^

Regression tests are a comprehensive set of build and run tests meant to ensure that new changes to the SCM do not break any existing capabilities. These tests are run on code changes before they are merged, and so ensure that the ``main`` branch is always free of major bugs in all facets of the system covered by the tests.

The latest set of Regression tests are run automatically for every new code change when a Pull Request is opened via GitHub's `Continuous Integration`_. Regression tests are also run manually on a wide variety of platforms in preparation for code release to ensure that all capabilities work as expected for a reasonable spectrum of possible machines a user might want to use.

Continuous Integration
^^^^^^^^^^^^^^^^^^^^^^

The CCPP SCM GitHub repository is set up with Continuous Integration (CI) testing for building the SCM and running some simple test cases. These tests are run automatically on code changes before they are merged, and so ensures that new changes to the SCM do not break basic capabilities. The latest set of tests use the following combinations of SCM prerequisites:

**Regression tests**
 - GNU compilers 11.4.0, Python 3.9.12, netCDF-c 4.7.3, netCDF-FORTRAN 4.5.3, bacio 2.4.1, sp 2.3.3, and w3emc 2.9.2

**Build tests**

All tests use the same versions of NCEP-supported libraries: bacio 2.4.1, sp 2.3.3, and w3emc 2.9.2. Detailed information on these tests can be found in the definition files for these tests, stored in the SCM repository under ``ccpp-scm/.github/workflows``.
