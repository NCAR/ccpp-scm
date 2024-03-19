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
| ``├── CODEOWNERS`` - list of code maintainers/developers who are automatically assigned to review Pull Requests on GitHub
| ``├── LICENSE``
| ``├── README.md``
| ``├── ccpp`` - contains the CCPP prebuild configuration file
| ``│   ├── config``
| ``│   ├── framework`` - Contains CCPP framework submodule. See https://github.com/NCAR/ccpp-framework for contents
| ``│   ├── physics`` - Contains CCPP physics submodule. See https://github.com/NCAR/ccpp-physics for contents
| ``│   ├── physics_namelists`` - contains physics namelist files associated with suites
| ``│   └── suites`` - contains suite definition files
| ``├── contrib``
| ``│   ├── get_all_static_data.sh`` - script for downloading/extracting the processed SCM case data
| ``│   ├── get_mg_inccn_data.sh`` - script for downloading/extracting the Morrison-Gettelman data
| ``│   └── get_thompson_tables.sh`` - script for downloading/extracting the Thompson lookup tables
| ``├── docker``
| ``│   └── Dockerfile`` - contains Docker instructions for building the CCPP SCM image
| ``├── environment-suite-sim.yml`` - Python environment dependency file for the CCPP Suite Simulator
| ``├── environment-ufsreplay.yml`` - Python environment dependency file for the UFS Replay capability
| ``├── environment.yml`` - Python environment dependency file for the SCM
| ``├── scm``
| ``│   ├── LICENSE.txt`` - Contains licensing information
| ``│   ├── data`` - Directory where data is staged by scripts in the ``ccpp/contrib/`` directory
| ``│   │   └── vert_coord_data`` - contains data to calculate vertical coordinates (from GSM-based GFS only)
| ``│   ├── doc``
| ``│   │   └── TechGuide`` - Contains source code and other files for this User’s/Technical Guide
| ``│   ├── etc`` - contains case configuration, machine setup scripts, and plotting scripts
| ``│   │   ├── CENTOS_docker_setup.sh`` - contains machine setup for Docker container
| ``│   │   ├── case_config`` - contains case configuration files
| ``│   │   ├── modules`` - contains module files for loading build environments on both pre-configured and custom platforms
| ``│   │   ├── scm_qsub_example.py`` - example ``qsub`` (LSF) run script
| ``│   │   ├── scm_slurm_example.py`` - example ``srun`` (SLURM) run script
| ``│   │   ├── scripts`` - Python scripts for setting up cases, plotting, and the CCPP Suite Simulator
| ``│   │   │   ├── ccpp_suite_sim`` - Python scripts for the CCPP Suite Simulator
| ``│   │   │   ├── plot_configs`` - plot configuration files
| ``│   │   └── tracer_config`` - tracer configuration files
| ``│   └── src`` - source code for SCM infrastructure, Python run script, CMakeLists.txt for the SCM, example multirun setup files, suite_info.py
| ``└── test`` - Contains scripts for regression testing, Continuous Integration tests

