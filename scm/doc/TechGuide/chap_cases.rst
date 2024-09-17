Cases
=====

How to run cases
----------------

Only two files are needed to set up and run a case with the SCM. The
first is a configuration namelist file found in ``ccpp-scm/scm/etc/case_config`` that contains parameters
for the SCM infrastructure. The second necessary file is a NetCDF file
containing data to initialize the column state and time-dependent data
to force the column state. The two files are described below.

.. _`case config`:

Case configuration namelist parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``case_config`` namelist expects the following parameters:

-  ``case_name``

   -  Identifier for which dataset (initialization and forcing) to load.
      This string must correspond to a dataset included in the directory
      ``ccpp-scm/scm/data/processed_case_input/`` (without the file extension).

-  ``sfc_roughness_length_cm``

   -  Surface roughness length in cm for calculating surface-related
      fields from specified surface fluxes (only used if surface fluxes
      are specified).

-  ``reference_profile_choice``

   -  An integer representing the choice of reference profile to use
      above the supplied initialization and forcing data (1 :math:`=`
      “McClatchey” profile, 2 :math:`=` mid-latitude summer standard
      atmosphere)

-  ``column_area``

   -  A list of floating point values representing the characteristic
      horizontal domain area of each atmospheric column in square meters
      (this could be analogous to a 3D model’s horizontal grid size or
      the characteristic horizontal scale of an observation array; these
      values are used in scale-aware schemes; if using multiple columns,
      you may specify an equal number of column areas)

-  ``model_ics``

   -  A boolean set to ``.true.`` if UFS atmosphere initial conditions are used
      rather than field campaign-based initial conditions

-  ``C_RES``

   -  An integer representing the grid size of the UFS atmosphere
      initial conditions; the integer represents the number of grid
      points in each horizontal direction of each cube tile

-  ``input_type``

   -  1 => DEPHY-SCM format.

Optional variables (that may be overridden via run script command line
arguments) are:

-  ``npz_type``

   - Changes the type of FV3 vertical grid to produce (see src/scm_vgrid.F90 for
     valid values), default=''.

-  ``vert_coord_file``

   -  File containing FV3 vertical grid coefficients, default=''.

-  ``n_levels``

   -  Specify the integer number of vertical levels, default=127.

-  ``dt``

   - Specify the timestep to use (if different than the default specified in
     ../../src/suite_info.py), default=600.

-  ``do_spinup``

   - Set to ``.true.`` when allowing the model to spin up before the "official"
     model integration starts. 

-  ``spinup_timesteps``

   - Number of timesteps to spin up when ``do_spinup`` is true

-  ``lsm_ics``

   - Set to ``.true.`` when LSM initial conditions are included (but not all ICs from
     another model)

.. _`case input`:

Case input data file (DEPHY format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Development and Evaluation of Physics in atmospheric models (DEPHY)
format is an internationally-adopted data format intended for use by SCM
and LESs. The initialization and forcing data for each case in the CCPP SCM
repository is stored in a NetCDF (version 4) file. Additional cases in DEPHY
format, not maintained by the DTC, can be cloned from the DEPHY-SCM repository,
and run by providing the DEPHY-SCM file location to the SCM. For example:

.. code:: bash

   cd [...]/ccpp-scm/scm/data
   git clone https://github.com/GdR-DEPHY/DEPHY-SCM DEPHY-SCM
   cd [...]/ccpp-scm/scm/bin
   ./run_scm.py -c MAGIC_LEG04A --case_data_dir [...]/ccpp-scm/scm/data/DEPHY-SCM/MAGIC/LEG04A -v

Each DEPHY file has three dimensions (``time``, ``t0``, ``levels``) and contains the initial
conditions (``t0``, ``levels``) and forcing data (``time``, ``levels``). Not all fields
are required for all cases. For example, the fields ``hfss`` and ``hfls`` are only needed if the
global attributes ``surface_forcing_temp`` or ``surface_forcing_moisture`` are set to ``surface_flux``
and state nudging variables are only required if the ``nudging_*`` terms in the global attributes
are turned on. Using an active LSM (Noah, NoahMP, RUC) requires many more variables than are listed
here. Example files for using with Noah and NoahMP LSMs are included in
ccpp-scm/scm/data/processed_case_input/fv3_model_point_noah[mp].nc. More information on the DEPHY
format requirements can be found at `DEPHY <https://github.com/GdR-DEPHY/DEPHY-SCM>`__.

.. _`case input dephy`:
.. literalinclude:: dephy_case_header.txt
    :name: lst_case_input_netcdf_header_dephy
    :caption: example NetCDF file (DEPHY format) header for case initialization and forcing data

Included Cases
--------------

Several cases are included in the repository to serve as examples for
users to create their own and for basic research. All case configuration
namelist files for included cases can be found in ``ccpp-scm/scm/etc/case_config`` and represent the
following observational field campaigns:

-  Tropical Warm Pool – International Cloud Experiment (TWP-ICE)
   maritime deep convection

-  Atmospheric Radiation Measurement (ARM) Southern Great Plains (SGP)
   Summer 1997 continental deep convection

-  Atlantic Stratocumulus Transition EXperiment (ASTEX) maritime
   stratocumulus-to-cumulus transition

-  Barbados Oceanographic and Meteorological EXperiment (BOMEX) maritime
   shallow convection

-  Large eddy simulation ARM Symbiotic Simulation and Observation
   (LASSO) for May 18, 2016 (with capability to run all LASSO dates -
   see :numref:`Section %s <lasso>`) continental shallow convection

-  GEWEX Atmospheric Boundary Layer Study (GABLS3) for July 1, 2006
   development of a nocturnal low-level jet

-  Multidisciplinary drifting Observatory for the Study of Arctic Climate
   expedition (MOSAiC)
   -  SS: Strongly stably stratified boundary layer (March 2-10 2020)
   -  AMPS: Arctic mixed-phase stratocumuluas cloud (Oct 31 - Nov 5 2019)

-  Cold-Air Outbreaks in the Marine Boundary Layer Experiment (COMBLE) for
   March 12, 2020 mixed phased clouds in the polar marine boundary layer

For the ARM SGP case, several case configuration files representing
different time periods of the observational dataset are included,
denoted by a trailing letter. The LASSO case may be run with different
forcing applied, so three case configuration files corresponding to
these different forcing are included.

In addition, cases can be generated from UFS initial conditions See
:numref:`Section %s <UFScasegen>` for information on how to generate these
files for other locations and dates, given appropriate UFS Atmosphere
initial conditions and output.

How to set up new cases
-----------------------

Setting up a new case involves preparing the two types of files listed
above. For the case initialization and forcing data file, this typically
involves writing a custom script or program to parse the data from its
original format to the DEPHY format, listed above. Formatting for DEPHY
is documented in the `DEPHY repository
<https://github.com/GdR-DEPHY/DEPHY-SCM/blob/master/DEPHY-SCM_CommonFormat_v1.0.pdf>`__.

For reference, the following formulas are used:

.. math:: \theta_{il} = \theta - \frac{\theta}{T}\left(\frac{L_v}{c_p}q_l + \frac{L_s}{c_p}q_i\right)

.. math:: q_t = q_v + q_l + q_i

where :math:`\theta_{il}` is the ice-liquid water potential temperature,
:math:`\theta` is the potential temperature, :math:`L_v` is the latent
heat of vaporization, :math:`L_s` is the latent heat of sublimation
:math:`c_p` is the specific heat capacity of air at constant pressure,
:math:`T` is absolute temperature, :math:`q_t` is the total water
specific humidity, :math:`q_v` is the water vapor specific humidity,
:math:`q_l` is the suspended liquid water specific humidity, and
:math:`q_i` is the suspended ice water specific humidity.

For the case configuration file, it is most efficient to copy an
existing file in ``ccpp-scm/scm/etc/case_config`` and edit it to suit one’s case. Recall from subsection
:numref:`Subsection %s <case config>` that this file is used to configure
the SCM framework parameters for a given case. Be sure to check that
model timing parameters such as the time step and output frequency are
appropriate for the physics suite being used. There is likely some
stability criterion that governs the maximum time step based on the
chosen parameterizations and number of vertical levels (grid spacing).
The ``case_name`` parameter should match the name of the case input data file that was
configured for the case (without the file extension). The parameter
should be less than or equal to the length of the forcing data unless
the desired behavior of the simulation is to proceed with the last
specified forcing values after the length of the forcing data has been
surpassed. If the case input data is
specified to a lower altitude than the vertical domain, the remainder of
the column will be filled in with values from a reference profile. There
is a tropical profile and mid-latitude summer profile provided, although
one may add more choices by adding a data file to ``ccpp-scm/scm/data/processed_case_input`` and adding a parser
section to the subroutine ``get_reference_profile`` in ``scm/src/scm_input.f90``. Surface fluxes can either be specified in
the case input data file or calculated using a surface scheme using
surface properties. In addition, one must specify a ``column_area`` for each column.

To control the forcing method, one must choose how the momentum and
scalar variable forcing are applied. The three methods of Randall and
Cripe (1999, JGR) have been implemented: “revealed forcing” where total
(horizontal :math:`+` vertical) advective tendencies are applied (type
1), “horizontal advective forcing” where horizontal advective tendencies
are applied and vertical advective tendencies are calculated from a
prescribed vertical velocity and the calculated (modeled) profiles (type
2), and “relaxation forcing” where nudging to observed profiles replaces
horizontal advective forcing combined with vertical advective forcing
from prescribed vertical velocity (type 3). If relaxation forcing is
chosen, a ``relax_time`` that represents the timescale over which the profile would
return to the nudging profiles must be specified.

.. _`lasso`:

Using other LASSO cases
-----------------------

In order to use other LASSO cases than the one provided, perform the
following steps:

#. Access http://archive.arm.gov/lassobrowser and use the navigation on
   the left to choose the dates for which you would like to run a SCM
   simulation. Pay attention to the “Large Scale Forcing” tab where you
   can choose how the large scale forcing was generated, with options
   for ECMWF, MSDA, and VARANAL. All are potentially valid, and it is
   likely worth exploring the differences among forcing methods. Click
   on Submit to view a list of simulations for the selected criteria.
   Choose from the simulations (higher skill scores are preferred) and
   check the “Config Obs Model Tar” box to download the data. Once the
   desired simulations have been checked, order the data (you may need
   to create an ARM account to do so).

#. Once the data is downloaded, decompress it. From the ``config`` directory, copy
   the files ``input_ls_forcing.nc``, ``input_sfc_forcing.nc``, and ``wrfinput_d01.nc`` into their own directory under ``ccpp-scm/scm/data/raw_case_input/``.

#. Modify ``ccpp-scm/scm/etc/scripts/lasso1_forcing_file_generator_gjf.py`` to point to the input files listed above. Execute the script
   in order to generate a case input file for the SCM (to be put in ``ccpp-scm/scm/data/processed_case_input/``):

   .. code:: bash

      ./lasso1_forcing_file_generator_gjf.py

#. Create a new case configuration file (or copy and modify an existing
   one) in ``ccpp-scm/scm/etc/case_config``. Be sure that the ``case_name`` variable points to the newly
   created/processed case input file from above.

.. _`UFScasegen`:

Using UFS Output to Create SCM Cases: UFS Case Generation
---------------------------------------------------------

.. _`pydepend_casegen`:

Python Dependencies
~~~~~~~~~~~~~~~~~~~

The scripts here require a few python packages that may not be found by
default in all python installations. There is a YAML file with the
python environment needed to run the script in ``ccpp-scm/environment-ufscasegen.yml``. To create and activate
this environment using conda:

Create environment (only once):

.. code:: bash

  > conda env create -f environment-ufscasegen.yml

This will create the conda environment ``env_ufscasegen``

Activate environment:

.. code:: bash

  > conda activate env_ufscasegen

.. _`ufscasegen`:

UFS_case_gen.py
~~~~~~~~~~~~~~~

A script exists in ``scm/etc/scripts/UFS_case_gen.py`` to read in UFS history (output) files and their
initial conditions to generate a SCM case input data file, in DEPHY
format.

.. code:: bash

   ./UFS_case_gen.py [-h] (-l LOCATION LOCATION | -ij INDEX INDEX) -d
   DATE -i IN_DIR -g GRID_DIR -f FORCING_DIR -n
   CASE_NAME [-t {1,2,3,4,5,6,7}] [-a AREA] [-oc]
   [-lam] [-sc] [-near] [-fm] [-vm] [-wn] [-geos]

Mandatory arguments:

#. ``--location (-l)`` OR ``--index (-ij)``: Either longitude and latitude in decimal degrees east and north
   of a location OR the UFS grid index with the tile number

   -  -l 261.51 38.2 (two floating point values separated by a space)

   -  -ij 8 49 (two integer values separated by a space; this option
      must also use the ``--tile (-t)`` argument to specify the tile number)

#. ``--date (-d)`` YYYYMMDDHHMMSS: date corresponding to the UFS initial conditions

#. ``--in_dir (-i)``: path to the directory containing the UFS initial conditions

#. ``--grid_dir (-g)``: path to the directory containing the UFS supergrid files (AKA "fix"
   directory)

#. ``--forcing_dir (-f)``: path to the directory containing the UFS history files

#. ``--case_name (-n)``: name of case

Optional arguments:

#. ``--tile (-t)``: if one already knows the correct tile for the given longitude and
   latitude OR one is specifying the UFS grid index (``--index`` argument)

#. ``--area (-a)``: area of grid cell in :math:`m^2` (if known or different than the
   value calculated from the supergrid file)

#. ``--old_chgres (-oc)``: flag if UFS initial conditions were generated using older version
   of chgres (global_chgres); might be the case for pre-2018 data

#. ``--lam (-lam)``: flag to signal that the ICs and forcing is from a limited-area
   model run

#. ``--save_comp (-sc)``: flag to create UFS reference file for comparison

#. ``--use_nearest (-near)``: flag to indicate using the nearest UFS history file gridpoint for calculation 
   of forcing; only valid for use with -fm=1 or -fm=3

#. ``--forcing_method (-fm)``: method used to calculate forcing (1=total tendencies from UFS dycore,
   2=advective terms calculated from UFS history files, 3=total time tendency terms calculated), default=2

#. ``--vertical_method (-vm)``: method used to calculate vertical advective forcing (1=vertical advective
   terms calculated from UFS history files and added to total, 2=smoothed vertical velocity provided), default=2

#. ``--wind_nudge (-wn)``: flag to turn on wind nudging to UFS profiles

#. ``--geostrophic (-geos)``: flag to turn on geostrophic wind forcing

Notes Regarding Implemented Forcing Methods

The ``--forcing_method`` option hides some complexity that should be understood when running this script since 
each method has a particular use case and produces potentially very different forcing terms. Forcing method 1 is 
designed to be used in concert with the three-dimensional UFS. I.e., the UFS must be run with diagnostic tendencies 
activated so that the `nophysics` term is calculated and output for all grid points 
(see https://ccpp-techdoc.readthedocs.io/en/latest/ParamSpecificOutput.html#tendencies). This diagnostic term 
represents the tendency produced for each state variable by the UFS between calls to the "slow" physics. This 
includes the tendency due to advection, but also any tendencies due to other non-physics processes, e.g. "fast"
physics, coupling to external components, data assimilation, etc. Within the SCM, this diagnostic is used as the 
forcing term for each state variable. Although one can achieve results as close as possible between a UFS column 
and the single column model using this method, it will NOT be bit-for-bit for many reasons. Some of these reasons
include: diagnostic output is not typically instantaneous for every timestep, the UFS' vertical coordinate is 
semi-Lagrangian and includes a remapping step as the surface pressure changes for each column, whereas the SCM
uses a Eulerian vertical coordinate without the vertical remapping step, and some interpolation happens in the
UFS_case_gen.py script due to the UFS initial conditions and history files using different grids. This method
can only be used when the UFS has been configured and run with the anticipation of running the SCM using this
forcing method afterward because it requires considerable extra disk space for the additional output.

The `--forcing_method` 2 option is the most general in the sense that the same method could apply to any three-dimensional 
model output. For a given column, it uses a configurable number of neighboring grid points to calculate the horizontal 
advective tendencies using the horizontal components of the three-dimensional wind and horizontal derivatives of the 
state variables. Note that the script performs some smoothing in the vertical profiles used to calculate the advective 
terms in order to eliminate small-scale noise in the forcing terms and the derivatives are calculated using a second- or
fourth-order centered difference scheme, depending on the number of neighboring points used. Vertical advective terms 
are calculated based on the specification of `--vertical_method` (-vm). For vertical_method 1, the vertical advective 
terms are calculated from the history files using UFS vertical velocities and the same modeled smoothed vertical profiles 
of the state variables using the upstream scheme. Note that while the horizontal terms use neighboring points, the vertical 
advective terms only use the central, chosen column. This method is sometimes referred to as "total advective forcing" and
tends to be less "responsive" to the SCM-modeled state. I.e., a SCM run using vertical method 2 has a greater chance of 
deviating from the UFS column state and not being able to "recover". For this reason, vertical method 2 is often used 
in the literature, whereby the vertical velocity profile from the three-dimensional model is provided as forcing to the SCM
and the vertical advective terms are calculated during the SCM integration using the SCM-modeled state variable profiles.

The final forcing method, 3, uses the three-dimensional history files to calculate profiles of the total time-rate of change
of the state variables to use as forcing for the SCM. Note that this is tantamount to strongly nudging the SCM state to the 
UFS state and already intrinsically includes both the physics and dynamics tendencies. While a simulation using this forcing 
is more-or-less guaranteed to produce a SCM simulation that closely matches the three-dimensional output of the state variables,
it strongly minimizes the contribution of physics in the SCM simulation. Indeed, an SCM simulation without running a physics suite 
at all would still be expected to closely track the mean state of the three-dimensional column, so this method will likely be of 
limited use for physics studies.

Forcing the horizontal components of the wind can be notoriously difficult in SCMs, and the most straightforward method is to 
simply nudge them to the three-dimensional modeled state. This method is achieved by using the `--wind_nudge` (-wn) option and uses a nudging 
timescale of one hour. It should be possible to calculate a nudging timescale based on the magnitude of the wind in the neighboring 
grid cells, although this is not implemented yet.

The second method to force the horizontal wind components is to calculate the geostrophic wind using the "large scale" pressure 
gradient from the three-dimensional model. This is achieved by using the `--geostrophic` (-geos) option. What qualifies as large 
enough of a horizontal gridscale to achieve geostrophic balance is when the Rossby number is much less than one. The script uses a 
configurable Rossby number (around 0.1) to expand the number of neighboring grid points such that geostrophic balance can be assumed 
given the particular UFS history file grid. The geostrophic winds are calculated using the horizontal geopotential gradient and the 
local latitude-dependent Coriolis parameter. From the PBL top downward, the geostrophic winds are assumed to go to zero. In testing 
with this method, if the initial horizontal winds have a significant ageostrophic component (the initial condition winds are 
appreciably different than the calculated geostrophic winds), this often leads to spurious clockwise turning of the mean modeled winds 
with time. An option exists within the script to assume that the mean three-dimensional winds are, in fact, identical to the 
geostrophic winds as well. Using this option eliminates any spurious turning.

.. _`ufsforcingensemblegenerator`:

UFS_forcing_ensemble_generator.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is an additional script in ``scm/etc/scripts/UFS_forcing_ensemble_generator.py`` to create UFS-caseGen case(s) starting
with output from UFS Weather Model (UWM) Regression Tests (RTs).

.. code:: bash

   UFS_forcing_ensemble_generator.py [-h] -d DIR -n CASE_NAME
   (-lonl LON_1 LON_2 -latl LAT_1 LAT_2 -nens NENSMEMBERS |
   -lons [LON_LIST] -lats [LAT_LIST] |
   -fxy [LON_LAT_FILE])
   [-dt TIMESTEP] [-cres C_RES] [-sdf SUITE] [-sc] [-near] [-fm] [-vm] [-wn] [-geos]

Mandatory arguments:

#. ``--dir (-d)``: path to UFS Regression Test output

#. ``--case_name (-n)``: name of cases

#. Either: (see examples below)

   -  ``--lon_limits (-lonl)`` AND ``--lat_limits (-latl)`` AND ``--nensmembers (-nens)``: longitude range, latitude range, and number of cases to
      create

   -  ``--lon_list (-lons)`` AND ``--lat_list (-lats)``: longitude and latitude of cases

   -  ``--lonlat_file (fxy)``: file containing longitudes and latitudes

Optional arguments:

#. ``--timestep (-dt)``: SCM timestep, in seconds

#. ``--C_res (-cres)``: UFS spatial resolution

#. ``--suite (-sdf)``: CCPP suite definition file to use for ensemble

#. ``--save_comp (-sc)``: flag to create UFS reference file for comparison

#. ``--use_nearest (-near)``: flag to indicate using the nearest UFS history file gridpoint

#. ``--forcing_method (-fm)``: method used to calculate forcing (1=total tendencies from UFS dycore,
   2=advective terms calculated from UFS history files, 3=total time tendency terms calculated), default=2

#. ``--vertical_method (-vm)``: method used to calculate vertical advective forcing (1=vertical advective
   terms calculated from UFS history files and added to total, 2=smoothed vertical velocity provided), default=2

#. ``--wind_nudge (-wn)``: flag to turn on wind nudging to UFS profiles

#. ``--geostrophic (-geos)``: flag to turn on geostrophic wind forcing

Examples to run from within the ``scm/etc/scripts`` directory to create SCM cases starting
with the output from a UFS Weather Model regression test(s):

On the supported platforms Derecho (NCAR) and Hera (NOAA), there are
staged UWM RTs located at:

-  Derecho ``/glade/scratch/epicufsrt/GMTB/CCPP-SCM/UFS_RTs``
-  Hera ``/scratch1/BMC/gmtb/CCPP-SCM/UFS_RTs``

.. _`example1`:

Example 1: UFS-caseGen for single point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UFS regression test, ``control_c192``, for single point.

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d [path_to_regression_tests_output]/control_c192_intel/ -sc --C_RES 192 -dt 360  -n control_c192 -lons 300 -lats 34

Upon successful completion of the script, the command to run the case(s)
will print to the screen. For example,

.. code:: bash

   ./run_scm.py --npz_type gfs --file scm_ufsens_control_c192.py --timestep 360

The file ``scm_ufsens_control_c192.py`` is created in ``ccpp-scm/scm/bin/``, where the SCM run script is to be exectued.

.. _`example2`:

Example 2: UFS-caseGen for list of points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UFS regression test, ``control_c384``, for multiple points.

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d /glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/NEMSfv3gfs/develop-20240607/control_c384_intel/ -sc --C_RES 384 -dt 225 -n control_c384 -lons 300 300 300 300 -lats 34 35 35 37

Upon successful completion of the script, the command to run the case(s)
will print to the screen. For example,

.. code:: bash

   ./run_scm.py --npz_type gfs --file scm_ufsens_control_c384.py --timestep 225

The file ``scm_ufsens_control_c384.py`` contains ALL of the cases created. Each case created will have the
naming convention ``case_name_nXXX``, where the suffix ``XXX`` is the case number from 0 to the
number of points provided. The contents of the file should look like:

.. code:: bash

   run_list = [{"case": "control_c384_n000", "suite": "SCM_GFS_v16"},
               {"case": "control_c384_n001", "suite": "SCM_GFS_v16"},
               {"case": "control_c384_n002", "suite": "SCM_GFS_v16"},
               {"case": "control_c384_n003", "suite": "SCM_GFS_v16"}]

.. _`example3`:

Example 3: UFS-caseGen for an ensemble of points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UFS regression test, ``control_p8``, for an ensemble (10) of randomly selected points
over a specified longitude (:math:`300-320^oW`) and latitude
(:math:`40-50^oN`) range

But first, to use the ``control_p8`` test we need to rerun the regression test to
generate UFS history files with a denser and constant output interval.
First, in ``control_p8/model_configure``, change ``--output_fh`` to ``"interval -1"``, where ``interval`` is the UFS history file output frequency
(in hours), see `UFS Weather Model Users
Guide <https://ufs-weather-model.readthedocs.io/en/latest/InputsOutputs.html>`__
for more details.

For the purposes of this example the ``control_p8`` test has already been rerun, but if
starting from your own UWM RTs, you can rerun the UWM regression test,
on Derecho for example, by running the following command in the RT
directory: ``qsub job_card``

Now the cases can be generated with the following command:

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d /glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/NEMSfv3gfs/develop-20240607/control_p8_intel -sc --C_RES 96 -dt 720 -n control_p8 -lonl 300 320 -latl 40 50 -nens 10 -sdf SCM_GFS_v17_p8

Upon successful completion of the script, the command to run the case(s)
will print to the screen. For example,

.. code:: bash

   ./run_scm.py --npz_type gfs --file scm_ufsens_control_p8.py --timestep 720

The file ``scm_ufsens_control_p8.py`` contains ten cases (n000-n009) to be run. The contents of the
file should look like:

.. code:: bash

   run_list = [{"case": "control_p8_n000", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n001", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n002", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n003", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n004", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n005", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n006", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n007", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n008", "suite": "SCM_GFS_v17_p8"},
               {"case": "control_p8_n009", "suite": "SCM_GFS_v17_p8"}]
