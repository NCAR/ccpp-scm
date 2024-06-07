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

-  ``runtime``

   -  Specify the model runtime in seconds (integer). This should
      correspond with the forcing dataset used. If a runtime is
      specified that is longer than the supplied forcing, the forcing is
      held constant at the last specified values.

-  ``thermo_forcing_type``

   -  An integer representing how forcing for temperature and moisture
      state variables is applied (1 :math:`=` total advective
      tendencies, 2 :math:`=` horizontal advective tendencies with
      prescribed vertical motion, 3 :math:`=` relaxation to observed
      profiles with vertical motion prescribed)

-  ``mom_forcing_type``

   -  An integer representing how forcing for horizontal momentum state
      variables is applied (1 :math:`=` total advective tendencies; not
      implemented yet, 2 :math:`=` horizontal advective tendencies with
      prescribed vertical motion, 3 :math:`=` relaxation to observed
      profiles with vertical motion prescribed)

-  ``relax_time``

   -  A floating point number representing the timescale in seconds for
      the relaxation forcing (only used if ``thermo_forcing_type = 3`` or ``mom_forcing_type = 3``)

-  ``sfc_flux_spec``

   -  A boolean set to ``.true.`` if surface flux are specified from the forcing
      data (there is no need to have surface schemes in a suite
      definition file if so)

-  ``sfc_roughness_length_cm``

   -  Surface roughness length in cm for calculating surface-related
      fields from specified surface fluxes (only used if ``sfc_flux_spec`` is True).

-  ``sfc_type``

   -  An integer representing the character of the surface (0 :math:`=`
      sea surface, 1 :math:`=` land surface, 2 :math:`=` sea-ice
      surface)

-  ``reference_profile_choice``

   -  An integer representing the choice of reference profile to use
      above the supplied initialization and forcing data (1 :math:`=`
      “McClatchey” profile, 2 :math:`=` mid-latitude summer standard
      atmosphere)

-  ``year``

   -  An integer representing the year of the initialization time

-  ``month``

   -  An integer representing the month of the initialization time

-  ``day``

   -  An integer representing the day of the initialization time

-  ``hour``

   -  An integer representing the hour of the initialization time

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

   -  0 => original DTC format, 1 => DEPHY-SCM format.

Optional variables (that may be overridden via run script command line
arguments) are:

-  ``vert_coord_file``

   -  File containing FV3 vertical grid coefficients.

-  ``n_levels``

   -  Specify the integer number of vertical levels.

.. _`case input`:

Case input data file (CCPP-SCM format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The initialization and forcing data for each case is stored in a NetCDF
(version 4) file within the ``ccpp-scm/scm/data/processed_case_input`` directory. Each file has at least two
dimensions (``time`` and ``levels``, potentially with additions for vertical snow and soil
levels) and is organized into 3 groups: scalars, initial, and forcing.
Not all fields are required for all cases. For example the fields ``sh_flux_sfc`` and ``lh_flux_sfc``
are only needed if the variable ``sfc_flx_spec = .true.`` in the case configuration file
and state nudging variables are only required if ``thermo_forcing_type = 3`` or ``mom_forcing_type = 3``
. Using an active LSM (Noah, NoahMP, RUC) requires many more variables
than are listed here. Example files for using with Noah and NoahMP LSMs
are included in ``ccpp-scm/scm/data/processed_case_input/fv3_model_point_noah[mp].nc``.

.. _`case input arm`:
.. literalinclude:: arm_case_header.txt
    :name: lst_case_input_netcdf_header_arm
    :caption: example NetCDF file (CCPP-SCM format) header for case initialization and forcing data
 
Case input data file (DEPHY format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Development and Evaluation of Physics in atmospheric models (DEPHY)
format is an internationally-adopted data format intended for use by SCM
and LESs. The initialization and forcing data for each case is stored in
a NetCDF (version 4) file, although these files are not by default
included in the CCPP SCM repository. To access these cases you need to
clone the DEPHY-SCM repository, and provide the DEPHY-SCM file location
to the SCM. For example:

.. code:: bash

   cd [...]/ccpp-scm/scm/data
   git clone https://github.com/GdR-DEPHY/DEPHY-SCM DEPHY-SCM
   cd [...]/ccpp-scm/scm/bin
   ./run_scm.py -c MAGIC_LEG04A --case_data_dir [...]/ccpp-scm/scm/data/DEPHY-SCM/MAGIC/LEG04A -v

Each DEPHY file has three dimensions (``time``, ``t0``, ``levels``) and contains the initial
conditions (``t0``, ``levels``) and forcing data (``time``, ``levels``). Just as when using the CCPP-SCM
formatted inputs, :numref:`Subsection %s <case input>`, not all fields
are required for all cases. More information on the DEPHY format
requirements can be found at
`DEPHY <https://github.com/GdR-DEPHY/DEPHY-SCM>`__.

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

For the ARM SGP case, several case configuration files representing
different time periods of the observational dataset are included,
denoted by a trailing letter. The LASSO case may be run with different
forcing applied, so three case configuration files corresponding to
these different forcing are included. In addition, two example cases are
included for using UFS Atmosphere initial conditions:

-  UFS initial conditions for 38.1 N, 98.5 W (central Kansas) for 00Z on
   Oct. 3, 2016 with Noah variables on the C96 FV3 grid (``fv3_model_point_noah.nc``)

-  UFS initial conditions for 38.1 N, 98.5 W (central Kansas) for 00Z on
   Oct. 3, 2016 with NoahMP variables on the C96 FV3 grid (``fv3_model_point_noahmp.nc``)

See :numref:`Section %s <UFSreplay>` for information on how to generate these
files for other locations and dates, given appropriate UFS Atmosphere
initial conditions and output.

How to set up new cases
-----------------------

Setting up a new case involves preparing the two types of files listed
above. For the case initialization and forcing data file, this typically
involves writing a custom script or program to parse the data from its
original format to the format that the SCM expects, listed above. An
example of this type of script written in Python is included in ``ccpp-scm/scm/etc/scripts/twpice_forcing_file_generator.py``. The
script reads in the data as supplied from its source, converts any
necessary variables, and writes a NetCDF (version 4) file in the format
described in subsections :numref:`Subsection %s <case input>` and
:numref:`Subsection %s <case input dephy>`. For reference, the following
formulas are used:

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

As shown in the example NetCDF header, the SCM expects that the vertical
dimension is pressure levels (index 1 is the surface) and the time
dimension is in seconds. The initial conditions expected are the height
of the pressure levels in meters, and arrays representing vertical
columns of :math:`\theta_{il}` in K, :math:`q_t`, :math:`q_l`, and
:math:`q_i` in kg kg\ :math:`^{-1}`, :math:`u` and :math:`v` in m
s\ :math:`^{-1}`, turbulence kinetic energy in m\ :math:`^2`
s\ :math:`^{-2}` and ozone mass mixing ratio in kg kg\ :math:`^{-1}`.

For forcing data, the SCM expects a time series of the following
variables: latitude and longitude in decimal degrees [in case the
column(s) is moving in time (e.g., Lagrangian column)], the surface
pressure (Pa) and surface temperature (K). If surface fluxes are
specified for the new case, one must also include a time series of the
kinematic surface sensible heat flux (K m s\ :math:`^{-1}`) and
kinematic surface latent heat flux (kg kg\ :math:`^{-1}` m
s\ :math:`^{-1}`). The following variables are expected as 2-dimensional
arrays (vertical levels first, time second): the geostrophic u (E-W) and
v (N-S) winds (m s\ :math:`^{-1}`), and the horizontal and vertical
advective tendencies of :math:`\theta_{il}` (K s\ :math:`^{-1}`) and
:math:`q_t` (kg kg\ :math:`^{-1}` s\ :math:`^{-1}`), the large scale
vertical velocity (m s\ :math:`^{-1}`), large scale pressure vertical
velocity (Pa s\ :math:`^{-1}`), the prescribed radiative heating rate (K
s\ :math:`^{-1}`), and profiles of u, v, T, :math:`\theta_{il}` and
:math:`q_t` to use for nudging.

Although it is expected that all variables are in the NetCDF file, only
those that are used with the chosen forcing method are required to be
nonzero. For example, the following variables are required depending on
the values of ``mom_forcing_type`` and ``thermo_forcing_type`` specified in the case configuration file:

-  ``mom_forcing_type = 1``

   -  Not implemented yet

-  ``mom_forcing_type = 2``

   -  geostrophic winds and large scale vertical velocity

-  ``mom_forcing_type = 3``

   -  u and v nudging profiles

-  ``thermo_forcing_type = 1``

   -  horizontal and vertical advective tendencies of
      :math:`\theta_{il}` and :math:`q_t` and prescribed radiative
      heating (can be zero if radiation scheme is active)

-  ``thermo_forcing_type = 2``

   -  horizontal advective tendencies of :math:`\theta_{il}` and
      :math:`q_t`, prescribed radiative heating (can be zero if
      radiation scheme is active), and the large scale vertical pressure
      velocity

-  ``thermo_forcing_type = 3``

   -  :math:`\theta_{il}` and :math:`q_t` nudging profiles and the large
      scale vertical pressure velocity

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
surpassed. The initial date and time should fall within the forcing
period specified in the case input data file. If the case input data is
specified to a lower altitude than the vertical domain, the remainder of
the column will be filled in with values from a reference profile. There
is a tropical profile and mid-latitude summer profile provided, although
one may add more choices by adding a data file to ``ccpp-scm/scm/data/processed_case_input`` and adding a parser
section to the subroutine ``get_reference_profile`` in ``scm/src/scm_input.f90``. Surface fluxes can either be specified in
the case input data file or calculated using a surface scheme using
surface properties. If surface fluxes are specified from data, set ``sfc_flux_spec`` to ``.true.``
and specify ``sfc_roughness_length_cm`` for the surface over which the column resides. Otherwise,`
specify a ``sfc_type``. In addition, one must specify a ``column_area`` for each column.

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

.. _`UFSreplay`:

Using UFS Output to Create SCM Cases: UFS-Replay
------------------------------------------------

.. _`pydepend_replay`:

Python Dependencies
~~~~~~~~~~~~~~~~~~~

The scripts here require a few python packages that may not be found by
default in all python installations. There is a YAML file with the
python environment needed to run the script in ``ccpp-scm/environment-ufsreplay.yml``. To create and activate
this environment using conda:

Create environment (only once):

.. code:: bash

  > conda env create -f environment-ufsreplay.yml

This will create the conda environment ``env_ufsreplay``

Activate environment:

.. code:: bash

  > conda activate env_ufsreplay

.. _`ufsicgenerator`:

UFS_case_gen.py
~~~~~~~~~~~~~~~~~~~

A script exists in ``scm/etc/scripts/UFS_case_gen.py`` to read in UFS history (output) files and their
initial conditions to generate a SCM case input data file, in DEPHY
format.

.. code:: bash

   ./UFS_case_gen.py [-h] (-l LOCATION LOCATION | -ij INDEX INDEX) -d
   DATE -i IN_DIR -g GRID_DIR -f FORCING_DIR -n
   CASE_NAME [-t {1,2,3,4,5,6,7}] [-a AREA] [-oc]
   [-lam] [-sc] [-near]

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

#. ``--use_nearest (-near)``: flag to indicate using the nearest UFS history file gridpoint

.. _`ufsforcingensemblegenerator`:

UFS_forcing_ensemble_generator.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is an additional script in ``scm/etc/scripts/UFS_forcing_ensemble_generator.py`` to create UFS-replay case(s) starting
with output from UFS Weather Model (UWM) Regression Tests (RTs).

.. code:: bash

   UFS_forcing_ensemble_generator.py [-h] -d DIR -n CASE_NAME
   (-lonl LON_1 LON_2 -latl LAT_1 LAT_2 -nens NENSMEMBERS |
   -lons [LON_LIST] -lats [LAT_LIST])
   [-dt TIMESTEP] [-cres C_RES] [-sdf SUITE] [-sc] [-near]

Mandatory arguments:

#. ``--dir (-d)``: path to UFS Regression Test output

#. ``--case_name (-n)``: name of cases

#. Either: (see examples below)

   -  ``--lon_limits (-lonl)`` AND ``--lat_limits (-latl)`` AND ``--nensmembers (-nens)``: longitude range, latitude range, and number of cases to
      create

   -  ``--lon_list (-lons)`` AND ``--lat_list (-lats)``: longitude and latitude of cases

Optional arguments:

#. ``--timestep (-dt)``: SCM timestep, in seconds

#. ``--C_res (-cres)``: UFS spatial resolution

#. ``--suite (-sdf)``: CCPP suite definition file to use for ensemble

#. ``--save_comp (-sc)``: flag to create UFS reference file for comparison

#. ``--use_nearest (-near)``: flag to indicate using the nearest UFS history file gridpoint

Examples to run from within the ``scm/etc/scripts`` directory to create SCM cases starting
with the output from a UFS Weather Model regression test(s):

On the supported platforms Cheyenne (NCAR) and Hera (NOAA), there are
staged UWM RTs located at:

-  Cheyenne ``/glade/scratch/epicufsrt/GMTB/CCPP-SCM/UFS_RTs``
-  Hera ``/scratch1/BMC/gmtb/CCPP-SCM/UFS_RTs``

.. _`example1`:

Example 1: UFS-replay for single point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UFS regression test, ``control_c192``, for single point.

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d /glade/scratch/epicufsrt/GMTB/CCPP-SCM/UFS_RTs/control_c192/ -sc --C_RES 192 -dt 360  -n control_c192 -lons 300 -lats 34

Upon successful completion of the script, the command to run the case(s)
will print to the screen. For example,

.. code:: bash

   ./run_scm.py --npz_type gfs --file scm_ufsens_control_c192.py --timestep 360

The file ``scm_ufsens_control_c192.py`` is created in ``ccpp-scm/scm/bin/``, where the SCM run script is to be exectued.

.. _`example2`:

Example 2: UFS-replay for list of points
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

UFS regression test, ``control_c384``, for multiple points.

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d /glade/scratch/epicufsrt/GMTB/CCPP-SCM/UFS_RTs/control_c384/ -sc --C_RES 384 -dt 225 -n control_c384 -lons 300 300 300 300 -lats 34 35 35 37

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

Example 3: UFS-replay for an ensemble of points
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
on Cheyenne for example, by running the following command in the RT
directory: ``qsub job_card``

Now the cases can be generated with the following command:

.. code:: bash

   ./UFS_forcing_ensemble_generator.py -d /glade/scratch/epicufsrt/GMTB/CCPP-SCM/UFS_RTs/control_p8/ -sc --C_RES 96 -dt 720 -n control_p8 -lonl 300 320 -latl 40 50 -nens 10 -sdf SCM_GFS_v17_p8

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
