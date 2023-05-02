!> \file scm_input.f90
!!  Contains input-related subroutines -- reading in model configuration from file or the command line and reading in the case
!!  initial conditions and forcing; also contains reference profile input (temporarily hard-coded).

module scm_input

use scm_kinds, only : sp, dp, qp
use netcdf
use scm_type_defs, only: character_length

implicit none

integer :: missing_snow_layers = 3
integer :: missing_soil_layers = 4
integer :: missing_ice_layers = 2

contains

!> \ingroup SCM
!! @{
!! \defgroup input scm_input
!! @{
!! Contains input-related subroutines -- reading in model configuration from file or the command line and reading in the case
!! initial conditions and forcing; also contains reference profile input (temporarily hard-coded).

!> Subroutine to get basic model configuration data from a namelist file (whose name is specified as the first argument on the command line) and from
!! data entered on the command line with the format: "var1='string' var2=d.d var3=i". Namelist variables listed on the command line
!! override those that are specified in the external namelist file. The case configuration namelist variables are also written out to
!! a namelist file placed in the output directory. Note: This routine uses GET_COMMAND which is an intrinsic routine in the Fortran 2003 standard. This
!! requires that the compiler supports this standard.
subroutine get_config_nml(scm_state)
  use scm_type_defs, only : scm_state_type

  type(scm_state_type), target, intent(inout) :: scm_state

  character(len=character_length)    :: experiment_name !< name of the experiment configuration file (usually case name)
  character(len=character_length)    :: npz_type !< used for generating different FV3 vertical grids
  character(len=character_length)    :: vert_coord_file !< file containing FV3 vertical grid coefficients
  character(len=character_length)    :: case_name !< name of case initialization and forcing dataset
  real(kind=dp)        :: dt !< time step in seconds
  real(kind=dp)        :: runtime !< total runtime in seconds
  integer              :: n_itt_out !< multiple of timestep for writing output
  integer              :: n_itt_diag !< multiple of timestep for resetting diagnostics (overwrites fhzero from physics namelist if present)
  integer              :: n_levels !< number of model levels (currently only 64 supported)
  integer              :: n_soil   !< number of model soil levels (currently only 4 supported)
  integer              :: n_snow   !< number of model snow levels (currently only 3 supported)
  integer              :: n_columns !< number of columns to use
  integer              :: n_time_levels
  integer              :: time_scheme !< 1 => forward Euler, 2 => filtered leapfrog
  character(len=character_length)    :: output_dir !< name of the output directory
  character(len=character_length)    :: output_file !< name of the output file (without the file extension)
  integer              :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer              :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer              :: C_RES            !< reference "C" resoltiion of FV3 grid (needed for GWD and mountain blocking)
  integer              :: spinup_timesteps
  real(kind=dp)        :: relax_time !< relaxation time scale (s)
  logical              :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
  real(kind=dp)        :: sfc_roughness_length_cm !< surface roughness length used for calculating surface layer parameters from specified fluxes
  integer              :: sfc_type !< 0: sea surface, 1: land surface, 2: sea-ice surface
  logical              :: model_ics !<  true means have land info too
  logical              :: lsm_ics !< true when LSM initial conditions are included (but not all ICs from another model)
  logical              :: do_spinup
  integer              :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere
  integer              :: year, month, day, hour, min
  real(kind=dp)        :: column_area
  integer              :: input_type !< 0 => original DTC format, 1 => DEPHY-SCM format

  character(len=character_length)    :: physics_suite !< name of the physics suite name (currently only GFS_operational supported)
  character(len=character_length)    :: physics_nml
  
  character(len=character_length), allocatable, dimension(:) :: tracer_names
  integer,                         allocatable, dimension(:) :: tracer_types

  integer                          :: ioerror

  CHARACTER(LEN=*), parameter :: experiment_namelist = 'input_experiment.nml'

  NAMELIST /case_config/ npz_type, vert_coord_file, case_name, dt, runtime, n_itt_out, n_itt_diag, &
    n_levels, output_dir, thermo_forcing_type, model_ics, &
    lsm_ics, do_spinup, C_RES, spinup_timesteps, mom_forcing_type, relax_time, sfc_type, sfc_flux_spec, &
    sfc_roughness_length_cm, reference_profile_choice, year, month, day, hour, min, &
    column_area, input_type
    
  NAMELIST /physics_config/ physics_suite, physics_nml

  !>  \section get_config_alg Algorithm
  !!  @{

  !> Define default values for experiment configuration (to be overridden by external namelist file or command line arguments)
  npz_type = ''
  vert_coord_file = ''
  n_columns = 1
  case_name = 'twpice'
  dt = 600.0
  time_scheme = 1
  runtime = 2138400.0
  n_itt_out = 1
  n_itt_diag = -999
  n_levels = 127
  n_soil   = 4
  n_snow   = 3
  output_dir = 'output'
  output_file = 'output'
  thermo_forcing_type = 2
  mom_forcing_type = 3
  C_RES = 384
  spinup_timesteps = 0
  relax_time = 7200.0
  sfc_flux_spec = .false.
  sfc_roughness_length_cm = 1.0
  sfc_type = 0
  model_ics = .false.
  lsm_ics = .false.
  do_spinup = .false.
  reference_profile_choice = 1
  year = 2006
  month = 1
  day = 19
  hour = 3
  min = 0
  input_type = 0
  
  open(unit=10, file=experiment_namelist, status='old', action='read', iostat=ioerror)
  if(ioerror /= 0) then
    write(*,'(a,i0)') 'There was an error opening the file ' // experiment_namelist // &
                      '; error code = ', ioerror
    STOP
  else
    read(10, NML=case_config, iostat=ioerror)
  end if

  if(ioerror /= 0) then
    write(*,'(a,i0)') 'There was an error reading the namelist case_config in the file '&
                      // experiment_namelist // '; error code = ',ioerror
    STOP
  end if

  !The current implementation of GFS physics does not support more than one column, since radiation sub schemes use
  !internal module variables. This means that one cannot specify different ways to treat O3, CO2 etc., and also that
  !the code crashes in GFS_initialize_scm_run and later in radiation_gases.f, because it tries to allocate module
  !variables that are already allocated. For now, throw an error and abort.
  if (n_columns>1) then
    write(*,'(a)') 'The current implementation does not allow to run more than one column at a time.'
    STOP
  end if
  
  !read in the physics suite and namelist
  read(10, NML=physics_config, iostat=ioerror)
  close(10)

  select case(time_scheme)
    case(1)
      n_time_levels = 1
    case(2)
      n_time_levels = 2
    case default
      n_time_levels = 2
  end select
  
  call get_tracers(tracer_names, tracer_types)
  
  call scm_state%create(n_columns, n_levels, n_soil, n_snow, n_time_levels, tracer_names, tracer_types)

  scm_state%experiment_name = experiment_name
  scm_state%npz_type = npz_type
  scm_state%vert_coord_file = vert_coord_file
  scm_state%output_dir = output_dir
  scm_state%output_file = output_file
  scm_state%case_name = case_name
  scm_state%physics_suite_name = physics_suite
  scm_state%physics_nml = physics_nml
  scm_state%area(:) = column_area

  scm_state%n_cols = n_columns
  scm_state%n_levels = n_levels
  scm_state%n_time_levels = n_time_levels
  scm_state%dt = dt
  scm_state%n_itt_out = n_itt_out
  scm_state%n_itt_diag = n_itt_diag
  scm_state%runtime = runtime
  scm_state%time_scheme = time_scheme
  scm_state%init_year = year
  scm_state%init_month = month
  scm_state%init_day = day
  scm_state%init_hour = hour
  scm_state%init_min = min

  scm_state%output_period = n_itt_out*dt
  scm_state%thermo_forcing_type = thermo_forcing_type
  scm_state%mom_forcing_type = mom_forcing_type
  scm_state%C_RES            = C_RES
  scm_state%spinup_timesteps = spinup_timesteps
  scm_state%sfc_flux_spec = sfc_flux_spec
  scm_state%sfc_roughness_length_cm(:) = sfc_roughness_length_cm
  scm_state%sfc_type = REAL(sfc_type, kind=dp)
  scm_state%model_ics = model_ics
  scm_state%lsm_ics = lsm_ics
  scm_state%do_spinup = do_spinup
  scm_state%reference_profile_choice = reference_profile_choice
  scm_state%relax_time = relax_time
  scm_state%input_type = input_type
  
  deallocate(tracer_names)
!> @}
end subroutine get_config_nml


!> Subroutine to read the netCDF file containing case initialization and forcing. The forcing files (netCDF4) should be located in the
!! "processed_case_input" directory.
subroutine get_case_init(scm_state, scm_input)
  use scm_type_defs, only : scm_state_type, scm_input_type
  use NetCDF_read, only: NetCDF_read_var, check, missing_value
  type(scm_state_type), intent(in) :: scm_state
  type(scm_input_type), target, intent(inout) :: scm_input
  
  integer                     :: input_nlev !< number of levels in the input file
  integer                     :: input_nsoil !< number of soil levels in the input file
  integer                     :: input_nsnow !< number of snow levels in the input file
  integer                     :: input_nice !< number of sea ice levels in the input file
  integer                     :: input_ntimes !< number of times represented in the input file
  integer                     :: input_nsoil_plus_nsnow !< number of combined snow and soil levels in the input file
  real(kind=dp)               :: input_lat !< column latitude (deg)
  real(kind=dp)               :: input_lon !< column longitude (deg)
  real(kind=dp)               :: input_area    !< surface area [m^2]
  
  integer                     :: input_vegsrc !< vegetation source
  
  real(kind=dp)               :: input_slmsk   !< sea land ice mask [0,1,2]
  real(kind=dp)               :: input_tsfco !< input sea surface temperature OR surface skin temperature over land OR surface skin temperature over ice (depending on slmsk) (K)
  real(kind=dp)               :: input_weasd !< water equivalent accumulated snow depth (mm)
  real(kind=dp)               :: input_tg3     !< deep soil temperature (K)
  real(kind=dp)               :: input_zorl  !< composite surface roughness length (cm)
  real(kind=dp)               :: input_alvsf !< 60 degree vis albedo with strong cosz dependency
  real(kind=dp)               :: input_alvwf !< 60 degree vis albedo with weak cosz dependency
  real(kind=dp)               :: input_alnsf !< 60 degree nir albedo with strong cosz dependency
  real(kind=dp)               :: input_alnwf !< 60 degree nir albedo with weak cosz dependency
  real(kind=dp)               :: input_facsf !< fractional coverage with strong cosz dependency
  real(kind=dp)               :: input_facwf !< fractional coverage with weak cosz dependency
  real(kind=dp)               :: input_vegfrac  !< vegetation fraction
  real(kind=dp)               :: input_canopy  !< amount of water stored in canopy (kg m-2)
  real(kind=dp)               :: input_f10m  !< ratio of sigma level 1 wind and 10m wind
  real(kind=dp)               :: input_t2m    !< 2-meter absolute temperature (K)
  real(kind=dp)               :: input_q2m    !< 2-meter specific humidity (kg kg-1)
  integer                     :: input_vegtyp !< vegetation type
  integer                     :: input_soiltyp!< soil type
  real(kind=dp)               :: input_uustar  !< surface friction velocity (m s-1)
  real(kind=dp)               :: input_ffmm    !< Monin-Obukhov similarity function for momentum
  real(kind=dp)               :: input_ffhh    !< Monin-Obukhov similarity function for heat
  real(kind=dp)               :: input_hice    !< sea ice thickness (m)
  real(kind=dp)               :: input_fice    !< ice fraction (frac)
  real(kind=dp)               :: input_tisfc   !< ice surface temperature (K)
  real(kind=dp)               :: input_tprcp   !< instantaneous total precipitation amount (m)
  real(kind=dp)               :: input_srflag  !< snow/rain flag for precipitation
  real(kind=dp)               :: input_snwdph  !< water equivalent snow depth (mm)
  real(kind=dp)               :: input_shdmin  !< minimun vegetation fraction
  real(kind=dp)               :: input_shdmax  !< maximun vegetation fraction
  integer                     :: input_slopetype !< slope type
  real(kind=dp)               :: input_snoalb  !< maximum snow albedo (frac)
  real(kind=dp)               :: input_sncovr  !< snow area fraction (frac)
  real(kind=dp)               :: input_snodl   !< snowd on land portion of cell
  real(kind=dp)               :: input_weasdl  !< weasd on land portion of cell
  real(kind=dp)               :: input_tsfc    !< tsfc composite
  real(kind=dp)               :: input_tsfcl   !< surface skin temperature over land (K)
  real(kind=dp)               :: input_zorlw    !< surfce roughness length over water [cm]
  real(kind=dp)               :: input_zorll   !< surface roughness length over land (cm)
  real(kind=dp)               :: input_zorli   !< surface roughness length over ice (cm)
  real(kind=dp)               :: input_albdirvis_lnd !<
  real(kind=dp)               :: input_albdirnir_lnd !<
  real(kind=dp)               :: input_albdifvis_lnd !<
  real(kind=dp)               :: input_albdifnir_lnd !<
  real(kind=dp)               :: input_emis_lnd   !<
  real(kind=dp)               :: input_albdirvis_ice !<
  real(kind=dp)               :: input_albdirnir_ice !<
  real(kind=dp)               :: input_albdifvis_ice !<
  real(kind=dp)               :: input_albdifnir_ice !<
  real(kind=dp)               :: input_zorlwav   !< surface roughness length from wave model (cm)
  
  real(kind=dp), allocatable  :: input_stc(:) !< soil temperature (K)
  real(kind=dp), allocatable  :: input_smc(:) !< total soil moisture content (fraction)  
  real(kind=dp), allocatable  :: input_slc(:) !< liquid soil moisture content (fraction)
  real(kind=dp), allocatable  :: input_tiice(:)   !< sea ice internal temperature (K)

  real(kind=dp)               :: input_stddev !< standard deviation of subgrid orography (m)
  real(kind=dp)               :: input_convexity !< convexity of subgrid orography 
  real(kind=dp)               :: input_ol1 !< fraction of grid box with subgrid orography higher than critical height 1
  real(kind=dp)               :: input_ol2 !< fraction of grid box with subgrid orography higher than critical height 2
  real(kind=dp)               :: input_ol3 !< fraction of grid box with subgrid orography higher than critical height 3
  real(kind=dp)               :: input_ol4 !< fraction of grid box with subgrid orography higher than critical height 4
  real(kind=dp)               :: input_oa1 !< assymetry of subgrid orography 1
  real(kind=dp)               :: input_oa2 !< assymetry of subgrid orography 2
  real(kind=dp)               :: input_oa3 !< assymetry of subgrid orography 3
  real(kind=dp)               :: input_oa4 !< assymetry of subgrid orography 4
  real(kind=dp)               :: input_sigma !< slope of subgrid orography
  real(kind=dp)               :: input_theta !< angle with respect to east of maximum subgrid orographic variations (deg)
  real(kind=dp)               :: input_gamma !< anisotropy of subgrid orography
  real(kind=dp)               :: input_elvmax!< maximum of subgrid orography (m)
  real(kind=dp)               :: input_oro !< orography (m)
  real(kind=dp)               :: input_oro_uf !< unfiltered orography (m)
  real(kind=dp)               :: input_landfrac !< fraction of horizontal grid area occupied by land
  real(kind=dp)               :: input_lakefrac !< fraction of horizontal grid area occupied by lake
  real(kind=dp)               :: input_lakedepth !< lake depth (m)

  real(kind=dp)               :: input_tvxy !< vegetation temperature (K)
  real(kind=dp)               :: input_tgxy !< ground temperature for Noahmp (K)
  real(kind=dp)               :: input_tahxy !< canopy air temperature (K)
  real(kind=dp)               :: input_canicexy !< canopy intercepted ice mass (mm)
  real(kind=dp)               :: input_canliqxy !< canopy intercepted liquid water (mm)
  real(kind=dp)               :: input_eahxy !< canopy air vapor pressure (Pa)
  real(kind=dp)               :: input_cmxy !< surface drag coefficient for momentum for noahmp
  real(kind=dp)               :: input_chxy !< surface exchange coeff heat & moisture for noahmp
  real(kind=dp)               :: input_fwetxy !< area fraction of canopy that is wetted/snowed
  real(kind=dp)               :: input_sneqvoxy !< snow mass at previous time step (mm)
  real(kind=dp)               :: input_alboldxy !< snow albedo at previous time step (frac)
  real(kind=dp)               :: input_qsnowxy !< snow precipitation rate at surface (mm s-1)
  real(kind=dp)               :: input_wslakexy !< lake water storage (mm)
  real(kind=dp)               :: input_taussxy !< non-dimensional snow age
  real(kind=dp)               :: input_waxy !< water storage in aquifer (mm)
  real(kind=dp)               :: input_wtxy !< water storage in aquifer and saturated soil (mm)
  real(kind=dp)               :: input_zwtxy !< water table depth (m)
  real(kind=dp)               :: input_xlaixy !< leaf area index
  real(kind=dp)               :: input_xsaixy !< stem area index
  real(kind=dp)               :: input_lfmassxy !< leaf mass (g m-2)
  real(kind=dp)               :: input_stmassxy !< stem mass (g m-2)
  real(kind=dp)               :: input_rtmassxy !< fine root mass (g m-2)
  real(kind=dp)               :: input_woodxy !< wood mass including woody roots (g m-2)
  real(kind=dp)               :: input_stblcpxy !< stable carbon in deep soil (g m-2)
  real(kind=dp)               :: input_fastcpxy !< short-lived carbon in shallow soil (g m-2)
  real(kind=dp)               :: input_smcwtdxy !< soil water content between the bottom of the soil and the water table (m3 m-3)
  real(kind=dp)               :: input_deeprechxy !< recharge to or from the water table when deep (m)
  real(kind=dp)               :: input_rechxy !< recharge to or from the water table when shallow (m)
  real(kind=dp)               :: input_snowxy !< number of snow layers
  
  real(kind=dp), allocatable  :: input_snicexy(:) !< snow layer ice (mm)
  real(kind=dp), allocatable  :: input_snliqxy(:) !< snow layer liquid (mm)
  real(kind=dp), allocatable  :: input_tsnoxy(:) !< snow temperature (K)
  real(kind=dp), allocatable  :: input_smoiseq(:) !< equilibrium soil water content (m3 m-3)
  real(kind=dp), allocatable  :: input_zsnsoxy(:) !< layer bottom depth from snow surface (m)
  
  real(kind=dp)               :: input_tref !< sea surface reference temperature for NSST (K)
  real(kind=dp)               :: input_z_c !< sub-layer cooling thickness for NSST (m)
  real(kind=dp)               :: input_c_0 !< coefficient 1 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp)               :: input_c_d !< coefficient 2 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp)               :: input_w_0 !< coefficient 3 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp)               :: input_w_d !< coefficient 4 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp)               :: input_xt !< heat content in diurnal thermocline layer for NSST (K m)
  real(kind=dp)               :: input_xs !< salinity content in diurnal thermocline layer for NSST (ppt m)
  real(kind=dp)               :: input_xu !< u-current in diurnal thermocline layer for NSST (m2 s-1)
  real(kind=dp)               :: input_xv !< v-current in diurnal thermocline layer for NSST (m2 s-1)
  real(kind=dp)               :: input_xz !< thickness of diurnal thermocline layer for NSST (m)
  real(kind=dp)               :: input_zm !< thickness of ocean mixed layer for NSST (m)
  real(kind=dp)               :: input_xtts !< sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST (m)
  real(kind=dp)               :: input_xzts !< sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST (m K-1)
  real(kind=dp)               :: input_d_conv !< thickness of free convection layer for NSST (m)
  real(kind=dp)               :: input_ifd !< index to start DTM run for NSST
  real(kind=dp)               :: input_dt_cool !< sub-layer cooling amount for NSST (K)
  real(kind=dp)               :: input_qrain !< sensible heat due to rainfall for NSST (W)
  
  real(kind=dp)               :: input_wetness !< normalized soil wetness for RUC LSM
  real(kind=dp)               :: input_clw_surf_land !< cloud condensed water mixing ratio at surface over land for RUC LSM (kg kg-1)
  real(kind=dp)               :: input_clw_surf_ice !< cloud condensed water mixing ratio at surface over ice for RUC LSM (kg kg-1)
  real(kind=dp)               :: input_qwv_surf_land !< water vapor mixing ratio at surface over land for RUC LSM (kg kg-1)
  real(kind=dp)               :: input_qwv_surf_ice !< water vapor mixing ratio at surface over ice for RUC LSM (kg kg-1)
  real(kind=dp)               :: input_tsnow_land !< snow temperature at the bottom of the first snow layer over land for RUC LSM (K)
  real(kind=dp)               :: input_tsnow_ice !< snow temperature at the bottom of the first snow layer over ice for RUC LSM (K)
  real(kind=dp)               :: input_snowfallac_land !< run-total snow accumulation on the ground over land for RUC LSM (kg m-2)
  real(kind=dp)               :: input_snowfallac_ice !< run-total snow accumulation on the ground over ice for RUC LSM (kg m-2)
  real(kind=dp)               :: input_sncovr_ice !<
  real(kind=dp)               :: input_sfalb_lnd
  real(kind=dp)               :: input_sfalb_lnd_bck
  real(kind=dp)               :: input_sfalb_ice
  real(kind=dp)               :: input_emis_ice
  real(kind=dp)               :: input_lai !< leaf area index for RUC LSM
  
  real(kind=dp), allocatable  :: input_tslb(:)    !< soil temperature for RUC LSM (K)
  real(kind=dp), allocatable  :: input_smois(:)   !< volume fraction of soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_sh2o(:)    !< volume fraction of unfrozen soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_smfr(:)    !< volume fraction of frozen soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_flfr(:)    !< flag for frozen soil physics
  
  ! dimension variables
  !real(kind=dp), allocatable  :: input_pres_i(:) !< interface pressures
  !real(kind=dp), allocatable  :: input_pres_l(:) !< layer pressures
  real(kind=dp), allocatable  :: input_pres(:) !< input file pressure levels (Pa)
  real(kind=dp), allocatable  :: input_time(:) !< input file times (seconds since the beginning of the case)

  !initial profile variables
  real(kind=dp), allocatable  :: input_thetail(:) !< ice-liquid water potential temperature profile (K)
  real(kind=dp), allocatable  :: input_temp(:) !< temperature profile (K)
  real(kind=dp), allocatable  :: input_qt(:) !< total water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_ql(:) !< liquid water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_qi(:) !< ice water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_u(:) !< east-west horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_v(:) !< north-south horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_tke(:) !< TKE profile (m^2 s^-2)
  real(kind=dp), allocatable  :: input_ozone(:) !< ozone profile (kg kg^-1)
  
  real(kind=dp), allocatable  :: input_pres_surf(:) !< time-series of surface pressure (Pa)
  real(kind=dp), allocatable  :: input_T_surf(:) !< time-series of surface temperature (K)

  real(kind=dp), allocatable  :: input_w_ls(:,:) !< 2D vertical velocity (m s^-1)
  real(kind=dp), allocatable  :: input_omega(:,:) !< 2D pressure vertical velocity (Pa s^-1)
  real(kind=dp), allocatable  :: input_u_g(:,:) !< 2D geostrophic east-west wind (m s^-1)
  real(kind=dp), allocatable  :: input_v_g(:,:) !< 2D geostrophic north-south wind (m s^-1)
  real(kind=dp), allocatable  :: input_u_nudge(:,:) !< 2D nudging east-west wind (m s^-1)
  real(kind=dp), allocatable  :: input_v_nudge(:,:) !< 2D nudging north-south wind (m s^-1)
  real(kind=dp), allocatable  :: input_T_nudge(:,:) !< 2D nudging abs. temperature (K)
  real(kind=dp), allocatable  :: input_thil_nudge(:,:) !< 2D nudging liq. pot. temperature (K)
  real(kind=dp), allocatable  :: input_qt_nudge(:,:) !< 2D nudging specific humidity (kg kg^-1)
  real(kind=dp), allocatable  :: input_dT_dt_rad(:,:) !< 2D radiative heating rate (K s^-1)
  real(kind=dp), allocatable  :: input_h_advec_thetail(:,:) !< 2D theta_il tendency due to large-scale horizontal advection (K s^-1)
  real(kind=dp), allocatable  :: input_h_advec_qt(:,:) !< 2D q_t tendency due to large-scale horizontal advection (kg kg^-1 s^-1)
  real(kind=dp), allocatable  :: input_v_advec_thetail(:,:) !< 2D theta_il tendency due to large-scale vertical advection (K s^-1)
  real(kind=dp), allocatable  :: input_v_advec_qt(:,:) !< 2D q_t tendency due to large-scale horizontal vertical (kg kg^-1 s^-1)

  real(kind=dp), allocatable  :: input_sh_flux_sfc(:) !< time-series of surface sensible heat flux (K m s^-1)
  real(kind=dp), allocatable  :: input_lh_flux_sfc(:) !< time-series of surface latent heat flux (kg kg^-1 m s^-1)
  
  CHARACTER(LEN=nf90_max_name)      :: tmpName
  integer                           :: ncid, varID, grp_ncid, allocate_status,ierr
  real(kind=dp)                     :: nc_missing_value

  !>  \section get_case_init_alg Algorithm
  !!  @{

  !> - Open the case input file found in the processed_case_input dir corresponding to the experiment name.
  call check(NF90_OPEN(trim(adjustl(scm_state%case_name))//'.nc',nf90_nowrite,ncid),"nf90_open()")
  
  !> - Read in missing value from file (replace module variable if present)
  ierr = NF90_GET_ATT(ncid, NF90_GLOBAL, 'missing_value', nc_missing_value)
  if(ierr == NF90_NOERR) then
    missing_value = nc_missing_value
  end if
  
  !> - Get the dimensions (global group).
  
  !required dimensions
  call check(NF90_INQ_DIMID(ncid,"levels",varID),"nf90_inq_dimid(levels)")
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nlev),"nf90_inq_dim(levels)")
  call check(NF90_INQ_DIMID(ncid,"time",varID),"inq_dimid(time)")
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_ntimes),"nf90_inq_dim(time)")
  
  !possible dimensions (if using model ICs)
  ierr = NF90_INQ_DIMID(ncid,"nsoil",varID)
  if(ierr /= NF90_NOERR) then
    input_nsoil = missing_soil_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nsoil),"nf90_inq_dim(nsoil)")
  end if
  ierr = NF90_INQ_DIMID(ncid,"nsnow",varID)
  if(ierr /= NF90_NOERR) then
    input_nsnow = missing_snow_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nsnow),"nf90_inq_dim(nsnow)")
  end if
  ierr = NF90_INQ_DIMID(ncid,"nice",varID)
  if(ierr /= NF90_NOERR) then
    input_nice = missing_ice_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nice),"nf90_inq_dim(nice)")
  end if  
  
  !> - Allocate the dimension variables.
  allocate(input_pres(input_nlev),input_time(input_ntimes), stat=allocate_status)

  !> - Read in the dimension variables (required).
  call NetCDF_read_var(ncid, "levels", .True., input_pres)
  call NetCDF_read_var(ncid, "time", .True., input_time)

  !> - Read in the initial conditions.

  !>  - Find group ncid for initial group.
  call check(NF90_INQ_GRP_NCID(ncid,"initial",grp_ncid),"nf90_inq_grp_ncid(initial)")

  !>  - Allocate the initial profiles (required). One of thetail or temp is required.
  allocate(input_thetail(input_nlev), input_temp(input_nlev), input_qt(input_nlev), input_ql(input_nlev), input_qi(input_nlev), &
    input_u(input_nlev), input_v(input_nlev), input_tke(input_nlev), input_ozone(input_nlev), stat=allocate_status)
  
  !>  - Read in the initial profiles. The variable names in all input files are expected to be identical.
  
  !Either thetail or T must be present
  call NetCDF_read_var(grp_ncid, "thetail", .False., input_thetail)
  call NetCDF_read_var(grp_ncid, "temp", .False., input_temp)
  if (maxval(input_thetail) < 0 .and. maxval(input_temp) < 0) then
    write(*,*) "One of thetail or temp variables must be present in ",trim(adjustl(scm_state%case_name))//'.nc',". Stopping..."
    STOP
  end if
  call NetCDF_read_var(grp_ncid, "qt",    .True., input_qt   )
  call NetCDF_read_var(grp_ncid, "ql",    .True., input_ql   )
  call NetCDF_read_var(grp_ncid, "qi",    .True., input_qi   )
  call NetCDF_read_var(grp_ncid, "u",     .True., input_u    )
  call NetCDF_read_var(grp_ncid, "v",     .True., input_v    )
  call NetCDF_read_var(grp_ncid, "tke",   .True., input_tke  )
  call NetCDF_read_var(grp_ncid, "ozone", .True., input_ozone)
  
  !possible initial profiles
  !needed for Noah LSM and others (when running with model ICs)
  allocate(input_stc(input_nsoil), input_smc(input_nsoil), input_slc(input_nsoil), &
           stat=allocate_status)
  call NetCDF_read_var(grp_ncid, "stc", .False., input_stc)
  call NetCDF_read_var(grp_ncid, "smc", .False., input_smc)
  call NetCDF_read_var(grp_ncid, "slc", .False., input_slc)
  
  !needed for NoahMP LSM (when running with model ICs)
  allocate(input_snicexy(input_nsnow), input_snliqxy(input_nsnow), input_tsnoxy(input_nsnow), &
     input_smoiseq(input_nsoil), input_zsnsoxy(input_nsnow + input_nsoil))
  call NetCDF_read_var(grp_ncid, "snicexy", .False., input_snicexy)
  call NetCDF_read_var(grp_ncid, "snliqxy", .False., input_snliqxy)
  call NetCDF_read_var(grp_ncid, "tsnoxy",  .False., input_tsnoxy )
  call NetCDF_read_var(grp_ncid, "smoiseq", .False., input_smoiseq)
  call NetCDF_read_var(grp_ncid, "zsnsoxy", .False., input_zsnsoxy)
  
  !needed for fractional grid (when running with model ICs)
  allocate(input_tiice(input_nice))
  call NetCDF_read_var(grp_ncid, "tiice", .False., input_tiice)
  
  !needed for RUC LSM (when running with model ICs)
  allocate(input_tslb(input_nsoil), input_smois(input_nsoil), input_sh2o(input_nsoil), &
      input_smfr(input_nsoil), input_flfr(input_nsoil))
  call NetCDF_read_var(grp_ncid, "tslb",  .False., input_tslb )
  call NetCDF_read_var(grp_ncid, "smois", .False., input_smois)
  call NetCDF_read_var(grp_ncid, "sh2o",  .False., input_sh2o )
  call NetCDF_read_var(grp_ncid, "smfr",  .False., input_smfr )
  call NetCDF_read_var(grp_ncid, "flfr",  .False., input_flfr )
  
  !>  - Find group ncid for scalar group.
  call check(NF90_INQ_GRP_NCID(ncid,"scalars",grp_ncid),"nf90_inq_grp_ncid(scalars)")
  
  !required
  call NetCDF_read_var(grp_ncid, "lat", .True., input_lat)
  call NetCDF_read_var(grp_ncid, "lon", .True., input_lon)
  !time data in file ignored?
  call NetCDF_read_var(grp_ncid, "area", .False., input_area)
  
  !possible scalars
  !Noah LSM parameters (when running with model ICs)
  call NetCDF_read_var(grp_ncid, "vegsrc",  .False., input_vegsrc   )
  call NetCDF_read_var(grp_ncid, "slmsk",   .False., input_slmsk)
  call NetCDF_read_var(grp_ncid, "tsfco",   .False., input_tsfco)
  call NetCDF_read_var(grp_ncid, "weasd",   .False., input_weasd)
  call NetCDF_read_var(grp_ncid, "tg3",     .False., input_tg3)
  call NetCDF_read_var(grp_ncid, "zorl",    .False., input_zorl)
  call NetCDF_read_var(grp_ncid, "alvsf",   .False., input_alvsf)
  call NetCDF_read_var(grp_ncid, "alvwf",   .False., input_alvwf)
  call NetCDF_read_var(grp_ncid, "alnsf",   .False., input_alnsf)
  call NetCDF_read_var(grp_ncid, "alnwf",   .False., input_alnwf)
  call NetCDF_read_var(grp_ncid, "facsf",   .False., input_facsf)
  call NetCDF_read_var(grp_ncid, "facwf",   .False., input_facwf)
  call NetCDF_read_var(grp_ncid, "vegfrac", .False., input_vegfrac)
  call NetCDF_read_var(grp_ncid, "canopy",  .False., input_canopy)
  call NetCDF_read_var(grp_ncid, "f10m",    .False., input_f10m)
  call NetCDF_read_var(grp_ncid, "t2m",     .False., input_t2m)
  call NetCDF_read_var(grp_ncid, "q2m",     .False., input_q2m)
  call NetCDF_read_var(grp_ncid, "vegtyp",  .False., input_vegtyp   )
  call NetCDF_read_var(grp_ncid, "soiltyp", .False., input_soiltyp  )
  call NetCDF_read_var(grp_ncid, "uustar",  .False., input_uustar)
  call NetCDF_read_var(grp_ncid, "ffmm",    .False., input_ffmm)
  call NetCDF_read_var(grp_ncid, "ffhh",    .False., input_ffhh)
  call NetCDF_read_var(grp_ncid, "hice",    .False., input_hice)
  call NetCDF_read_var(grp_ncid, "fice",    .False., input_fice)
  call NetCDF_read_var(grp_ncid, "tisfc",   .False., input_tisfc)
  call NetCDF_read_var(grp_ncid, "tprcp",   .False., input_tprcp)
  call NetCDF_read_var(grp_ncid, "srflag",  .False., input_srflag)
  call NetCDF_read_var(grp_ncid, "snwdph",  .False., input_snwdph)
  call NetCDF_read_var(grp_ncid, "shdmin",  .False., input_shdmin)
  call NetCDF_read_var(grp_ncid, "shdmax",  .False., input_shdmax)
  call NetCDF_read_var(grp_ncid, "slopetyp",.False., input_slopetype)
  call NetCDF_read_var(grp_ncid, "snoalb",  .False., input_snoalb)
  call NetCDF_read_var(grp_ncid, "sncovr",  .False., input_sncovr)
  call NetCDF_read_var(grp_ncid, "snodl",   .False., input_snodl)
  call NetCDF_read_var(grp_ncid, "weasdl",  .False., input_weasdl)
  call NetCDF_read_var(grp_ncid, "tsfc",    .False., input_tsfc)
  call NetCDF_read_var(grp_ncid, "tsfcl",   .False., input_tsfcl)
  call NetCDF_read_var(grp_ncid, "zorlw",   .False., input_zorlw)
  call NetCDF_read_var(grp_ncid, "zorll",   .False., input_zorll)
  call NetCDF_read_var(grp_ncid, "zorli",   .False., input_zorli)
  call NetCDF_read_var(grp_ncid, "albdirvis_lnd", .False., input_albdirvis_lnd)
  call NetCDF_read_var(grp_ncid, "albdirnir_lnd", .False., input_albdirnir_lnd)
  call NetCDF_read_var(grp_ncid, "albdifvis_lnd", .False., input_albdifvis_lnd)
  call NetCDF_read_var(grp_ncid, "albdifnir_lnd", .False., input_albdifnir_lnd)
  call NetCDF_read_var(grp_ncid, "emis_lnd",      .False., input_emis_lnd)
  call NetCDF_read_var(grp_ncid, "albdirvis_ice", .False., input_albdirvis_ice)
  call NetCDF_read_var(grp_ncid, "albdirnir_ice", .False., input_albdirnir_ice)
  call NetCDF_read_var(grp_ncid, "albdifvis_ice", .False., input_albdifvis_ice)
  call NetCDF_read_var(grp_ncid, "albdifnir_ice", .False., input_albdifnir_ice)
  call NetCDF_read_var(grp_ncid, "zorlwav", .False., input_zorlwav)
  
  !orographic parameters
  call NetCDF_read_var(grp_ncid, "stddev",    .False., input_stddev)
  call NetCDF_read_var(grp_ncid, "convexity", .False., input_convexity)
  call NetCDF_read_var(grp_ncid, "oa1",       .False., input_oa1)
  call NetCDF_read_var(grp_ncid, "oa2",       .False., input_oa2)
  call NetCDF_read_var(grp_ncid, "oa3",       .False., input_oa3)
  call NetCDF_read_var(grp_ncid, "oa4",       .False., input_oa4)
  call NetCDF_read_var(grp_ncid, "ol1",       .False., input_ol1)
  call NetCDF_read_var(grp_ncid, "ol2",       .False., input_ol2)
  call NetCDF_read_var(grp_ncid, "ol3",       .False., input_ol3)
  call NetCDF_read_var(grp_ncid, "ol4",       .False., input_ol4)
  call NetCDF_read_var(grp_ncid, "theta",     .False., input_theta)
  call NetCDF_read_var(grp_ncid, "gamma",     .False., input_gamma)
  call NetCDF_read_var(grp_ncid, "sigma",     .False., input_sigma)
  call NetCDF_read_var(grp_ncid, "elvmax",    .False., input_elvmax)
  call NetCDF_read_var(grp_ncid, "oro",       .False., input_oro)
  call NetCDF_read_var(grp_ncid, "oro_uf",    .False., input_oro_uf)
  call NetCDF_read_var(grp_ncid, "landfrac",  .False., input_landfrac)
  call NetCDF_read_var(grp_ncid, "lakefrac",  .False., input_lakefrac)
  call NetCDF_read_var(grp_ncid, "lakedepth", .False., input_lakedepth)
  
  !NoahMP parameters
  call NetCDF_read_var(grp_ncid, "tvxy",      .False., input_tvxy)
  call NetCDF_read_var(grp_ncid, "tgxy",      .False., input_tgxy)
  call NetCDF_read_var(grp_ncid, "tahxy",     .False., input_tahxy)
  call NetCDF_read_var(grp_ncid, "canicexy",  .False., input_canicexy)
  call NetCDF_read_var(grp_ncid, "canliqxy",  .False., input_canliqxy)
  call NetCDF_read_var(grp_ncid, "eahxy",     .False., input_eahxy)
  call NetCDF_read_var(grp_ncid, "cmxy",      .False., input_cmxy)
  call NetCDF_read_var(grp_ncid, "chxy",      .False., input_chxy)
  call NetCDF_read_var(grp_ncid, "fwetxy",    .False., input_fwetxy)
  call NetCDF_read_var(grp_ncid, "sneqvoxy",  .False., input_sneqvoxy)
  call NetCDF_read_var(grp_ncid, "alboldxy",  .False., input_alboldxy)
  call NetCDF_read_var(grp_ncid, "qsnowxy",   .False., input_qsnowxy)
  call NetCDF_read_var(grp_ncid, "wslakexy",  .False., input_wslakexy)
  call NetCDF_read_var(grp_ncid, "taussxy",   .False., input_taussxy)
  call NetCDF_read_var(grp_ncid, "waxy",      .False., input_waxy)
  call NetCDF_read_var(grp_ncid, "wtxy",      .False., input_wtxy)
  call NetCDF_read_var(grp_ncid, "zwtxy",     .False., input_zwtxy)
  call NetCDF_read_var(grp_ncid, "xlaixy",    .False., input_xlaixy)
  call NetCDF_read_var(grp_ncid, "xsaixy",    .False., input_xsaixy)
  call NetCDF_read_var(grp_ncid, "lfmassxy",  .False., input_lfmassxy)
  call NetCDF_read_var(grp_ncid, "stmassxy",  .False., input_stmassxy)
  call NetCDF_read_var(grp_ncid, "rtmassxy",  .False., input_rtmassxy)
  call NetCDF_read_var(grp_ncid, "woodxy",    .False., input_woodxy)
  call NetCDF_read_var(grp_ncid, "stblcpxy",  .False., input_stblcpxy)
  call NetCDF_read_var(grp_ncid, "fastcpxy",  .False., input_fastcpxy)
  call NetCDF_read_var(grp_ncid, "smcwtdxy",  .False., input_smcwtdxy)
  call NetCDF_read_var(grp_ncid, "deeprechxy",.False., input_deeprechxy)
  call NetCDF_read_var(grp_ncid, "rechxy",    .False., input_rechxy)
  call NetCDF_read_var(grp_ncid, "snowxy",    .False., input_snowxy)
  
  !NSST variables
  call NetCDF_read_var(grp_ncid, "tref",    .False., input_tref)
  call NetCDF_read_var(grp_ncid, "z_c",     .False., input_z_c)
  call NetCDF_read_var(grp_ncid, "c_0",     .False., input_c_0)
  call NetCDF_read_var(grp_ncid, "c_d",     .False., input_c_d)
  call NetCDF_read_var(grp_ncid, "w_0",     .False., input_w_0)
  call NetCDF_read_var(grp_ncid, "w_d",     .False., input_w_d)
  call NetCDF_read_var(grp_ncid, "xt",      .False., input_xt)
  call NetCDF_read_var(grp_ncid, "xs",      .False., input_xs)
  call NetCDF_read_var(grp_ncid, "xu",      .False., input_xu)
  call NetCDF_read_var(grp_ncid, "xv",      .False., input_xv)
  call NetCDF_read_var(grp_ncid, "xz",      .False., input_xz)
  call NetCDF_read_var(grp_ncid, "zm",      .False., input_zm)
  call NetCDF_read_var(grp_ncid, "xtts",    .False., input_xtts)
  call NetCDF_read_var(grp_ncid, "xzts",    .False., input_xzts)
  call NetCDF_read_var(grp_ncid, "d_conv",  .False., input_d_conv)
  call NetCDF_read_var(grp_ncid, "ifd",     .False., input_ifd)
  call NetCDF_read_var(grp_ncid, "dt_cool", .False., input_dt_cool)
  call NetCDF_read_var(grp_ncid, "qrain",   .False., input_qrain)
  
  !RUC LSM variables
  call NetCDF_read_var(grp_ncid, "wetness",          .False., input_wetness)
  call NetCDF_read_var(grp_ncid, "clw_surf_land",    .False., input_clw_surf_land)
  call NetCDF_read_var(grp_ncid, "clw_surf_ice",     .False., input_clw_surf_ice)
  call NetCDF_read_var(grp_ncid, "qwv_surf_land",    .False., input_qwv_surf_land)
  call NetCDF_read_var(grp_ncid, "qwv_surf_ice",     .False., input_qwv_surf_ice)
  call NetCDF_read_var(grp_ncid, "tsnow_land",       .False., input_tsnow_land)
  call NetCDF_read_var(grp_ncid, "tsnow_ice",        .False., input_tsnow_ice)
  call NetCDF_read_var(grp_ncid, "snowfallac_land",  .False., input_snowfallac_land)
  call NetCDF_read_var(grp_ncid, "snowfallac_ice",   .False., input_snowfallac_ice)
  call NetCDF_read_var(grp_ncid, "sncovr_ice",       .False., input_sncovr_ice)
  call NetCDF_read_var(grp_ncid, "sfalb_lnd",        .False., input_sfalb_lnd)
  call NetCDF_read_var(grp_ncid, "sfalb_lnd_bck",    .False., input_sfalb_lnd_bck)
  call NetCDF_read_var(grp_ncid, "emis_ice",         .False., input_emis_ice)
  call NetCDF_read_var(grp_ncid, "lai",              .False., input_lai)
  
  !> - Read in the forcing data.

  !>  - Find group ncid for forcing group.
  call check(NF90_INQ_GRP_NCID(ncid,"forcing",grp_ncid),"nf90_inq_grp_ncid(forcing)")

  !>  - (Recall that multidimensional arrays need to be read in with the order of dimensions reversed from the netCDF file).

  !>  - Allocate the time-series and 2D forcing data.
  allocate(input_pres_surf(input_ntimes), input_T_surf(input_ntimes),            &
    input_sh_flux_sfc(input_ntimes), input_lh_flux_sfc(input_ntimes), input_w_ls(input_ntimes, input_nlev), &
    input_omega(input_ntimes, input_nlev), input_u_g(input_ntimes, input_nlev), input_v_g(input_ntimes, input_nlev), &
    input_dT_dt_rad(input_ntimes, input_nlev), input_h_advec_thetail(input_ntimes, input_nlev), &
    input_h_advec_qt(input_ntimes, input_nlev), input_v_advec_thetail(input_ntimes, input_nlev), &
    input_v_advec_qt(input_ntimes, input_nlev), input_u_nudge(input_ntimes, input_nlev), input_v_nudge(input_ntimes, input_nlev),  &
    input_T_nudge(input_ntimes, input_nlev), input_thil_nudge(input_ntimes, input_nlev), input_qt_nudge(input_ntimes, input_nlev), &
    stat=allocate_status)

  !>  - Read in the time-series and 2D forcing data.
  call NetCDF_read_var(grp_ncid, "p_surf", .True., input_pres_surf)
  call NetCDF_read_var(grp_ncid, "T_surf", .True., input_T_surf)
  call NetCDF_read_var(grp_ncid, "sh_flux_sfc", .False., input_sh_flux_sfc)
  call NetCDF_read_var(grp_ncid, "lh_flux_sfc", .False., input_lh_flux_sfc)
  
  call NetCDF_read_var(grp_ncid, "w_ls", .True., input_w_ls)
  call NetCDF_read_var(grp_ncid, "omega", .True., input_omega)
  call NetCDF_read_var(grp_ncid, "u_g", .True., input_u_g)
  call NetCDF_read_var(grp_ncid, "v_g", .True., input_v_g)
  call NetCDF_read_var(grp_ncid, "u_nudge", .True., input_u_nudge)
  call NetCDF_read_var(grp_ncid, "v_nudge", .True., input_v_nudge)
  call NetCDF_read_var(grp_ncid, "T_nudge", .True., input_T_nudge)
  call NetCDF_read_var(grp_ncid, "thil_nudge", .True., input_thil_nudge)
  call NetCDF_read_var(grp_ncid, "qt_nudge", .True., input_qt_nudge)
  call NetCDF_read_var(grp_ncid, "dT_dt_rad", .True., input_dT_dt_rad)
  call NetCDF_read_var(grp_ncid, "h_advec_thetail", .True., input_h_advec_thetail)
  call NetCDF_read_var(grp_ncid, "h_advec_qt", .True., input_h_advec_qt)
  call NetCDF_read_var(grp_ncid, "v_advec_thetail", .True., input_v_advec_thetail)
  call NetCDF_read_var(grp_ncid, "v_advec_qt", .True., input_v_advec_qt)

  call check(NF90_CLOSE(NCID=ncid),"nf90_close()")

  call scm_input%create(input_ntimes, input_nlev, input_nsoil, input_nsnow, input_nice)
    
  ! GJF already done in scm_input%create routine
  !scm_input%input_nlev = input_nlev
  !scm_input%input_ntimes = input_ntimes

  scm_input%input_pres = input_pres
  scm_input%input_time = input_time
  scm_input%input_temp = input_temp
  scm_input%input_thetail = input_thetail
  scm_input%input_qt = input_qt
  scm_input%input_ql = input_ql
  scm_input%input_qi = input_qi
  scm_input%input_u = input_u
  scm_input%input_v = input_v
  scm_input%input_tke = input_tke
  scm_input%input_ozone = input_ozone
  scm_input%input_lat = input_lat
  scm_input%input_lon = input_lon
  scm_input%input_pres_surf = input_pres_surf
  scm_input%input_T_surf = input_T_surf
  scm_input%input_sh_flux_sfc_kin = input_sh_flux_sfc
  scm_input%input_lh_flux_sfc_kin = input_lh_flux_sfc
  scm_input%input_w_ls = input_w_ls
  scm_input%input_omega = input_omega
  scm_input%input_u_g = input_u_g
  scm_input%input_v_g = input_v_g
  scm_input%input_dT_dt_rad = input_dT_dt_rad
  scm_input%input_h_advec_thetail = input_h_advec_thetail
  scm_input%input_h_advec_qt = input_h_advec_qt
  scm_input%input_v_advec_thetail = input_v_advec_thetail
  scm_input%input_v_advec_qt = input_v_advec_qt
  scm_input%input_u_nudge = input_u_nudge
  scm_input%input_v_nudge = input_v_nudge
  scm_input%input_T_nudge = input_T_nudge
  scm_input%input_thil_nudge = input_thil_nudge
  scm_input%input_qt_nudge = input_qt_nudge
  
  scm_input%input_stc   = input_stc  
  scm_input%input_smc   = input_smc  
  scm_input%input_slc   = input_slc  
  
  scm_input%input_snicexy    = input_snicexy
  scm_input%input_snliqxy    = input_snliqxy
  scm_input%input_tsnoxy     = input_tsnoxy
  scm_input%input_smoiseq    = input_smoiseq
  scm_input%input_zsnsoxy    = input_zsnsoxy
  
  scm_input%input_tiice      = input_tiice
  scm_input%input_tslb       = input_tslb
  scm_input%input_smois      = input_smois
  scm_input%input_sh2o       = input_sh2o
  scm_input%input_smfr       = input_smfr
  scm_input%input_flfr       = input_flfr
  
  scm_input%input_vegsrc   = input_vegsrc
  
  scm_input%input_slmsk    = input_slmsk
  scm_input%input_canopy   = input_canopy
  scm_input%input_hice     = input_hice
  scm_input%input_fice     = input_fice
  scm_input%input_tisfc    = input_tisfc
  scm_input%input_snwdph   = input_snwdph
  scm_input%input_snoalb   = input_snoalb
  scm_input%input_sncovr   = input_sncovr
  if (input_area > missing_value) then
    scm_input%input_area   = input_area
  end if
  scm_input%input_tsfco    = input_tsfco
  scm_input%input_weasd    = input_weasd
  scm_input%input_tg3      = input_tg3
  scm_input%input_zorl     = input_zorl
  scm_input%input_alvsf    = input_alvsf
  scm_input%input_alvwf    = input_alvwf
  scm_input%input_alnsf    = input_alnsf
  scm_input%input_alnwf    = input_alnwf
  scm_input%input_facsf    = input_facsf
  scm_input%input_facwf    = input_facwf
  scm_input%input_vegfrac  = input_vegfrac
  scm_input%input_canopy   = input_canopy
  scm_input%input_f10m     = input_f10m
  scm_input%input_t2m      = input_t2m
  scm_input%input_q2m      = input_q2m
  scm_input%input_vegtyp   = input_vegtyp
  scm_input%input_soiltyp  = input_soiltyp
  scm_input%input_uustar   = input_uustar
  scm_input%input_ffmm     = input_ffmm
  scm_input%input_ffhh     = input_ffhh
  scm_input%input_hice     = input_hice
  scm_input%input_fice     = input_fice
  scm_input%input_tisfc    = input_tisfc
  scm_input%input_tprcp    = input_tprcp
  scm_input%input_srflag   = input_srflag
  scm_input%input_snwdph   = input_snwdph
  scm_input%input_shdmin   = input_shdmin
  scm_input%input_shdmax   = input_shdmax
  scm_input%input_slopetype = input_slopetype
  scm_input%input_snoalb   = input_snoalb
  scm_input%input_sncovr   = input_sncovr
  scm_input%input_snodl    = input_snodl
  scm_input%input_weasdl   = input_weasdl
  scm_input%input_tsfc     = input_tsfc
  scm_input%input_tsfcl    = input_tsfcl
  scm_input%input_zorlw    = input_zorlw
  scm_input%input_zorll    = input_zorll
  scm_input%input_zorli    = input_zorli
  scm_input%input_albdirvis_lnd = input_albdirvis_lnd
  scm_input%input_albdirnir_lnd = input_albdirnir_lnd
  scm_input%input_albdifvis_lnd = input_albdifvis_lnd
  scm_input%input_albdifnir_lnd = input_albdifnir_lnd
  scm_input%input_emis_lnd = input_emis_lnd
  scm_input%input_albdirvis_ice = input_albdirvis_ice
  scm_input%input_albdirnir_ice = input_albdirnir_ice
  scm_input%input_albdifvis_ice = input_albdifvis_ice
  scm_input%input_albdifnir_ice = input_albdifnir_ice
  scm_input%input_zorlwav  = input_zorlwav
  
  scm_input%input_stddev   = input_stddev
  scm_input%input_convexity= input_convexity
  scm_input%input_oa1      = input_oa1
  scm_input%input_oa2      = input_oa2
  scm_input%input_oa3      = input_oa3
  scm_input%input_oa4      = input_oa4
  scm_input%input_ol1      = input_ol1
  scm_input%input_ol2      = input_ol2
  scm_input%input_ol3      = input_ol3
  scm_input%input_ol4      = input_ol4
  scm_input%input_sigma    = input_sigma
  scm_input%input_theta    = input_theta
  scm_input%input_gamma    = input_gamma
  scm_input%input_elvmax   = input_elvmax
  scm_input%input_oro      = input_oro
  scm_input%input_oro_uf   = input_oro_uf
  scm_input%input_landfrac = input_landfrac
  scm_input%input_lakefrac = input_lakefrac
  scm_input%input_lakedepth= input_lakedepth
  
  scm_input%input_tvxy = input_tvxy
  scm_input%input_tgxy = input_tgxy
  scm_input%input_tahxy = input_tahxy
  scm_input%input_canicexy = input_canicexy
  scm_input%input_canliqxy = input_canliqxy
  scm_input%input_eahxy = input_eahxy
  scm_input%input_cmxy = input_cmxy
  scm_input%input_chxy = input_chxy
  scm_input%input_fwetxy = input_fwetxy
  scm_input%input_sneqvoxy = input_sneqvoxy
  scm_input%input_alboldxy = input_alboldxy
  scm_input%input_qsnowxy = input_qsnowxy
  scm_input%input_wslakexy = input_wslakexy
  scm_input%input_taussxy = input_taussxy
  scm_input%input_waxy = input_waxy
  scm_input%input_wtxy = input_wtxy
  scm_input%input_zwtxy = input_zwtxy
  scm_input%input_xlaixy = input_xlaixy
  scm_input%input_xsaixy = input_xsaixy
  scm_input%input_lfmassxy = input_lfmassxy
  scm_input%input_stmassxy = input_stmassxy
  scm_input%input_rtmassxy = input_rtmassxy
  scm_input%input_woodxy = input_woodxy
  scm_input%input_stblcpxy = input_stblcpxy
  scm_input%input_fastcpxy = input_fastcpxy
  scm_input%input_smcwtdxy = input_smcwtdxy
  scm_input%input_deeprechxy = input_deeprechxy
  scm_input%input_rechxy = input_rechxy
  scm_input%input_snowxy = input_snowxy
  
  scm_input%input_tref    = input_tref
  scm_input%input_z_c     = input_z_c
  scm_input%input_c_0     = input_c_0
  scm_input%input_c_d     = input_c_d
  scm_input%input_w_0     = input_w_0
  scm_input%input_w_d     = input_w_d
  scm_input%input_xt      = input_xt
  scm_input%input_xs      = input_xs
  scm_input%input_xu      = input_xu
  scm_input%input_xv      = input_xv
  scm_input%input_xz      = input_xz
  scm_input%input_zm      = input_zm
  scm_input%input_xtts    = input_xtts
  scm_input%input_xzts    = input_xzts
  scm_input%input_d_conv  = input_d_conv
  scm_input%input_ifd     = input_ifd
  scm_input%input_dt_cool = input_dt_cool
  scm_input%input_qrain   = input_qrain
  
  scm_input%input_area     = input_area
  
  scm_input%input_wetness         = input_wetness
  scm_input%input_clw_surf_land   = input_clw_surf_land
  scm_input%input_clw_surf_ice    = input_clw_surf_ice
  scm_input%input_qwv_surf_land   = input_qwv_surf_land
  scm_input%input_qwv_surf_ice    = input_qwv_surf_ice
  scm_input%input_tsnow_land      = input_tsnow_land
  scm_input%input_tsnow_ice       = input_tsnow_ice
  scm_input%input_snowfallac_land = input_snowfallac_land
  scm_input%input_snowfallac_ice  = input_snowfallac_ice
  scm_input%input_sncovr_ice      = input_sncovr_ice
  scm_input%input_sfalb_lnd       = input_sfalb_lnd
  scm_input%input_sfalb_lnd_bck   = input_sfalb_lnd_bck
  scm_input%input_emis_ice        = input_emis_ice
  scm_input%input_lai             = input_lai
  
!> @}
end subroutine get_case_init

subroutine get_case_init_DEPHY(scm_state, scm_input)
  !corresponds to the DEPHY-SCM specs, version 1
  
  use scm_type_defs, only : scm_state_type, scm_input_type
  use NetCDF_read, only: NetCDF_read_var, NetCDF_read_att, NetCDF_conditionally_read_var, check, missing_value, missing_value_int
  use scm_physical_constants, only: con_hvap, con_hfus, con_cp, con_rocp, con_rd
  use scm_utils, only: find_vertical_index_pressure, find_vertical_index_height
  
  type(scm_state_type), intent(inout) :: scm_state
  type(scm_input_type), target, intent(inout) :: scm_input
  
  ! dimension variables
  real(kind=dp), allocatable  :: input_t0(:) !< input initialization times (seconds since global attribute "startdate")
  real(kind=dp), allocatable  :: input_time(:) !< input forcing times (seconds since the beginning of the case)
  real(kind=dp), allocatable  :: input_lev(:) !< corresponds to either pressure or height (depending on attribute) - why is this needed when both pressure and height also provided in ICs?
  
  !non-standard dimensions (may or may not exist in the file)
  real(kind=dp), allocatable :: input_soil(:) !< soil depth
  
  ! global attributes
  character(len=19)    :: char_startDate, char_endDate !format YYYY-MM-DD HH:MM:SS
  integer :: init_year, init_month, init_day, init_hour, init_min, init_sec, end_year, end_month, end_day, end_hour, end_min, end_sec
  integer :: adv_u, adv_v, adv_temp, adv_theta, adv_thetal, adv_qv, adv_qt, adv_rv, adv_rt, forc_w, forc_omega, forc_geo
  integer :: rad_temp, rad_theta, rad_thetal
  character(len=3) :: char_rad_temp, char_rad_theta, char_rad_thetal
  integer :: nudging_temp, nudging_theta, nudging_thetal, nudging_qv, nudging_qt, nudging_rv, nudging_rt, nudging_u, nudging_v
  real(kind=sp) :: z_nudging_temp, z_nudging_theta, z_nudging_thetal, z_nudging_qv, z_nudging_qt, z_nudging_rv, z_nudging_rt, z_nudging_u, z_nudging_v
  real(kind=sp) :: p_nudging_temp, p_nudging_theta, p_nudging_thetal, p_nudging_qv, p_nudging_qt, p_nudging_rv, p_nudging_rt, p_nudging_u, p_nudging_v
  character(len=5) :: input_surfaceType
  character(len=11) :: input_surfaceForcingWind='',input_surfaceForcingMoist='',input_surfaceForcingLSM='',input_surfaceForcingTemp=''
  
  ! initial variables (IC = Initial Condition)
  real(kind=dp), allocatable :: input_lat(:) !< column latitude (deg)
  real(kind=dp), allocatable :: input_lon(:) !< column longitude (deg)
  real(kind=dp), allocatable :: input_z0 (:)    ! surfce_forcing: surface_roughness_length_for_momentum_in_air
  real(kind=sp), allocatable :: input_pres(:,:) !< IC pressure levels (Pa)
  real(kind=sp), allocatable :: input_height(:,:) !< IC height levels (m)
  real(kind=sp), allocatable :: input_pres_surf(:) !< IC surface pressure (Pa)
  real(kind=sp), allocatable :: input_u(:,:) !< IC east-west horizontal wind profile (m s^-1)
  real(kind=sp), allocatable :: input_v(:,:) !< IC north-south horizontal wind profile (m s^-1)
  real(kind=sp), allocatable :: input_temp(:,:) !< IC temperature profile (K)
  real(kind=sp), allocatable :: input_theta(:,:) !< IC potential temperature profile (K)
  real(kind=sp), allocatable :: input_thetal(:,:) !< IC liquid potential temperature profile (K)
  real(kind=sp), allocatable :: input_qv(:,:) !< IC specific humidity profile (kg kg^-1)
  real(kind=sp), allocatable :: input_qt(:,:) !< IC total water specific humidity profile (kg kg^-1)
  real(kind=sp), allocatable :: input_ql(:,:) !< IC liquid water specific humidity profile (kg kg^-1)
  real(kind=sp), allocatable :: input_qi(:,:) !< IC ice water specific humidity profile (kg kg^-1)
  real(kind=sp), allocatable :: input_rv(:,:) !< IC water vapor mixing ratio profile (kg kg^-1)
  real(kind=sp), allocatable :: input_rt(:,:) !< IC total water mixing ratio profile (kg kg^-1)
  real(kind=sp), allocatable :: input_rl(:,:) !< IC liquid water mixing ratio profile (kg kg^-1)
  real(kind=sp), allocatable :: input_ri(:,:) !< IC ice water mixing ratio profile (kg kg^-1)
  real(kind=sp), allocatable :: input_rh(:,:) !< IC relative humidity profile (%)
  real(kind=sp), allocatable :: input_tke(:,:) !< IC TKE profile (m^2 s^-2)
  
  ! Model ICs (extension of DEPHY format)
  real(kind=dp), allocatable  :: input_ozone(:,:)   !< ozone profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_stc(:,:)     !< soil temperature (K)
  real(kind=dp), allocatable  :: input_smc(:,:)     !< total soil moisture content (fraction)  
  real(kind=dp), allocatable  :: input_slc(:,:)     !< liquid soil moisture content (fraction)
  real(kind=dp), allocatable  :: input_snicexy(:,:) !< snow layer ice (mm)
  real(kind=dp), allocatable  :: input_snliqxy(:,:) !< snow layer liquid (mm)
  real(kind=dp), allocatable  :: input_tsnoxy(:,:)  !< snow temperature (K)
  real(kind=dp), allocatable  :: input_smoiseq(:,:) !< equilibrium soil water content (m3 m-3)
  real(kind=dp), allocatable  :: input_zsnsoxy(:,:) !< layer bottom depth from snow surface (m)
  real(kind=dp), allocatable  :: input_tiice(:,:)   !< sea ice internal temperature (K)
  real(kind=dp), allocatable  :: input_tslb(:,:)    !< soil temperature for RUC LSM (K)
  real(kind=dp), allocatable  :: input_smois(:,:)   !< volume fraction of soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_sh2o(:,:)    !< volume fraction of unfrozen soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_smfr(:,:)    !< volume fraction of frozen soil moisture for RUC LSM (frac)
  real(kind=dp), allocatable  :: input_flfr(:,:)    !< flag for frozen soil physics
  
  real(kind=dp), allocatable  :: input_area(:)      !< surface area [m^2]
  real(kind=dp), allocatable  :: input_tsfco(:) !< input sea surface temperature OR surface skin temperature over land OR surface skin temperature over ice (depending on slmsk) (K)
  integer      , allocatable  :: input_vegsrc(:) !< vegetation source
  integer      , allocatable  :: input_vegtyp(:) !< vegetation type
  integer      , allocatable  :: input_soiltyp(:)!< soil type
  integer      , allocatable  :: input_slopetype(:) !< slope type
  real(kind=dp), allocatable  :: input_vegfrac(:)  !< vegetation fraction
  real(kind=dp), allocatable  :: input_shdmin(:)  !< minimun vegetation fraction
  real(kind=dp), allocatable  :: input_shdmax(:)  !< maximun vegetation fraction
  real(kind=dp), allocatable  :: input_slmsk(:)   !< sea land ice mask [0,1,2]
  real(kind=dp), allocatable  :: input_canopy(:)  !< amount of water stored in canopy (kg m-2)
  real(kind=dp), allocatable  :: input_hice(:)    !< sea ice thickness (m)
  real(kind=dp), allocatable  :: input_fice(:)    !< ice fraction (frac)
  real(kind=dp), allocatable  :: input_tisfc(:)   !< ice surface temperature (K)
  real(kind=dp), allocatable  :: input_snwdph(:)  !< water equivalent snow depth (mm)
  real(kind=dp), allocatable  :: input_snoalb(:)  !< maximum snow albedo (frac)
  real(kind=dp), allocatable  :: input_sncovr(:)  !< snow area fraction (frac)
  real(kind=dp), allocatable  :: input_tg3(:)     !< deep soil temperature (K)
  real(kind=dp), allocatable  :: input_uustar(:)  !< surface friction velocity (m s-1)
  real(kind=dp), allocatable  :: input_alvsf(:) !< 60 degree vis albedo with strong cosz dependency
  real(kind=dp), allocatable  :: input_alnsf(:) !< 60 degree nir albedo with strong cosz dependency
  real(kind=dp), allocatable  :: input_alvwf(:) !< 60 degree vis albedo with weak cosz dependency
  real(kind=dp), allocatable  :: input_alnwf(:) !< 60 degree nir albedo with weak cosz dependency
  real(kind=dp), allocatable  :: input_facsf(:) !< fractional coverage with strong cosz dependency
  real(kind=dp), allocatable  :: input_facwf(:) !< fractional coverage with weak cosz dependency
  real(kind=dp), allocatable  :: input_weasd(:) !< water equivalent accumulated snow depth (mm)
  real(kind=dp), allocatable  :: input_f10m(:)  !< ratio of sigma level 1 wind and 10m wind
  real(kind=dp), allocatable  :: input_t2m(:)    !< 2-meter absolute temperature (K)
  real(kind=dp), allocatable  :: input_q2m(:)    !< 2-meter specific humidity (kg kg-1)
  real(kind=dp), allocatable  :: input_ffmm(:)    !< Monin-Obukhov similarity function for momentum
  real(kind=dp), allocatable  :: input_ffhh(:)    !< Monin-Obukhov similarity function for heat
  real(kind=dp), allocatable  :: input_tprcp(:)   !< instantaneous total precipitation amount (m)
  real(kind=dp), allocatable  :: input_srflag(:)  !< snow/rain flag for precipitation
  real(kind=dp), allocatable  :: input_tsfcl(:)   !< surface skin temperature over land (K)
  real(kind=dp), allocatable  :: input_zorll(:)   !< surface roughness length over land (cm)
  real(kind=dp), allocatable  :: input_zorli(:)   !< surface roughness length over ice (cm)
  real(kind=dp), allocatable  :: input_zorlw(:)   !< surface roughness length from wave model (cm)
  
  real(kind=dp), allocatable  :: input_stddev(:) !< standard deviation of subgrid orography (m)
  real(kind=dp), allocatable  :: input_convexity(:) !< convexity of subgrid orography 
  real(kind=dp), allocatable  :: input_ol1(:) !< fraction of grid box with subgrid orography higher than critical height 1
  real(kind=dp), allocatable  :: input_ol2(:) !< fraction of grid box with subgrid orography higher than critical height 2
  real(kind=dp), allocatable  :: input_ol3(:) !< fraction of grid box with subgrid orography higher than critical height 3
  real(kind=dp), allocatable  :: input_ol4(:) !< fraction of grid box with subgrid orography higher than critical height 4
  real(kind=dp), allocatable  :: input_oa1(:) !< assymetry of subgrid orography 1
  real(kind=dp), allocatable  :: input_oa2(:) !< assymetry of subgrid orography 2
  real(kind=dp), allocatable  :: input_oa3(:) !< assymetry of subgrid orography 3
  real(kind=dp), allocatable  :: input_oa4(:) !< assymetry of subgrid orography 4
  real(kind=dp), allocatable  :: input_sigma(:) !< slope of subgrid orography
  real(kind=dp), allocatable  :: input_theta_oro(:) !< angle with respect to east of maximum subgrid orographic variations (deg)
  real(kind=dp), allocatable  :: input_gamma(:) !< anisotropy of subgrid orography
  real(kind=dp), allocatable  :: input_elvmax(:)!< maximum of subgrid orography (m)
  real(kind=dp), allocatable  :: input_oro(:) !< orography (m)
  real(kind=dp), allocatable  :: input_oro_uf(:) !< unfiltered orography (m)
  real(kind=dp), allocatable  :: input_landfrac(:) !< fraction of horizontal grid area occupied by land
  real(kind=dp), allocatable  :: input_lakefrac(:) !< fraction of horizontal grid area occupied by lake
  real(kind=dp), allocatable  :: input_lakedepth(:) !< lake depth (m)
  
  real(kind=dp), allocatable  :: input_tvxy(:) !< vegetation temperature (K)
  real(kind=dp), allocatable  :: input_tgxy(:) !< ground temperature for Noahmp (K)
  real(kind=dp), allocatable  :: input_tahxy(:) !< canopy air temperature (K)
  real(kind=dp), allocatable  :: input_canicexy(:) !< canopy intercepted ice mass (mm)
  real(kind=dp), allocatable  :: input_canliqxy(:) !< canopy intercepted liquid water (mm)
  real(kind=dp), allocatable  :: input_eahxy(:) !< canopy air vapor pressure (Pa)
  real(kind=dp), allocatable  :: input_cmxy(:) !< surface drag coefficient for momentum for noahmp
  real(kind=dp), allocatable  :: input_chxy(:) !< surface exchange coeff heat & moisture for noahmp
  real(kind=dp), allocatable  :: input_fwetxy(:) !< area fraction of canopy that is wetted/snowed
  real(kind=dp), allocatable  :: input_sneqvoxy(:) !< snow mass at previous time step (mm)
  real(kind=dp), allocatable  :: input_alboldxy(:) !< snow albedo at previous time step (frac)
  real(kind=dp), allocatable  :: input_qsnowxy(:) !< snow precipitation rate at surface (mm s-1)
  real(kind=dp), allocatable  :: input_wslakexy(:) !< lake water storage (mm)
  real(kind=dp), allocatable  :: input_taussxy(:) !< non-dimensional snow age
  real(kind=dp), allocatable  :: input_waxy(:) !< water storage in aquifer (mm)
  real(kind=dp), allocatable  :: input_wtxy(:) !< water storage in aquifer and saturated soil (mm)
  real(kind=dp), allocatable  :: input_zwtxy(:) !< water table depth (m)
  real(kind=dp), allocatable  :: input_xlaixy(:) !< leaf area index
  real(kind=dp), allocatable  :: input_xsaixy(:) !< stem area index
  real(kind=dp), allocatable  :: input_lfmassxy(:) !< leaf mass (g m-2)
  real(kind=dp), allocatable  :: input_stmassxy(:) !< stem mass (g m-2)
  real(kind=dp), allocatable  :: input_rtmassxy(:) !< fine root mass (g m-2)
  real(kind=dp), allocatable  :: input_woodxy(:) !< wood mass including woody roots (g m-2)
  real(kind=dp), allocatable  :: input_stblcpxy(:) !< stable carbon in deep soil (g m-2)
  real(kind=dp), allocatable  :: input_fastcpxy(:) !< short-lived carbon in shallow soil (g m-2)
  real(kind=dp), allocatable  :: input_smcwtdxy(:) !< soil water content between the bottom of the soil and the water table (m3 m-3)
  real(kind=dp), allocatable  :: input_deeprechxy(:) !< recharge to or from the water table when deep (m)
  real(kind=dp), allocatable  :: input_rechxy(:) !< recharge to or from the water table when shallow (m)
  real(kind=dp), allocatable  :: input_snowxy(:) !< number of snow layers
  
  real(kind=dp), allocatable  :: input_tref(:) !< sea surface reference temperature for NSST (K)
  real(kind=dp), allocatable  :: input_z_c(:) !< sub-layer cooling thickness for NSST (m)
  real(kind=dp), allocatable  :: input_c_0(:) !< coefficient 1 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp), allocatable  :: input_c_d(:) !< coefficient 2 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp), allocatable  :: input_w_0(:) !< coefficient 3 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp), allocatable  :: input_w_d(:) !< coefficient 4 to calculate d(Tz)/d(Ts) for NSST
  real(kind=dp), allocatable  :: input_xt(:) !< heat content in diurnal thermocline layer for NSST (K m)
  real(kind=dp), allocatable  :: input_xs(:) !< salinity content in diurnal thermocline layer for NSST (ppt m)
  real(kind=dp), allocatable  :: input_xu(:) !< u-current in diurnal thermocline layer for NSST (m2 s-1)
  real(kind=dp), allocatable  :: input_xv(:) !< v-current in diurnal thermocline layer for NSST (m2 s-1)
  real(kind=dp), allocatable  :: input_xz(:) !< thickness of diurnal thermocline layer for NSST (m)
  real(kind=dp), allocatable  :: input_zm(:) !< thickness of ocean mixed layer for NSST (m)
  real(kind=dp), allocatable  :: input_xtts(:) !< sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST (m)
  real(kind=dp), allocatable  :: input_xzts(:) !< sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST (m K-1)
  real(kind=dp), allocatable  :: input_d_conv(:) !< thickness of free convection layer for NSST (m)
  real(kind=dp), allocatable  :: input_ifd(:) !< index to start DTM run for NSST
  real(kind=dp), allocatable  :: input_dt_cool(:) !< sub-layer cooling amount for NSST (K)
  real(kind=dp), allocatable  :: input_qrain(:) !< sensible heat due to rainfall for NSST (W)
  
  real(kind=dp), allocatable  :: input_wetness(:)         !< normalized soil wetness for RUC LSM
  real(kind=dp), allocatable  :: input_lai(:)             !< leaf area index for RUC LSM
  real(kind=dp), allocatable  :: input_clw_surf_land(:)   !< cloud condensed water mixing ratio at surface over land for RUC LSM 
  real(kind=dp), allocatable  :: input_clw_surf_ice(:)    !< cloud condensed water mixing ratio at surface over ice for RUC LSM
  real(kind=dp), allocatable  :: input_qwv_surf_land(:)   !< water vapor mixing ratio at surface over land for RUC LSM 
  real(kind=dp), allocatable  :: input_qwv_surf_ice(:)    !< water vapor mixing ratio at surface over ice for RUC LSM 
  real(kind=dp), allocatable  :: input_tsnow_land(:)      !< snow temperature at the bottom of the first snow layer over land for RUC LSM 
  real(kind=dp), allocatable  :: input_tsnow_ice(:)       !< snow temperature at the bottom of the first snow layer over ice for RUC LSM 
  real(kind=dp), allocatable  :: input_snowfallac_land(:) !< run-total snow accumulation on the ground over land for RUC LSM 
  real(kind=dp), allocatable  :: input_snowfallac_ice(:)  !< run-total snow accumulation on the ground over ice for RUC LSM 
  real(kind=dp), allocatable  :: input_sncovr_ice(:)      !<
  real(kind=dp), allocatable  :: input_sfalb_lnd(:)       !<
  real(kind=dp), allocatable  :: input_sfalb_lnd_bck(:)   !<
  real(kind=dp), allocatable  :: input_sfalb_ice(:)       !<
  real(kind=dp), allocatable  :: input_emis_ice(:)        !<
  
  ! forcing variables
  real(kind=sp), allocatable :: input_force_pres_surf(:) !< forcing surface pressure (Pa)
  real(kind=sp), allocatable :: input_force_height(:,:) !< forcing height levels (m)
  real(kind=sp), allocatable :: input_force_pres(:,:) !< forcing pressure levels (Pa)
  real(kind=sp), allocatable :: input_force_u_g(:,:) !< forcing 2D geostrophic east-west wind (m s^-1)
  real(kind=sp), allocatable :: input_force_v_g(:,:) !< forcing 2D geostrophic north-south wind (m s^-1)
  real(kind=sp), allocatable :: input_force_w(:,:)
  real(kind=sp), allocatable :: input_force_omega(:,:)
  real(kind=sp), allocatable :: input_force_u_adv(:,:)
  real(kind=sp), allocatable :: input_force_v_adv(:,:)
  real(kind=sp), allocatable :: input_force_temp_adv(:,:)
  real(kind=sp), allocatable :: input_force_theta_adv(:,:)
  real(kind=sp), allocatable :: input_force_thetal_adv(:,:)
  real(kind=sp), allocatable :: input_force_qt_adv(:,:)
  real(kind=sp), allocatable :: input_force_qv_adv(:,:)
  real(kind=sp), allocatable :: input_force_rt_adv(:,:)
  real(kind=sp), allocatable :: input_force_rv_adv(:,:)
  real(kind=sp), allocatable :: input_force_temp_rad(:,:)
  real(kind=sp), allocatable :: input_force_theta_rad(:,:)
  real(kind=sp), allocatable :: input_force_thetal_rad(:,:)
  real(kind=sp), allocatable :: input_force_sfc_sens_flx(:)
  real(kind=sp), allocatable :: input_force_sfc_lat_flx(:)
  real(kind=sp), allocatable :: input_force_wpthetap(:)
  real(kind=sp), allocatable :: input_force_wpqvp(:)
  real(kind=sp), allocatable :: input_force_wpqtp(:)
  real(kind=sp), allocatable :: input_force_wprvp(:)
  real(kind=sp), allocatable :: input_force_wprtp(:)
  real(kind=sp), allocatable :: input_force_ts(:)
  real(kind=sp), allocatable :: input_force_ustar(:)
  real(kind=sp), allocatable :: input_force_temp_nudging(:,:)
  real(kind=sp), allocatable :: input_force_theta_nudging(:,:)
  real(kind=sp), allocatable :: input_force_thetal_nudging(:,:)
  real(kind=sp), allocatable :: input_force_qv_nudging(:,:)
  real(kind=sp), allocatable :: input_force_qt_nudging(:,:)
  real(kind=sp), allocatable :: input_force_rv_nudging(:,:)
  real(kind=sp), allocatable :: input_force_rt_nudging(:,:)
  real(kind=sp), allocatable :: input_force_u_nudging(:,:)
  real(kind=sp), allocatable :: input_force_v_nudging(:,:)
  
  integer :: ncid, varID, allocate_status, ierr, i, k
  integer :: active_lon, active_lat, active_init_time
  CHARACTER(LEN=nf90_max_name) :: tmpName
  real(kind=sp), parameter :: p0 = 100000.0
  real(kind=sp) :: exner, exner_inv, rho, elapsed_sec, missing_value_eps
  real(kind=dp) :: rinc(5)
  integer :: jdat(1:8), idat(1:8) !(yr, mon, day, t-zone, hr, min, sec, mil-sec)

  integer :: input_n_init_times, input_n_forcing_times, input_n_lev, input_n_snow, input_n_ice, input_n_soil
  
  missing_value_eps = missing_value + 0.01
  
  !> - Open the case input file found in the processed_case_input dir corresponding to the experiment name.
  call check(NF90_OPEN(trim(adjustl(scm_state%case_name))//'_SCM_driver.nc',nf90_nowrite,ncid),"nf90_open()")
  
  !> - Get the dimensions.
  
  !required dimensions
  call check(NF90_INQ_DIMID(ncid,"t0",varID),"nf90_inq_dimid(t0)")
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_init_times),"nf90_inq_dim(t0)")
  call check(NF90_INQ_DIMID(ncid,"time",varID),"nf90_inq_dimid(time)")
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_forcing_times),"nf90_inq_dim(time)")
  call check(NF90_INQ_DIMID(ncid,"lev",varID),"nf90_inq_dimid(lev)")
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_lev),"nf90_inq_dim(lev)")
  !Check whether long_name = 'altitude', units='m' OR long_name = 'pressure', units='Pa'?
  !It may not matter, because 'lev' may not be needed when the IC pressure and height are BOTH already provided

  !### TO BE USED IF DEPHY-SCM can be extended to include model ICs ###
  !possible dimensions (if using model ICs)
  ierr = NF90_INQ_DIMID(ncid,"nsoil",varID)
  if(ierr /= NF90_NOERR) then
    input_n_soil = missing_soil_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_soil),"nf90_inq_dim(nsoil)")
  end if
  ierr = NF90_INQ_DIMID(ncid,"nsnow",varID)
  if(ierr /= NF90_NOERR) then
    input_n_snow = missing_snow_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_snow),"nf90_inq_dim(nsnow)")
  end if
  ierr = NF90_INQ_DIMID(ncid,"nice",varID)
  if(ierr /= NF90_NOERR) then
    input_n_ice = missing_ice_layers
  else
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_n_ice),"nf90_inq_dim(nice)")
  end if  
  
  !> - Allocate the dimension variables.
  allocate(input_t0    (input_n_init_times),                                        &
           input_time  (input_n_forcing_times),                                     &
           input_lev   (input_n_lev),                                               &
    stat=allocate_status)

  !> - Read in the dimension variables (required).
  call NetCDF_read_var(ncid, "t0", .True., input_t0)
  call NetCDF_read_var(ncid, "time", .True., input_time)
  call NetCDF_read_var(ncid, "lev", .True., input_lev)

  !> - Read in global attributes

  call NetCDF_read_att(ncid, NF90_GLOBAL, 'start_date', .True., char_startDate)
    
  read(char_startDate(1:4),'(i4)')   init_year
  read(char_startDate(6:7),'(i2)')   init_month
  read(char_startDate(9:10),'(i2)')  init_day
  read(char_startDate(12:13),'(i2)') init_hour
  read(char_startDate(15:16),'(i2)') init_min
  read(char_startDate(18:19),'(i2)') init_sec
  
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'end_date', .True., char_endDate)
  
  read(char_endDate(1:4),'(i4)')   end_year
  read(char_endDate(6:7),'(i2)')   end_month
  read(char_endDate(9:10),'(i2)')  end_day
  read(char_endDate(12:13),'(i2)') end_hour
  read(char_endDate(15:16),'(i2)') end_min
  read(char_endDate(18:19),'(i2)') end_sec
  
  !compare init time to what was in case config file? replace?
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_ua',             .False., adv_u)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_va',             .False., adv_v)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_ta',             .False., adv_temp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_theta',          .False., adv_theta)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_thetal',         .False., adv_thetal)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'rad_ta',             .False., rad_temp, char_rad_temp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'rad_theta',          .False., rad_theta, char_rad_theta)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'rad_thetal',         .False., rad_thetal, char_rad_thetal)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_qv',             .False., adv_qv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_qt',             .False., adv_qt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_rv',             .False., adv_rv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'adv_rt',             .False., adv_rt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'forc_wa',            .False., forc_w)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'forc_wap',           .False., forc_omega)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'forc_geo',           .False., forc_geo)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_ta',         .False., nudging_temp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_theta',      .False., nudging_theta)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_thetal',     .False., nudging_thetal)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_qv',         .False., nudging_qv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_qt',         .False., nudging_qt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_rv',         .False., nudging_rv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_rt',         .False., nudging_rt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_ua',         .False., nudging_u)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'nudging_va',         .False., nudging_v)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_ta',      .False., z_nudging_temp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_theta',   .False., z_nudging_theta)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_thetal',  .False., z_nudging_thetal)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_qv',      .False., z_nudging_qv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_qt',      .False., z_nudging_qt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_rv',      .False., z_nudging_rv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_rt',      .False., z_nudging_rt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_ua',      .False., z_nudging_u)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'zh_nudging_va',      .False., z_nudging_v)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_ta',      .False., p_nudging_temp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_theta',   .False., p_nudging_theta)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_thetal',  .False., p_nudging_thetal)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_qv',      .False., p_nudging_qv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_qt',      .False., p_nudging_qt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_rv',      .False., p_nudging_rv)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_rt',      .False., p_nudging_rt)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_ua',      .False., p_nudging_u)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'pa_nudging_va',      .False., p_nudging_v)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'surface_type',       .False., input_surfaceType)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'surface_forcing_temp',     .False., input_surfaceForcingTemp)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'surface_forcing_moisture', .False., input_surfaceForcingMoist)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'surface_forcing_wind',     .False., input_surfaceForcingWind)
  call NetCDF_read_att(ncid, NF90_GLOBAL, 'surface_forcing_lsm',      .False., input_surfaceForcingLSM)

  !> - Allocate the initial variables.
  allocate(input_pres     (input_n_lev, input_n_init_times), &
           input_height   (input_n_lev, input_n_init_times), &
           input_pres_surf(input_n_init_times),              &
           input_u        (input_n_lev, input_n_init_times), &
           input_v        (input_n_lev, input_n_init_times), &
           input_temp     (input_n_lev, input_n_init_times), &
           input_theta    (input_n_lev, input_n_init_times), &
           input_thetal   (input_n_lev, input_n_init_times), &
           input_qv       (input_n_lev, input_n_init_times), &
           input_qt       (input_n_lev, input_n_init_times), &
           input_ql       (input_n_lev, input_n_init_times), &
           input_qi       (input_n_lev, input_n_init_times), &
           input_rv       (input_n_lev, input_n_init_times), &
           input_rt       (input_n_lev, input_n_init_times), &
           input_rl       (input_n_lev, input_n_init_times), &
           input_ri       (input_n_lev, input_n_init_times), &
           input_rh       (input_n_lev, input_n_init_times), &
           input_tke      (input_n_lev, input_n_init_times), &
           stat=allocate_status)

  if (trim(input_surfaceForcingLSM) == "lsm") then
    !if model ICs are included in the file
    scm_state%lsm_ics = .true.

    !variables with vertical extent
    allocate(input_ozone   (input_n_lev,  input_n_init_times), &
             input_stc     (input_n_soil, input_n_init_times), &
             input_smc     (input_n_soil, input_n_init_times), &
             input_slc     (input_n_soil, input_n_init_times), &
             input_snicexy (input_n_snow, input_n_init_times), &
             input_snliqxy (input_n_snow, input_n_init_times), &
             input_tsnoxy  (input_n_snow, input_n_init_times), &
             input_smoiseq (input_n_soil, input_n_init_times), &
             input_zsnsoxy (input_n_soil + input_n_snow, input_n_init_times), &
             input_tiice   (input_n_ice,  input_n_init_times), &
             input_tslb    (input_n_soil, input_n_init_times), &
             input_smois   (input_n_soil, input_n_init_times), &
             input_sh2o    (input_n_soil, input_n_init_times), &
             input_smfr    (input_n_soil, input_n_init_times), &
             input_flfr    (input_n_soil, input_n_init_times), &
             stat=allocate_status)
             
             
    !variables without vertical extent
    allocate(input_area      (          input_n_init_times), &
             input_tsfco     (          input_n_init_times), &
             input_vegsrc    (          input_n_init_times), &
             input_vegtyp    (          input_n_init_times), &
             input_soiltyp   (          input_n_init_times), &
             input_slopetype (          input_n_init_times), &
             input_vegfrac   (          input_n_init_times), &
             input_shdmin    (          input_n_init_times), &
             input_shdmax    (          input_n_init_times), &
             input_slmsk     (          input_n_init_times), &
             input_canopy    (          input_n_init_times), &
             input_hice      (          input_n_init_times), &
             input_fice      (          input_n_init_times), &
             input_tisfc     (          input_n_init_times), &
             input_snwdph    (          input_n_init_times), &
             input_snoalb    (          input_n_init_times), &
             input_sncovr    (          input_n_init_times), &
             input_tg3       (          input_n_init_times), &
             input_uustar    (          input_n_init_times), &
             input_alvsf     (          input_n_init_times), &
             input_alnsf     (          input_n_init_times), &
             input_alvwf     (          input_n_init_times), &
             input_alnwf     (          input_n_init_times), &
             input_facsf     (          input_n_init_times), &
             input_facwf     (          input_n_init_times), &
             input_weasd     (          input_n_init_times), &
             input_f10m      (          input_n_init_times), &
             input_t2m       (          input_n_init_times), &
             input_q2m       (          input_n_init_times), &
             input_ffmm      (          input_n_init_times), &
             input_ffhh      (          input_n_init_times), &
             input_tprcp     (          input_n_init_times), &
             input_srflag    (          input_n_init_times), &
             input_tsfcl     (          input_n_init_times), &
             input_zorll     (          input_n_init_times), &
             input_zorli     (          input_n_init_times), &
             input_zorlw     (          input_n_init_times), &
             stat=allocate_status)
    allocate(input_stddev    (          input_n_init_times), &
             input_convexity (          input_n_init_times), &
             input_ol1       (          input_n_init_times), &
             input_ol2       (          input_n_init_times), &
             input_ol3       (          input_n_init_times), &
             input_ol4       (          input_n_init_times), &
             input_oa1       (          input_n_init_times), &
             input_oa2       (          input_n_init_times), &
             input_oa3       (          input_n_init_times), &
             input_oa4       (          input_n_init_times), &
             input_sigma     (          input_n_init_times), &
             input_theta_oro (          input_n_init_times), &
             input_gamma     (          input_n_init_times), &
             input_elvmax    (          input_n_init_times), &
             input_oro       (          input_n_init_times), &
             input_oro_uf    (          input_n_init_times), &
             input_landfrac  (          input_n_init_times), &
             input_lakefrac  (          input_n_init_times), &
             input_lakedepth (          input_n_init_times), &
             stat=allocate_status)
    allocate(input_tvxy      (          input_n_init_times), &
             input_tgxy      (          input_n_init_times), & 
             input_tahxy     (          input_n_init_times), &
             input_canicexy  (          input_n_init_times), &
             input_canliqxy  (          input_n_init_times), &
             input_eahxy     (          input_n_init_times), &
             input_cmxy      (          input_n_init_times), &
             input_chxy      (          input_n_init_times), &
             input_fwetxy    (          input_n_init_times), &
             input_sneqvoxy  (          input_n_init_times), &
             input_alboldxy  (          input_n_init_times), &
             input_qsnowxy   (          input_n_init_times), &
             input_wslakexy  (          input_n_init_times), &
             input_taussxy   (          input_n_init_times), &
             input_waxy      (          input_n_init_times), &
             input_wtxy      (          input_n_init_times), &
             input_zwtxy     (          input_n_init_times), &
             input_xlaixy    (          input_n_init_times), &
             input_xsaixy    (          input_n_init_times), &
             input_lfmassxy  (          input_n_init_times), &
             input_stmassxy  (          input_n_init_times), &
             input_rtmassxy  (          input_n_init_times), &
             input_woodxy    (          input_n_init_times), &
             input_stblcpxy  (          input_n_init_times), &
             input_fastcpxy  (          input_n_init_times), &
             input_smcwtdxy  (          input_n_init_times), &
             input_deeprechxy(          input_n_init_times), &
             input_rechxy    (          input_n_init_times), & 
             input_snowxy    (          input_n_init_times), & 
             stat=allocate_status)
    allocate(input_tref      (          input_n_init_times), &
             input_z_c       (          input_n_init_times), &
             input_c_0       (          input_n_init_times), &
             input_c_d       (          input_n_init_times), &
             input_w_0       (          input_n_init_times), &
             input_w_d       (          input_n_init_times), &
             input_xt        (          input_n_init_times), &
             input_xs        (          input_n_init_times), &
             input_xu        (          input_n_init_times), &
             input_xv        (          input_n_init_times), &
             input_xz        (          input_n_init_times), &
             input_zm        (          input_n_init_times), &
             input_xtts      (          input_n_init_times), &
             input_xzts      (          input_n_init_times), &
             input_d_conv    (          input_n_init_times), &
             input_ifd       (          input_n_init_times), &
             input_dt_cool   (          input_n_init_times), &
             input_qrain     (          input_n_init_times), &
             stat=allocate_status)
    allocate(input_wetness         (          input_n_init_times), &
             input_lai             (          input_n_init_times), &
             input_clw_surf_land   (          input_n_init_times), &
             input_clw_surf_ice    (          input_n_init_times), &
             input_qwv_surf_land   (          input_n_init_times), &
             input_qwv_surf_ice    (          input_n_init_times), &
             input_tsnow_land      (          input_n_init_times), &
             input_tsnow_ice       (          input_n_init_times), &
             input_snowfallac_land (          input_n_init_times), &
             input_snowfallac_ice  (          input_n_init_times), &
             input_sncovr_ice      (          input_n_init_times), &
             input_sfalb_lnd       (          input_n_init_times), &
             input_sfalb_lnd_bck   (          input_n_init_times), &
             input_sfalb_ice       (          input_n_init_times), &
             input_emis_ice        (          input_n_init_times), &
             stat=allocate_status)
  end if

  !>  - Read in the initial profiles.
  call NetCDF_read_var(ncid, "pa", .True., input_pres)
  call NetCDF_read_var(ncid, "zh", .True., input_height)
  call NetCDF_read_var(ncid, "ps", .True., input_pres_surf)
  call NetCDF_read_var(ncid, "ua", .True., input_u)
  call NetCDF_read_var(ncid, "va", .True., input_v)
  
  !one of the following should be present, but not all, hence they are not requried
  call NetCDF_read_var(ncid, "ta", .False., input_temp)
  call NetCDF_read_var(ncid, "theta", .False., input_theta)
  call NetCDF_read_var(ncid, "thetal", .False., input_thetal)
  
  !one or more of the following should be present, but not all, hence they are not requried
  call NetCDF_read_var(ncid, "qv",  .False., input_qv)
  call NetCDF_read_var(ncid, "qt",  .False., input_qt)
  call NetCDF_read_var(ncid, "ql",  .False., input_ql)
  call NetCDF_read_var(ncid, "qi",  .False., input_qi)
  call NetCDF_read_var(ncid, "rv",  .False., input_rv)
  call NetCDF_read_var(ncid, "rt",  .False., input_rt)
  call NetCDF_read_var(ncid, "rl",  .False., input_rl)
  call NetCDF_read_var(ncid, "ri",  .False., input_ri)
  call NetCDF_read_var(ncid, "hur", .False., input_rh)
 
  call NetCDF_read_var(ncid, "tke", .True., input_tke)
  
  if (trim(input_surfaceForcingLSM) == "lsm") then
    call NetCDF_read_var(ncid, "o3",      .True.,  input_ozone)
    call NetCDF_read_var(ncid, "area",    .True.,  input_area)
    
    !orographic parameters
    call NetCDF_read_var(ncid, "stddev",    .True., input_stddev)
    call NetCDF_read_var(ncid, "convexity", .True., input_convexity)
    call NetCDF_read_var(ncid, "oa1",       .True., input_oa1)
    call NetCDF_read_var(ncid, "oa2",       .True., input_oa2)
    call NetCDF_read_var(ncid, "oa3",       .True., input_oa3)
    call NetCDF_read_var(ncid, "oa4",       .True., input_oa4)
    call NetCDF_read_var(ncid, "ol1",       .True., input_ol1)
    call NetCDF_read_var(ncid, "ol2",       .True., input_ol2)
    call NetCDF_read_var(ncid, "ol3",       .True., input_ol3)
    call NetCDF_read_var(ncid, "ol4",       .True., input_ol4)
    call NetCDF_read_var(ncid, "theta_oro", .True., input_theta_oro)
    call NetCDF_read_var(ncid, "gamma",     .True., input_gamma)
    call NetCDF_read_var(ncid, "sigma",     .True., input_sigma)
    call NetCDF_read_var(ncid, "elvmax",    .True., input_elvmax)
    call NetCDF_read_var(ncid, "oro",       .True., input_oro)
    call NetCDF_read_var(ncid, "oro_uf",    .True., input_oro_uf)
    call NetCDF_read_var(ncid, "landfrac",  .True., input_landfrac)
    call NetCDF_read_var(ncid, "lakefrac",  .True., input_lakefrac)
    call NetCDF_read_var(ncid, "lakedepth", .True., input_lakedepth)
    
    !NSST variables
    call NetCDF_read_var(ncid, "tref",    .True., input_tref)
    call NetCDF_read_var(ncid, "z_c",     .True., input_z_c)
    call NetCDF_read_var(ncid, "c_0",     .True., input_c_0)
    call NetCDF_read_var(ncid, "c_d",     .True., input_c_d)
    call NetCDF_read_var(ncid, "w_0",     .True., input_w_0)
    call NetCDF_read_var(ncid, "w_d",     .True., input_w_d)
    call NetCDF_read_var(ncid, "xt",      .True., input_xt)
    call NetCDF_read_var(ncid, "xs",      .True., input_xs)
    call NetCDF_read_var(ncid, "xu",      .True., input_xu)
    call NetCDF_read_var(ncid, "xv",      .True., input_xv)
    call NetCDF_read_var(ncid, "xz",      .True., input_xz)
    call NetCDF_read_var(ncid, "zm",      .True., input_zm)
    call NetCDF_read_var(ncid, "xtts",    .True., input_xtts)
    call NetCDF_read_var(ncid, "xzts",    .True., input_xzts)
    call NetCDF_read_var(ncid, "d_conv",  .True., input_d_conv)
    call NetCDF_read_var(ncid, "ifd",     .True., input_ifd)
    call NetCDF_read_var(ncid, "dt_cool", .True., input_dt_cool)
    call NetCDF_read_var(ncid, "qrain",   .True., input_qrain)
  end if
  
  !> - Allocate the forcing variables.
  
  !allocate all, but conditionally read forcing variables given global atts; set unused forcing variables to missing

  allocate(input_lat                  (input_n_forcing_times),              &
           input_lon                  (input_n_forcing_times),              &
           input_z0                   (input_n_forcing_times),              &
           input_force_pres_surf      (input_n_forcing_times),              &
           input_force_pres           (input_n_lev, input_n_forcing_times), &
           input_force_height         (input_n_lev, input_n_forcing_times), &
           input_force_u_g            (input_n_lev, input_n_forcing_times), &
           input_force_v_g            (input_n_lev, input_n_forcing_times), &
           input_force_w              (input_n_lev, input_n_forcing_times), &
           input_force_omega          (input_n_lev, input_n_forcing_times), &
           input_force_u_adv          (input_n_lev, input_n_forcing_times), &
           input_force_v_adv          (input_n_lev, input_n_forcing_times), &
           input_force_temp_adv       (input_n_lev, input_n_forcing_times), &
           input_force_theta_adv      (input_n_lev, input_n_forcing_times), &
           input_force_thetal_adv     (input_n_lev, input_n_forcing_times), &
           input_force_qt_adv         (input_n_lev, input_n_forcing_times), &
           input_force_qv_adv         (input_n_lev, input_n_forcing_times), &
           input_force_rt_adv         (input_n_lev, input_n_forcing_times), &
           input_force_rv_adv         (input_n_lev, input_n_forcing_times), &
           input_force_temp_rad       (input_n_lev, input_n_forcing_times), &
           input_force_theta_rad      (input_n_lev, input_n_forcing_times), &
           input_force_thetal_rad     (input_n_lev, input_n_forcing_times), &
           input_force_sfc_sens_flx   (input_n_forcing_times),              &
           input_force_sfc_lat_flx    (input_n_forcing_times),              &
           input_force_wpthetap       (input_n_forcing_times),              &
           input_force_wpqvp          (input_n_forcing_times),              &
           input_force_wpqtp          (input_n_forcing_times),              &
           input_force_wprvp          (input_n_forcing_times),              &
           input_force_wprtp          (input_n_forcing_times),              &
           input_force_ts             (input_n_forcing_times),              &
           input_force_ustar          (input_n_forcing_times),              &
           input_force_u_nudging      (input_n_lev, input_n_forcing_times), &
           input_force_v_nudging      (input_n_lev, input_n_forcing_times), &
           input_force_temp_nudging   (input_n_lev, input_n_forcing_times), &
           input_force_theta_nudging  (input_n_lev, input_n_forcing_times), &
           input_force_thetal_nudging (input_n_lev, input_n_forcing_times), &
           input_force_qt_nudging     (input_n_lev, input_n_forcing_times), &
           input_force_qv_nudging     (input_n_lev, input_n_forcing_times), &
           input_force_rt_nudging     (input_n_lev, input_n_forcing_times), &
           input_force_rv_nudging     (input_n_lev, input_n_forcing_times), &
    stat=allocate_status)

  call NetCDF_read_var(ncid, "lat",     .True., input_lat)
  call NetCDF_read_var(ncid, "lon",     .True., input_lon)  
  call NetCDF_read_var(ncid, "ps_forc", .True., input_force_pres_surf)
  call NetCDF_read_var(ncid, "zh_forc", .True., input_force_height)
  call NetCDF_read_var(ncid, "pa_forc", .True., input_force_pres)
  
  !conditionally read forcing vars (or set to missing); if the global attribute is set to expect a variable and it doesn't exist, stop the model
  call NetCDF_conditionally_read_var(adv_u,      "adv_ua",     "tnua_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_u_adv)
  call NetCDF_conditionally_read_var(adv_v,      "adv_va",     "tnva_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_v_adv)
  call NetCDF_conditionally_read_var(adv_temp,   "adv_ta",     "tnta_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_temp_adv)
  call NetCDF_conditionally_read_var(adv_theta,  "adv_theta",  "tntheta_adv",  trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_theta_adv)
  call NetCDF_conditionally_read_var(adv_thetal, "adv_thetal", "tnthetal_adv", trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_thetal_adv)
  call NetCDF_conditionally_read_var(adv_qt,     "adv_qt",     "tnqt_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_qt_adv)
  call NetCDF_conditionally_read_var(adv_qv,     "adv_qv",     "tnqv_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_qv_adv)
  call NetCDF_conditionally_read_var(adv_rt,     "adv_rt",     "tnrt_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_rt_adv)
  call NetCDF_conditionally_read_var(adv_rv,     "adv_rv",     "tnrv_adv",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_rv_adv)
  call NetCDF_conditionally_read_var(rad_temp,   "rad_ta",     "tnta_rad",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_temp_rad)
  call NetCDF_conditionally_read_var(rad_theta,  "rad_theta",  "tntheta_rad",  trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_theta_rad)
  call NetCDF_conditionally_read_var(rad_thetal, "rad_thetal", "tnthetal_rad", trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_thetal_rad)
  !need to also handle the case when rad_[temp,theta,thetal]_char = 'adv' (make sure [temp,theta,thetal]_adv is not missing)
  !need to also turn off radiation when radiation is being forced (put in a warning that this is not supported for now?)
  
  call NetCDF_conditionally_read_var(forc_w,     "forc_w",   "wa",    trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_w)
  call NetCDF_conditionally_read_var(forc_omega, "forc_wap", "wap",   trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_omega)
  call NetCDF_conditionally_read_var(forc_geo,   "forc_geo", "ug",    trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_u_g)
  call NetCDF_conditionally_read_var(forc_geo,   "forc_geo", "vg",    trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_v_g)
  
  call NetCDF_conditionally_read_var(nudging_u,      "nudging_u",      "ua_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_u_nudging)
  call NetCDF_conditionally_read_var(nudging_v,      "nudging_v",      "va_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_v_nudging)
  call NetCDF_conditionally_read_var(nudging_temp,   "nudging_temp",   "ta_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_temp_nudging)
  call NetCDF_conditionally_read_var(nudging_theta,  "nudging_theta",  "theta_nud",  trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_theta_nudging)
  call NetCDF_conditionally_read_var(nudging_thetal, "nudging_thetal", "thetal_nud", trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_thetal_nudging)
  call NetCDF_conditionally_read_var(nudging_qv,     "nudging_qv",     "qv_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_qv_nudging)
  call NetCDF_conditionally_read_var(nudging_qt,     "nudging_qt",     "qt_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_qt_nudging)
  call NetCDF_conditionally_read_var(nudging_rv,     "nudging_rv",     "rv_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_rv_nudging)
  call NetCDF_conditionally_read_var(nudging_rt,     "nudging_rt",     "rt_nud",     trim(adjustl(scm_state%case_name))//'.nc', ncid, input_force_rt_nudging)

  !
  ! Surface forcing: Temperature
  !
  if (trim(input_surfaceForcingTemp) == 'kinematic') then
     call NetCDF_read_var(ncid, "wpthetap_s", .False., input_force_wpthetap)
     call NetCDF_read_var(ncid, "ts_forc",    .False., input_force_ts)
  else if (trim(input_surfaceForcingTemp) == 'surface_flux') then
     call NetCDF_read_var(ncid, "hfss",       .False., input_force_sfc_sens_flx)
     call NetCDF_read_var(ncid, "ts_forc",    .False., input_force_ts)
  else if (trim(input_surfaceForcingTemp) == 'ts') then
     call NetCDF_read_var(ncid, "ts_forc",    .False., input_force_ts)
  endif

  !
  ! Surface forcing: Moisture
  !
  if (trim(input_surfaceForcingMoist) == 'kinematic') then
     call NetCDF_read_var(ncid, "wpqvp_s",    .False., input_force_wpqvp)
     call NetCDF_read_var(ncid, "wpqtp_s",    .False., input_force_wpqtp)
     call NetCDF_read_var(ncid, "wprvp_s",    .False., input_force_wprvp)
     call NetCDF_read_var(ncid, "wprtp_s",    .False., input_force_wprtp)
  else if (trim(input_surfaceForcingMoist) == 'surface_flux') then
     call NetCDF_read_var(ncid, "hfls",       .False., input_force_sfc_sens_flx)
  endif

  !
  ! Surface forcing: Wind
  !
  if (trim(input_surfaceForcingWind) == 'z0') then
     call NetCDF_read_var(ncid, "z0",    .False., input_z0)
  else if (trim(input_surfaceForcingWind) == 'ustar') then
     call NetCDF_read_var(ncid, "ustar", .False., input_force_ustar)
  end if

  !
  ! Surface forcing Model LSM ICs
  !
  if (trim(input_surfaceForcingLSM) == "lsm") then
     call NetCDF_read_var(ncid, "stc",     .True., input_stc)
     call NetCDF_read_var(ncid, "smc",     .True., input_smc)
     call NetCDF_read_var(ncid, "slc",     .True., input_slc)
     call NetCDF_read_var(ncid, "snicexy", .True., input_snicexy)
     call NetCDF_read_var(ncid, "snliqxy", .True., input_snliqxy)
     call NetCDF_read_var(ncid, "tsnoxy",  .True., input_tsnoxy )
     call NetCDF_read_var(ncid, "smoiseq", .True., input_smoiseq)
     call NetCDF_read_var(ncid, "zsnsoxy", .True., input_zsnsoxy)
     call NetCDF_read_var(ncid, "tiice",   .True., input_tiice)
     call NetCDF_read_var(ncid, "tslb",    .True., input_tslb )
     call NetCDF_read_var(ncid, "smois",   .True., input_smois)
     call NetCDF_read_var(ncid, "sh2o",    .True., input_sh2o )
     call NetCDF_read_var(ncid, "smfr",    .True., input_smfr )
     call NetCDF_read_var(ncid, "flfr",    .True., input_flfr )
     
     call NetCDF_read_var(ncid, "vegsrc",   .True., input_vegsrc   )
     call NetCDF_read_var(ncid, "vegtyp",   .True., input_vegtyp   )
     call NetCDF_read_var(ncid, "soiltyp",  .True., input_soiltyp  )
     call NetCDF_read_var(ncid, "slopetyp", .True., input_slopetype)
     call NetCDF_read_var(ncid, "tsfco",    .True., input_tsfco)
     call NetCDF_read_var(ncid, "vegfrac",  .True., input_vegfrac)
     call NetCDF_read_var(ncid, "shdmin",   .True., input_shdmin)
     call NetCDF_read_var(ncid, "shdmax",   .True., input_shdmax)
     call NetCDF_read_var(ncid, "slmsk",    .True., input_slmsk)
     call NetCDF_read_var(ncid, "canopy",   .True., input_canopy)
     call NetCDF_read_var(ncid, "hice",     .True., input_hice)
     call NetCDF_read_var(ncid, "fice",     .True., input_fice)
     call NetCDF_read_var(ncid, "tisfc",    .True., input_tisfc)
     call NetCDF_read_var(ncid, "snowd",    .True., input_snwdph)
     call NetCDF_read_var(ncid, "snoalb",   .True., input_snoalb)
     call NetCDF_read_var(ncid, "tg3",      .True., input_tg3)
     call NetCDF_read_var(ncid, "uustar",   .True., input_uustar)
     call NetCDF_read_var(ncid, "alvsf",    .True., input_alvsf)
     call NetCDF_read_var(ncid, "alnsf",    .True., input_alnsf)
     call NetCDF_read_var(ncid, "alvwf",    .True., input_alvwf)
     call NetCDF_read_var(ncid, "alnwf",    .True., input_alnwf)
     call NetCDF_read_var(ncid, "facsf",    .True., input_facsf)
     call NetCDF_read_var(ncid, "facwf",    .True., input_facwf)
     call NetCDF_read_var(ncid, "weasd",    .True., input_weasd)
     call NetCDF_read_var(ncid, "f10m",     .True., input_f10m)
     call NetCDF_read_var(ncid, "t2m",      .True., input_t2m)
     call NetCDF_read_var(ncid, "q2m",      .True., input_q2m)
     call NetCDF_read_var(ncid, "ffmm",     .True., input_ffmm)
     call NetCDF_read_var(ncid, "ffhh",     .True., input_ffhh)
     call NetCDF_read_var(ncid, "tprcp",    .True., input_tprcp)
     call NetCDF_read_var(ncid, "srflag",   .True., input_srflag)
     call NetCDF_read_var(ncid, "sncovr",   .True., input_sncovr)
     call NetCDF_read_var(ncid, "tsfcl",    .True., input_tsfcl)
     call NetCDF_read_var(ncid, "zorll",    .True., input_zorll)
     call NetCDF_read_var(ncid, "zorli",    .True., input_zorli)
     call NetCDF_read_var(ncid, "zorlw",    .True., input_zorlw)
 
     !NoahMP parameters
     call NetCDF_read_var(ncid, "tvxy",      .False., input_tvxy)
     call NetCDF_read_var(ncid, "tgxy",      .False., input_tgxy)
     call NetCDF_read_var(ncid, "tahxy",     .False., input_tahxy)
     call NetCDF_read_var(ncid, "canicexy",  .False., input_canicexy)
     call NetCDF_read_var(ncid, "canliqxy",  .False., input_canliqxy)
     call NetCDF_read_var(ncid, "eahxy",     .False., input_eahxy)
     call NetCDF_read_var(ncid, "cmxy",      .False., input_cmxy)
     call NetCDF_read_var(ncid, "chxy",      .False., input_chxy)
     call NetCDF_read_var(ncid, "fwetxy",    .False., input_fwetxy)
     call NetCDF_read_var(ncid, "sneqvoxy",  .False., input_sneqvoxy)
     call NetCDF_read_var(ncid, "alboldxy",  .False., input_alboldxy)
     call NetCDF_read_var(ncid, "qsnowxy",   .False., input_qsnowxy)
     call NetCDF_read_var(ncid, "wslakexy",  .False., input_wslakexy)
     call NetCDF_read_var(ncid, "taussxy",   .False., input_taussxy)
     call NetCDF_read_var(ncid, "waxy",      .False., input_waxy)
     call NetCDF_read_var(ncid, "wtxy",      .False., input_wtxy)
     call NetCDF_read_var(ncid, "zwtxy",     .False., input_zwtxy)
     call NetCDF_read_var(ncid, "xlaixy",    .False., input_xlaixy)
     call NetCDF_read_var(ncid, "xsaixy",    .False., input_xsaixy)
     call NetCDF_read_var(ncid, "lfmassxy",  .False., input_lfmassxy)
     call NetCDF_read_var(ncid, "stmassxy",  .False., input_stmassxy)
     call NetCDF_read_var(ncid, "rtmassxy",  .False., input_rtmassxy)
     call NetCDF_read_var(ncid, "woodxy",    .False., input_woodxy)
     call NetCDF_read_var(ncid, "stblcpxy",  .False., input_stblcpxy)
     call NetCDF_read_var(ncid, "fastcpxy",  .False., input_fastcpxy)
     call NetCDF_read_var(ncid, "smcwtdxy",  .False., input_smcwtdxy)
     call NetCDF_read_var(ncid, "deeprechxy",.False., input_deeprechxy)
     call NetCDF_read_var(ncid, "rechxy",    .False., input_rechxy)
     call NetCDF_read_var(ncid, "snowxy",    .False., input_snowxy)
     !RUC LSM variables
     call NetCDF_read_var(ncid, "wetness",          .False., input_wetness)
     call NetCDF_read_var(ncid, "clw_surf_land",    .False., input_clw_surf_land)
     call NetCDF_read_var(ncid, "clw_surf_ice",     .False., input_clw_surf_ice)
     call NetCDF_read_var(ncid, "qwv_surf_land",    .False., input_qwv_surf_land)
     call NetCDF_read_var(ncid, "qwv_surf_ice",     .False., input_qwv_surf_ice)
     call NetCDF_read_var(ncid, "tsnow_land",       .False., input_tsnow_land)
     call NetCDF_read_var(ncid, "tsnow_ice",        .False., input_tsnow_ice)
     call NetCDF_read_var(ncid, "snowfallac_land",  .False., input_snowfallac_land)
     call NetCDF_read_var(ncid, "snowfallac_ice",   .False., input_snowfallac_ice)
     call NetCDF_read_var(ncid, "sncovr_ice",       .False., input_sncovr_ice)
     call NetCDF_read_var(ncid, "sfalb_lnd",        .False., input_sfalb_lnd)
     call NetCDF_read_var(ncid, "sfalb_lnd_bck",    .False., input_sfalb_lnd_bck)
     call NetCDF_read_var(ncid, "emis_ice",         .False., input_emis_ice)
     call NetCDF_read_var(ncid, "lai",              .False., input_lai)
  else
     write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that an LSM should be used, but the required initial conditions are missing. Stopping ...'
     stop
  end if
  
  call check(NF90_CLOSE(NCID=ncid),"nf90_close()")
  
  call scm_input%create(input_n_forcing_times, input_n_lev, input_n_soil, input_n_snow, input_n_ice)
  
  !fill the scm_input DDT
  
  !There may need to be logic to control which of the lon, lat, and init_times to use in the future, but for now, just take the first
  active_lon = 1
  active_lat = 1
  active_init_time = 1
  
  rinc(1:5)   = 0
  idat = 0
  jdat = 0
  idat(1) = init_year
  idat(2) = init_month
  idat(3) = init_day
  idat(5) = init_hour
  idat(6) = init_min
  idat(7) = init_sec
  jdat(1) = end_year
  jdat(2) = end_month
  jdat(3) = end_day
  jdat(5) = end_hour
  jdat(6) = end_min
  jdat(7) = end_sec
  call w3difdat(jdat,idat,4,rinc)
  elapsed_sec = rinc(4)
  
  !the following variables replace what is in the case configuration file
  scm_state%init_year = init_year
  scm_state%init_month = init_month
  scm_state%init_day = init_day
  scm_state%init_hour = init_hour
  scm_state%init_min = init_min
  scm_state%runtime = elapsed_sec
  
  scm_input%input_time = input_time
  scm_input%input_pres_surf(1) = input_pres_surf(active_init_time) !perhaps input_pres_surf should only be equal to input_force_pres_surf?
  scm_input%input_pres = input_pres(:,active_init_time)
  scm_input%input_u = input_u(:,active_init_time)
  scm_input%input_v = input_v(:,active_init_time)
  scm_input%input_tke = input_tke(:,active_init_time)
  
  !if mixing ratios are present, and not specific humidities, convert from mixing ratio to specific humidities
  if ((maxval(input_qv(:,active_init_time)) < 0 .and. &
       maxval(input_qt(:,active_init_time)) < 0) .and. &
      (maxval(input_rv(:,active_init_time)) > 0 .or. &
       maxval(input_rt(:,active_init_time)) > 0)) then
     if (maxval(input_rv(:,active_init_time)) > 0) then
       do k=1, input_n_lev
         input_qv(k,active_init_time) = input_rv(k,active_init_time)/&
            (1.0 + input_rv(k,active_init_time))
       end do
     end if
     if (maxval(input_rt(:,active_init_time)) > 0) then
       do k=1, input_n_lev
         input_qt(k,active_init_time) = input_rt(k,active_init_time)/&
            (1.0 + input_rt(k,active_init_time))
       end do
     end if
     if (maxval(input_rl(:,active_init_time)) > 0) then
       do k=1, input_n_lev
         input_ql(k,active_init_time) = input_rl(k,active_init_time)/&
            (1.0 + input_rl(k,active_init_time))
       end do
     end if
     if (maxval(input_ri(:,active_init_time)) > 0) then
       do k=1, input_n_lev
         input_qi(k,active_init_time) = input_ri(k,active_init_time)/&
            (1.0 + input_ri(k,active_init_time))
       end do
     end if
  end if
  
  !make sure that one of qv or qt (and rv or rt due to above conversion) is present (add support for rh later)
  if (maxval(input_qv(:,active_init_time)) >= 0) then
    if (maxval(input_qt(:,active_init_time)) >= 0) then
      if (maxval(input_ql(:,active_init_time)) >= 0) then
        if (maxval(input_qi(:,active_init_time)) >= 0) then
          !all of qv, qt, ql, qi (need to check for consistency that they add up correctly?)
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qt = input_qt(:,active_init_time)
          scm_input%input_ql = input_ql(:,active_init_time)
          scm_input%input_qi = input_qi(:,active_init_time)
        else !qv, qt, ql, but not qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qt = input_qt(:,active_init_time)
          scm_input%input_ql = input_ql(:,active_init_time)
          !derive qi
          do k=1, input_n_lev
            scm_input%input_qi(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_qv(k) - scm_input%input_ql(k))
          end do
        end if !qi test
      else 
        if (maxval(input_qi(:,active_init_time)) >= 0) then !qv, qt, qi, but no ql
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qt = input_qt(:,active_init_time)
          scm_input%input_qi = input_qi(:,active_init_time)
          !derive ql
          do k=1, input_n_lev
            scm_input%input_ql(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_qv(k) - scm_input%input_qi(k))
          end do
        else !qv, qt, no ql or qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qt = input_qt(:,active_init_time)
          !assume that all cloud is liquid for now (could implement partitioning later)
          do k=1, input_n_lev
            scm_input%input_ql(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_qv(k))
          end do
          scm_input%input_qi = 0.0
        end if !qi test
      end if !ql test
    else !qv, but not qt
      if (maxval(input_ql(:,active_init_time)) >= 0) then
        if (maxval(input_qi(:,active_init_time)) >= 0) then !qv, no qt, ql, qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_ql = input_ql(:,active_init_time)
          scm_input%input_qi = input_qi(:,active_init_time)
          !derive qt
          do k=1, input_n_lev
            scm_input%input_qt(k) = max(0.0, scm_input%input_qv(k) + scm_input%input_ql(k) + scm_input%input_qi(k))
          end do
        else !qv, no qt, ql, no qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_ql = input_ql(:,active_init_time)
          !derive qt
          do k=1, input_n_lev
            scm_input%input_qt(k) = max(0.0, scm_input%input_qv(k) + scm_input%input_ql(k))
          end do
          scm_input%input_qi = 0.0          
        end if ! qi test
      else
        if (maxval(input_qi(:,active_init_time)) >= 0) then !qv, no qt, no ql, qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qi = input_qi(:,active_init_time)
          !derive qt
          do k=1, input_n_lev
            scm_input%input_qt(k) = max(0.0, scm_input%input_qv(k) + scm_input%input_qi(k))
          end do
          scm_input%input_ql = 0.0
        else !qv, no qt, no ql, no qi
          scm_input%input_qv = input_qv(:,active_init_time)
          scm_input%input_qt = scm_input%input_qv
          scm_input%input_ql = 0.0
          scm_input%input_qi = 0.0
        end if ! qi test
      end if ! ql test
    end if !qt test
  else if (maxval(input_qt(:,active_init_time)) >= 0) then !qt, but not qv
    if (maxval(input_ql(:,active_init_time)) >= 0) then
      if (maxval(input_qi(:,active_init_time)) >= 0) then !no qv, qt, ql, qi
        scm_input%input_qt = input_qt(:,active_init_time)
        scm_input%input_ql = input_ql(:,active_init_time)
        scm_input%input_qi = input_qi(:,active_init_time)
        !derive qv
        do k=1, input_n_lev
          scm_input%input_qv(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_ql(k) - scm_input%input_qi(k))
        end do
      else !no qv, qt, ql, no qi
        scm_input%input_qt = input_qt(:,active_init_time)
        scm_input%input_ql = input_ql(:,active_init_time)
        !derive qv
        do k=1, input_n_lev
          scm_input%input_qv(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_ql(k))
        end do
        scm_input%input_qi = 0.0
      end if
    else
      if (maxval(input_qi(:,active_init_time)) >= 0) then !no qv, qt, no ql, qi
        scm_input%input_qt = input_qt(:,active_init_time)
        scm_input%input_qi = input_qi(:,active_init_time)
        !derive qv
        do k=1, input_n_lev
          scm_input%input_qv(k) = max(0.0, scm_input%input_qt(k) - scm_input%input_qi(k))
        end do
        scm_input%input_ql = 0.0
      else !no qv, qt, no ql, no qi
        scm_input%input_qt = input_qt(:,active_init_time)
        scm_input%input_qv = scm_input%input_qt
        scm_input%input_ql = 0.0
        scm_input%input_qi = 0.0
      end if
    end if !ql test
  else
    !no qv or qt
    write(*,*) 'When reading '//trim(adjustl(scm_state%case_name))//'.nc, all of the supported moisture variables (qv, qt, rv, rt) were missing. Stopping...'
    stop
  end if
  
  !make sure that at least one of temp, theta, thetal is present;
  !the priority for use is temp, thetal, theta
  if (maxval(input_temp(:,active_init_time)) > 0) then
    scm_input%input_temp = input_temp(:,active_init_time)
    !since temperature was present (and is ultimately needed in the physics), choose to use it, and set the alternative to missing, even if it is also present in the file
    scm_input%input_thetail = missing_value
  else if (maxval(input_thetal(:,active_init_time)) > 0) then
    !convert thetal to thetail
    do k=1, input_n_lev
      exner_inv = (p0/scm_input%input_pres(k))**con_rocp
      scm_input%input_thetail(k) = input_thetal(k,active_init_time) - &
        (con_hfus/con_cp)*exner_inv*(scm_input%input_qi(k)/(1.0 - scm_input%input_qi(k)))
    end do
    !since thetail is present, choose to use it, and set the alternative temperature to missing, even if it is also present in the file
    scm_input%input_temp = missing_value
  else if (maxval(input_theta(:,active_init_time)) > 0) then
    !convert theta to thetail
    do k=1, input_n_lev
      exner_inv = (p0/scm_input%input_pres(k))**con_rocp
      scm_input%input_thetail(k) = input_theta(k,active_init_time) - &
        (con_hvap/con_cp)*exner_inv*(scm_input%input_ql(k)/(1.0 - scm_input%input_ql(k))) - &
        (con_hfus/con_cp)*exner_inv*(scm_input%input_qi(k)/(1.0 - scm_input%input_qi(k)))
    end do
    !since thetail is present, choose to use it, and set the alternative temperature to missing, even if it is also present in the file
    scm_input%input_temp = missing_value
  else
    write(*,*) 'When reading '//trim(adjustl(scm_state%case_name))//'.nc, all of the supported temperature variables (temp, theta, thetal) were missing. Stopping...'
    stop
  end if

  if (trim(input_surfaceForcingLSM) == "lsm") then
    scm_input%input_ozone = input_ozone(:,active_init_time)
    scm_input%input_area = input_area(active_init_time)
    
    scm_input%input_stddev   = input_stddev(active_init_time)
    scm_input%input_convexity= input_convexity(active_init_time)
    scm_input%input_oa1      = input_oa1(active_init_time)
    scm_input%input_oa2      = input_oa2(active_init_time)
    scm_input%input_oa3      = input_oa3(active_init_time)
    scm_input%input_oa4      = input_oa4(active_init_time)
    scm_input%input_ol1      = input_ol1(active_init_time)
    scm_input%input_ol2      = input_ol2(active_init_time)
    scm_input%input_ol3      = input_ol3(active_init_time)
    scm_input%input_ol4      = input_ol4(active_init_time)
    scm_input%input_sigma    = input_sigma(active_init_time)
    scm_input%input_theta    = input_theta_oro(active_init_time)
    scm_input%input_gamma    = input_gamma(active_init_time)
    scm_input%input_elvmax   = input_elvmax(active_init_time)
    scm_input%input_oro      = input_oro(active_init_time)
    scm_input%input_oro_uf   = input_oro_uf(active_init_time)
    scm_input%input_landfrac = input_landfrac(active_init_time)
    scm_input%input_lakefrac = input_lakefrac(active_init_time)
    scm_input%input_lakedepth= input_lakedepth(active_init_time)
    
    scm_input%input_tref    = input_tref(active_init_time)
    scm_input%input_z_c     = input_z_c(active_init_time)
    scm_input%input_c_0     = input_c_0(active_init_time)
    scm_input%input_c_d     = input_c_d(active_init_time)
    scm_input%input_w_0     = input_w_0(active_init_time)
    scm_input%input_w_d     = input_w_d(active_init_time)
    scm_input%input_xt      = input_xt(active_init_time)
    scm_input%input_xs      = input_xs(active_init_time)
    scm_input%input_xu      = input_xu(active_init_time)
    scm_input%input_xv      = input_xv(active_init_time)
    scm_input%input_xz      = input_xz(active_init_time)
    scm_input%input_zm      = input_zm(active_init_time)
    scm_input%input_xtts    = input_xtts(active_init_time)
    scm_input%input_xzts    = input_xzts(active_init_time)
    scm_input%input_d_conv  = input_d_conv(active_init_time)
    scm_input%input_ifd     = input_ifd(active_init_time)
    scm_input%input_dt_cool = input_dt_cool(active_init_time)
    scm_input%input_qrain   = input_qrain(active_init_time)
  else
    !### what to do about ozone??? ### read in standard profile if not included in DEPHY file as part of model ICs?
    scm_input%input_ozone = 0.0
  end if
  scm_input%input_lat = input_lat(active_lat)
  scm_input%input_lon = input_lon(active_lon)
  
  scm_input%input_pres_surf = input_force_pres_surf(:)
  
  do i=1, input_n_forcing_times
    scm_input%input_pres_forcing(i,:) = input_force_pres(:,i)
  end do
    
  if (input_SurfaceType == 'ocean') then
    scm_state%sfc_type = 0.0
  else if (input_SurfaceType == 'land') then
    scm_state%sfc_type = 1.0
  end if
  !no sea ice type?
  
  if (input_surfaceForcingTemp == 'ts') then
    if (maxval(input_force_ts) < 0) then
      write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that the variable ts should be present, but it is missing. Stopping ...'
      stop
    else
      !overwrite sfc_flux_spec
      scm_state%sfc_flux_spec = .false.
      scm_input%input_T_surf = input_force_ts(:)
      scm_state%surface_thermo_control = 2
    end if
  else if (input_surfaceForcingTemp == 'surface_flux') then
    !overwrite sfc_flux_spec
    scm_state%sfc_flux_spec = .true.
    scm_state%surface_thermo_control = 0
    
    if (maxval(input_force_ts) < 0) then
      !since no surface temperature is given, assume that the surface temperature is equivalent to the static, surface-adjacent temperature in the initial profile
      if (maxval(scm_input%input_temp) > 0) then
        !temperature profile is available
        scm_input%input_T_surf = scm_input%input_temp(1)
      else
        !ice-liquid potential temperature profile is available
        exner = (scm_input%input_pres(1)/p0)**con_rocp
        exner_inv = (p0/scm_input%input_pres(1))**con_rocp
        scm_input%input_T_surf = exner*(scm_input%input_thetail(1) + &
          (con_hvap/con_cp)*exner_inv*(scm_input%input_ql(1)/(1.0 - scm_input%input_ql(1))) + &
          (con_hfus/con_cp)*exner_inv*(scm_input%input_qi(1)/(1.0 - scm_input%input_qi(1))))
      end if
    else
      scm_input%input_T_surf = input_force_ts(:)
    end if
    
    !kinematic surface fluxes are specified (but may need to be converted)
    if (maxval(input_force_wpthetap(:)) < missing_value_eps) then
      write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that the variable wpthetap should be present, but it is missing. Stopping ...'
      stop
    else
      !convert from theta to T
      do i=1, input_n_forcing_times
        exner = (scm_input%input_pres_surf(i)/p0)**con_rocp
        scm_input%input_sh_flux_sfc_kin(i) = exner*input_force_wpthetap(i)
      end do
    end if
    
    !if mixing ratios are present, and not specific humidities, convert from mixing ratio to specific humidities
    if ((maxval(input_force_wpqvp(:)) < missing_value_eps .and. &
         maxval(input_force_wpqtp(:)) < missing_value_eps) .and. &
        (maxval(input_force_wprvp(:)) > missing_value_eps .or. &
         maxval(input_force_wprtp(:)) > missing_value_eps)) then
       if (maxval(input_force_wprvp(:)) > missing_value_eps) then
         do i=1, input_n_forcing_times
           input_force_wpqvp(i) = input_force_wprvp(i)/&
              (1.0 + input_force_wprvp(i))
         end do
       end if
       if (maxval(input_force_wprtp(:)) > missing_value_eps) then
         do i=1, input_n_forcing_times
           input_force_wpqtp(i) = input_force_wprtp(i)/&
              (1.0 + input_force_wprtp(i))
         end do
       end if
    end if
    
    if (maxval(input_force_wpqvp(:)) < missing_value_eps .and. maxval(input_force_wpqtp(:)) < missing_value_eps) then
      write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that the variable wpqvp, wpqtp, wprvp, or wprtp should be present, but all are missing. Stopping ...'
      stop
    else
      if (maxval(input_force_wpqvp(:)) > missing_value_eps) then !use wpqvp if available
        scm_input%input_lh_flux_sfc_kin = input_force_wpqvp(:)
      else
        !surface total flux of water should just be vapor
        scm_input%input_lh_flux_sfc_kin = input_force_wpqtp(:)
      end if
    end if
  else if (input_surfaceForcingTemp == 'surface_flux') then
    !overwrite sfc_flux_spec
    scm_state%sfc_flux_spec = .true.
    scm_state%surface_thermo_control = 1
    
    if (maxval(input_force_ts) < 0) then
      !since no surface temperature is given, assume that the surface temperature is equivalent to the static, surface-adjacent temperature in the initial profile
      if (maxval(scm_input%input_temp) > 0) then
        !temperature profile is available
        scm_input%input_T_surf = scm_input%input_temp(1)
      else
        !ice-liquid potential temperature profile is available
        exner = (scm_input%input_pres(1)/p0)**con_rocp
        exner_inv = (p0/scm_input%input_pres(1))**con_rocp
        scm_input%input_T_surf = exner*(scm_input%input_thetail(1) + &
          (con_hvap/con_cp)*exner_inv*(scm_input%input_ql(1)/(1.0 - scm_input%input_ql(1))) + &
          (con_hfus/con_cp)*exner_inv*(scm_input%input_qi(1)/(1.0 - scm_input%input_qi(1))))
      end if
    else
      scm_input%input_T_surf = input_force_ts(:)
    end if
    
    
    if (maxval(input_force_sfc_sens_flx(:)) < missing_value_eps) then
      write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that the variable sfc_sens_flx should be present, but it is missing. Stopping ...'
      stop
    else
      scm_input%input_sh_flux_sfc = input_force_sfc_sens_flx(:)
    end if
    
    if (maxval(input_force_sfc_lat_flx(:)) < missing_value_eps) then
      write(*,*) 'The global attribute surfaceForcing in '//trim(adjustl(scm_state%case_name))//'.nc indicates that the variable sfc_lat_flx should be present, but it is missing. Stopping ...'
      stop
    else
      scm_input%input_lh_flux_sfc = input_force_sfc_lat_flx(:)
    end if
  else if (trim(input_surfaceForcingLSM) == 'lsm') then
    !these were considered required variables above, so they should not need to be checked for missing
    scm_input%input_stc   = input_stc(:,active_init_time)
    scm_input%input_smc   = input_smc(:,active_init_time)  
    scm_input%input_slc   = input_slc(:,active_init_time)  
    
    scm_input%input_snicexy    = input_snicexy(:,active_init_time)
    scm_input%input_snliqxy    = input_snliqxy(:,active_init_time)
    scm_input%input_tsnoxy     = input_tsnoxy(:,active_init_time)
    scm_input%input_smoiseq    = input_smoiseq(:,active_init_time)
    scm_input%input_zsnsoxy    = input_zsnsoxy(:,active_init_time)
    
    scm_input%input_tiice      = input_tiice(:,active_init_time)
    scm_input%input_tslb       = input_tslb(:,active_init_time)
    scm_input%input_smois      = input_smois(:,active_init_time)
    scm_input%input_sh2o       = input_sh2o(:,active_init_time)
    scm_input%input_smfr       = input_smfr(:,active_init_time)
    scm_input%input_flfr       = input_flfr(:,active_init_time)
    
    scm_input%input_vegsrc   = input_vegsrc(active_init_time)
    scm_input%input_vegtyp   = REAL(input_vegtyp(active_init_time), kind=dp)
    scm_input%input_soiltyp  = REAL(input_soiltyp(active_init_time), kind=dp)
    scm_input%input_slopetype = REAL(input_slopetype(active_init_time), kind=dp)
    scm_input%input_tsfco    = input_tsfco(active_init_time)
    scm_input%input_vegfrac  = input_vegfrac(active_init_time)
    scm_input%input_shdmin   = input_shdmin(active_init_time)
    scm_input%input_shdmax   = input_shdmax(active_init_time)
    scm_input%input_slmsk    = input_slmsk(active_init_time)
    scm_input%input_canopy   = input_canopy(active_init_time)
    scm_input%input_hice     = input_hice(active_init_time)
    scm_input%input_fice     = input_fice(active_init_time)
    scm_input%input_tisfc    = input_tisfc(active_init_time)
    scm_input%input_snwdph   = input_snwdph(active_init_time)
    scm_input%input_snoalb   = input_snoalb(active_init_time)
    scm_input%input_sncovr   = input_sncovr(active_init_time)
    scm_input%input_tg3      = input_tg3(active_init_time)
    scm_input%input_uustar   = input_uustar(active_init_time)
    scm_input%input_alvsf    = input_alvsf(active_init_time)
    scm_input%input_alnsf    = input_alnsf(active_init_time)
    scm_input%input_alvwf    = input_alvwf(active_init_time)
    scm_input%input_alnwf    = input_alnwf(active_init_time)
    scm_input%input_facsf    = input_facsf(active_init_time)
    scm_input%input_facwf    = input_facwf(active_init_time)
    scm_input%input_weasd    = input_weasd(active_init_time)
    scm_input%input_f10m     = input_f10m(active_init_time)
    scm_input%input_t2m      = input_t2m(active_init_time)
    scm_input%input_q2m      = input_q2m(active_init_time)
    scm_input%input_ffmm     = input_ffmm(active_init_time)
    scm_input%input_ffhh     = input_ffhh(active_init_time)
    scm_input%input_tprcp    = input_tprcp(active_init_time)
    scm_input%input_srflag   = input_srflag(active_init_time)
    scm_input%input_tsfcl    = input_tsfcl(active_init_time)
    scm_input%input_zorll    = input_zorll(active_init_time)
    scm_input%input_zorli    = input_zorli(active_init_time)
    scm_input%input_zorlw    = input_zorlw(active_init_time)
    
    scm_input%input_tvxy     = input_tvxy(active_init_time)
    scm_input%input_tgxy     = input_tgxy(active_init_time)
    scm_input%input_tahxy    = input_tahxy(active_init_time)
    scm_input%input_canicexy = input_canicexy(active_init_time)
    scm_input%input_canliqxy = input_canliqxy(active_init_time)
    scm_input%input_eahxy    = input_eahxy(active_init_time)
    scm_input%input_cmxy     = input_cmxy(active_init_time)
    scm_input%input_chxy     = input_chxy(active_init_time)
    scm_input%input_fwetxy   = input_fwetxy(active_init_time)
    scm_input%input_sneqvoxy = input_sneqvoxy(active_init_time)
    scm_input%input_alboldxy = input_alboldxy(active_init_time)
    scm_input%input_qsnowxy  = input_qsnowxy(active_init_time)
    scm_input%input_wslakexy = input_wslakexy(active_init_time)
    scm_input%input_taussxy  = input_taussxy(active_init_time)
    scm_input%input_waxy     = input_waxy(active_init_time)
    scm_input%input_wtxy     = input_wtxy(active_init_time)
    scm_input%input_zwtxy    = input_zwtxy(active_init_time)
    scm_input%input_xlaixy   = input_xlaixy(active_init_time)
    scm_input%input_xsaixy   = input_xsaixy(active_init_time)
    scm_input%input_lfmassxy = input_lfmassxy(active_init_time)
    scm_input%input_stmassxy = input_stmassxy(active_init_time)
    scm_input%input_rtmassxy = input_rtmassxy(active_init_time)
    scm_input%input_woodxy   = input_woodxy(active_init_time)
    scm_input%input_stblcpxy = input_stblcpxy(active_init_time)
    scm_input%input_fastcpxy = input_fastcpxy(active_init_time)
    scm_input%input_smcwtdxy = input_smcwtdxy(active_init_time)
    scm_input%input_deeprechxy = input_deeprechxy(active_init_time)
    scm_input%input_rechxy   = input_rechxy(active_init_time)
    scm_input%input_snowxy   = input_snowxy(active_init_time)
    scm_input%input_wetness    = input_wetness(active_init_time)
    scm_input%input_lai        = input_lai(active_init_time)
    scm_input%input_clw_surf_land   = input_clw_surf_land(active_init_time)
    scm_input%input_clw_surf_ice    = input_clw_surf_ice(active_init_time)
    scm_input%input_qwv_surf_land   = input_qwv_surf_land(active_init_time)
    scm_input%input_qwv_surf_ice    = input_qwv_surf_ice(active_init_time)
    scm_input%input_tsnow_land      = input_tsnow_land(active_init_time)
    scm_input%input_tsnow_ice       = input_tsnow_ice(active_init_time)
    scm_input%input_snowfallac_land = input_snowfallac_land(active_init_time)
    scm_input%input_snowfallac_ice  = input_snowfallac_ice(active_init_time)
    scm_input%input_sncovr_ice      = input_sncovr_ice(active_init_time)
    scm_input%input_sfalb_lnd       = input_sfalb_lnd(active_init_time)
    scm_input%input_sfalb_lnd_bck   = input_sfalb_lnd_bck(active_init_time)
    scm_input%input_emis_ice        = input_emis_ice(active_init_time)
  end if
  
  if (input_surfaceForcingWind == 'z0') then
    scm_state%surface_momentum_control = 0
    scm_state%sfc_roughness_length_cm = input_z0*100.0 !convert from m to cm
  else if (input_surfaceForcingWind == 'ustar') then
    !not supported
    scm_state%surface_momentum_control = 1
    write(*,*) 'The global attribute surfaceForcingWind in '//trim(adjustl(scm_state%case_name))//'.nc indicates that surface wind is controlled by a specified time-series of ustar. This is currently not supported. Stopping ...'
    stop
  end if
  
  if (forc_omega > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_omega(i,:) = input_force_omega(:,i)
    end do
    scm_state%force_omega = .true.
    scm_state%force_w = .false. ! only one of forc_w, forc_omega should be true, with forc_omega having higher priority
    !set all individual w forcing controls to .true. until finer control is available from the input file
    scm_state%force_sub_for_T = .true.
    scm_state%force_sub_for_qv = .true.
    scm_state%force_sub_for_u = .true.
    scm_state%force_sub_for_v = .true.
  else if (forc_w > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_w_ls(i,:) = input_force_w(:,i)
    end do
    scm_state%force_w = .true.
    scm_state%force_omega = .false.
    !set all individual w forcing controls to .true. until finer control is available from the input file
    scm_state%force_sub_for_T = .true.
    scm_state%force_sub_for_qv = .true.
    scm_state%force_sub_for_u = .true.
    scm_state%force_sub_for_v = .true.
  end if

  if (forc_geo > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_u_g(i,:) = input_force_u_g(:,i)
      scm_input%input_v_g(i,:) = input_force_v_g(:,i)
    end do
    scm_state%force_geo = .true.
  end if
  
  if (adv_temp > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_T(i,:) = input_force_temp_adv(:,i)
    end do
    scm_state%force_adv_T = 1
  else if (adv_theta > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_theta(i,:) = input_force_theta_adv(:,i)
    end do
    scm_state%force_adv_T = 2
  else if (adv_thetal > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_thetal(i,:) = input_force_thetal_adv(:,i)
    end do
    scm_state%force_adv_T = 3
  end if
  
  if (adv_qv > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_qv(i,:) = input_force_qv_adv(:,i)
    end do
    scm_state%force_adv_qv = .true.
  else if (adv_qt > 0) then
    !since there is no information about individual advected species, assume it is all vapor
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_qv(i,:) = input_force_qt_adv(:,i)
    end do
    scm_state%force_adv_qv = .true.
  else if (adv_rv > 0 .or. adv_rt > 0) then
    !convert to specific humidity
    if (adv_rv > 0) then
      do i=1, input_n_forcing_times
        do k=1, input_n_lev
          scm_input%input_tot_advec_qv(i,k) = input_force_rv_adv(k,i)/&
            (1.0 + input_force_rv_adv(k,i))
        end do
      end do
    else if (adv_rt > 0) then
      !since there is no information about individual advected species, assume it is all vapor
      do i=1, input_n_forcing_times
        do k=1, input_n_lev
          scm_input%input_tot_advec_qv(i,k) = input_force_rt_adv(k,i)/&
            (1.0 + input_force_rt_adv(k,i))
        end do
      end do
    end if
    scm_state%force_adv_qv = .true.
  end if
  
  if (adv_u > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_u(i,:) = input_force_u_adv(:,i)
    end do
    scm_state%force_adv_u = .true.
  end if
  
  if (adv_v > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_tot_advec_v(i,:) = input_force_v_adv(:,i)
    end do
    scm_state%force_adv_v = .true.
  end if
  
  if (char_rad_temp == 'adv' .or. char_rad_theta == 'adv' .or. char_rad_thetal == 'adv') then
    scm_state%force_rad_T = 4
    if (scm_state%force_adv_T == 0) then
      write(*,*) 'The global attribute rad_temp, rad_theta, or rad_thetal in '//trim(adjustl(scm_state%case_name))//'.nc indicates that radiative forcing is included in the advection term, but there is no advection term. Stopping ...'
      stop
    end if
  else if (rad_temp > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_dT_dt_rad(i,:) = input_force_temp_rad(:,i)
    end do
    scm_state%force_rad_T = 1
  else if (rad_theta > 0) then
    scm_state%force_rad_T = 2
    do i=1, input_n_forcing_times
      do k=1, input_n_lev
        exner = (scm_input%input_pres(k)/p0)**con_rocp
        scm_input%input_dT_dt_rad(i,k) = exner*input_force_theta_rad(k,i)
      end do
    end do
  else if (rad_thetal > 0) then
    scm_state%force_rad_T = 3
    do i=1, input_n_forcing_times
      do k=1, input_n_lev
        exner = (scm_input%input_pres(k)/p0)**con_rocp
        scm_input%input_dT_dt_rad(i,k) = exner*input_force_thetal_rad(k,i)
      end do
    end do
  else
    scm_state%force_rad_T = 0
  end if
  
  if (nudging_temp > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_T_nudge(i,:) = input_force_temp_nudging(:,i)
    end do
    scm_state%force_nudging_T = 1
    scm_state%force_nudging_T_time = nudging_temp
    if (p_nudging_temp > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_temp, input_force_pres(:,i), scm_input%input_k_T_nudge(i))
        if (scm_input%input_k_T_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_T_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_temp > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_temp, input_force_height(:,i), scm_input%input_k_T_nudge(i))
        if (scm_input%input_k_T_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_T_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_T_nudge = 1
    end if
  else if (nudging_theta > 0) then
    !assume no cloud water since there is no associate [ql,qi]_nudge in the input?
    do i=1, input_n_forcing_times
      scm_input%input_thil_nudge(i,:) = input_force_theta_nudging(:,i)
    end do
    scm_state%force_nudging_T = 2
    scm_state%force_nudging_T_time = nudging_theta
    if (p_nudging_theta > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_theta, input_force_pres(:,i), scm_input%input_k_thil_nudge(i))
        if (scm_input%input_k_thil_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_thil_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_theta > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_theta, input_force_height(:,i), scm_input%input_k_thil_nudge(i))
        if (scm_input%input_k_thil_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_thil_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_thil_nudge = 1
    end if
  else if (nudging_thetal > 0) then
    !assume no cloud water since there is no associate [ql,qi]_nudge in the input?
    do i=1, input_n_forcing_times
      scm_input%input_thil_nudge(i,:) = input_force_thetal_nudging(:,i)
    end do
    scm_state%force_nudging_T = 3
    scm_state%force_nudging_T_time = nudging_thetal
    if (p_nudging_thetal > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_thetal, input_force_height(:,i), scm_input%input_k_thil_nudge(i))
        if (scm_input%input_k_thil_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_thil_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_thetal > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_thetal, input_force_height(:,i), scm_input%input_k_thil_nudge(i))
        if (scm_input%input_k_thil_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_thil_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_thil_nudge = 1
    end if
  end if
  
  if (nudging_qv > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_qt_nudge(i,:) = input_force_qv_nudging(:,i)
    end do
    scm_state%force_nudging_qv = .true.
    scm_state%force_nudging_qv_time = nudging_qv
    if (p_nudging_qv > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_qv, input_force_pres(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_qv > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_qv, input_force_height(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_qt_nudge = 1
    end if
  else if (nudging_qt > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_qt_nudge(i,:) = input_force_qt_nudging(:,i)
    end do
    scm_state%force_nudging_qv = .true.
    scm_state%force_nudging_qv_time = nudging_qt
    if (p_nudging_qt > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_qt, input_force_pres(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_qt > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_qt, input_force_height(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_qt_nudge = 1
    end if
  else if (nudging_rv > 0) then
    do i=1, input_n_forcing_times
      do k=1, input_n_lev
        scm_input%input_qt_nudge(i,:) = input_force_rv_nudging(:,i)/&
          (1.0 + input_force_rv_nudging(:,i))
      end do
    end do
    scm_state%force_nudging_qv = .true.
    scm_state%force_nudging_qv_time = nudging_rv
    if (p_nudging_rv > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_rv, input_force_pres(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_rv > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_rv, input_force_height(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_qt_nudge = 1
    end if
  else if (nudging_rt > 0) then
    do i=1, input_n_forcing_times
      do k=1, input_n_lev
        scm_input%input_qt_nudge(i,:) = input_force_rt_nudging(:,i)/&
          (1.0 + input_force_rt_nudging(:,i))
      end do
    end do
    scm_state%force_nudging_qv = .true.
    scm_state%force_nudging_qv_time = nudging_rt
    if (p_nudging_rt > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_rt, input_force_pres(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_rt > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_rt, input_force_height(:,i), scm_input%input_k_qt_nudge(i))
        if (scm_input%input_k_qt_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_qt_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_qt_nudge = 1
    end if
  end if
  
  if (nudging_u > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_u_nudge(i,:) = input_force_u_nudging(:,i)
    end do
    scm_state%force_nudging_u = .true.
    scm_state%force_nudging_u_time = nudging_u
    if (p_nudging_u > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_u, input_force_pres(:,i), scm_input%input_k_u_nudge(i))
        if (scm_input%input_k_u_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_u_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_u > 0) then 
      do i=1, input_n_forcing_times 
        call find_vertical_index_height(z_nudging_u, input_force_height(:,i), scm_input%input_k_u_nudge(i))
        if (scm_input%input_k_u_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_u_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_u_nudge = 1
    end if
  end if
  
  if (nudging_v > 0) then
    do i=1, input_n_forcing_times
      scm_input%input_v_nudge(i,:) = input_force_v_nudging(:,i)
    end do
    scm_state%force_nudging_v = .true.
    scm_state%force_nudging_v_time = nudging_v
    if (p_nudging_v > 0) then
      do i=1, input_n_forcing_times
        call find_vertical_index_pressure(p_nudging_v, input_force_pres(:,i), scm_input%input_k_v_nudge(i))
        if (scm_input%input_k_v_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_v_nudge(i) = input_n_lev
        end if
      end do
    else if (z_nudging_v > 0) then 
      do i=1, input_n_forcing_times
        call find_vertical_index_height(z_nudging_v, input_force_height(:,i), scm_input%input_k_v_nudge(i))
        if (scm_input%input_k_v_nudge(i) < 0) then
          !if the vertical index is not found (when it is less than 0), set the nudging index to the top of the input profile so that nudging is turned off
          scm_input%input_k_v_nudge(i) = input_n_lev
        end if
      end do
    else
      scm_input%input_k_v_nudge = 1
    end if
  end if
  
end subroutine get_case_init_DEPHY

!> Subroutine to get reference profile to use above the case data (temporarily hard-coded profile)
subroutine get_reference_profile(scm_state, scm_reference)
  use scm_type_defs, only : scm_state_type, scm_reference_type
  use NetCDF_read, only: check
  
  type(scm_state_type), target, intent(in) :: scm_state
  type(scm_reference_type), target, intent(inout) :: scm_reference

  integer  :: nlev !< number of pressure levels in the reference profile
  real(kind=dp), allocatable :: pres(:)  !< reference profile pressure (Pa)
  real(kind=dp), allocatable :: T(:) !< reference profile temperature (K)
  real(kind=dp), allocatable :: qv(:) !< reference profile specific humidity (kg kg^-1)
  real(kind=dp), allocatable :: ozone(:) !< reference profile ozone concentration (kg kg^-1)

  integer :: i, ioerror
  character(len=120)                 :: line
  real :: dummy

  integer                           :: ncid, varID, allocate_status
  CHARACTER(LEN=nf90_max_name)      :: tmpName

  select case (scm_state%reference_profile_choice)
    case (1)
      open(unit=1, file='McCProfiles.dat', status='old', action='read', iostat=ioerror)
      if(ioerror /= 0) then
        write(*,*) 'There was an error opening the file McCprofiles.dat in the processed_case_input directory. &
          Error code = ',ioerror
        stop
      endif

      ! find number of records
      read(1, '(a)', iostat=ioerror) line !first line is header
      nlev = 0
      do
        read(1, '(a)', iostat=ioerror) line
        if (ioerror /= 0) exit
        nlev = nlev + 1
      end do

      allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev))

      rewind(1)
      read(1, '(a)', iostat=ioerror) line !first line is header
      do i=1, nlev
        read(1,'(5ES16.4)', iostat=ioerror) dummy, pres(i), T(i), qv(i), ozone(i)
      END DO
      close(1)
    case (2)
      call check(NF90_OPEN('mid_lat_summer_std.nc',nf90_nowrite,ncid),"nf90_open()")

      call check(NF90_INQ_DIMID(ncid,"height",varID),"nf90_inq_dimid(height)")
      call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, nlev),"nf90_inq_dim(height)")

      !> - Allocate the dimension variables.
      allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev), stat=allocate_status)

      call check(NF90_INQ_VARID(ncid,"pressure",varID),"nf90_inq_varid(pressure)")
      call check(NF90_GET_VAR(ncid,varID,pres),"nf90_get_var(pressure)")
      call check(NF90_INQ_VARID(ncid,"temperature",varID),"nf90_inq_varid(temperature)")
      call check(NF90_GET_VAR(ncid,varID,T),"nf90_get_var(temperature)")
      call check(NF90_INQ_VARID(ncid,"q_v",varID),"nf90_inq_varid(q_v)")
      call check(NF90_GET_VAR(ncid,varID,qv),"nf90_get_var(q_v)")
      call check(NF90_INQ_VARID(ncid,"o3",varID),"nf90_inq_varid(o3)")
      call check(NF90_GET_VAR(ncid,varID,ozone),"nf90_get_var(o3)")

      call check(NF90_CLOSE(NCID=ncid),"nf90_close()")

  end select

  call scm_reference%create(nlev)

  scm_reference%ref_nlev = nlev
  scm_reference%ref_pres = pres
  scm_reference%ref_T = T
  scm_reference%ref_qv = qv
  scm_reference%ref_ozone = ozone

end subroutine get_reference_profile

!> Subroutine to get reference profile to use above the case data (temporarily hard-coded profile)
subroutine get_reference_profile_old(nlev, pres, T, qv, ozone)
  integer, intent(out)  :: nlev !< number of pressure levels in the reference profile
  real(kind=dp), allocatable, intent(out) :: pres(:)  !< reference profile pressure (Pa)
  real(kind=dp), allocatable, intent(out) :: T(:) !< reference profile temperature (K)
  real(kind=dp), allocatable, intent(out) :: qv(:) !< reference profile specific humidity (kg kg^-1)
  real(kind=dp), allocatable, intent(out) :: ozone(:) !< reference profile ozone concentration (kg kg^-1)

  !> \todo write a more sophisticated reference profile subroutine (can choose between reference profiles)

  !> - For the prototype, the 'McClatchey' sounding used in Jennifer Fletcher's code is hardcoded.
  nlev = 20

  allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev))

  pres = (/ 1030.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 281.0, 209.0, 130.0, 59.5, 27.7, 13.2, 6.52, 3.33, &
    0.951, 0.0671, 0.000300 /)
  pres = pres*100.0

  T = (/ 294., 290., 285., 279., 273., 267., 261., 255., 248., 235., 222., 216., 218., 224., 234., 245., 258., 276., 218., 210. /)

  qv = (/ 11.75, 8.611, 6.047, 3.877, 2.363, 1.387, 0.9388, 0.6364, 0.4019, 0.1546, 0.01976, 0.004002, 0.003999, 0.004011, &
    0.004002, 0.004004, 0.003994, 0.003995, 0.003996, 0.004000 /)
  qv = qv*1.0E-3

  ozone = (/ 6., 6., 6., 6.2, 6.4, 6.6, 6.9, 7.5, 7.9, 9., 12., 19., 34., 30., 20., 9.2, 4.1, 0.43, 0.0086, 0.0000043 /)
  ozone = ozone*1.0E-5

end subroutine get_reference_profile_old

subroutine get_tracers(tracer_names, tracer_types)
  character(len=character_length), allocatable, intent(inout), dimension(:) :: tracer_names
  integer,                         allocatable, intent(inout), dimension(:) :: tracer_types

  character(len=*), parameter :: file_name = 'tracers.txt'

  character(len=100) :: name!, std_name, units
  integer            :: i, fu, rc, n_lines

  open (action='read', file=FILE_NAME, iostat=rc, newunit=fu)
    if (rc == 0) then
        n_lines = 0
        do
            read (fu, *, iostat=rc) name!, std_name, units
            if (rc /= 0) exit
            n_lines = n_lines + 1
        end do
        allocate(tracer_names(n_lines),tracer_types(n_lines))
        rewind(fu)
        do i=1,n_lines
            read (fu, *, iostat=rc) name!, std_name, units
            if (rc /= 0) exit
            tracer_names(i) = trim(name)
            tracer_types(i) = 0 ! temporary until SCM is configured to work with GOCART
        end do
    else
        write(*,'(a,i0)') 'There was an error opening the file ' // FILE_NAME // &
                          '; error code = ', rc
        stop
    end if

    close (fu)
end subroutine get_tracers

!> @}
!> @}
end module scm_input
