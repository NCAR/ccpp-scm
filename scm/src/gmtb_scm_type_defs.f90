!> \file gmtb_scm_type_defs.f90
!!  Contains type definitions for SCM-related variables and physics-related variables

module gmtb_scm_type_defs

!> \section arg_table_gmtb_scm_type_defs
!! \htmlinclude gmtb_scm_type_defs.html
!!

  use gmtb_scm_kinds, only : sp, dp, qp
  use GFS_typedefs, only: GFS_control_type, GFS_statein_type, GFS_stateout_type, GFS_sfcprop_type, GFS_coupling_type, &
    GFS_grid_type, GFS_tbd_type, GFS_cldprop_type, GFS_radtend_type, GFS_diag_type, GFS_interstitial_type, &
    GFS_init_type
  use machine, only: kind_phys

  implicit none

  integer, parameter :: character_length = 80
  integer, parameter :: int_zero = 0
  integer, parameter :: int_one = 1
  integer, parameter :: int_neg_one = -1
  real(kind=dp), parameter :: real_zero = 0.0
  real(kind=dp), parameter :: real_one = 1.0

  character(len = character_length) :: clear_char = ''


  type scm_state_type

    character(len=character_length)                 :: experiment_name !> name of model configuration file
    character(len=character_length)                 :: model_name !< name of "host" model (must be "GFS" for prototype)
    character(len=character_length)                 :: output_dir !< name of output directory to place netCDF file
    character(len=character_length)                 :: case_data_dir !< location of the case initialization and forcing data files (relative to the executable path)
    character(len=character_length)                 :: vert_coord_data_dir !< location of the vertical coordinate data files (relative to the executable path)
    character(len=character_length)                 :: output_file !< name of output file (without the file extension)
    character(len=character_length)                 :: case_name !< name of case initialization and forcing to use (different than experiment name, which names the model run (as a control, experiment_1, etc.))
    character(len=character_length), allocatable    :: physics_suite_name(:) !< name of physics suite (must be "GFS_operational" for prototype)
    character(len=64), allocatable                  :: physics_nml(:)

    integer                           :: n_levels !< number of model levels (must be 64 for prototype)
    integer                           :: n_nsoil  !< number of model levels (must be 4 for prototype)
    integer                           :: itt !< current model iteration
    integer                           :: itt_out  !< output iteration counter
    integer                           :: time_scheme !< 1=> forward Euler, 2=> filtered leapfrog
    integer                           :: n_cols !< number of columns
    integer                           :: n_timesteps !< number of timesteps needed to integrate over runtime
    integer                           :: n_time_levels !< number of time levels to keep track of for time-integration scheme (2 for leapfrog)
    integer                           :: n_itt_swrad !< number of iterations between calls to SW rad
    integer                           :: n_itt_lwrad !< number of iterations between calls to LW rad
    integer                           :: n_itt_out !< number of iterations between calls to write the output
    integer                           :: n_levels_smooth !< the number of levels over which the input profiles are smoothed into the reference profiles
    integer                           :: n_tracers !< number of tracers
    integer                           :: water_vapor_index
    integer                           :: ozone_index  !< index for ozone in the tracer array
    integer                           :: cloud_water_index  !< index for cloud water in the tracer array
    integer                           :: cloud_ice_index  !< index for ice cloud water in the tracer array
    integer                           :: rain_index  !< index for rain water in the tracer array
    integer                           :: snow_index   !< index for snow water in the tracer array
    integer                           :: graupel_index    !< index for graupel water in the tracer array
    integer                           :: cloud_amount_index   !< index for cloud amount in the tracer array
    integer                           :: cloud_droplet_nc_index !< index for liquid cloud droplet number concentration in the tracer array
    integer                           :: cloud_ice_nc_index !< index for ice cloud particle number concentration in the tracer array
    integer                           :: rain_nc_index !< index for rain number concentration in the tracer array
    integer                           :: snow_nc_index !< index for snow number concentration in the tracer array
    integer                           :: graupel_nc_index !< index for graupel number concentration in the tracer array
    integer                           :: tke_index !< index for TKE in the tracer array
    integer                           :: water_friendly_aerosol_index !< index for water-friendly aerosols in the tracer array
    integer                           :: ice_friendly_aerosol_index !< index for ice-friendly aerosols in the tracer array
    integer                           :: init_year, init_month, init_day, init_hour
    character(len=32), allocatable    :: tracer_names(:) !< name of physics suite (must be "GFS_operational" for prototype)
    integer, allocatable              :: blksz(:)

    logical                           :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
    integer                           :: sfc_type !< 0: sea surface, 1: land surface, 2: sea-ice surface
    real(kind=dp)                     :: sfc_type_real(1)
    integer                           :: veg_type !
    real(kind=dp)                     :: veg_type_real(1)
    integer                           :: veg_frac !< 0: sea surface, 1: land surface, 2: sea-ice surface
    real(kind=dp)                     :: veg_frac_real(1)
    real(kind=dp)                     :: slopetype(1)
    real(kind=dp)                     :: shdmin(1)
    real(kind=dp)                     :: shdmax(1)
    real(kind=dp)                     :: tg3(1)
    real(kind=dp)                     :: slmsk(1)
    real(kind=dp)                     :: canopy(1)
    real(kind=dp)                     :: hice(1)
    real(kind=dp)                     :: fice(1)
    real(kind=dp)                     :: tisfc(1)
    real(kind=dp)                     :: snwdph(1)
    real(kind=dp)                     :: snoalb(1)
    real(kind=dp)                     :: sncovr(1)
    real(kind=dp)                     :: uustar(1)
    integer                           :: soil_type !< 1-14?
    real(kind=dp)                     :: soil_type_real(1)
    integer                           :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
    integer                           :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
    integer                           :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere
    integer                           :: C_RES !< effective model resolution in cubed sphere resolution (needed for GWD)
    logical                           :: model_ics !<  true means have land info too

    real(kind=dp)                           :: model_time !< elapsed model time (s)
    real(kind=dp)                           :: dt !< physics time step (s)
    real(kind=dp)                           :: dt_now !< time step currently being used (if it changes due to time-stepping scheme)
    real(kind=dp)                           :: runtime !< total runtime (s)
    real(kind=dp)                           :: output_frequency !< how often output is written (s)
    real(kind=dp)                           :: relax_time !< time scale for hor. wind nudging (s)
    real(kind=dp)                           :: deg_to_rad_const !< conversion constant from degrees to radians
    real(kind=dp)                           :: c_filter !< parameter that controls the amount of damping in the leapfrog filter

    !> - Define the SCM state variables; variables with appended "i" are interface; variables with appended "l" are layer-centered.
    !!  - index order for grid is (horizontal, vertical);
    !!  - index order for state variables is (horizontal, vertical, timesteps);
    !!  - index order for tracer is (horizontal, vertical, tracer_index, timesteps)
    real(kind=dp), allocatable              :: pres_i(:,:,:), pres_l(:,:,:) !< pressure on grid interfaces, centers (Pa)
    real(kind=dp), allocatable              :: si(:,:,:), sl(:,:,:) !< sigma on grid interfaces, centers
    real(kind=dp), allocatable              :: exner_i(:,:,:), exner_l(:,:,:) !< exner function on grid interfaces, centers
    real(kind=dp), allocatable              :: geopotential_i(:,:,:), geopotential_l(:,:,:) !< geopotential on grid interfaces, centers
    real(kind=dp), allocatable              :: a_k(:,:), b_k(:,:) !< used to determine grid sigma and pressure levels

    real(kind=dp), allocatable              :: lat(:,:), lon(:,:) !< latitude and longitude (radians)
    real(kind=dp), allocatable              :: area(:,:) !< area over which the column represents a mean (analogous to grid size or observational array area)

    real(kind=dp), allocatable              :: state_T(:,:,:,:) !< model state absolute temperature at grid centers (K)
    real(kind=dp), allocatable              :: state_u(:,:,:,:), state_v(:,:,:,:) !< model state horizontal winds at grid centers (m/s)
    real(kind=dp), allocatable              :: state_tracer(:,:,:,:,:) !< model state tracer at grid centers
    real(kind=dp), allocatable              :: temp_T(:,:,:,:), temp_u(:,:,:,:), temp_v(:,:,:,:), temp_tracer(:,:,:,:,:) !< used for time-filtering

    !> - Define forcing-related variables (indexing is (horizontal, vertical)).
    real(kind=dp), allocatable              :: u_force_tend(:,:), v_force_tend(:,:), T_force_tend(:,:), qv_force_tend(:,:) !< total u, v, T, q forcing (units/s) (horizontal, vertical)
    real(kind=dp), allocatable              :: w_ls(:,:), omega(:,:,:), u_g(:,:), v_g(:,:), dT_dt_rad(:,:), h_advec_thil(:,:), &
      h_advec_qt(:,:), v_advec_thil(:,:), v_advec_qt(:,:), u_nudge(:,:), v_nudge(:,:), T_nudge(:,:), thil_nudge(:,:), qt_nudge(:,:) !< forcing terms interpolated to the model time and grid
    real(kind=dp), allocatable              :: T_surf(:,:), pres_surf(:,:) !< surface temperature and pressure interpolated to the model time
    real(kind=dp), allocatable              :: sh_flux(:), lh_flux(:) !< surface sensible and latent heat fluxes interpolated to the model time
    real(kind=dp), allocatable              :: sfc_roughness_length_cm(:) !< surface roughness length used for calculating surface layer parameters from specified fluxes
    real(kind=dp), allocatable              :: alvsf(:), alnsf(:),alvwf(:),alnwf(:) !< surface  albedos
    real(kind=dp), allocatable              :: facsf(:), facwf(:),stddev(:),hprime(:,:)          !< other surface stuff
    real(kind=dp), allocatable              :: stc(:,:,:,:) !< soil temperature 
    real(kind=dp), allocatable              :: smc(:,:,:,:) !< soil moisture
    real(kind=dp), allocatable              :: slc(:,:,:,:) !< soil liquid content

    contains
      procedure :: create  => scm_state_create

  end type scm_state_type

  type scm_input_type
    !> - Define the case-specific initialization and forcing variables.
    integer                           :: input_nlev !< number of levels in the input file
    integer                           :: input_nsoil !< number of soil levels in the input file
    integer                           :: input_ntimes !< number of times in the input file where forcing is available
    integer                           :: input_vegsrc !<
    integer                           :: input_vegtyp !<
    integer                           :: input_soiltyp !<
    integer                           :: input_slopetype !<
    real(kind=dp)                     :: input_vegfrac  !<
    real(kind=dp)                     :: input_shdmin   !<
    real(kind=dp)                     :: input_shdmax   !<
    real(kind=dp)                     :: input_zorl     !<
    real(kind=dp)                     :: input_slmsk   !< sea land ice mask [0,1,2]
    real(kind=dp)                     :: input_canopy  !< amount of water stored in canopy
    real(kind=dp)                     :: input_hice    !< ice thickness
    real(kind=dp)                     :: input_fice    !< ice fraction
    real(kind=dp)                     :: input_tisfc   !< ice surface temperature
    real(kind=dp)                     :: input_snwdph  !< snow depth
    real(kind=dp)                     :: input_snoalb  !< snow albedo
    real(kind=dp)                     :: input_sncovr  !< snow cover
    real(kind=dp)                     :: input_area     !<
    real(kind=dp)                     :: input_tg3      !<
    real(kind=dp)                     :: input_uustar   !<
    real(kind=dp)                     :: input_alvsf !< uv+visible black sky albedo (z=60 degree)
    real(kind=dp)                     :: input_alnsf !< near IR black sky albedo (z=60 degree)
    real(kind=dp)                     :: input_alvwf !< uv+visible white sky albedo
    real(kind=dp)                     :: input_alnwf !< near IR white sky albedo
    real(kind=dp)                     :: input_convexity !< surface orograpy standard deviation
    real(kind=dp)                     :: input_stddev    !< surface orograpy standard deviation
    real(kind=dp)                     :: input_oa1       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_oa2       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_oa3       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_oa4       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_ol1       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_ol2       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_ol3       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_ol4       !< surface orograpy standard deviation
    real(kind=dp)                     :: input_theta     !< surface orograpy standard deviation
    real(kind=dp)                     :: input_gamma     !< surface orograpy standard deviation
    real(kind=dp)                     :: input_sigma     !< surface orograpy standard deviation
    real(kind=dp)                     :: input_elvmax    !< surface orograpy standard deviation
    real(kind=dp)                     :: input_facsf !< fraction of surface depedent on sun angle
    real(kind=dp)                     :: input_facwf !< fraction of surface not depedent on sun angle
    real(kind=dp), allocatable        :: input_pres_i(:) !< pressure (Pa) of input interface
    real(kind=dp), allocatable        :: input_pres_l(:) !< pressure (Pa) of input levels
    real(kind=dp), allocatable              :: input_pres(:) !< pressure (Pa) of input levels
    real(kind=dp), allocatable              :: input_time(:) !< time (s) since beginning of forcing file
    real(kind=dp), allocatable              :: input_height(:) !< height of input levels (m) (initial)
    real(kind=dp), allocatable              :: input_thetail(:) !< ice-liquid water potential temperature(K) (initial)
    real(kind=dp), allocatable              :: input_temp(:) !< temperature(K) (initial)
    real(kind=dp), allocatable              :: input_qt(:) !< total water specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_ql(:) !< suspended liquid specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_qi(:) !< suspended ice specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_u(:) !< zonal wind (m/s) (initial)
    real(kind=dp), allocatable              :: input_v(:) !< meridional wind (m/s) (initial)
    real(kind=dp), allocatable              :: input_tke(:) !< turbulence kinetic energy (m^2/s^2) (initial)
    real(kind=dp), allocatable              :: input_ozone(:) !< ozone mass mixing ratio (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_stc(:) !< soil temperature (k) (initial)
    real(kind=dp), allocatable              :: input_smc(:) !< soil moisture conteng (g/g) (initial)
    real(kind=dp), allocatable              :: input_slc(:) !< soil liquid content (g/g) (initial)
    real(kind=dp), allocatable              :: input_lat(:) !< time-series of latitude of column center
    real(kind=dp), allocatable              :: input_lon(:) !< time-series of longitude of column center
    real(kind=dp), allocatable              :: input_pres_surf(:) !< time-series of surface pressure (Pa)
    real(kind=dp), allocatable              :: input_T_surf(:) !< time-series of surface temperture
    real(kind=dp), allocatable              :: input_w_ls(:,:) !< large-scale vertical velocity (m/s) (time, levels)
    real(kind=dp), allocatable              :: input_omega(:,:) !< large-scale pressure vertical velocity (Pa/s) (time, levels)
    real(kind=dp), allocatable              :: input_u_g(:,:) !< geostrophic zonal wind (m/s) (time, levels)
    real(kind=dp), allocatable              :: input_v_g(:,:) !< geostrophic meridional wind (m/s) (time, levels)
    real(kind=dp), allocatable              :: input_u_nudge(:,:) !< E-W wind profile to nudge towards (m/s) (time, levels)
    real(kind=dp), allocatable              :: input_v_nudge(:,:) !< N-S wind profile to nudge towards (m/s) (time, levels)
    real(kind=dp), allocatable              :: input_T_nudge(:,:) !< absolute temperature profile to nudge towards (K) (time, levels)
    real(kind=dp), allocatable              :: input_thil_nudge(:,:) !< liquid potential temperature profile to nudge towards (K) (time, levels)
    real(kind=dp), allocatable              :: input_qt_nudge(:,:) !< specific humidity profile to nudge towards (kg/kg) (time, levels)
    real(kind=dp), allocatable              :: input_dT_dt_rad(:,:) !< large-scale T-tendency (K/s) (time, levels)
    real(kind=dp), allocatable              :: input_h_advec_thetail(:,:) !< large-scale tendency of ice-liquid potential temperature due to horizontal advection (K/s) (time, levels)
    real(kind=dp), allocatable              :: input_h_advec_qt(:,:) !< large-scale tendency of total water specific humidity due to horizontal advection (kg/kg /s) (time, levels)
    real(kind=dp), allocatable              :: input_v_advec_thetail(:,:) !< large-scale tendency of ice-liquid potential temperature due to vertical advection (K/s) (time, levels)
    real(kind=dp), allocatable              :: input_v_advec_qt(:,:) !< large-scale tendency of total water specific humidity due to vertical advection (kg/kg /s) (time, levels)
    real(kind=dp), allocatable              :: input_sh_flux_sfc(:) !< time-series of surface sensible heat flux (K m s^-1)
    real(kind=dp), allocatable              :: input_lh_flux_sfc(:) !< time-series of surface latent heat flux (kg kg^-1 m s^-1)

    contains
      procedure :: create  => scm_input_create
      procedure :: create_modelics  => scm_input_create_model_ics

  end type scm_input_type

  type scm_reference_type
    !> - Define the reference profile variables.
    integer                                 :: ref_nlev !< number of levels in the reference profile
    real(kind=dp), allocatable              :: ref_pres(:) !< pressure (Pa) of the reference profile levels
    real(kind=dp), allocatable              :: ref_T(:) !< absolute T (K) of the reference profile levels
    real(kind=dp), allocatable              :: ref_qv(:) !< water vapor specific humidity (kg/kg) of the reference profile levels
    real(kind=dp), allocatable              :: ref_ozone(:) !< ozone mass mixing ratio (kg/kg) of the reference profile level

    contains
      procedure :: create => scm_reference_create

  end type scm_reference_type

!> \section arg_table_physics_type
!! \htmlinclude physics_type.html
!!
  type physics_type

    type(GFS_control_type), allocatable      :: Model(:)
    type(GFS_statein_type), allocatable      :: Statein(:)
    type(GFS_stateout_type), allocatable     :: Stateout(:)
    type(GFS_sfcprop_type), allocatable      :: Sfcprop(:)
    type(GFS_coupling_type), allocatable     :: Coupling(:)
    type(GFS_grid_type), allocatable         :: Grid(:)
    type(GFS_tbd_type), allocatable          :: Tbd(:)
    type(GFS_cldprop_type), allocatable      :: Cldprop(:)
    type(GFS_radtend_type), allocatable      :: Radtend(:)
    type(GFS_diag_type), allocatable         :: Diag(:)
    type(GFS_interstitial_type), allocatable :: Interstitial(:)
    type(GFS_init_type), allocatable         :: Init_parm(:)

    contains
      procedure :: create => physics_create
      procedure :: associate => physics_associate
  end type physics_type

  contains

  subroutine scm_state_create(scm_state, n_columns, n_levels, n_nsoil, n_time_levels)
    class(scm_state_type)             :: scm_state
    integer, intent(in)               :: n_columns, n_levels, n_nsoil, n_time_levels

    scm_state%experiment_name = clear_char
    scm_state%model_name = clear_char
    scm_state%output_dir = clear_char
    scm_state%case_data_dir = clear_char
    scm_state%vert_coord_data_dir = clear_char
    scm_state%output_file = clear_char
    scm_state%case_name = clear_char

    allocate(scm_state%physics_suite_name(n_columns), scm_state%physics_nml(n_columns))
    scm_state%physics_suite_name = clear_char
    scm_state%physics_nml = clear_char

    scm_state%n_levels = n_levels
    scm_state%itt = int_zero
    scm_state%itt_out = int_zero
    scm_state%time_scheme = int_zero
    scm_state%n_cols = n_columns
    scm_state%n_timesteps = int_zero
    scm_state%n_time_levels = n_time_levels
    scm_state%n_tracers = 16
    allocate(scm_state%tracer_names(scm_state%n_tracers))
    scm_state%water_vapor_index = 1
    scm_state%ozone_index = 2
    scm_state%cloud_water_index = 3
    scm_state%cloud_ice_index = 4
    scm_state%rain_index = 5
    scm_state%snow_index = 6
    scm_state%graupel_index = 7
    scm_state%cloud_amount_index = 8
    scm_state%cloud_droplet_nc_index = 9
    scm_state%cloud_ice_nc_index = 10
    scm_state%rain_nc_index = 11
    scm_state%snow_nc_index = 12
    scm_state%graupel_nc_index = 13
    scm_state%tke_index = 14
    scm_state%water_friendly_aerosol_index = 15
    scm_state%ice_friendly_aerosol_index = 16
    scm_state%tracer_names(1) = 'vap_wat'
    scm_state%tracer_names(2) = 'o3mr'
    scm_state%tracer_names(3) = 'liq_wat'
    scm_state%tracer_names(4) = 'ice_wat'
    scm_state%tracer_names(5) = 'rainwat'
    scm_state%tracer_names(6) = 'snowwat'
    scm_state%tracer_names(7) = 'graupel'
    scm_state%tracer_names(8) = 'cld_amt'
    scm_state%tracer_names(9) = 'water_nc'
    scm_state%tracer_names(10)= 'ice_nc'
    scm_state%tracer_names(11)= 'rain_nc'
    scm_state%tracer_names(12)= 'snow_nc'
    scm_state%tracer_names(13)= 'graupel_nc'
    scm_state%tracer_names(14)= 'sgs_tke'
    scm_state%tracer_names(15)= 'liq_aero'
    scm_state%tracer_names(16)= 'ice_aero'
    scm_state%n_itt_swrad = int_zero
    scm_state%n_itt_lwrad = int_zero
    scm_state%n_itt_out = int_zero
    scm_state%n_levels_smooth = 5
    allocate(scm_state%blksz(n_columns))
    scm_state%blksz = int_one

    scm_state%sfc_flux_spec = .false.
    scm_state%sfc_type = int_zero
    scm_state%sfc_type_real = real_zero
    scm_state%veg_type = int_zero
    scm_state%veg_type_real = real_zero
    scm_state%veg_frac = int_zero
    scm_state%veg_frac_real = real_zero
    scm_state%shdmin = real_zero
    scm_state%shdmax = real_zero
    scm_state%tg3 = real_zero
    scm_state%uustar = real_zero
    scm_state%soil_type = int_zero
    scm_state%soil_type_real = real_zero
    scm_state%C_RES            = int_zero
    scm_state%mom_forcing_type = int_zero
    scm_state%thermo_forcing_type = int_zero
    scm_state%reference_profile_choice = int_zero

    scm_state%model_time = real_zero
    scm_state%dt = real_zero
    scm_state%dt_now = real_zero
    scm_state%runtime = real_zero
    scm_state%output_frequency = real_zero
    scm_state%relax_time = real_zero
    scm_state%deg_to_rad_const = real_zero
    scm_state%c_filter = 0.15

    scm_state%init_year = int_zero
    scm_state%init_month = int_zero
    scm_state%init_day = int_zero
    scm_state%init_hour = int_zero

    allocate(scm_state%pres_l(n_columns, 1, n_levels), scm_state%pres_i(n_columns, 1, n_levels+1), &
      scm_state%exner_l(n_columns, 1, n_levels), scm_state%exner_i(n_columns, 1, n_levels+1), &
      scm_state%geopotential_l(n_columns, 1, n_levels), scm_state%geopotential_i(n_columns, 1, n_levels+1))
    scm_state%pres_l = real_zero
    scm_state%pres_i = real_zero
    scm_state%exner_l = real_zero
    scm_state%exner_i = real_zero
    scm_state%geopotential_l = real_zero
    scm_state%geopotential_i = real_zero

    allocate(scm_state%a_k(1, n_levels+1), scm_state%b_k(1, n_levels+1), scm_state%si(n_columns, 1, n_levels+1), &
      scm_state%sl(n_columns, 1, n_levels))
    scm_state%a_k = real_zero
    scm_state%b_k = real_zero
    scm_state%si = real_zero
    scm_state%sl = real_zero

    allocate(scm_state%lat(n_columns,1), scm_state%lon(n_columns,1), scm_state%area(n_columns,1))
    scm_state%lat = real_zero
    scm_state%lon = real_zero
    scm_state%area = real_one

    allocate(scm_state%state_T(n_columns, 1, n_levels, n_time_levels), &
      scm_state%state_u(n_columns, 1, n_levels, n_time_levels), scm_state%state_v(n_columns, 1, n_levels, n_time_levels), &
      scm_state%state_tracer(n_columns, 1, n_levels, scm_state%n_tracers, n_time_levels))
    scm_state%state_T = real_zero
    scm_state%state_u = real_zero
    scm_state%state_v = real_zero
    scm_state%state_tracer = real_zero

    allocate(scm_state%temp_tracer(n_columns, 1, n_levels, scm_state%n_tracers, n_time_levels), &
      scm_state%temp_T(n_columns, 1, n_levels, n_time_levels), &
      scm_state%temp_u(n_columns, 1, n_levels, n_time_levels), scm_state%temp_v(n_columns, 1, n_levels, n_time_levels))
    allocate(scm_state%stc(n_columns, 1, n_nsoil, n_time_levels))
    allocate(scm_state%smc(n_columns, 1, n_nsoil, n_time_levels))
    allocate(scm_state%slc(n_columns, 1, n_nsoil, n_time_levels))
    scm_state%temp_tracer = real_zero
    scm_state%temp_T = real_zero
    scm_state%temp_u = real_zero
    scm_state%temp_v = real_zero

    allocate(scm_state%w_ls(n_columns, n_levels), scm_state%omega(n_columns, 1, n_levels), scm_state%u_g(n_columns, n_levels), &
      scm_state%v_g(n_columns, n_levels), scm_state%dT_dt_rad(n_columns, n_levels), scm_state%h_advec_thil(n_columns, n_levels), &
      scm_state%h_advec_qt(n_columns, n_levels), scm_state%v_advec_thil(n_columns, n_levels), &
      scm_state%v_advec_qt(n_columns, n_levels), scm_state%pres_surf(n_columns,1), scm_state%T_surf(n_columns,1), &
      scm_state%alvsf(n_columns), scm_state%alnsf(n_columns), scm_state%alvwf(n_columns), scm_state%alnwf(n_columns), &
      scm_state%facsf(n_columns), scm_state%facwf(n_columns), &
      scm_state%hprime(n_columns,14), scm_state%stddev(n_columns), &
      scm_state%u_nudge(n_columns, n_levels), scm_state%v_nudge(n_columns, n_levels), &
      scm_state%T_nudge(n_columns, n_levels), scm_state%thil_nudge(n_columns, n_levels), &
      scm_state%qt_nudge(n_columns, n_levels), scm_state%sh_flux(n_columns), scm_state%lh_flux(n_columns), &
      scm_state%u_force_tend(n_columns,n_levels), scm_state%v_force_tend(n_columns,n_levels), &
      scm_state%T_force_tend(n_columns,n_levels), scm_state%qv_force_tend(n_columns,n_levels), &
      scm_state%sfc_roughness_length_cm(n_columns))
    scm_state%w_ls = real_zero
    scm_state%omega = real_zero
    scm_state%u_g = real_zero
    scm_state%v_g = real_zero
    scm_state%dT_dt_rad = real_zero
    scm_state%h_advec_thil = real_zero
    scm_state%h_advec_qt = real_zero
    scm_state%v_advec_thil = real_zero
    scm_state%v_advec_qt = real_zero
    scm_state%pres_surf = real_zero
    scm_state%T_surf = real_zero
    scm_state%alvsf = real_zero
    scm_state%alnsf = real_zero
    scm_state%alvwf = real_zero
    scm_state%alnwf = real_zero
    scm_state%facsf = real_zero
    scm_state%facwf = real_zero
    scm_state%hprime = real_zero
    scm_state%stddev = real_zero
    scm_state%u_nudge = real_zero
    scm_state%v_nudge = real_zero
    scm_state%T_nudge = real_zero
    scm_state%thil_nudge = real_zero
    scm_state%qt_nudge = real_zero
    scm_state%sh_flux = real_zero
    scm_state%lh_flux = real_zero
    scm_state%u_force_tend = real_zero
    scm_state%v_force_tend = real_zero
    scm_state%T_force_tend = real_zero
    scm_state%qv_force_tend = real_zero
    scm_state%sfc_roughness_length_cm = real_one

  end subroutine scm_state_create

  subroutine scm_input_create_model_ics(scm_input, nlev_soil,nlev)
    class(scm_input_type)             :: scm_input
    integer, intent(in)               :: nlev,nlev_soil

    scm_input%input_nsoil= nlev_soil
    allocate(scm_input%input_stc(nlev_soil), scm_input%input_smc(nlev_soil), scm_input%input_slc(nlev_soil))
    allocate(scm_input%input_temp(nlev), scm_input%input_pres_i(nlev+1),scm_input%input_pres_l(nlev))
    scm_input%input_stc    = real_zero
    scm_input%input_smc    = real_zero
    scm_input%input_slc    = real_zero
    scm_input%input_temp   = real_zero
    scm_input%input_pres_i = real_zero
    scm_input%input_pres_l = real_zero

  end subroutine scm_input_create_model_ics

  subroutine scm_input_create(scm_input, ntimes, nlev)
    class(scm_input_type)             :: scm_input
    integer, intent(in)               :: ntimes, nlev

    scm_input%input_nlev = nlev
    scm_input%input_ntimes = ntimes

    allocate(scm_input%input_pres(nlev),scm_input%input_time(ntimes))
    scm_input%input_pres = real_zero
    scm_input%input_time = real_zero

    allocate(scm_input%input_height(nlev), scm_input%input_thetail(nlev), scm_input%input_qt(nlev), scm_input%input_ql(nlev), &
      scm_input%input_qi(nlev), scm_input%input_u(nlev), scm_input%input_v(nlev), scm_input%input_tke(nlev), &
      scm_input%input_ozone(nlev))
    scm_input%input_height = real_zero
    scm_input%input_thetail = real_zero
    scm_input%input_qt = real_zero
    scm_input%input_ql = real_zero
    scm_input%input_qi = real_zero
    scm_input%input_u = real_zero
    scm_input%input_v = real_zero
    scm_input%input_tke = real_zero
    scm_input%input_ozone = real_zero

    allocate(scm_input%input_lat(ntimes), scm_input%input_lon(ntimes), scm_input%input_pres_surf(ntimes), &
      scm_input%input_T_surf(ntimes), scm_input%input_sh_flux_sfc(ntimes), scm_input%input_lh_flux_sfc(ntimes), &
      scm_input%input_w_ls(ntimes, nlev), scm_input%input_omega(ntimes, nlev), scm_input%input_u_g(ntimes, nlev), &
      scm_input%input_v_g(ntimes, nlev), scm_input%input_dT_dt_rad(ntimes, nlev), scm_input%input_h_advec_thetail(ntimes, nlev), &
      scm_input%input_h_advec_qt(ntimes, nlev), scm_input%input_v_advec_thetail(ntimes, nlev), &
      scm_input%input_v_advec_qt(ntimes, nlev), scm_input%input_u_nudge(ntimes, nlev), scm_input%input_v_nudge(ntimes, nlev),    &
      scm_input%input_T_nudge(ntimes, nlev), scm_input%input_thil_nudge(ntimes, nlev), scm_input%input_qt_nudge(ntimes, nlev))
    scm_input%input_lat = real_zero
    scm_input%input_lon = real_zero
    scm_input%input_pres_surf = real_zero
    scm_input%input_T_surf = real_zero
    scm_input%input_pres_surf = real_zero
    scm_input%input_T_surf = real_zero
    scm_input%input_alvsf        = real_zero
    scm_input%input_alnsf        = real_zero
    scm_input%input_alvwf        = real_zero
    scm_input%input_alnwf        = real_zero
    scm_input%input_stddev       = real_zero
    scm_input%input_convexity    = real_zero
    scm_input%input_oa1          = real_zero
    scm_input%input_oa2          = real_zero
    scm_input%input_oa3          = real_zero
    scm_input%input_oa4          = real_zero
    scm_input%input_ol1          = real_zero
    scm_input%input_ol2          = real_zero
    scm_input%input_ol3          = real_zero
    scm_input%input_ol4          = real_zero
    scm_input%input_theta        = real_zero
    scm_input%input_gamma        = real_zero
    scm_input%input_sigma        = real_zero
    scm_input%input_elvmax       = real_zero
    scm_input%input_facsf        = real_zero
    scm_input%input_sh_flux_sfc = real_zero
    scm_input%input_lh_flux_sfc = real_zero
    scm_input%input_w_ls = real_zero
    scm_input%input_omega = real_zero
    scm_input%input_u_g = real_zero
    scm_input%input_v_g = real_zero
    scm_input%input_dT_dt_rad = real_zero
    scm_input%input_h_advec_thetail = real_zero
    scm_input%input_h_advec_qt = real_zero
    scm_input%input_v_advec_thetail = real_zero
    scm_input%input_v_advec_qt = real_zero
    scm_input%input_sh_flux_sfc = real_zero
    scm_input%input_lh_flux_sfc = real_zero
    scm_input%input_w_ls = real_zero
    scm_input%input_omega = real_zero
    scm_input%input_u_g = real_zero
    scm_input%input_v_g = real_zero
    scm_input%input_dT_dt_rad = real_zero
    scm_input%input_h_advec_thetail = real_zero
    scm_input%input_h_advec_qt = real_zero
    scm_input%input_v_advec_thetail = real_zero
    scm_input%input_v_advec_qt = real_zero
    scm_input%input_u_nudge = real_zero
    scm_input%input_v_nudge = real_zero
    scm_input%input_T_nudge = real_zero
    scm_input%input_thil_nudge = real_zero
    scm_input%input_qt_nudge = real_zero

  end subroutine scm_input_create

  subroutine scm_reference_create(scm_reference, nlev)
    class(scm_reference_type)             :: scm_reference
    integer, intent(in)               :: nlev

    scm_reference%ref_nlev = nlev

    allocate(scm_reference%ref_pres(nlev), scm_reference%ref_T(nlev), scm_reference%ref_qv(nlev), scm_reference%ref_ozone(nlev))
    scm_reference%ref_pres = real_zero
    scm_reference%ref_T = real_zero
    scm_reference%ref_qv = real_zero
    scm_reference%ref_ozone = real_zero

  end subroutine scm_reference_create

  subroutine physics_create(physics, n_columns)
    class(physics_type) :: physics
    integer, intent(in) :: n_columns
    
    real(kind=kind_phys) :: kind_phys_zero

    integer :: i
    integer, dimension(8) :: zeroes_8

    zeroes_8(:) = int_zero
    kind_phys_zero = real_zero

    allocate(physics%Model(n_columns), physics%Statein(n_columns), physics%Stateout(n_columns), physics%Sfcprop(n_columns), &
      physics%Coupling(n_columns), physics%Grid(n_columns), physics%Tbd(n_columns), physics%Cldprop(n_columns), &
      physics%Radtend(n_columns), physics%Diag(n_columns), physics%Interstitial(n_columns), &
      physics%Init_parm(n_columns))

    do i=1, n_columns
      physics%Init_parm(i)%me = int_zero
      physics%Init_parm(i)%master = int_zero
      physics%Init_parm(i)%isc = int_one
      physics%Init_parm(i)%jsc = int_one
      physics%Init_parm(i)%nx = int_one
      physics%Init_parm(i)%ny = int_one
      physics%Init_parm(i)%levs = int_one
      physics%Init_parm(i)%cnx = int_one
      physics%Init_parm(i)%cny = int_one
      physics%Init_parm(i)%gnx = int_one
      physics%Init_parm(i)%gny = int_one
      physics%Init_parm(i)%nlunit = int_one
      physics%Init_parm(i)%logunit= 10 + i
      physics%Init_parm(i)%bdat(:) = zeroes_8(:)
      physics%Init_parm(i)%cdat(:) = zeroes_8(:)
      physics%Init_parm(i)%dt_dycore = kind_phys_zero
      physics%Init_parm(i)%dt_phys = kind_phys_zero
      physics%Init_parm(i)%iau_offset = int_zero
      physics%Init_parm(i)%ak => null()
      physics%Init_parm(i)%bk => null()
      physics%Init_parm(i)%xlon => null()
      physics%Init_parm(i)%xlat => null()
      physics%Init_parm(i)%area => null()
      physics%Init_parm(i)%tracer_names => null()
      physics%Init_parm(i)%blksz => null()
      physics%Init_parm(i)%restart = .false.
      physics%Init_parm(i)%hydrostatic = .true.
      physics%Init_parm(i)%tile_num = int_one
    end do

  end subroutine physics_create

  subroutine physics_associate(physics, scm_state, col)
    class(physics_type) :: physics
    type(scm_state_type), target, intent(in) :: scm_state
    integer, intent(in) :: col

    physics%Statein(col)%phii => scm_state%geopotential_i(col,:,:)
    physics%Statein(col)%prsi => scm_state%pres_i(col,:,:)
    physics%Statein(col)%prsik => scm_state%exner_i(col,:,:)
    physics%Statein(col)%phil => scm_state%geopotential_l(col,:,:)
    physics%Statein(col)%prsl => scm_state%pres_l(col,:,:)
    physics%Statein(col)%prslk => scm_state%exner_l(col,:,:)

    physics%Statein(col)%pgr => scm_state%pres_surf(col,:)
    physics%Statein(col)%ugrs => scm_state%state_u(col,:,:,1)
    physics%Statein(col)%vgrs => scm_state%state_v(col,:,:,1)
    physics%Statein(col)%vvl => scm_state%omega(col,:,:)
    physics%Statein(col)%tgrs => scm_state%state_T(col,:,:,1)
    physics%Statein(col)%qgrs => scm_state%state_tracer(col,:,:,:,1)

    physics%Sfcprop(col)%tsfc => scm_state%T_surf(col,:)
    physics%Sfcprop(col)%tref => scm_state%T_surf(col,:)
    physics%Sfcprop(col)%slmsk => scm_state%sfc_type_real
    physics%Sfcprop(col)%vtype => scm_state%veg_type_real
    physics%Sfcprop(col)%vfrac => scm_state%veg_frac_real
    physics%Sfcprop(col)%slope => scm_state%slopetype
    physics%Interstitial(col)%sigmaf = min(scm_state%veg_frac_real,0.01)
    physics%Sfcprop(col)%shdmax => scm_state%shdmax        
    physics%Sfcprop(col)%shdmin => scm_state%shdmin        
    physics%Sfcprop(col)%tg3 => scm_state%tg3        
    physics%Sfcprop(col)%uustar => scm_state%uustar        
    physics%Sfcprop(col)%stype => scm_state%soil_type_real
    physics%Sfcprop(col)%alvsf => scm_state%alvsf
    physics%Sfcprop(col)%alnsf => scm_state%alnsf
    physics%Sfcprop(col)%hprim => scm_state%stddev
    physics%Sfcprop(col)%hprime=> scm_state%hprime 
    physics%Sfcprop(col)%alvwf => scm_state%alvwf
    physics%Sfcprop(col)%alnwf => scm_state%alnwf
    physics%Sfcprop(col)%facsf => scm_state%facsf
    physics%Sfcprop(col)%facwf => scm_state%facwf

    physics%Sfcprop(col)%stc   => scm_state%stc(col,:,:,1)
    physics%Sfcprop(col)%smc   => scm_state%smc(col,:,:,1)
    physics%Sfcprop(col)%slc   => scm_state%slc(col,:,:,1)
    
    !GJF : the following logic was introduced into FV3GFS_io.F90 as part of the fractional landmask update (additional logic exists in the same file if the fractional landmask is actually used!)
    if (physics%Sfcprop(col)%slmsk(1) > 1.9) then
      physics%Sfcprop(col)%landfrac(1)=0.
    else
      physics%Sfcprop(col)%landfrac(1) = physics%Sfcprop(col)%slmsk(1) 
    end if
    if (physics%Sfcprop(col)%lakefrac(1) > 0.) then
      physics%Sfcprop(col)%oceanfrac(1) = 0. ! lake & ocean don't coexist in a cell, lake dominates
    else
      physics%Sfcprop(col)%oceanfrac(1) = 1.-physics%Sfcprop(col)%landfrac(1) !LHS:ocean frac [0:1]
    end if
    !GJF

    if(scm_state%time_scheme == 2) then
      physics%Stateout(col)%gu0 => scm_state%state_u(col,:,:,2)
      physics%Stateout(col)%gv0 => scm_state%state_v(col,:,:,2)
      physics%Stateout(col)%gt0 => scm_state%state_T(col,:,:,2)
      physics%Stateout(col)%gq0 => scm_state%state_tracer(col,:,:,:,2)
    else
      physics%Stateout(col)%gu0 => scm_state%state_u(col,:,:,1)
      physics%Stateout(col)%gv0 => scm_state%state_v(col,:,:,1)
      physics%Stateout(col)%gt0 => scm_state%state_T(col,:,:,1)
      physics%Stateout(col)%gq0 => scm_state%state_tracer(col,:,:,:,1)
    endif

    if(scm_state%sfc_flux_spec) then
      physics%Sfcprop(col)%spec_sh_flux => scm_state%sh_flux
      physics%Sfcprop(col)%spec_lh_flux => scm_state%lh_flux
      physics%Sfcprop(col)%zorl => scm_state%sfc_roughness_length_cm
    endif
    if(scm_state%model_ics) then
      physics%Sfcprop(col)%zorl => scm_state%sfc_roughness_length_cm
      physics%Sfcprop(col)%canopy => scm_state%canopy        
      physics%Sfcprop(col)%hice => scm_state%hice        
      physics%Sfcprop(col)%fice => scm_state%fice        
      physics%Sfcprop(col)%tisfc => scm_state%tisfc        
      physics%Sfcprop(col)%snowd  => scm_state%snwdph        
      physics%Sfcprop(col)%snoalb => scm_state%snoalb        
      physics%Sfcprop(col)%sncovr => scm_state%sncovr        
    endif


  end subroutine physics_associate

end module gmtb_scm_type_defs
