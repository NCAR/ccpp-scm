!> \file gmtb_scm_type_defs.f90
!!  Contains type definitions for SCM-related variables and physics-related variables

module gmtb_scm_type_defs

!> \section arg_table_gmtb_scm_type_defs
!! \htmlinclude gmtb_scm_type_defs.html
!!

  use gmtb_scm_kinds, only: sp, dp, qp
  use GFS_typedefs,   only: GFS_control_type,      &
                            GFS_statein_type,      &
                            GFS_stateout_type,     &
                            GFS_sfcprop_type,      &
                            GFS_coupling_type,     &
                            GFS_grid_type,         &
                            GFS_tbd_type,          &
                            GFS_cldprop_type,      &
                            GFS_radtend_type,      &
                            GFS_diag_type,         &
                            GFS_interstitial_type, &
                            GFS_init_type
  use machine,        only: kind_phys
  use ccpp_api,       only: ccpp_t

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
    character(len=character_length)                 :: physics_suite_name !< name of physics suite (must be "GFS_operational" for prototype)
    character(len=character_length)                 :: physics_nml

    integer                           :: n_levels !< number of model levels (must be 64 for prototype)
    integer                           :: n_soil  !< number of model levels (must be 4 for prototype)
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
    integer                           :: mass_weighted_rime_factor_index !< index for mass-weighted rime factor
    integer                           :: init_year, init_month, init_day, init_hour
    character(len=32), allocatable    :: tracer_names(:) !< name of physics suite (must be "GFS_operational" for prototype)
    integer, allocatable              :: blksz(:)

    logical                           :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
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
    real(kind=dp), allocatable              :: pres_i(:,:), pres_l(:,:) !< pressure on grid interfaces, centers (Pa)
    real(kind=dp), allocatable              :: si(:,:), sl(:,:) !< sigma on grid interfaces, centers
    real(kind=dp), allocatable              :: exner_i(:,:), exner_l(:,:) !< exner function on grid interfaces, centers
    real(kind=dp), allocatable              :: geopotential_i(:,:), geopotential_l(:,:) !< geopotential on grid interfaces, centers
    real(kind=dp), allocatable              :: a_k(:), b_k(:) !< used to determine grid sigma and pressure levels

    real(kind=dp), allocatable              :: lat(:), lon(:) !< latitude and longitude (radians)
    real(kind=dp), allocatable              :: area(:) !< area over which the column represents a mean (analogous to grid size or observational array area)

    real(kind=dp), allocatable              :: state_T(:,:,:) !< model state absolute temperature at grid centers (K)
    real(kind=dp), allocatable              :: state_u(:,:,:), state_v(:,:,:) !< model state horizontal winds at grid centers (m/s)
    real(kind=dp), allocatable              :: state_tracer(:,:,:,:) !< model state tracer at grid centers
    real(kind=dp), allocatable              :: temp_T(:,:,:), temp_u(:,:,:), temp_v(:,:,:), temp_tracer(:,:,:,:) !< used for time-filtering

    !> - Define forcing-related variables (indexing is (horizontal, vertical)).
    real(kind=dp), allocatable              :: u_force_tend(:,:), v_force_tend(:,:), T_force_tend(:,:), qv_force_tend(:,:) !< total u, v, T, q forcing (units/s) (horizontal, vertical)
    real(kind=dp), allocatable              :: w_ls(:,:), omega(:,:), u_g(:,:), v_g(:,:), dT_dt_rad(:,:), h_advec_thil(:,:), &
      h_advec_qt(:,:), v_advec_thil(:,:), v_advec_qt(:,:), u_nudge(:,:), v_nudge(:,:), T_nudge(:,:), thil_nudge(:,:), qt_nudge(:,:) !< forcing terms interpolated to the model time and grid
    real(kind=dp), allocatable              :: T_surf(:), pres_surf(:) !< surface temperature and pressure interpolated to the model time
    real(kind=dp), allocatable              :: sh_flux(:), lh_flux(:) !< surface sensible and latent heat fluxes interpolated to the model time
    real(kind=dp), allocatable              :: sfc_roughness_length_cm(:) !< surface roughness length used for calculating surface layer parameters from specified fluxes
    
    real(kind=dp), allocatable              :: sfc_type(:) !< 0: sea surface, 1: land surface, 2: sea-ice surface
        
    contains
      procedure :: create  => scm_state_create

  end type scm_state_type

  type scm_input_type
    !> - Define the case-specific initialization and forcing variables.
    integer                           :: input_nlev !< number of levels in the input file
    integer                           :: input_nsoil !< number of soil levels in the input file
    integer                           :: input_nsnow !< number of snow layers in the input file
    integer                           :: input_nice !< number of sea ice layers in the input file
    integer                           :: input_ntimes !< number of times in the input file where forcing is available
    integer                           :: input_vegsrc !< vegetation source
    real(kind=dp)                     :: input_vegtyp !< vegetation type classification
    real(kind=dp)                     :: input_soiltyp !<
    real(kind=dp)                     :: input_slopetype !< surface slope classification
    real(kind=dp)                     :: input_lat !< latitude of column center
    real(kind=dp)                     :: input_lon !< longitude of column center
    real(kind=dp)                     :: input_tsfco !< input sea surface temperature OR surface skin temperature over land OR surface skin temperature over ice (depending on slmsk) (K)
    real(kind=dp)                     :: input_vegfrac  !< vegetation area fraction
    real(kind=dp)                     :: input_shdmin   !< minimun vegetation fraction
    real(kind=dp)                     :: input_shdmax   !< maximun vegetation fraction
    real(kind=dp)                     :: input_zorlo    !< surfce roughness length over ocean [cm]
    real(kind=dp)                     :: input_slmsk   !< sea land ice mask [0,1,2]
    real(kind=dp)                     :: input_canopy  !< amount of water stored in canopy (kg m-2)
    real(kind=dp)                     :: input_hice    !< sea ice thickness (m)
    real(kind=dp)                     :: input_fice    !< ice fraction (frac)
    real(kind=dp)                     :: input_tisfc   !< ice surface temperature (K)
    real(kind=dp)                     :: input_snwdph  !< water equivalent snow depth (mm)
    real(kind=dp)                     :: input_snoalb  !< maximum snow albedo (frac)
    real(kind=dp)                     :: input_sncovr  !< snow area fraction (frac)
    real(kind=dp)                     :: input_area    !< surface area [m^2]
    real(kind=dp)                     :: input_tg3     !< deep soil temperature (K)
    real(kind=dp)                     :: input_uustar  !< surface friction velocity (m s-1)
    real(kind=dp)                     :: input_alvsf   !< 60 degree vis albedo with strong cosz dependency
    real(kind=dp)                     :: input_alnsf     !< 60 degree nir albedo with strong cosz dependency
    real(kind=dp)                     :: input_alvwf     !< 60 degree vis albedo with weak cosz dependency
    real(kind=dp)                     :: input_alnwf     !< 60 degree nir albedo with weak cosz dependency
    real(kind=dp)                     :: input_facsf     !< fractional coverage with strong cosz dependency
    real(kind=dp)                     :: input_facwf     !< fractional coverage with weak cosz dependency
    real(kind=dp)                     :: input_weasd     !< water equivalent accumulated snow depth (mm)
    real(kind=dp)                     :: input_f10m      !< ratio of sigma level 1 wind and 10m wind
    real(kind=dp)                     :: input_t2m     !< 2-meter absolute temperature (K)
    real(kind=dp)                     :: input_q2m     !< 2-meter specific humidity (kg kg-1)
    real(kind=dp)                     :: input_ffmm    !< Monin-Obukhov similarity function for momentum
    real(kind=dp)                     :: input_ffhh    !< Monin-Obukhov similarity function for heat
    real(kind=dp)                     :: input_tprcp   !< instantaneous total precipitation amount (m)
    real(kind=dp)                     :: input_srflag  !< snow/rain flag for precipitation
    real(kind=dp)                     :: input_tsfcl   !< surface skin temperature over land (K)
    real(kind=dp)                     :: input_zorll   !< surface roughness length over land (cm)
    real(kind=dp)                     :: input_zorli   !< surface roughness length over ice (cm)
    real(kind=dp)                     :: input_zorlw   !< surface roughness length from wave model (cm)
    
    real(kind=dp)                     :: input_stddev !< standard deviation of subgrid orography (m)
    real(kind=dp)                     :: input_convexity !< convexity of subgrid orography 
    real(kind=dp)                     :: input_ol1 !< fraction of grid box with subgrid orography higher than critical height 1
    real(kind=dp)                     :: input_ol2 !< fraction of grid box with subgrid orography higher than critical height 2
    real(kind=dp)                     :: input_ol3 !< fraction of grid box with subgrid orography higher than critical height 3
    real(kind=dp)                     :: input_ol4 !< fraction of grid box with subgrid orography higher than critical height 4
    real(kind=dp)                     :: input_oa1 !< assymetry of subgrid orography 1
    real(kind=dp)                     :: input_oa2 !< assymetry of subgrid orography 2
    real(kind=dp)                     :: input_oa3 !< assymetry of subgrid orography 3
    real(kind=dp)                     :: input_oa4 !< assymetry of subgrid orography 4
    real(kind=dp)                     :: input_sigma !< slope of subgrid orography
    real(kind=dp)                     :: input_theta !< angle with respect to east of maximum subgrid orographic variations (deg)
    real(kind=dp)                     :: input_gamma !< anisotropy of subgrid orography
    real(kind=dp)                     :: input_elvmax!< maximum of subgrid orography (m)
    real(kind=dp)                     :: input_oro !< orography (m)
    real(kind=dp)                     :: input_oro_uf !< unfiltered orography (m)
    real(kind=dp)                     :: input_landfrac !< fraction of horizontal grid area occupied by land
    real(kind=dp)                     :: input_lakefrac !< fraction of horizontal grid area occupied by lake
    real(kind=dp)                     :: input_lakedepth !< lake depth (m)
    
    real(kind=dp)                     :: input_tvxy !< vegetation temperature (K)
    real(kind=dp)                     :: input_tgxy !< ground temperature for Noahmp (K)
    real(kind=dp)                     :: input_tahxy !< canopy air temperature (K)
    real(kind=dp)                     :: input_canicexy !< canopy intercepted ice mass (mm)
    real(kind=dp)                     :: input_canliqxy !< canopy intercepted liquid water (mm)
    real(kind=dp)                     :: input_eahxy !< canopy air vapor pressure (Pa)
    real(kind=dp)                     :: input_cmxy !< surface drag coefficient for momentum for noahmp
    real(kind=dp)                     :: input_chxy !< surface exchange coeff heat & moisture for noahmp
    real(kind=dp)                     :: input_fwetxy !< area fraction of canopy that is wetted/snowed
    real(kind=dp)                     :: input_sneqvoxy !< snow mass at previous time step (mm)
    real(kind=dp)                     :: input_alboldxy !< snow albedo at previous time step (frac)
    real(kind=dp)                     :: input_qsnowxy !< snow precipitation rate at surface (mm s-1)
    real(kind=dp)                     :: input_wslakexy !< lake water storage (mm)
    real(kind=dp)                     :: input_taussxy !< non-dimensional snow age
    real(kind=dp)                     :: input_waxy !< water storage in aquifer (mm)
    real(kind=dp)                     :: input_wtxy !< water storage in aquifer and saturated soil (mm)
    real(kind=dp)                     :: input_zwtxy !< water table depth (m)
    real(kind=dp)                     :: input_xlaixy !< leaf area index
    real(kind=dp)                     :: input_xsaixy !< stem area index
    real(kind=dp)                     :: input_lfmassxy !< leaf mass (g m-2)
    real(kind=dp)                     :: input_stmassxy !< stem mass (g m-2)
    real(kind=dp)                     :: input_rtmassxy !< fine root mass (g m-2)
    real(kind=dp)                     :: input_woodxy !< wood mass including woody roots (g m-2)
    real(kind=dp)                     :: input_stblcpxy !< stable carbon in deep soil (g m-2)
    real(kind=dp)                     :: input_fastcpxy !< short-lived carbon in shallow soil (g m-2)
    real(kind=dp)                     :: input_smcwtdxy !< soil water content between the bottom of the soil and the water table (m3 m-3)
    real(kind=dp)                     :: input_deeprechxy !< recharge to or from the water table when deep (m)
    real(kind=dp)                     :: input_rechxy !< recharge to or from the water table when shallow (m)
    real(kind=dp)                     :: input_snowxy !< number of snow layers
    
    real(kind=dp)                     :: input_tref !< sea surface reference temperature for NSST (K)
    real(kind=dp)                     :: input_z_c !< sub-layer cooling thickness for NSST (m)
    real(kind=dp)                     :: input_c_0 !< coefficient 1 to calculate d(Tz)/d(Ts) for NSST
    real(kind=dp)                     :: input_c_d !< coefficient 2 to calculate d(Tz)/d(Ts) for NSST
    real(kind=dp)                     :: input_w_0 !< coefficient 3 to calculate d(Tz)/d(Ts) for NSST
    real(kind=dp)                     :: input_w_d !< coefficient 4 to calculate d(Tz)/d(Ts) for NSST
    real(kind=dp)                     :: input_xt !< heat content in diurnal thermocline layer for NSST (K m)
    real(kind=dp)                     :: input_xs !< salinity content in diurnal thermocline layer for NSST (ppt m)
    real(kind=dp)                     :: input_xu !< u-current in diurnal thermocline layer for NSST (m2 s-1)
    real(kind=dp)                     :: input_xv !< v-current in diurnal thermocline layer for NSST (m2 s-1)
    real(kind=dp)                     :: input_xz !< thickness of diurnal thermocline layer for NSST (m)
    real(kind=dp)                     :: input_zm !< thickness of ocean mixed layer for NSST (m)
    real(kind=dp)                     :: input_xtts !< sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST (m)
    real(kind=dp)                     :: input_xzts !< sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST (m K-1)
    real(kind=dp)                     :: input_d_conv !< thickness of free convection layer for NSST (m)
    real(kind=dp)                     :: input_ifd !< index to start DTM run for NSST
    real(kind=dp)                     :: input_dt_cool !< sub-layer cooling amount for NSST (K)
    real(kind=dp)                     :: input_qrain !< sensible heat due to rainfall for NSST (W)
    
    real(kind=dp)                     :: input_wetness !< normalized soil wetness for RUC LSM
    real(kind=dp)                     :: input_clw_surf !< cloud condensed water mixing ratio at surface for RUC LSM (kg kg-1)
    real(kind=dp)                     :: input_qwv_surf !< water vapor mixing ratio at surface for RUC LSM (kg kg-1)
    real(kind=dp)                     :: input_tsnow !< snow temperature at the bottom of the first snow layer for RUC LSM (K)
    real(kind=dp)                     :: input_snowfallac !< run-total snow accumulation on the ground for RUC LSM (kg m-2)
    real(kind=dp)                     :: input_acsnow !< snow water equivalent of run-total frozen precip for RUC LSM (kg m-2)
    real(kind=dp)                     :: input_lai !< leaf area index for RUC LSM
    
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
    real(kind=dp), allocatable              :: input_snicexy(:) !< snow layer ice (mm)
    real(kind=dp), allocatable              :: input_snliqxy(:) !< snow layer liquid (mm)
    real(kind=dp), allocatable              :: input_tsnoxy(:) !< snow temperature (K)
    real(kind=dp), allocatable              :: input_smoiseq(:) !< equilibrium soil water content (m3 m-3)
    real(kind=dp), allocatable              :: input_zsnsoxy(:) !< layer bottom depth from snow surface (m)
    real(kind=dp), allocatable              :: input_tiice(:) !< sea ice internal temperature (K)
    real(kind=dp), allocatable              :: input_tslb(:) !< soil temperature for RUC LSM (K)    
    real(kind=dp), allocatable              :: input_smois(:) !< volume fraction of soil moisture for RUC LSM (frac)
    real(kind=dp), allocatable              :: input_sh2o(:) !< volume fraction of unfrozen soil moisture for RUC LSM (frac)
    real(kind=dp), allocatable              :: input_smfr(:) !< volume fraction of frozen soil moisture for RUC LSM (frac)
    real(kind=dp), allocatable              :: input_flfr(:) !< flag for frozen soil physics
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

    type(GFS_control_type)                   :: Model
    type(GFS_statein_type)                   :: Statein
    type(GFS_stateout_type)                  :: Stateout
    type(GFS_sfcprop_type)                   :: Sfcprop
    type(GFS_coupling_type)                  :: Coupling
    type(GFS_grid_type)                      :: Grid
    type(GFS_tbd_type)                       :: Tbd
    type(GFS_cldprop_type)                   :: Cldprop
    type(GFS_radtend_type)                   :: Radtend
    type(GFS_diag_type)                      :: Diag
    type(GFS_interstitial_type)              :: Interstitial
    type(GFS_init_type)                      :: Init_parm

    contains
      procedure :: create => physics_create
      procedure :: associate => physics_associate
      procedure :: set => physics_set
  end type physics_type

  type(physics_type), target :: physics

  type(ccpp_t),       target :: cdata

  contains

  subroutine scm_state_create(scm_state, n_columns, n_levels, n_soil, n_snow, n_time_levels, tracers)
    class(scm_state_type)             :: scm_state
    integer, intent(in)               :: n_columns, n_levels, n_soil, n_snow, n_time_levels
    character(len=character_length), intent(in), dimension(:) :: tracers
    
    integer :: i

    scm_state%experiment_name = clear_char
    scm_state%model_name = clear_char
    scm_state%output_dir = clear_char
    scm_state%case_data_dir = clear_char
    scm_state%vert_coord_data_dir = clear_char
    scm_state%output_file = clear_char
    scm_state%case_name = clear_char

    scm_state%physics_suite_name = clear_char
    scm_state%physics_nml = clear_char

    scm_state%n_levels = n_levels
    scm_state%itt = int_zero
    scm_state%itt_out = int_zero
    scm_state%time_scheme = int_zero
    scm_state%n_cols = n_columns
    scm_state%n_timesteps = int_zero
    scm_state%n_time_levels = n_time_levels
    
    scm_state%n_tracers = size(tracers)
    allocate(scm_state%tracer_names(scm_state%n_tracers))
    scm_state%tracer_names = tracers
    scm_state%water_vapor_index               = get_tracer_index(scm_state%tracer_names,"sphum")
    scm_state%ozone_index                     = get_tracer_index(scm_state%tracer_names,"o3mr")
    scm_state%cloud_water_index               = get_tracer_index(scm_state%tracer_names,"liq_wat")
    scm_state%cloud_ice_index                 = get_tracer_index(scm_state%tracer_names,"ice_wat")
    scm_state%rain_index                      = get_tracer_index(scm_state%tracer_names,"rainwat")
    scm_state%snow_index                      = get_tracer_index(scm_state%tracer_names,"snowwat")
    scm_state%graupel_index                   = get_tracer_index(scm_state%tracer_names,"graupel")
    scm_state%cloud_amount_index              = get_tracer_index(scm_state%tracer_names,"cld_amt")
    scm_state%cloud_droplet_nc_index          = get_tracer_index(scm_state%tracer_names,"water_nc")
    scm_state%cloud_ice_nc_index              = get_tracer_index(scm_state%tracer_names,"ice_nc")
    scm_state%rain_nc_index                   = get_tracer_index(scm_state%tracer_names,"rain_nc")
    scm_state%snow_nc_index                   = get_tracer_index(scm_state%tracer_names,"snow_nc")
    scm_state%graupel_nc_index                = get_tracer_index(scm_state%tracer_names,"graupel_nc")
    scm_state%tke_index                       = get_tracer_index(scm_state%tracer_names,"sgs_tke")
    scm_state%water_friendly_aerosol_index    = get_tracer_index(scm_state%tracer_names,"liq_aero")
    scm_state%ice_friendly_aerosol_index      = get_tracer_index(scm_state%tracer_names,"ice_aero")
    scm_state%mass_weighted_rime_factor_index = get_tracer_index(scm_state%tracer_names,"q_rimef")
    
    scm_state%n_itt_swrad = int_zero
    scm_state%n_itt_lwrad = int_zero
    scm_state%n_itt_out = int_zero
    scm_state%n_levels_smooth = 5
    allocate(scm_state%blksz(n_columns))
    scm_state%blksz = int_one

    scm_state%sfc_flux_spec = .false.
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

    allocate(scm_state%pres_l(n_columns, n_levels), scm_state%pres_i(n_columns, n_levels+1), &
      scm_state%exner_l(n_columns, n_levels), scm_state%exner_i(n_columns, n_levels+1), &
      scm_state%geopotential_l(n_columns, n_levels), scm_state%geopotential_i(n_columns, n_levels+1))
    scm_state%pres_l = real_zero
    scm_state%pres_i = real_zero
    scm_state%exner_l = real_zero
    scm_state%exner_i = real_zero
    scm_state%geopotential_l = real_zero
    scm_state%geopotential_i = real_zero

    allocate(scm_state%a_k(n_levels+1), scm_state%b_k(n_levels+1), scm_state%si(n_columns, n_levels+1), &
      scm_state%sl(n_columns, n_levels))
    scm_state%a_k = real_zero
    scm_state%b_k = real_zero
    scm_state%si = real_zero
    scm_state%sl = real_zero

    allocate(scm_state%lat(n_columns), scm_state%lon(n_columns), scm_state%area(n_columns))
    scm_state%lat = real_zero
    scm_state%lon = real_zero
    scm_state%area = real_one

    allocate(scm_state%state_T(n_columns, n_levels, n_time_levels), &
      scm_state%state_u(n_columns, n_levels, n_time_levels), scm_state%state_v(n_columns, n_levels, n_time_levels), &
      scm_state%state_tracer(n_columns, n_levels, scm_state%n_tracers, n_time_levels))
    scm_state%state_T = real_zero
    scm_state%state_u = real_zero
    scm_state%state_v = real_zero
    scm_state%state_tracer = real_zero

    allocate(scm_state%temp_tracer(n_columns, n_levels, scm_state%n_tracers, n_time_levels), &
      scm_state%temp_T(n_columns, n_levels, n_time_levels), &
      scm_state%temp_u(n_columns, n_levels, n_time_levels), scm_state%temp_v(n_columns, n_levels, n_time_levels))
    scm_state%temp_tracer = real_zero
    scm_state%temp_T = real_zero
    scm_state%temp_u = real_zero
    scm_state%temp_v = real_zero

    allocate(scm_state%w_ls(n_columns, n_levels), scm_state%omega(n_columns, n_levels), scm_state%u_g(n_columns, n_levels), &
      scm_state%v_g(n_columns, n_levels), scm_state%dT_dt_rad(n_columns, n_levels), scm_state%h_advec_thil(n_columns, n_levels), &
      scm_state%h_advec_qt(n_columns, n_levels), scm_state%v_advec_thil(n_columns, n_levels), &
      scm_state%v_advec_qt(n_columns, n_levels), scm_state%pres_surf(n_columns), scm_state%T_surf(n_columns), &
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
    
    allocate(scm_state%sfc_type(n_columns))
    
    scm_state%sfc_type = real_zero
    
  end subroutine scm_state_create

  subroutine scm_input_create(scm_input, ntimes, nlev, nsoil, nsnow, nice)
    class(scm_input_type)             :: scm_input
    integer, intent(in)               :: ntimes, nlev, nsoil, nsnow, nice

    scm_input%input_nlev = nlev
    scm_input%input_ntimes = ntimes
    scm_input%input_nsoil = nsoil
    scm_input%input_nsnow = nsnow
    scm_input%input_nice = nice

    allocate(scm_input%input_pres(nlev),scm_input%input_time(ntimes))
    scm_input%input_pres = real_zero
    scm_input%input_time = real_zero

    allocate(scm_input%input_height(nlev), scm_input%input_thetail(nlev), scm_input%input_qt(nlev), scm_input%input_ql(nlev), &
      scm_input%input_qi(nlev), scm_input%input_u(nlev), scm_input%input_v(nlev), scm_input%input_tke(nlev), &
      scm_input%input_ozone(nlev), scm_input%input_temp(nlev))
    scm_input%input_height = real_zero
    scm_input%input_thetail = real_zero
    scm_input%input_temp = real_zero
    scm_input%input_qt = real_zero
    scm_input%input_ql = real_zero
    scm_input%input_qi = real_zero
    scm_input%input_u = real_zero
    scm_input%input_v = real_zero
    scm_input%input_tke = real_zero
    scm_input%input_ozone = real_zero
    
    allocate(scm_input%input_pres_i(nlev+1),scm_input%input_pres_l(nlev))
    scm_input%input_pres_i = real_zero
    scm_input%input_pres_l = real_zero

    allocate(scm_input%input_pres_surf(ntimes), &
      scm_input%input_T_surf(ntimes), scm_input%input_sh_flux_sfc(ntimes), scm_input%input_lh_flux_sfc(ntimes), &
      scm_input%input_w_ls(ntimes, nlev), scm_input%input_omega(ntimes, nlev), scm_input%input_u_g(ntimes, nlev), &
      scm_input%input_v_g(ntimes, nlev), scm_input%input_dT_dt_rad(ntimes, nlev), scm_input%input_h_advec_thetail(ntimes, nlev), &
      scm_input%input_h_advec_qt(ntimes, nlev), scm_input%input_v_advec_thetail(ntimes, nlev), &
      scm_input%input_v_advec_qt(ntimes, nlev), scm_input%input_u_nudge(ntimes, nlev), scm_input%input_v_nudge(ntimes, nlev),    &
      scm_input%input_T_nudge(ntimes, nlev), scm_input%input_thil_nudge(ntimes, nlev), scm_input%input_qt_nudge(ntimes, nlev))
    scm_input%input_pres_surf = real_zero
    scm_input%input_T_surf = real_zero
    scm_input%input_lat = real_zero
    scm_input%input_lon = real_zero
    
    allocate(scm_input%input_stc(nsoil), scm_input%input_smc(nsoil), scm_input%input_slc(nsoil))
    scm_input%input_stc = real_zero
    scm_input%input_smc = real_zero
    scm_input%input_slc = real_zero
    
    allocate(scm_input%input_snicexy(nsnow),scm_input%input_snliqxy(nsnow), scm_input%input_tsnoxy(nsnow), &
      scm_input%input_smoiseq(nsoil), scm_input%input_zsnsoxy(nsnow + nsoil))
    scm_input%input_snicexy = real_zero
    scm_input%input_snliqxy = real_zero
    scm_input%input_tsnoxy  = real_zero
    scm_input%input_smoiseq = real_zero
    scm_input%input_zsnsoxy = real_zero
    
    allocate(scm_input%input_tiice(nice))
    scm_input%input_tiice = real_zero
    
    allocate(scm_input%input_tslb(nsoil), scm_input%input_smois(nsoil), scm_input%input_sh2o(nsoil), &
        scm_input%input_smfr(nsoil), scm_input%input_flfr(nsoil))
    scm_input%input_tslb = real_zero
    scm_input%input_smois = real_zero
    scm_input%input_sh2o = real_zero
    scm_input%input_smfr = real_zero
    scm_input%input_flfr = real_zero
    
    scm_input%input_tsfco = real_zero
    scm_input%input_vegsrc = int_zero
    scm_input%input_vegtyp = real_zero
    scm_input%input_soiltyp = real_zero
    scm_input%input_slopetype = real_zero
    scm_input%input_vegfrac = real_zero
    scm_input%input_shdmin = real_zero
    scm_input%input_shdmax = real_zero
    scm_input%input_zorlo = real_zero
    scm_input%input_slmsk = real_zero
    scm_input%input_canopy = real_zero
    scm_input%input_hice = real_zero
    scm_input%input_fice = real_zero
    scm_input%input_tisfc = real_zero
    scm_input%input_snwdph = real_zero
    scm_input%input_snoalb = real_zero
    scm_input%input_sncovr = real_zero
    scm_input%input_area = real_zero
    scm_input%input_tg3 = real_zero
    scm_input%input_uustar = real_zero
    scm_input%input_alvsf        = real_zero
    scm_input%input_alnsf        = real_zero
    scm_input%input_alvwf        = real_zero
    scm_input%input_alnwf        = real_zero
    scm_input%input_facsf        = real_zero
    scm_input%input_facwf        = real_zero
    scm_input%input_weasd        = real_zero
    scm_input%input_f10m         = real_zero
    scm_input%input_t2m         = real_zero
    scm_input%input_q2m         = real_zero
    scm_input%input_ffmm         = real_zero
    scm_input%input_ffhh         = real_zero
    scm_input%input_tprcp         = real_zero
    scm_input%input_srflag         = real_zero
    scm_input%input_tsfcl         = real_zero
    scm_input%input_zorll         = real_zero
    scm_input%input_zorli         = real_zero
    scm_input%input_zorlw         = real_zero
    
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
    scm_input%input_oro          = real_zero
    scm_input%input_oro_uf       = real_zero
    scm_input%input_landfrac     = real_zero
    scm_input%input_lakefrac     = real_zero
    scm_input%input_lakedepth    = real_zero
    
    scm_input%input_tvxy = real_zero
    scm_input%input_tgxy = real_zero
    scm_input%input_tahxy = real_zero
    scm_input%input_canicexy = real_zero
    scm_input%input_canliqxy = real_zero
    scm_input%input_eahxy = real_zero
    scm_input%input_cmxy = real_zero
    scm_input%input_chxy = real_zero
    scm_input%input_fwetxy = real_zero
    scm_input%input_sneqvoxy = real_zero
    scm_input%input_alboldxy = real_zero
    scm_input%input_qsnowxy = real_zero
    scm_input%input_wslakexy = real_zero
    scm_input%input_taussxy = real_zero
    scm_input%input_waxy = real_zero
    scm_input%input_wtxy = real_zero
    scm_input%input_zwtxy = real_zero
    scm_input%input_xlaixy = real_zero
    scm_input%input_xsaixy = real_zero
    scm_input%input_lfmassxy = real_zero
    scm_input%input_stmassxy = real_zero
    scm_input%input_rtmassxy = real_zero
    scm_input%input_woodxy = real_zero
    scm_input%input_stblcpxy = real_zero
    scm_input%input_fastcpxy = real_zero
    scm_input%input_smcwtdxy = real_zero
    scm_input%input_deeprechxy = real_zero
    scm_input%input_rechxy = real_zero
    scm_input%input_snowxy = real_zero
    
    scm_input%input_tref       = real_zero  !< sea surface reference temperature for NSST (K)
    scm_input%input_z_c        = real_zero  !< sub-layer cooling thickness for NSST (m)
    scm_input%input_c_0        = real_zero  !< coefficient 1 to calculate d(Tz)/d(Ts) for NSST
    scm_input%input_c_d        = real_zero  !< coefficient 2 to calculate d(Tz)/d(Ts) for NSST
    scm_input%input_w_0        = real_zero  !< coefficient 3 to calculate d(Tz)/d(Ts) for NSST
    scm_input%input_w_d        = real_zero  !< coefficient 4 to calculate d(Tz)/d(Ts) for NSST
    scm_input%input_xt         = real_zero  !< heat content in diurnal thermocline layer for NSST (K m)
    scm_input%input_xs         = real_zero  !< salinity content in diurnal thermocline layer for NSST (ppt m)
    scm_input%input_xu         = real_zero  !< u-current in diurnal thermocline layer for NSST (m2 s-1)
    scm_input%input_xv         = real_zero  !< v-current in diurnal thermocline layer for NSST (m2 s-1)
    scm_input%input_xz         = real_zero  !< thickness of diurnal thermocline layer for NSST (m)
    scm_input%input_zm         = real_zero  !< thickness of ocean mixed layer for NSST (m)
    scm_input%input_xtts       = real_zero  !< sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST (m)
    scm_input%input_xzts       = real_zero  !< sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST (m K-1)
    scm_input%input_d_conv     = real_zero  !< thickness of free convection layer for NSST (m)
    scm_input%input_ifd        = real_zero  !< index to start DTM run for NSST
    scm_input%input_dt_cool    = real_zero  !< sub-layer cooling amount for NSST (K)
    scm_input%input_qrain      = real_zero  !< sensible heat due to rainfall for NSST (W)
    
    scm_input%input_wetness    = real_zero  !< normalized soil wetness for RUC LSM
    scm_input%input_clw_surf   = real_zero  !< cloud condensed water mixing ratio at surface for RUC LSM (kg kg-1)
    scm_input%input_qwv_surf   = real_zero  !< water vapor mixing ratio at surface for RUC LSM (kg kg-1)
    scm_input%input_tsnow      = real_zero  !< snow temperature at the bottom of the first snow layer for RUC LSM (K)
    scm_input%input_snowfallac = real_zero  !< run-total snow accumulation on the ground for RUC LSM (kg m-2)
    scm_input%input_acsnow     = real_zero  !< snow water equivalent of run-total frozen precip for RUC LSM (kg m-2)
    scm_input%input_lai        = real_zero  !< leaf area index for RUC LSM
    
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

    physics%Init_parm%me = int_zero
    physics%Init_parm%master = int_zero
    physics%Init_parm%isc = int_one
    physics%Init_parm%jsc = int_one
    physics%Init_parm%nx = int_one
    physics%Init_parm%ny = int_one
    physics%Init_parm%levs = int_one
    physics%Init_parm%cnx = int_one
    physics%Init_parm%cny = int_one
    physics%Init_parm%gnx = int_one
    physics%Init_parm%gny = int_one
    physics%Init_parm%nlunit = int_one
    physics%Init_parm%logunit= 10
    physics%Init_parm%bdat(:) = zeroes_8(:)
    physics%Init_parm%cdat(:) = zeroes_8(:)
    physics%Init_parm%dt_dycore = kind_phys_zero
    physics%Init_parm%dt_phys = kind_phys_zero
    physics%Init_parm%iau_offset = int_zero
    physics%Init_parm%ak => null()
    physics%Init_parm%bk => null()
    physics%Init_parm%xlon => null()
    physics%Init_parm%xlat => null()
    physics%Init_parm%area => null()
    physics%Init_parm%tracer_names => null()
    physics%Init_parm%blksz => null()
    physics%Init_parm%restart = .false.
    physics%Init_parm%hydrostatic = .true.
    physics%Init_parm%tile_num = int_one

  end subroutine physics_create

  subroutine physics_associate(physics, scm_state)
    class(physics_type) :: physics
    type(scm_state_type), target, intent(in) :: scm_state
    
    integer :: i

    physics%Statein%phii => scm_state%geopotential_i
    physics%Statein%prsi => scm_state%pres_i
    physics%Statein%prsik => scm_state%exner_i
    physics%Statein%phil => scm_state%geopotential_l
    physics%Statein%prsl => scm_state%pres_l
    physics%Statein%prslk => scm_state%exner_l

    physics%Statein%pgr => scm_state%pres_surf
    physics%Statein%ugrs => scm_state%state_u(:,:,1)
    physics%Statein%vgrs => scm_state%state_v(:,:,1)
    physics%Statein%vvl => scm_state%omega
    physics%Statein%tgrs => scm_state%state_T(:,:,1)
    physics%Statein%qgrs => scm_state%state_tracer(:,:,:,1)
    
    if (.not. scm_state%model_ics) then
      do i =1, physics%Model%ncols
        if (scm_state%sfc_type(i) == 0) then
          physics%Sfcprop%tsfco => scm_state%T_surf
        elseif (scm_state%sfc_type(i) == 1) then
          physics%Sfcprop%tsfcl => scm_state%T_surf
        elseif (scm_state%sfc_type(i) == 2) then
          physics%Sfcprop%tisfc => scm_state%T_surf
        end if
      end do
    end if
    
    if(scm_state%time_scheme == 2) then
      physics%Stateout%gu0 => scm_state%state_u(:,:,2)
      physics%Stateout%gv0 => scm_state%state_v(:,:,2)
      physics%Stateout%gt0 => scm_state%state_T(:,:,2)
      physics%Stateout%gq0 => scm_state%state_tracer(:,:,:,2)
    else
      physics%Stateout%gu0 => scm_state%state_u(:,:,1)
      physics%Stateout%gv0 => scm_state%state_v(:,:,1)
      physics%Stateout%gt0 => scm_state%state_T(:,:,1)
      physics%Stateout%gq0 => scm_state%state_tracer(:,:,:,1)
    endif

    if(scm_state%sfc_flux_spec) then
      physics%Sfcprop%spec_sh_flux => scm_state%sh_flux
      physics%Sfcprop%spec_lh_flux => scm_state%lh_flux
    endif
    
  end subroutine physics_associate
  
  subroutine physics_set(physics, scm_input, scm_state)
    !used for initializing variables in the physics DDT (but not pointer association); 
    !this should be utilized for variables that cannot be modified or "forced" by the SCM;
    !most of this routine follows what is in FV3/io/FV3GFS_io.F90/sfc_prop_restart_read
    use gmtb_scm_physical_constants, only: con_tice
    
    class(physics_type) :: physics
    type(scm_input_type), intent(in) :: scm_input
    type(scm_state_type), intent(in) :: scm_state
    
    integer :: i
    real(kind=dp) :: tem1, tem
    
    !double check under what circumstances these should actually be set from input!!! (these overwrite the initialzation in GFS_typedefs)
    do i = 1, physics%Model%ncols
      if (scm_state%model_ics) then
        !check for non-missing values
        physics%Sfcprop%hprime(i,1)  = scm_input%input_stddev
        physics%Sfcprop%hprime(i,2)  = scm_input%input_convexity
        physics%Sfcprop%hprime(i,3)  = scm_input%input_oa1
        physics%Sfcprop%hprime(i,4)  = scm_input%input_oa2
        physics%Sfcprop%hprime(i,5)  = scm_input%input_oa3
        physics%Sfcprop%hprime(i,6)  = scm_input%input_oa4
        physics%Sfcprop%hprime(i,7)  = scm_input%input_ol1
        physics%Sfcprop%hprime(i,8)  = scm_input%input_ol2
        physics%Sfcprop%hprime(i,9)  = scm_input%input_ol3
        physics%Sfcprop%hprime(i,10) = scm_input%input_ol4
        physics%Sfcprop%hprime(i,11) = scm_input%input_theta
        physics%Sfcprop%hprime(i,12) = scm_input%input_gamma
        physics%Sfcprop%hprime(i,13) = scm_input%input_sigma
        physics%Sfcprop%hprime(i,14) = scm_input%input_elvmax
        physics%Sfcprop%oro(i)       = scm_input%input_oro
        physics%Sfcprop%oro_uf(i)    = scm_input%input_oro_uf
        physics%Sfcprop%landfrac(i)  = scm_input%input_landfrac
        physics%Sfcprop%lakefrac(i)  = scm_input%input_lakefrac
        physics%Sfcprop%lakedepth(i) = scm_input%input_lakedepth
        
        physics%Sfcprop%slmsk(i)  = scm_input%input_slmsk
        physics%Sfcprop%tsfco(i)  = scm_input%input_tsfco
        physics%Sfcprop%weasd(i)  = scm_input%input_weasd
        physics%Sfcprop%tg3(i)    = scm_input%input_tg3
        physics%Sfcprop%zorlo(i)  = scm_input%input_zorlo
        physics%Sfcprop%alvsf(i)  = scm_input%input_alvsf
        physics%Sfcprop%alnsf(i)  = scm_input%input_alnsf
        physics%Sfcprop%alvwf(i)  = scm_input%input_alvwf
        physics%Sfcprop%alnwf(i)  = scm_input%input_alnwf
        physics%Sfcprop%facsf(i)  = scm_input%input_facsf
        physics%Sfcprop%facwf(i)  = scm_input%input_facwf
        physics%Sfcprop%vfrac(i)  = scm_input%input_vegfrac
        !GJF: is this needed anymore (not in FV3GFS_io)?
        physics%Interstitial%sigmaf(i) = min(scm_input%input_vegfrac,0.01)
        physics%Sfcprop%canopy(i) = scm_input%input_canopy
        physics%Sfcprop%f10m(i)   = scm_input%input_f10m
        physics%Sfcprop%t2m(i)    = scm_input%input_t2m
        physics%Sfcprop%q2m(i)    = scm_input%input_q2m
        physics%Sfcprop%vtype(i)  = scm_input%input_vegtyp
        physics%Sfcprop%stype(i)  = scm_input%input_soiltyp
        physics%Sfcprop%uustar(i) = scm_input%input_uustar
        physics%Sfcprop%ffmm(i)   = scm_input%input_ffmm
        physics%Sfcprop%ffhh(i)   = scm_input%input_ffhh
        physics%Sfcprop%hice(i)   = scm_input%input_hice
        physics%Sfcprop%fice(i)   = scm_input%input_fice
        physics%Sfcprop%tisfc(i)  = scm_input%input_tisfc
        physics%Sfcprop%tprcp(i)  = scm_input%input_tprcp
        physics%Sfcprop%srflag(i) = scm_input%input_srflag
        physics%Sfcprop%snowd(i)  = scm_input%input_snwdph
        physics%Sfcprop%shdmin(i) = scm_input%input_shdmin
        physics%Sfcprop%shdmax(i) = scm_input%input_shdmax
        physics%Sfcprop%slope(i)  = scm_input%input_slopetype
        physics%Sfcprop%snoalb(i) = scm_input%input_snoalb
        physics%Sfcprop%sncovr(i) = scm_input%input_sncovr
        physics%Sfcprop%tsfcl(i)  = scm_input%input_tsfcl
        physics%Sfcprop%zorll(i)  = scm_input%input_zorll
        physics%Sfcprop%zorli(i)  = scm_input%input_zorli
        if (physics%Model%cplwav) then
          physics%Sfcprop%zorlw(i)  = scm_input%input_zorlw
        else
          physics%Sfcprop%zorlw(i)  = physics%Sfcprop%zorlo(i)
        end if
      else
        physics%Sfcprop%slmsk(i) = scm_state%sfc_type(i)
        ! tsfco is already pointing to T_surf forcing in physics_associate
        ! physics%Sfcprop%tsfco(i) => scm_state%T_surf
        physics%Sfcprop%zorlo(i) = scm_state%sfc_roughness_length_cm(i)
        ! tisfc is already pointing to T_surf forcing in physics_associate
        ! physics%Sfcprop%tisfc(i) => scm_state%T_surf
        ! tsfcl is already pointing to T_surf forcing in physics_associate
        ! physics%Sfcprop%tsfcl(i) => scm_state%T_surf
        if (physics%Sfcprop%slmsk(i) > 1.9_dp) physics%Sfcprop%fice(i) = 1.0 !needed to calculate tsfc and zorl below when model_ics == .false.
      end if
      
      !check for nonmissing values (this overwrites what is in the suite namelist file -- is that desirable?)
      if (scm_state%model_ics) then
        physics%Model%ivegsrc = scm_input%input_vegsrc
      end if
      
      if(scm_state%model_ics .and. physics%Model%frac_grid) then ! obtain slmsk from landfrac
        !! next 5 lines are temporary till lake model is available
        if (physics%Sfcprop%lakefrac(i) > real_zero) then
          !physics%Sfcprop%lakefrac(i) = nint(physics%Sfcprop%lakefrac(i))
          physics%Sfcprop%landfrac(i) = real_one - physics%Sfcprop%lakefrac(i)
          if (physics%Sfcprop%lakefrac(i) == real_zero) physics%Sfcprop%fice(i) = real_zero
        endif 
        physics%Sfcprop%slmsk(i) = ceiling(physics%Sfcprop%landfrac(i))
        if (physics%Sfcprop%fice(i) > physics%Model%min_lakeice .and. physics%Sfcprop%landfrac(i) == real_zero) physics%Sfcprop%slmsk(i) = 2 ! land dominates ice if co-exist
      else ! obtain landfrac from slmsk
        if (physics%Sfcprop%slmsk(i) > 1.9_dp) then
          physics%Sfcprop%landfrac(i) = real_zero
        else
          physics%Sfcprop%landfrac(i) = physics%Sfcprop%slmsk(i)
        endif
      endif
      
      if (physics%Sfcprop%lakefrac(i) > real_zero) then
        physics%Sfcprop%oceanfrac(i) = real_zero ! lake & ocean don't coexist in a cell
        if (physics%Sfcprop%fice(i) < physics%Model%min_lakeice) then
           physics%Sfcprop%fice(i) = real_zero
           if (physics%Sfcprop%slmsk(i) == 2) physics%Sfcprop%slmsk(i) = 0
        endif
      else
        physics%Sfcprop%oceanfrac(i) = real_one - physics%Sfcprop%landfrac(i)
        if (physics%Sfcprop%fice(i) < physics%Model%min_seaice) then
           physics%Sfcprop%fice(i) = real_zero
           if (physics%Sfcprop%slmsk(i) == 2) physics%Sfcprop%slmsk(i) = 0
        endif
      endif
      
      !--- NSSTM variables
      if (physics%Model%nstf_name(1) > 0) then
        if (physics%Model%nstf_name(2) == 1 .or. .not. scm_state%model_ics) then             ! nsst spinup
          !--- nsstm tref
          physics%Sfcprop%tref(i)    = physics%Sfcprop%tsfco(i)
          physics%Sfcprop%z_c(i)     = real_zero
          physics%Sfcprop%c_0(i)     = real_zero
          physics%Sfcprop%c_d(i)     = real_zero
          physics%Sfcprop%w_0(i)     = real_zero
          physics%Sfcprop%w_d(i)     = real_zero
          physics%Sfcprop%xt(i)      = real_zero
          physics%Sfcprop%xs(i)      = real_zero
          physics%Sfcprop%xu(i)      = real_zero
          physics%Sfcprop%xv(i)      = real_zero
          physics%Sfcprop%xz(i)      = 30.0_dp
          physics%Sfcprop%zm(i)      = real_zero
          physics%Sfcprop%xtts(i)    = real_zero
          physics%Sfcprop%xzts(i)    = real_zero
          physics%Sfcprop%d_conv(i)  = real_zero
          physics%Sfcprop%ifd(i)     = real_zero
          physics%Sfcprop%dt_cool(i) = real_zero
          physics%Sfcprop%qrain(i)   = real_zero
        elseif (physics%Model%nstf_name(2) == 0) then         ! nsst restart
          !check for non-missing values
          physics%Sfcprop%tref(i)    = scm_input%input_tref
          physics%Sfcprop%z_c(i)     = scm_input%input_z_c
          physics%Sfcprop%c_0(i)     = scm_input%input_c_0
          physics%Sfcprop%c_d(i)     = scm_input%input_c_d
          physics%Sfcprop%w_0(i)     = scm_input%input_w_0
          physics%Sfcprop%w_d(i)     = scm_input%input_w_d
          physics%Sfcprop%xt(i)      = scm_input%input_xt
          physics%Sfcprop%xs(i)      = scm_input%input_xs
          physics%Sfcprop%xu(i)      = scm_input%input_xu
          physics%Sfcprop%xv(i)      = scm_input%input_xv
          physics%Sfcprop%xz(i)      = scm_input%input_xz
          physics%Sfcprop%zm(i)      = scm_input%input_zm
          physics%Sfcprop%xtts(i)    = scm_input%input_xtts
          physics%Sfcprop%xzts(i)    = scm_input%input_xzts
          physics%Sfcprop%d_conv(i)  = scm_input%input_d_conv
          physics%Sfcprop%ifd(i)     = scm_input%input_ifd
          physics%Sfcprop%dt_cool(i) = scm_input%input_dt_cool
          physics%Sfcprop%qrain(i)   = scm_input%input_qrain
        endif
      endif
      
      if (scm_state%model_ics .and. physics%Model%lsm == physics%Model%lsm_ruc) then !.and. warm_start (not implemented here -- assuming that RUC LSM has warm start data from file)
        !--- Extra RUC LSM variables
        !check for nonmissing values
        physics%Sfcprop%wetness(i)    = scm_input%input_wetness
        physics%Sfcprop%clw_surf(i)   = scm_input%input_clw_surf
        physics%Sfcprop%qwv_surf(i)   = scm_input%input_qwv_surf
        physics%Sfcprop%tsnow(i)      = scm_input%input_tsnow
        physics%Sfcprop%snowfallac(i) = scm_input%input_snowfallac
        physics%Sfcprop%acsnow(i)     = scm_input%input_acsnow
        if (physics%Model%lsm == physics%Model%lsm_ruc .and. physics%Model%rdlai) then
          physics%Sfcprop%xlaixy(i) = scm_input%input_lai
        end if
      elseif (scm_state%model_ics .and. physics%Model%lsm == physics%Model%lsm_noahmp) then
        !check for nonmissing values
        physics%Sfcprop%snowxy    (i) = scm_input%input_snowxy
        physics%Sfcprop%tvxy      (i) = scm_input%input_tvxy
        physics%Sfcprop%tgxy      (i) = scm_input%input_tgxy
        physics%Sfcprop%canicexy  (i) = scm_input%input_canicexy
        physics%Sfcprop%canliqxy  (i) = scm_input%input_canliqxy
        physics%Sfcprop%eahxy     (i) = scm_input%input_eahxy
        physics%Sfcprop%tahxy     (i) = scm_input%input_tahxy
        physics%Sfcprop%cmxy      (i) = scm_input%input_cmxy
        physics%Sfcprop%chxy      (i) = scm_input%input_chxy
        physics%Sfcprop%fwetxy    (i) = scm_input%input_fwetxy
        physics%Sfcprop%sneqvoxy  (i) = scm_input%input_sneqvoxy
        physics%Sfcprop%alboldxy  (i) = scm_input%input_alboldxy
        physics%Sfcprop%qsnowxy   (i) = scm_input%input_qsnowxy
        physics%Sfcprop%wslakexy  (i) = scm_input%input_wslakexy
        physics%Sfcprop%zwtxy     (i) = scm_input%input_zwtxy
        physics%Sfcprop%waxy      (i) = scm_input%input_waxy
        physics%Sfcprop%wtxy      (i) = scm_input%input_wtxy
        physics%Sfcprop%lfmassxy  (i) = scm_input%input_lfmassxy
        physics%Sfcprop%rtmassxy  (i) = scm_input%input_rtmassxy
        physics%Sfcprop%stmassxy  (i) = scm_input%input_stmassxy
        physics%Sfcprop%woodxy    (i) = scm_input%input_woodxy
        physics%Sfcprop%stblcpxy  (i) = scm_input%input_stblcpxy
        physics%Sfcprop%fastcpxy  (i) = scm_input%input_fastcpxy
        physics%Sfcprop%xsaixy    (i) = scm_input%input_xsaixy
        physics%Sfcprop%xlaixy    (i) = scm_input%input_xlaixy
        physics%Sfcprop%taussxy   (i) = scm_input%input_taussxy
        physics%Sfcprop%smcwtdxy  (i) = scm_input%input_smcwtdxy
        physics%Sfcprop%deeprechxy(i) = scm_input%input_deeprechxy
        physics%Sfcprop%rechxy    (i) = scm_input%input_rechxy
      end if
      
      if (scm_state%model_ics .and. physics%Model%lsm == physics%Model%lsm_noah .or. &
          physics%Model%lsm == physics%Model%lsm_noahmp) then    !.or. (.not. warm_start) from FV3GFS_io is not implemented
        !check for nonmissing values
        !--- 3D variables
        ! do lsoil = 1,physics%Model%lsoil
        !   physics%Sfcprop%stc(i,lsoil) = scm_input%input_stc(lsoil)    !--- stc
        !   physics%Sfcprop%smc(i,lsoil) = scm_input%input_smc(lsoil)    !--- smc
        !   physics%Sfcprop%slc(i,lsoil) = scm_input%input_slc(lsoil)    !--- slc
        ! enddo
        physics%Sfcprop%stc(i,:) = scm_input%input_stc(:)    !--- stc
        physics%Sfcprop%smc(i,:) = scm_input%input_smc(:)    !--- smc
        physics%Sfcprop%slc(i,:) = scm_input%input_slc(:)    !--- slc

        if (physics%Model%lsm == physics%Model%lsm_noahmp) then
          ! do lsoil = -2, 0
          !   physics%Sfcprop%snicexy(i,lsoil) = scm_input%input_snicexy(lsoil)
          !   physics%Sfcprop%snliqxy(i,lsoil) = scm_input%input_snliqxy(lsoil)
          !   physics%Sfcprop%tsnoxy(i,lsoil)  = scm_input%input_tsnoxy(lsoil)
          ! enddo
          physics%Sfcprop%snicexy(i,:) = scm_input%input_snicexy(:)
          physics%Sfcprop%snliqxy(i,:) = scm_input%input_snliqxy(:)
          physics%Sfcprop%tsnoxy(i,:)  = scm_input%input_tsnoxy(:)

          ! do lsoil = 1, 4
          !   physics%Sfcprop%smoiseq(i,lsoil)  = scm_input%input_smoiseq(lsoil) 
          ! enddo
          physics%Sfcprop%smoiseq(i,:)  = scm_input%input_smoiseq(:)

          ! do lsoil = -2, 4
          !   physics%Sfcprop%zsnsoxy(i,lsoil)  = scm_input%input_zsnsoxy(lsoil)
          ! enddo
          physics%Sfcprop%zsnsoxy(i,:)  = scm_input%input_zsnsoxy(:)
        endif

      else if (scm_state%model_ics .and. physics%Model%lsm == physics%Model%lsm_ruc) then
        !--- 3D variables
        ! do lsoil = 1,Model%lsoil_lsm
        !   physics%Sfcprop%tslb(i,lsoil) = scm_input%input_tslb(lsoil)
        !   physics%Sfcprop%smois(i,lsoil) = scm_input%input_smois(lsoil)
        !   physics%Sfcprop%sh2o(i,lsoil) = scm_input%input_sh2o(lsoil)
        !   physics%Sfcprop%keepsmfr(i,lsoil) = scm_input%input_smfr(lsoil)
        !   physics%Sfcprop%flag_frsoil(i,lsoil) = scm_input%input_flfr(lsoil)
        ! enddo
        physics%Sfcprop%tslb(i,:) = scm_input%input_tslb(:)
        physics%Sfcprop%smois(i,:) = scm_input%input_smois(:)
        physics%Sfcprop%sh2o(i,:) = scm_input%input_sh2o(:)
        physics%Sfcprop%keepsmfr(i,:) = scm_input%input_smfr(:)
        physics%Sfcprop%flag_frsoil(i,:) = scm_input%input_flfr(:)
      end if

      if (scm_state%model_ics) then
        !check for nonmissing values
        ! do k = 1,Model%kice
        !   physics%Sfcprop%tiice(i,k) = scm_input%input_tiice(k)   !--- internal ice temp
        ! enddo
        physics%Sfcprop%tiice(i,:) = scm_input%input_tiice(:)
      end if
      
      !GJF: computing sncovr if model_ics and sncovr is missing:
      ! Sfcprop(nb)%sncovr(ix) = zero
      ! if (Sfcprop(nb)%landfrac(ix) >= drythresh .or. Sfcprop(nb)%fice(ix) >= Model%min_seaice) then
      !   vegtyp = Sfcprop(nb)%vtype(ix)
      !   if (vegtyp == 0) vegtyp = 7
      !   rsnow  = 0.001_r8*Sfcprop(nb)%weasd(ix)/snupx(vegtyp)
      !   if (0.001_r8*Sfcprop(nb)%weasd(ix) < snupx(vegtyp)) then
      !     Sfcprop(nb)%sncovr(ix) = one - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
      !   else
      !     Sfcprop(nb)%sncovr(ix) = one
      !   endif
      ! endif
      
      if (scm_input%input_tsfcl < real_zero) then
        physics%Sfcprop%tsfcl(i) = physics%Sfcprop%tsfco(i) !--- compute tsfcl from existing variables
      end if
      
      if (scm_input%input_zorll < real_zero) then
        physics%Sfcprop%zorll(i) = physics%Sfcprop%zorlo(i) !--- compute zorll from existing variables
      end if
      
      if (scm_input%input_zorli < real_zero) then
        physics%Sfcprop%zorli(i) = physics%Sfcprop%zorlo(i) !--- compute zorli from existing variables
      end if
      
      if (scm_input%input_zorlw < real_zero) then
        physics%Sfcprop%zorlw(i) = physics%Sfcprop%zorlo(i) !--- compute zorlw from existing variables
      end if
      
      if(physics%Model%frac_grid .and. scm_state%model_ics) then ! 3-way composite
            physics%Sfcprop%tsfco(i) = max(con_tice, physics%Sfcprop%tsfco(i))
            tem1 = real_one - physics%Sfcprop%landfrac(i)
            tem  = tem1 * physics%Sfcprop%fice(i) ! tem = ice fraction wrt whole cell
            physics%Sfcprop%zorl(i) = physics%Sfcprop%zorll(i) * physics%Sfcprop%landfrac(i) &
                                 + physics%Sfcprop%zorli(i) * tem                      &
                                 + physics%Sfcprop%zorlo(i) * (tem1-tem)

            physics%Sfcprop%tsfc(i) = physics%Sfcprop%tsfcl(i) * physics%Sfcprop%landfrac(i) &
                                 + physics%Sfcprop%tisfc(i) * tem                      &
                                 + physics%Sfcprop%tsfco(i) * (tem1-tem)
      else
          !--- specify tsfcl/zorll/zorli from existing variable tsfco/zorlo
  !         physics%Sfcprop%tsfcl(i) = physics%Sfcprop%tsfco(i)
  !         physics%Sfcprop%zorll(i) = physics%Sfcprop%zorlo(i)
  !         physics%Sfcprop%zorli(i) = physics%Sfcprop%zorlo(i)
  !         physics%Sfcprop%zorl(i)  = physics%Sfcprop%zorlo(i)
  !         physics%Sfcprop%tsfc(i)  = physics%Sfcprop%tsfco(i)
            if (physics%Sfcprop%slmsk(i) == 1) then
              physics%Sfcprop%zorl(i) = physics%Sfcprop%zorll(i) 
              physics%Sfcprop%tsfc(i) = physics%Sfcprop%tsfcl(i)
            else
              tem = real_one - physics%Sfcprop%fice(i)
              physics%Sfcprop%zorl(i) = physics%Sfcprop%zorli(i) * physics%Sfcprop%fice(i) &
                                   + physics%Sfcprop%zorlo(i) * tem

              physics%Sfcprop%tsfc(i) = physics%Sfcprop%tisfc(i) * physics%Sfcprop%fice(i) &
                                   + physics%Sfcprop%tsfco(i) * tem
            endif
      endif ! if (Model%frac_grid)
      
      if (scm_state%model_ics .and. MAXVAL(scm_input%input_tiice) < real_zero) then
        physics%Sfcprop%tiice(i,1) = physics%Sfcprop%stc(i,1) !--- initialize internal ice temp from soil temp at layer 1
        physics%Sfcprop%tiice(i,2) = physics%Sfcprop%stc(i,2) !--- initialize internal ice temp from soil temp at layer 2
      end if
      
      !!!!NoahMP cold start code
      
    end do
    
  end subroutine physics_set
  
  function get_tracer_index (tracer_names, name)

    character(len=32), intent(in) :: tracer_names(:)
    character(len=*),  intent(in) :: name
    
    !--- local variables
    integer :: get_tracer_index
    integer :: i
    integer, parameter :: no_tracer = -99

    get_tracer_index = no_tracer

    do i=1, size(tracer_names)
       if (trim(name) == trim(tracer_names(i))) then
           get_tracer_index = i
           exit
       endif
    enddo

    if (get_tracer_index == no_tracer) then
      print *,'tracer with name '//trim(name)//' not found'
    else
      print *,'tracer FOUND:',trim(name)
    endif

    return
  end function get_tracer_index

end module gmtb_scm_type_defs
