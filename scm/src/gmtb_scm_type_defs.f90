!> \file gmtb_scm_type_defs.f90
!!  Contains type definitions for SCM-related variables and physics-related variables

module gmtb_scm_type_defs

  use gmtb_scm_kinds, only : sp, dp, qp
  use GFS_typedefs, only: GFS_control_type, GFS_statein_type, GFS_stateout_type, GFS_sfcprop_type, GFS_coupling_type, &
    GFS_grid_type, GFS_tbd_type, GFS_cldprop_type, GFS_radtend_type, GFS_diag_type, GFS_sfccycle_type, GFS_interstitial_type, &
    GFS_init_type, LTP
  use machine, only: kind_phys

  implicit none

  integer, parameter :: character_length = 80
  integer, parameter :: int_zero = 0
  integer, parameter :: int_one = 1
  real(kind=dp), parameter :: real_zero = 0.0
  real(kind=dp), parameter :: real_one = 1.0

  character(len = character_length) :: clear_char = ''

  type scm_state_type

    character(len=character_length)                 :: experiment_name !> name of model configuration file
    character(len=character_length)                 :: model_name !< name of "host" model (must be "GFS" for prototype)
    character(len=character_length)                 :: output_dir !< name of output directory to place netCDF file
    character(len=character_length)                 :: physics_suite_dir !< location of the physics suite XML files for the IPD (relative to the executable path)
    character(len=character_length)                 :: case_data_dir !< location of the case initialization and forcing data files (relative to the executable path)
    character(len=character_length)                 :: vert_coord_data_dir !< location of the vertical coordinate data files (relative to the executable path)
    character(len=character_length)                 :: output_file !< name of output file (without the file extension)
    character(len=character_length)                 :: case_name !< name of case initialization and forcing to use (different than experiment name, which names the model run (as a control, experiment_1, etc.))
    character(len=character_length), allocatable    :: physics_suite_name(:) !< name of physics suite (must be "GFS_operational" for prototype)
    character(len=65), allocatable    :: physics_nml(:)

    integer                           :: n_levels !< number of model levels (must be 64 for prototype)
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
    integer                           :: init_year, init_month, init_day, init_hour
    character(len=32), allocatable    :: tracer_names(:) !< name of physics suite (must be "GFS_operational" for prototype)
    integer, allocatable              :: blksz(:)

    logical                           :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
    integer                           :: sfc_type !< 0: sea surface, 1: land surface, 2: sea-ice surface
    real(kind=dp)                     :: sfc_type_real(1)
    integer                           :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
    integer                           :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
    integer                           :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere

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

    contains
      procedure :: create  => scm_state_create

  end type scm_state_type

  type scm_input_type
    !> - Define the case-specific initialization and forcing variables.
    integer                           :: input_nlev !< number of levels in the input file
    integer                           :: input_ntimes !< number of times in the input file where forcing is available
    real(kind=dp), allocatable              :: input_pres(:) !< pressure (Pa) of input levels
    real(kind=dp), allocatable              :: input_time(:) !< time (s) since beginning of forcing file
    real(kind=dp), allocatable              :: input_height(:) !< height of input levels (m) (initial)
    real(kind=dp), allocatable              :: input_thetail(:) !< ice-liquid water potential temperature(K) (initial)
    real(kind=dp), allocatable              :: input_qt(:) !< total water specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_ql(:) !< suspended liquid specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_qi(:) !< suspended ice specific humidity (kg/kg) (initial)
    real(kind=dp), allocatable              :: input_u(:) !< zonal wind (m/s) (initial)
    real(kind=dp), allocatable              :: input_v(:) !< meridional wind (m/s) (initial)
    real(kind=dp), allocatable              :: input_tke(:) !< turbulence kinetic energy (m^2/s^2) (initial)
    real(kind=dp), allocatable              :: input_ozone(:) !< ozone mass mixing ratio (kg/kg) (initial)
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

  end type scm_input_type

  type scm_reference_type
    !> - Define the reference profile variables.
    integer                           :: ref_nlev !< number of levels in the reference profile
    real(kind=dp), allocatable              :: ref_pres(:) !< pressure (Pa) of the reference profile levels
    real(kind=dp), allocatable              :: ref_T(:) !< absolute T (K) of the reference profile levels
    real(kind=dp), allocatable              :: ref_qv(:) !< water vapor specific humidity (kg/kg) of the reference profile levels
    real(kind=dp), allocatable              :: ref_ozone(:) !< ozone mass mixing ratio (kg/kg) of the reference profile level

    contains
      procedure :: create => scm_reference_create

  end type scm_reference_type

!> \section arg_table_physics_type
!! | local_name                       | standard_name                                               | long_name                                             | units         | rank | type                  |    kind   | intent | optional |
!! |--------------------------------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | physics%Model(i)                     | FV3-GFS_Control_type                                   | derived type GFS_control_type in FV3                    | DDT           |    0 | GFS_control_type      |           | none   | F        |
!! | physics%Cldprop(i)                   | FV3-GFS_Cldprop_type                                   | derived type GFS_cldprop_type in FV3                    | DDT           |    0 | GFS_cldprop_type      |           | none   | F        |
!! | physics%Coupling(i)                  | FV3-GFS_Coupling_type                                  | derived type GFS_coupling_type in FV3                   | DDT           |    0 | GFS_coupling_type     |           | none   | F        |
!! | physics%Diag(i)                      | FV3-GFS_Diag_type                                      | derived type GFS_diag_type in FV3                       | DDT           |    0 | GFS_diag_type         |           | none   | F        |
!! | physics%Grid(i)                      | FV3-GFS_Grid_type                                      | derived type GFS_grid_type in FV3                       | DDT           |    0 | GFS_grid_type         |           | none   | F        |
!! | physics%Radtend(i)                   | FV3-GFS_Radtend_type                                   | derived type GFS_radtend_type in FV3                    | DDT           |    0 | GFS_radtend_type      |           | none   | F        |
!! | physics%Sfccycle(i)                  | FV3-GFS_Sfccycle_type                                  | derived type GFS_sfccycle_type in FV3                   | DDT           |    0 | GFS_sfccycle_type     |           | none   | F        |
!! | physics%Sfcprop(i)                   | FV3-GFS_Sfcprop_type                                   | derived type GFS_sfcprop_type in FV3                    | DDT           |    0 | GFS_sfcprop_type      |           | none   | F        |
!! | physics%Statein(i)                   | FV3-GFS_Statein_type                                   | derived type GFS_statein_type in FV3                    | DDT           |    0 | GFS_statein_type      |           | none   | F        |
!! | physics%Stateout(i)                  | FV3-GFS_Stateout_type                                  | derived type GFS_stateout_type in FV3                   | DDT           |    0 | GFS_stateout_type     |           | none   | F        |
!! | physics%Tbd(i)                       | FV3-GFS_Tbd_type                                       | derived type GFS_tbd_type in FV3                        | DDT           |    0 | GFS_tbd_type          |           | none   | F        |
!! | physics%Interstitial(i)              | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | none   | F        |
!! | physics%Init_parm(i)                 | FV3-GFS_Init_type                                      | dervied type GFS_init_type in FV3                       | DDT           |    0 | GFS_init_type         |           | none   | F        |
!! | physics%Interstitial(i)%adjnirbmd                     | surface_downwelling_direct_near_infrared_shortwave_flux                                        | surface downwelling beam near-infrared shortwave flux at current time               | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjnirbmu                     | surface_upwelling_direct_near_infrared_shortwave_flux                                          | surface upwelling beam near-infrared shortwave flux at current time                 | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjnirdfd                     | surface_downwelling_diffuse_near_infrared_shortwave_flux                                       | surface downwelling diffuse near-infrared shortwave flux at current time            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjnirdfu                     | surface_upwelling_diffuse_near_infrared_shortwave_flux                                         | surface upwelling diffuse near-infrared shortwave flux at current time              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjsfcdlw                     | surface_downwelling_longwave_flux                                                              | surface downwelling longwave flux at current time                                   | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjsfcdsw                     | surface_downwelling_shortwave_flux                                                             | surface downwelling shortwave flux at current time                                  | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjsfcnsw                     | surface_net_downwelling_shortwave_flux                                                         | surface net downwelling shortwave flux at current time                              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjsfculw                     | surface_upwelling_longwave_flux                                                                | surface upwelling longwave flux at current time                                     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjvisbmd                     | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux                              | surface downwelling beam ultraviolet plus visible shortwave flux at current time    | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjvisbmu                     | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux                                | surface upwelling beam ultraviolet plus visible shortwave flux at current time      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjvisdfu                     | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux                               | surface upwelling diffuse ultraviolet plus visible shortwave flux at current time   | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%adjvisdfd                     | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux                             | surface downwelling diffuse ultraviolet plus visible shortwave flux at current time | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%aerodp                        | atmosphere_optical_thickness_due_to_ambient_aerosol_particles                                  | vertical integrated optical depth for various aerosol species                       | none          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cd                            | surface_drag_coefficient_for_momentum_in_air                                                   | surface exchange coeff for momentum                                                 | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cdq                           | surface_drag_coefficient_for_heat_and_moisture_in_air                                          | surface exchange coeff heat & moisture                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cice                          | sea_ice_concentration_for_physics                                                              | sea-ice concentration [0,1]                                                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cldf                          | cloud_area_fraction                                                                            | fraction of grid box area in which updrafts occur                                   | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cldsa                         | cloud_area_fraction_for_radiation                                                              | fraction of clouds for low, middle, high, total and BL                              | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cld1d                         | cloud_work_function                                                                            | cloud work function                                                                 | m2 s-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds                        |                                                                                                |                                                                                     | various       |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,1)                 | total_cloud_fraction                                                                           | layer total cloud fraction                                                          | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,2)                 | cloud_liquid_water_path                                                                        | layer cloud liquid water path                                                       | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,3)                 | mean_effective_radius_for_liquid_cloud                                                         | mean effective radius for liquid cloud                                              | micron        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,4)                 | cloud_ice_water_path                                                                           | layer cloud ice water path                                                          | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,5)                 | mean_effective_radius_for_ice_cloud                                                            | mean effective radius for ice cloud                                                 | micron        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,6)                 | cloud_rain_water_path                                                                          | cloud rain water path                                                               | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,7)                 | mean_effective_radius_for_rain_drop                                                            | mean effective radius for rain drop                                                 | micron        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,8)                 | cloud_snow_water_path                                                                          | cloud snow water path                                                               | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clouds(:,:,9)                 | mean_effective_radius_for_snow_flake                                                           | mean effective radius for snow flake                                                | micron        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clw                           | convective_transportable_tracers                                                               | array to contain cloud water and other convective trans. tracers                    | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clw(:,:,1)                    | cloud_ice_specific_humidity                                                                    | cloud ice specific humidity                                                         | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clw(:,:,2)                    | cloud_liquid_water_specific_humidity                                                           | cloud water specific humidity                                                       | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%clx                           | fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height                        | frac. of grid box with by subgrid orography higher than critical height             | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cnvc                          | convective_cloud_cover                                                                         | convective cloud cover                                                              | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cnvw                          | convective_cloud_water_specific_humidity                                                       | convective cloud water specific humidity                                            | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%cumabs                        | maximum_column_heating_rate                                                                    | maximum heating rate in column                                                      | K s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dd_mf                         | instantaneous_atmosphere_downdraft_convective_mass_flux                                        | (downdraft mass flux) * delt                                                        | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%del                           | air_pressure_difference_between_midlayers                                                      | air pressure difference between midlayers                                           | Pa            |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%del_gz                        | geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature              | difference between mid-layer geopotentials divided by mid-layer virtual temperature | m2 s-2 K-1    |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dkt                           | atmosphere_heat_diffusivity                                                                    | diffusivity for heat                                                                | m2 s-1        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dlength                       | characteristic_grid_length_scale                                                               | representative horizontal length scale of grid box                                  | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%domip                         | dominant_sleet_type                                                                            | dominant sleet type                                                                 | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%domr                          | dominant_rain_type                                                                             | dominant rain type                                                                  | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%doms                          | dominant_snow_type                                                                             | dominant snow type                                                                  | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%domzr                         | dominant_freezing_rain_type                                                                    | dominant freezing rain type                                                         | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dqdt                          | tendency_of_tracers_due_to_model_physics                                                       | updated tendency of the tracers                                                     | kg kg-1 s-1   |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dqsfc1                        | instantaneous_surface_upward_latent_heat_flux                                                  | surface upward latent heat flux                                                     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dq3dt_loc                     |                                                                                                |                                                                                     |               |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dq3dt_loc(:,:,6:6+physics%Interstitial(i)%oz_coeff-1) | change_in_ozone_concentration                                             | change in ozone concentration                                                       | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%drain                         | subsurface_runoff_flux                                                                         | subsurface runoff flux                                                              | g m-2 s-1     |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dtdt                          | tendency_of_air_temperature_due_to_model_physics                                               | air temperature tendency due to model physics                                       | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dtdtc                         | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky                        | clear sky radiative (shortwave + longwave) heating rate at current time             | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dtsfc1                        | instantaneous_surface_upward_sensible_heat_flux                                                | surface upward sensible heat flux                                                   | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dtzm                          | mean_change_over_depth_in_sea_water_temperature                                                | mean of dT(z)  (zsea1 to zsea2)                                                     | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dt_mf                         | instantaneous_atmosphere_detrainment_convective_mass_flux                                      | (detrainment mass flux) * delt                                                      | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dudt                          | tendency_of_x_wind_due_to_model_physics                                                        | zonal wind tendency due to model physics                                            | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dusfcg                        | instantaneous_x_stress_due_to_gravity_wave_drag                                                | zonal surface stress due to orographic gravity wave drag                            | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dusfc1                        | instantaneous_surface_x_momentum_flux                                                          | x momentum flux                                                                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dvdt                          | tendency_of_y_wind_due_to_model_physics                                                        | meridional wind tendency due to model physics                                       | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dvsfcg                        | instantaneous_y_stress_due_to_gravity_wave_drag                                                | meridional surface stress due to orographic gravity wave drag                       | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%dvsfc1                        | instantaneous_surface_y_momentum_flux                                                          | y momentum flux                                                                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%elvmax                        | maximum_subgrid_orography                                                                      | maximum of subgrid orography                                                        | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%ep1d                          | surface_upward_potential_latent_heat_flux                                                      | surface upward potential latent heat flux                                           | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%errmsg                        | error_message                                                                                  | error message for error handling in CCPP                                            | none          |    0 | character   | len=512   | none   | F        |
!! | physics%Interstitial(i)%errflg                        | error_flag                                                                                     | error flag for error handling in CCPP                                               | flag          |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%evap                          | kinematic_surface_upward_latent_heat_flux                                                      | kinematic surface upward latent heat flux                                           | kg kg-1 m s-1 |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%evbs                          | soil_upward_latent_heat_flux                                                                   | soil upward latent heat flux                                                        | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%evcw                          | canopy_upward_latent_heat_flux                                                                 | canopy upward latent heat flux                                                      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faerlw                        |                                                                                                | optical properties for longwave bands 01-16                                         | various       |    4 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faerlw(:,:,:,1)               | aerosol_optical_depth_for_longwave_bands_01-16                                                 | aerosol optical depth for longwave bands 01-16                                      | none          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faerlw(:,:,:,2)               | aerosol_single_scattering_albedo_for_longwave_bands_01-16                                      | aerosol single scattering albedo for longwave bands 01-16                           | frac          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faerlw(:,:,:,3)               | aerosol_asymmetry_parameter_for_longwave_bands_01-16                                           | aerosol asymmetry parameter for longwave bands 01-16                                | none          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faersw                        |                                                                                                | optical properties for shortwave bands 01-16                                        | various       |    4 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faersw(:,:,:,1)               | aerosol_optical_depth_for_shortwave_bands_01-16                                                | aerosol optical depth for shortwave bands 01-16                                     | none          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faersw(:,:,:,2)               | aerosol_single_scattering_albedo_for_shortwave_bands_01-16                                     | aerosol single scattering albedo for shortwave bands 01-16                          | frac          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%faersw(:,:,:,3)               | aerosol_asymmetry_parameter_for_shortwave_bands_01-16                                          | aerosol asymmetry parameter for shortwave bands 01-16                               | none          |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%fh2                           | Monin-Obukhov_similarity_function_for_heat_at_2m                                               | Monin-Obukhov similarity parameter for heat at 2m                                   | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%flag_guess                    | flag_for_guess_run                                                                             | flag for guess run                                                                  | flag          |    1 | logical     |           | none   | F        |
!! | physics%Interstitial(i)%flag_iter                     | flag_for_iteration                                                                             | flag for iteration                                                                  | flag          |    1 | logical     |           | none   | F        |
!! | physics%Interstitial(i)%fm10                          | Monin-Obukhov_similarity_function_for_momentum_at_10m                                          | Monin-Obukhov similarity parameter for momentum at 10m                              | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%frain                         | dynamics_to_physics_timestep_ratio                                                             | ratio of dynamics timestep to physics timestep                                      | none          |    0 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gabsbdlw                      | surface_downwelling_longwave_flux_absorbed_by_ground                                           | total sky surface downward longwave flux absorbed by the ground                     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gamma                         | anisotropy_of_subgrid_orography                                                                | anisotropy of subgrid orography                                                     | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gamq                          | countergradient_mixing_term_for_water_vapor                                                    | countergradient mixing term for water vapor                                         | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gamt                          | countergradient_mixing_term_for_temperature                                                    | countergradient mixing term for temperature                                         | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr                        |                                                                                                | gas volume mixing ratios                                                            | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,1)                 | volume_mixing_ratio_co2                                                                        | volume mixing ratio co2                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,2)                 | volume_mixing_ratio_n2o                                                                        | volume mixing ratio no2                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,3)                 | volume_mixing_ratio_ch4                                                                        | volume mixing ratio ch4                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,4)                 | volume_mixing_ratio_o2                                                                         | volume mixing ratio o2                                                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,5)                 | volume_mixing_ratio_co                                                                         | volume mixing ratio co                                                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,6)                 | volume_mixing_ratio_cfc11                                                                      | volume mixing ratio cfc11                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,7)                 | volume_mixing_ratio_cfc12                                                                      | volume mixing ratio cfc12                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,8)                 | volume_mixing_ratio_cfc22                                                                      | volume mixing ratio cfc22                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,9)                 | volume_mixing_ratio_ccl4                                                                       | volume mixing ratio ccl4                                                            | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gasvmr(:,:,10)                | volume_mixing_ratio_cfc113                                                                     | volume mixing ratio cfc113                                                          | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gflx                          | upward_heat_flux_in_soil                                                                       | soil heat flux                                                                      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gwdcu                         | tendency_of_x_wind_due_to_convective_gravity_wave_drag                                         | zonal wind tendency due to convective gravity wave drag                             | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%gwdcv                         | tendency_of_y_wind_due_to_convective_gravity_wave_drag                                         | meridional wind tendency due to convective gravity wave drag                        | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%hflx                          | kinematic_surface_upward_sensible_heat_flux                                                    | kinematic surface upward sensible heat flux                                         | K m s-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%hprime1                       | standard_deviation_of_subgrid_orography                                                        | standard deviation of subgrid orography                                             | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%idxday                        | daytime_points                                                                                 | daytime points                                                                      | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%im                            | horizontal_loop_extent                                                                         | horizontal loop extent                                                              | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%ipr                           | horizontal_index_of_printed_column                                                             | horizontal index of printed column                                                  | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%islmsk                        | sea_land_ice_mask                                                                              | sea/land/ice mask (=0/1/2)                                                          | flag          |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%iter                          | iteration_number                                                                               | number of iteration                                                                 | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%ix                            | horizontal_dimension                                                                           | horizontal dimension                                                                | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kb                            | vertical_index_difference_between_layer_and_lower_bound                                        | vertical index difference between layer and lower bound                             | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kbot                          | vertical_index_at_cloud_base                                                                   | vertical index at cloud base                                                        | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kcnv                          | flag_deep_convection                                                                           | flag indicating whether convection occurs in column (0 or 1)                        | flag          |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kd                            | vertical_index_difference_between_inout_and_local                                              | vertical index difference between in/out and local                                  | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kinver                        | index_of_highest_temperature_inversion                                                         | index of highest temperature inversion                                              | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kpbl                          | vertical_index_at_top_of_atmosphere_boundary_layer                                             | vertical index at top atmospheric boundary layer                                    | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%kt                            | vertical_index_difference_between_layer_and_upper_bound                                        | vertical index difference between layer and upper bound                             | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%ktop                          | vertical_index_at_cloud_top                                                                    | vertical index at cloud top                                                         | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%latidxprnt                    | latitude_index_in_debug_printouts                                                              | latitude index in debug printouts                                                   | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%levi                          | vertical_interface_dimension                                                                   | vertical interface dimension                                                        | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%levozp                        | vertical_dimension_of_ozone_forcing_data                                                       | number of vertical layers in ozone forcing data                                     | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%lm                            | vertical_layer_dimension_for_radiation                                                         | number of vertical layers for radiation                                             | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%lmk                           | adjusted_vertical_layer_dimension_for_radiation                                                | adjusted number of vertical layers for radiation                                    | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%lmp                           | adjusted_vertical_level_dimension_for_radiation                                                | adjusted number of vertical levels for radiation                                    | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%mbota                         | model_layer_number_at_cloud_base                                                               | vertical indices for low, middle and high cloud bases                               | index         |    2 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%mtopa                         | model_layer_number_at_cloud_top                                                                | vertical indices for low, middle and high cloud tops                                | index         |    2 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%nday                          | daytime_points_dimension                                                                       | daytime points dimension                                                            | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%ntk                           | index_of_TKE                                                                                   | index of TKE in the tracer array                                                    | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%nvdiff                        | number_of_vertical_diffusion_tracers                                                           | number of tracers to diffuse vertically                                             | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%oa4                           | asymmetry_of_subgrid_orography                                                                 | asymmetry of subgrid orography                                                      | none          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%oc                            | convexity_of_subgrid_orography                                                                 | convexity of subgrid orography                                                      | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%olyr                          | ozone_concentration_at_layer_for_radiation                                                     | ozone concentration layer                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%oz_coeff                      | number_of_coefficients_in_ozone_forcing_data                                                   | number of coefficients in ozone forcing data                                        | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%oz_pres                       | natural_log_of_ozone_forcing_data_pressure_levels                                              | natural log of ozone forcing data pressure levels                                   | log(Pa)       |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%plvl                          | air_pressure_at_interface_for_radiation_in_hPa                                                 | air pressure at vertical interface for radiation calculation                        | hPa           |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%plyr                          | air_pressure_at_layer_for_radiation_in_hPa                                                     | air pressure at vertical layer for radiation calculation                            | hPa           |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%qlyr                          | water_vapor_specific_humidity_at_layer_for_radiation                                           | specific humidity layer                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%qss                           | surface_specific_humidity                                                                      | surface air saturation specific humidity                                            | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%raddt                         | time_step_for_radiation                                                                        | radiation time step                                                                 | s             |    0 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%raincd                        | lwe_thickness_of_deep_convective_precipitation_amount                                          | deep convective rainfall amount on physics timestep                                 | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%raincs                        | lwe_thickness_of_shallow_convective_precipitation_amount                                       | shallow convective rainfall amount on physics timestep                              | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rainmcadj                     | lwe_thickness_of_moist_convective_adj_precipitation_amount                                     | adjusted moist convective rainfall amount on physics timestep                       | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rainp                         | tendency_of_rain_water_mixing_ratio_due_to_model_physics                                       | tendency of rain water mixing ratio due to model physics                            | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rainst                        | lwe_thickness_of_stratiform_precipitation_amount                                               | stratiform rainfall amount on physics timestep                                      | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rb                            | bulk_richardson_number_at_lowest_model_level                                                   | bulk Richardson number at the surface                                               | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rhc                           | critical_relative_humidity                                                                     | critical relative humidity                                                          | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rhcbot                        | critical_relative_humidity_at_surface                                                          | critical relative humidity at the surface                                           | frac          |    0 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rhcpbl                        | critical_relative_humidity_at_PBL_top                                                          | critical relative humidity at the PBL top                                           | frac          |    0 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%rhctop                        | critical_relative_humidity_at_top_of_atmosphere                                                | critical relative humidity at the top of atmosphere                                 | frac          |    0 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%runoff                        | surface_runoff_flux                                                                            | surface runoff flux                                                                 | g m-2 s-1     |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%save_qcw                      | cloud_condensed_water_specific_humidity_save                                                   | cloud condensed water specific humidity before entering a physics scheme            | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%save_qv                       | water_vapor_specific_humidity_save                                                             | water vapor specific humidity before entering a physics scheme                      | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%save_t                        | air_temperature_save                                                                           | air temperature before entering a physics scheme                                    | K             |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%save_u                        | x_wind_save                                                                                    | x-wind before entering a physics scheme                                             | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%save_v                        | y_wind_save                                                                                    | y-wind before entering a physics scheme                                             | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sbsno                         | snow_deposition_sublimation_upward_latent_heat_flux                                            | latent heat flux from snow depo/subl                                                | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%scmpsw                        | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes            | W m-2         |    1 | cmpfsw_type |           | none   | F        |
!! | physics%Interstitial(i)%sfcalb                        |                                                                                                |                                                                                     | frac          |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sfcalb(:,1)                   | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                           | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sfcalb(:,2)                   | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sfcalb(:,3)                   | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                            | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sfcalb(:,4)                   | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                                          | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sigma                         | slope_of_subgrid_orography                                                                     | slope of subgrid orography                                                          | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%sigmaf                        | vegetation_area_fraction                                                                       | areal fractional cover of green vegetation                                          | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%skip_macro                    | flag_skip_macro                                                                                | flag to skip cloud macrophysics in Morrison scheme                                  | flag          |    1 | logical     |           | none   | F        |
!! | physics%Interstitial(i)%slopetype                     | surface_slope_classification                                                                   | class of sfc slope                                                                  | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%snowc                         | surface_snow_area_fraction                                                                     | surface snow area fraction                                                          | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%snohf                         | snow_freezing_rain_upward_latent_heat_flux                                                     | latent heat flux due to snow and frz rain                                           | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%snowmt                        | surface_snow_melt                                                                              | snow melt during timestep                                                           | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%soiltype                      | cell_soil_type                                                                                 | soil type at each grid cell                                                         | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%stress                        | surface_wind_stress                                                                            | surface wind stress                                                                 | m2 s-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%theta                         | angle_from_east_of_maximum_subgrid_orographic_variations                                       | angle with_respect to east of maximum subgrid orographic variations                 | degrees       |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tice                          | sea_ice_temperature_for_physics                                                                | sea-ice surface temperature                                                         | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tlvl                          | air_temperature_at_interface_for_radiation                                                     | air temperature at vertical interface for radiation calculation                     | K             |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tlyr                          | air_temperature_at_layer_for_radiation                                                         | air temperature at vertical layer for radiation calculation                         | K             |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tracers_start_index           | start_index_of_other_tracers                                                                   | beginning index of the non-water tracer species                                     | index         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%tracers_total                 | number_of_total_tracers                                                                        | total number of tracers                                                             | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%trans                         | transpiration_flux                                                                             | total plant transpiration rate                                                      | kg m-2 s-1    |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tseal                         | surface_skin_temperature_for_nsst                                                              | ocean surface skin temperature                                                      | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tsfa                          | surface_air_temperature_for_radiation                                                          | lowest model layer air temperature for radiation                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tsfg                          | surface_ground_temperature_for_radiation                                                       | surface ground temperature for radiation                                            | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tsurf                         | surface_skin_temperature_after_iteration                                                       | surface skin temperature after iteration                                            | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%tracers_water                 | number_of_water_tracers                                                                        | number of water-related tracers                                                     | count         |    0 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%ud_mf                         | instantaneous_atmosphere_updraft_convective_mass_flux                                          | (updraft mass flux) * delt                                                          | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%vegtype                       | cell_vegetation_type                                                                           | vegetation type at each grid cell                                                   | index         |    1 | integer     |           | none   | F        |
!! | physics%Interstitial(i)%wind                          | wind_speed_at_lowest_model_layer                                                               | wind speed at lowest model level                                                    | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%work1                         | grid_size_related_coefficient_used_in_scale-sensitive_schemes                                  | grid size related coefficient used in scale-sensitive schemes                       | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%work2                         | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement                       | complement to work1                                                                 | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%work3                         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer                   | Exner function ratio bt midlayer and interface at 1st layer                         | ratio         |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%xcosz                         | instantaneous_cosine_of_zenith_angle                                                           | cosine of zenith angle at current time                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%xmu                           | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                                   | zenith angle temporal adjustment factor for shortwave                               | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Interstitial(i)%zice                          | sea_ice_thickness_for_physics                                                                  | sea-ice thickness                                                                   | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%LTP                             | extra_top_layer                                        | extra top layer for radiation                           | none          |    0 | integer               |           | none   | F        |
!! | physics%Statein(i)%phii        | geopotential_at_interface                              | geopotential at model layer interfaces                 | m2 s-2        |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prsi        | air_pressure_at_interface                              | air pressure at model layer interfaces                 | Pa            |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prsik       | dimensionless_exner_function_at_model_interfaces       | dimensionless Exner function at model layer interfaces | none          |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prsik(:,1)  | dimensionless_exner_function_at_lowest_model_interface | dimensionless Exner function at lowest model interface | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%phil        | geopotential                                           | geopotential at model layer centers                    | m2 s-2        |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prsl        | air_pressure                                           | mean layer pressure                                    | Pa            |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prsl(:,1)   | air_pressure_at_lowest_model_layer                     | mean pressure at lowest model layer                    | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prslk       | dimensionless_exner_function_at_model_layers           | dimensionless Exner function at model layer centers    | none          |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%prslk(:,1)  | dimensionless_exner_function_at_lowest_model_layer     | dimensionless Exner function at lowest model layer     | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%pgr         | surface_air_pressure                                   | surface pressure                                       | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%ugrs        | x_wind                                                 | zonal wind                                             | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%ugrs(:,1)   | x_wind_at_lowest_model_layer                           | zonal wind at lowest model layer                       | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%vgrs        | y_wind                                                 | meridional wind                                        | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%vgrs(:,1)   | y_wind_at_lowest_model_layer                           | meridional wind at lowest model layer                  | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%vvl         | omega                                                  | layer mean vertical velocity                           | Pa s-1        |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%tgrs        | air_temperature                                        | model layer mean temperature                           | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%tgrs(:,1)   | air_temperature_at_lowest_model_layer                  | mean temperature at lowest model layer                 | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%qgrs        | tracer_concentration                                   | model layer mean tracer concentration                  | kg kg-1       |    3 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%qgrs(:,:,scm_state%water_vapor_index) | water_vapor_specific_humidity                          | water vapor specific humidity                          | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Statein(i)%qgrs(:,1,scm_state%water_vapor_index) | specific_humidity_at_lowest_model_layer                | specific humidity at lowest model layer                | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gu0                       | x_wind_updated_by_physics                                  | zonal wind updated by physics                              | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gv0                       | y_wind_updated_by_physics                                  | meridional wind updated by physics                         | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gt0                       | air_temperature_updated_by_physics                         | temperature updated by physics                             | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gq0                       | tracer_concentration_updated_by_physics                    | tracer concentration updated by physics                    | kg kg-1       |    3 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gq0(:,:,scm_state%water_vapor_index)                | water_vapor_specific_humidity_updated_by_physics           | water vapor specific humidity updated by physics           | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gq0(:,:,scm_state%cloud_water_index) | cloud_condensed_water_specific_humidity_updated_by_physics | cloud condensed water specific humidity updated by physics | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Stateout(i)%gq0(:,:,scm_state%ozone_index) | ozone_concentration_updated_by_physics                     | ozone concentration updated by physics                     | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Model(i)%me                       | mpi_rank                                                                      | current MPI-rank                                        | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%master                   |                                                                               | master MPI-rank                                         | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nlunit                   |                                                                               | fortran unit number for file opens                      | none          |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%fn_nml                   |                                                                               | namelist filename                                       | none          |    0 | charater  |           | none   | F        |
!! | physics%Model(i)%fhzero                   |                                                                               | seconds between clearing of diagnostic buckets          | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%ldiag3d                  | flag_diagnostics_3D                                                           | flag for 3d diagnostic fields                           | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lssav                    | flag_diagnostics                                                              | logical flag for storing diagnostics                    | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%fhcyc                    |                                                                               | frequency for surface data cycling (secs)               | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%lgocart                  |                                                                               | flag for 3d diagnostic fields for gocart 1              | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%fhgoc3d                  |                                                                               | seconds between calls to gocart                         | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%thermodyn_id             |                                                                               | valid for GFS only for get_prs/phi                      | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%sfcpress_id              |                                                                               | valid for GFS only for get_prs/phi                      | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%gen_coord_hybrid         |                                                                               | flag for Henry's gen coord                              | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%isc                      |                                                                               | starting i-index for this MPI-domain                    | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%jsc                      |                                                                               | starting j-index for this MPI-domain                    | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nx                       |                                                                               | number of points in i-dir for this MPI rank             | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ny                       |                                                                               | number of points in j-dir for this MPI rank             | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%levs                     | vertical_dimension                                                            | number of vertical levels                               | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%cnx                      |                                                                               | number of points in i-dir for this cubed-sphere face    | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%cny                      |                                                                               | number of points in j-dir for this cubed-sphere face    | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%lonr                     | number_of_equatorial_longitude_points                                         | number of global points in x-dir (i) along the equator  | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%latr                     |                                                                               | number of global points in y-dir (j) along the meridian | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%cplflx                   |                                                                               | flag controlloing cplflx collection (default off)       | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%cplwav                   |                                                                               | flag controlloing cplwav collection (default off)       | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lsidea                   | flag_idealized_physics                                                        | flag for idealized physics                              | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%dtp                      | time_step_for_physics                                                         | physics timestep                                        | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%dtf                      | time_step_for_dynamics                                                        | dynamics timestep                                       | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%nscyc                    |                                                                               | trigger for surface data cycling                        |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nszero                   |                                                                               | trigger for zeroing diagnostic buckets                  |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%idat                     |                                                                               | initialization date and time                            |               |    1 | integer   |           | none   | F        |
!! | physics%Model(i)%idate                    |                                                                               | initial date with different size and ordering           |               |    1 | integer   |           | none   | F        |
!! | physics%Model(i)%fhswr                    |                                                                               | frequency for shortwave radiation                       | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%fhlwr                    |                                                                               | frequency for longwave radiation                        | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%nsswr                    |                                                                               | integer trigger for shortwave radiation                 |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nslwr                    |                                                                               | integer trigger for longwave  radiation                 |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%levr                     | number_of_vertical_layers_for_radiation_calculations                          | number of vertical levels for radiation calculations    | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nfxr                     |                                                                               | second dimension for fluxr diagnostic variable          |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%aero_in                  |                                                                               | aerosol flag for gbphys                                 |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lmfshal                  |                                                                               | parameter for radiation                                 |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lmfdeep2                 |                                                                               | parameter for radiation                                 |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%nrcm                     | array_dimension_of_random_number                                              | second dimension of random number stream for RAS        | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%iflip                    |                                                                               | iflip - is not the same as flipv                        |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%isol                     | flag_for_solar_constant                                                       | use prescribed solar constant                           | flag          |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ico2                     |                                                                               | prescribed global mean value (old opernl)               |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ialb                     |                                                                               | flag for using climatology alb, based on sfc type       |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%iems                     |                                                                               | use fixed value of 1.0                                  |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%iaer                     |                                                                               | default aerosol effect in sw only                       |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%iovr_sw                  |                                                                               | sw: max-random overlap clouds                           |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%iovr_lw                  |                                                                               | lw: max-random overlap clouds                           |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ictm                     | flag_for_initial_time-date_control                                            | flag for initial conditions and forcing                 | flag          |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%isubc_sw                 |                                                                               | flag for sw clouds without sub-grid approximation       |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%isubc_lw                 |                                                                               | flag for lw clouds without sub-grid approximation       |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%crick_proof              |                                                                               | CRICK-Proof cloud water                                 |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%ccnorm                   |                                                                               | Cloud condensate normalized by cloud cover              |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%norad_precip             |                                                                               | radiation precip flag for Ferrier/Moorthi               |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lwhtr                    |                                                                               | flag to output lw heating rate (Radtend%lwhc)           |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%swhtr                    |                                                                               | flag to output sw heating rate (Radtend%swhc)           |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%ncld                     | number_of_hydrometeors                                                        | choice of cloud scheme / number of hydrometeors         | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%zhao_mic                 |                                                                               | flag for Zhao-Carr microphysics                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%psautco                  | coefficient_from_cloud_ice_to_snow                                            | auto conversion coeff from ice to snow                  | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%prautco                  | coefficient_from_cloud_water_to_rain                                          | auto conversion coeff from cloud to rain                | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%evpco                    | coefficient_for_evaporation_of_rainfall                                       | coeff for evaporation of largescale rain                | none          |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%wminco                   | cloud_condensed_water_conversion_threshold                                    | water and ice minimum threshold for Zhao                | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%fprcp                    |                                                                               | no prognostic rain and snow (MG)                        |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%mg_dcs                   |                                                                               | Morrison-Gettleman microphysics parameters              |               |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%mg_qcvar                 |                                                                               |                                                         |               |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%mg_ts_auto_ice           |                                                                               | ice auto conversion time scale                          |               |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%lsm                      | flag_for_land_surface_scheme                                                  | flag for land surface model lsm=1 for noah lsm          | flag          |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%lsoil                    | soil_vertical_dimension                                                       | number of soil layers                                   | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ivegsrc                  | vegetation_type                                                               | land use classification                                 | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%isot                     | soil_type                                                                     | soil type classification                                | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%mom4ice                  | flag_for_mom4_coupling                                                        | flag controls mom4 sea ice                              | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%use_ufo                  |                                                                               | flag for gcycle surface option                          |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%ras                      |                                                                               | flag for ras convection scheme                          |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%flipv                    |                                                                               | flag for vertical direction flip (ras)                  |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%trans_trac               |                                                                               | flag for convective transport of tracers                |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%old_monin                |                                                                               | flag for diff monin schemes                             |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%cnvgwd                   |                                                                               | flag for conv gravity wave drag                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%mstrat                   |                                                                               | flag for moorthi approach for stratus                   |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%moist_adj                |                                                                               | flag for moist convective adjustment                    |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%cscnv                    |                                                                               | flag for Chikira-Sugiyama convection                    |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%cal_pre                  | flag_for_precipitation_type_algorithm                                         | flag controls precip type algorithm                     | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%do_aw                    |                                                                               | AW scale-aware option in cs convection                  |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%do_shoc                  |                                                                               | flag for SHOC                                           |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%shocaftcnv               |                                                                               | flag for SHOC                                           |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%shoc_cld                 |                                                                               | flag for clouds                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%uni_cld                  |                                                                               | flag for clouds in grrad                                |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%h2o_phys                 |                                                                               | flag for stratosphere h2o                               |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%pdfcld                   |                                                                               | flag for pdfcld                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%shcnvcw                  |                                                                               | flag for shallow convective cloud                       |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%redrag                   | flag_for_reduced_drag_coefficient_over_sea                                    | flag for reduced drag coeff. over sea                   | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%hybedmf                  |                                                                               | flag for hybrid edmf pbl scheme                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%dspheat                  | flag_TKE_dissipation_heating                                                  | flag for tke dissipative heating                        | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%cnvcld                   |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%random_clds              |                                                                               | flag controls whether clouds are random                 |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%shal_cnv                 |                                                                               | flag for calling shallow convection                     |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%imfshalcnv               |                                                                               | flag for mass-flux shallow convection scheme            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%imfdeepcnv               |                                                                               | flag for mass-flux deep convection scheme               |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nmtvr                    | number_of_statistical_measures_of_subgrid_orography                           | number of topographic variables in GWD                  | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%jcap                     |                                                                               | number of spectral wave trancation                      |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%cs_parm                  |                                                                               | tunable parameters for Chikira-Sugiyama convection      |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%flgmin                   |                                                                               | [in] ice fraction bounds                                |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%cgwf                     | multiplication_factors_for_convective_gravity_wave_drag                       | multiplication factor for convective GWD                | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%ccwf                     |                                                                               | multiplication factor for critical cloud workfunction   | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%cdmbgwd                  | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag | multiplication factors for cdmb and gwd                 | none          |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%sup                      |                                                                               | supersaturation in pdf cloud when t is very low         |               |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%ctei_rm                  |                                                                               | critical cloud top entrainment instability criteria     |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%crtrh                    |                                                                               | critical relative humidity at SFC, PBL top and TOA      |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%dlqf                     |                                                                               | factor for cloud condensate detrainment from cloud edges|               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%seed0                    |                                                                               | random seed for radiation                               |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%prslrd0                  | pressure_cutoff_for_rayleigh_damping                                          | pressure level from which Rayleigh Damping is applied   | Pa            |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%ral_ts                   | time_scale_for_rayleigh_damping                                               | time scale for Rayleigh damping in days                 | d             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%nst_anl                  |                                                                               | flag for NSSTM analysis in gcycle/sfcsub                |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lsea                     |                                                                               |                                                         |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%xkzm_m                   | atmosphere_momentum_diffusivity_background                                    | background vertical diffusion for momentum              | m2 s-1        |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%xkzm_h                   | atmosphere_heat_diffusivity_background                                        | background vertical diffusion for heat q                | m2 s-1        |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%xkzm_s                   | diffusivity_background_sigma_level                                            | sigma threshold for background mom. diffusion           | none          |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%nstf_name                |                                                                               |                                                         |               |    1 | integer   |           | none   | F        |
!! | physics%Model(i)%nstf_name(1)             | flag_for_nsstm_run                                                            | NSSTM flag: off/uncoupled/coupled=0/1/2                 | flag          |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nstf_name(4)             | vertical_temperature_average_range_lower_bound                                | zsea1 in mm                                             | mm            |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nstf_name(5)             | vertical_temperature_average_range_upper_bound                                | zsea2 in mm                                             | mm            |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%do_sppt                  |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%do_shum                  |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%do_skeb                  |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%do_vc                    |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%sppt                     |                                                                               | stochastic physics tendency amplitude                   |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%shum                     |                                                                               | stochastic boundary layer spf hum amp                   |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%skeb                     |                                                                               | stochastic KE backscatter amplitude                     |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%vcamp                    |                                                                               | stochastic vorticity confinment amp                     |               |    1 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%vc                       |                                                                               | deterministic vorticity confinement parameter           |               |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%tracer_names             |                                                                               | array of initialized tracers from dynamic core          |               |    1 | character |           | none   | F        |
!! | physics%Model(i)%ntrac                    |                                                                               | number of tracers                                       |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntoz                     |                                                                               | tracer index for ozone mixing ratio                     |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntcw                     | index_for_liquid_cloud_condensate                                             | tracer index for cloud condensate (or liquid water)     | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntiw                     |                                                                               | tracer index for  ice water                             |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntrw                     |                                                                               | tracer index for rain water                             |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntsw                     |                                                                               | tracer index for snow water                             |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntgl                     |                                                                               | tracer index for graupel                                |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntlnc                    |                                                                               | tracer index for liquid number concentration            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntinc                    |                                                                               | tracer index for ice    number concentration            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntrnc                    |                                                                               | tracer index for rain   number concentration            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntsnc                    |                                                                               | tracer index for snow   number concentration            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntke                     |                                                                               | tracer index for kinetic energy                         |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nto                      |                                                                               | tracer index for oxygen ion                             |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nto2                     |                                                                               | tracer index for oxygen                                 |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntot2d                   |                                                                               | total number of variables for phyf2d                    |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ntot3d                   |                                                                               | total number of variables for phyf3d                    |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%num_p2d                  |                                                                               | number of 2D arrays needed for microphysics             |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%num_p3d                  | array_dimension_of_microphysics                                               | number of 3D arrays needed for microphysics             | count         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nshoc_2d                 |                                                                               | number of 2d fields for SHOC                            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nshoc_3d                 |                                                                               | number of 3d fields for SHOC                            |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%ncnvcld3d                |                                                                               | number of convective 3d clouds fields                   |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%npdf3d                   |                                                                               | number of 3d arrays associated with pdf based clouds/mp |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%nctp                     |                                                                               | number of cloud types in Chikira-Sugiyama scheme        |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%debug                    |                                                                               | debug flag                                              |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%pre_rad                  |                                                                               | flag for testing purpose                                |               |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%ipt                      |                                                                               | index for diagnostic printout point                     |               |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%lprnt                    | flag_print                                                                    | control flag for diagnostic print out                   | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lsswr                    | flag_to_calc_sw                                                               | logical flags for sw radiation calls                    | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%lslwr                    | flag_to_calc_lw                                                               | logical flags for lw radiation calls                    | flag          |    0 | logical   |           | none   | F        |
!! | physics%Model(i)%solhr                    | forecast_hour                                                                 | hour time after 00z at the t-step                       | h             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%solcon                   | solar_constant                                                                | solar constant (sun-earth distant adjusted)             | W m-2         |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%slag                     | equation_of_time                                                              | equation of time (radian)                               | radians       |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%sdec                     | sine_of_solar_declination_angle                                               | sin of the solar declination angle                      | none          |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%cdec                     | cosine_of_solar_declination_angle                                             | cos of the solar declination angle                      | none          |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%clstp                    | convective_cloud_switch                                                       | index used by cnvc90 (for convective clouds)            | none          |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%phour                    |                                                                               | previous forecast time                                  | h             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%fhour                    | forecast_time                                                                 | curent forecast time                                    | h             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%zhour                    |                                                                               | previous hour diagnostic buckets emptied                | h             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%kdt                      | index_of_time_step                                                            | current forecast iteration                              | index         |    0 | integer   |           | none   | F        |
!! | physics%Model(i)%jdat                     |                                                                               | current forecast date and time                          |               |    1 | integer   |           | none   | F        |
!! | physics%Model(i)%sec                      | seconds_elapsed_since_model_initialization                                    | seconds elapsed since model initialization              | s             |    0 | real      | kind_phys | none   | F        |
!! | physics%Model(i)%blksz                    | horizontal_block_size                                                         | for explicit data blocking: block sizes of all blocks   | count         |    1 | integer   |           | none   | F        |
!! | physics%Tbd(i)%icsdsw                         | seed_random_numbers_sw                                                                         | random seeds for sub-column cloud generators sw         | none          |    1 | integer |           | none   | F        |
!! | physics%Tbd(i)%icsdlw                         | seed_random_numbers_lw                                                                         | random seeds for sub-column cloud generators lw         | none          |    1 | integer |           | none   | F        |
!! | physics%Tbd(i)%ozpl                           | ozone_forcing                                                                                  | ozone forcing data                                      | various       |    3 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%h2opl                          |                                                                                                | water forcing data                                      |               |    3 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%rann                           | random_number_array                                                                            | random number array (0-1)                               | none          |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%acv                            | accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90                            | accumulated convective rainfall amount for cnvc90 only  | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%acvb                           | smallest_cloud_base_vertical_index_encountered_thus_far                                        | smallest cloud base vertical index encountered thus far | index         |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%acvt                           | largest_cloud_top_vertical_index_encountered_thus_far                                          | largest cloud top vertical index encountered thus far   | index         |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%dtdtr                          |                                                                                                | temp. change due to radiative heating per time step     | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%dtotprcp                       |                                                                                                | change in totprcp  (diag_type)                          |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%dcnvprcp                       |                                                                                                | change in cnvprcp  (diag_type)                          |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%drain_cpl                      |                                                                                                | change in rain_cpl (coupling_type)                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%dsnow_cpl                      |                                                                                                | change in show_cpl (coupling_type)                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_fctd                       |                                                                                                | for CS convection                                       |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f2d                        |                                                                                                | 2d arrays saved for restart                             |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f2d(:,1)                   | surface_air_pressure_two_time_steps_back                                                       | surface air pressure two time steps back                | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f2d(:,2)                   | surface_air_pressure_at_previous_time_step                                                     | surface air pressure at previous time step              | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f2d(:,physics%Model(i)%num_p2d) | surface_wind_enhancement_due_to_convection                                                     | surface wind enhancement due to convection              | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f3d                        |                                                                                                | 3d arrays saved for restart                             |               |    3 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f3d(:,:,1)                 | air_temperature_two_time_steps_back                                                            | air temperature two time steps back                     | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f3d(:,:,2)                 | water_vapor_specific_humidity_two_time_steps_back                                              | water vapor specific humidity two time steps back       | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f3d(:,:,3)                 | air_temperature_at_previous_time_step                                                          | air temperature at previous time step                   | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%phy_f3d(:,:,4)                 | water_vapor_specific_humidity_at_previous_time_step                                            | water vapor specific humidity at previous time step     | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%blkno                          | block_number                                                                                   | for explicit data blocking: block number of this block  | index         |    0 | integer |           | none   | F        |
!! | physics%Tbd(i)%htlwc                          | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                     | total sky heating rate due to longwave radiation        | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%htlw0                          | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step  | clear sky heating rate due to longwave radiation        | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%htswc                          | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | total sky heating rate due to shortwave radiation       | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Tbd(i)%htsw0                          | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky heating rates due to shortwave radiation      | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | physics%Radtend(i)%sfcfsw               | sw_fluxes_sfc                                                                                 | sw radiation fluxes at sfc                              | W m-2         |    1 | sfcfsw_type |           | none   | F        |
!! | physics%Radtend(i)%sfcflw               | lw_fluxes_sfc                                                                                 | lw radiation fluxes at sfc                              | W m-2         |    1 | sfcflw_type |           | none   | F        |
!! | physics%Radtend(i)%htrsw                | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep                    | total sky sw heating rate                               | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%htrlw                | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep                     | total sky lw heating rate                               | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%sfalb                | surface_diffused_shortwave_albedo                                                             | mean surface diffused sw albedo                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%coszen               | cosine_of_zenith_angle                                                                        | mean cos of zenith angle over rad call period           | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%tsflw                | surface_midlayer_air_temperature_in_longwave_radiation                                        | surface air temp during lw calculation                  | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%semis                | surface_longwave_emissivity                                                                   | surface lw emissivity in fraction                       | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%coszdg               |                                                                                               | daytime mean cosz over rad call period                  | none          |    1 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%swhc                 | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_timestep | clear sky sw heating rates                              | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%lwhc                 | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_timestep  | clear sky lw heating rates                              | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Radtend(i)%lwhd                 |                                                                                               | idea sky lw heating rates                               | K s-1         |    3 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%fluxr                |                                                                         | accumulated 2-d fields, opt. includes aerosols                  |               |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%topfsw               | sw_fluxes_top_atmosphere                                                | sw radiation fluxes at toa                                      | W m-2         |    1 | topfsw_type |           | none   | F        |
!! | physics%Diag(i)%topflw               | lw_fluxes_top_atmosphere                                                | lw radiation fluxes at top                                      | W m-2         |    1 | topflw_type |           | none   | F        |
!! | physics%Diag(i)%srunoff              | surface_runoff                                                          | surface water runoff (from lsm)                                 | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%evbsa                |                                                                         | noah lsm diagnostics                                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%evcwa                |                                                                         | noah lsm diagnostics                                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%snohfa               |                                                                         | noah lsm diagnostics                                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%transa               |                                                                         | noah lsm diagnostics                                            | kg m-2 s-1    |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%sbsnoa               |                                                                         | noah lsm diagnostics                                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%snowca               |                                                                         | noah lsm diagnostics                                            |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%soilm                | soil_moisture_content                                                   | soil moisture                                                   | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%tmpmin               |                                                                         | min temperature at 2m height                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%tmpmax               |                                                                         | max temperature at 2m height                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dusfc                |                                                                         | u component of surface stress                                   |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dvsfc                |                                                                         | v component of surface stress                                   |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dtsfc                |                                                                         | sensible heat flux                                              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dqsfc                |                                                                         | latent heat flux                                                | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%totprcp              | accumulated_lwe_thickness_of_precipitation_amount                       | accumulated total precipitation                                 | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%gflux                |                                                                         | groud conductive heat flux                                      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dlwsfc               |                                                                         | time accumulated sfc dn lw flux                                 | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%ulwsfc               |                                                                         | time accumulated sfc up lw flux                                 | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%suntim               |                                                                         | sunshine duration time                                          | s             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%runoff               | total_runoff                                                            | total water runoff                                              | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%ep                   |                                                                         | potential evaporation                                           | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%cldwrk               |                                                                         | cloud workfunction (valid only with sas)                        |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dugwd                | time_integral_of_x_stress_due_to_gravity_wave_drag                      | vertically integrated u change by OGWD                          | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dvgwd                | time_integral_of_y_stress_due_to_gravity_wave_drag                      | vertically integrated v change by OGWD                          | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%psmean               |                                                                         | surface pressure                                                | kPa           |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%cnvprcp              |                                                                         | accumulated convective precipitation                            | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%spfhmin              |                                                                         | minimum specific humidity                                       | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%spfhmax              |                                                                         | maximum specific humidity                                       | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%rain                 | lwe_thickness_of_precipitation_amount_on_dynamics_timestep              | total rain at this time step                                    | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%rainc                | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep   | convective rain at this time step                               | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%ice                  |                                                                         | ice fall at this time step                                      |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%snow                 |                                                                         | snow fall at this time step                                     |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%graupel              |                                                                         | graupel fall at this time step                                  |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%totice               |                                                                         | accumulated ice precipitation                                   | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%totsnw               |                                                                         | accumulated snow precipitation                                  | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%totgrp               |                                                                         | accumulated graupel precipitation                               | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%u10m                 | x_wind_at_10m                                                           | 10 meter u wind speed                                           | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%v10m                 | y_wind_at_10m                                                           | 10 meter v wind speed                                           | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%zlvl                 | height_above_mean_sea_level_at_lowest_model_layer                       | layer 1 height                                                  | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%psurf                |                                                                         | surface pressure                                                | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%hpbl                 | atmosphere_boundary_layer_thickness                                     | pbl height                                                      | m             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%pwat                 | column_precipitable_water                                               | precipitable water                                              | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%t1                   |                                                                         | layer 1 temperature                                             | K             |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%q1                   |                                                                         | layer 1 specific humidity                                       | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%u1                   |                                                                         | layer 1 zonal wind                                              | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%v1                   |                                                                         | layer 1 merdional wind                                          | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%chh                  | surface_drag_mass_flux_for_heat_and_moisture_in_air                     | thermal exchange coefficient                                    | kg m-2 s-1    |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%cmm                  | surface_drag_wind_speed_for_momentum_in_air                             | momentum exchange coefficient                                   | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dlwsfci              |                                                                         | instantaneous sfc dnwd lw flux                                  | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%ulwsfci              |                                                                         | instantaneous sfc upwd lw flux                                  | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dswsfci              |                                                                         | instantaneous sfc dnwd sw flux                                  | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%uswsfci              |                                                                         | instantaneous sfc upwd sw flux                                  | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dusfci               |                                                                         | instantaneous u component of surface stress                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dvsfci               |                                                                         | instantaneous v component of surface stress                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dtsfci               |                                                                         | instantaneous sfc sensible heat flux                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dqsfci               |                                                                         | instantaneous sfc latent heat flux                              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%gfluxi               |                                                                         | instantaneous sfc ground heat flux                              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%epi                  |                                                                         | instantaneous sfc potential evaporation                         |               |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%smcwlt2              | volume_fraction_of_condensed_water_in_soil_at_wilting_point             | wilting point (volumetric)                                      | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%smcref2              | threshold_volume_fraction_of_condensed_water_in_soil                    | soil moisture threshold (volumetric)                            | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%wet1                 | normalized_soil_wetness                                                 | normalized soil wetness                                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%sr                   | ratio_of_snowfall_to_rainfall                                           | snow ratio: ratio of snow to total precipitation                | frac          |    1 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%du3dt                |                                                                         | u momentum change due to physics                                |               |    3 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%du3dt(:,:,2)         | cumulative_change_in_x_wind_due_to_surface_processes                    | cumulative change in x wind due to surface processes            | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%du3dt(:,:,4)         | cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag         | cumulative change in x wind due to convective gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dv3dt                |                                                                         | v momentum change due to physics                                |               |    3 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dv3dt(:,:,2)         | cumulative_change_in_y_wind_due_to_surface_processes                    | cumulative change in y wind due to surface processes            | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dv3dt(:,:,4)         | cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag         | cumulative change in y wind due to convective gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dt3dt                |                                                                         | temperature change due to physics                               |               |    3 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dt3dt(:,:,2)         | cumulative_change_in_temperature_due_to_surface_processes               | cumulative change in temperature due to surface processes       | K             |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dt3dt(:,:,6)         | large_scale_condensate_heating_rate_at_model_layers                     | large scale condensate heating rate at model layers             | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dq3dt                |                                                                         | moisture change due to physics                                  |               |    3 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dq3dt(:,:,4)         | large_scale_condensate_moistening_rate_at_model_layers                  | large scale condensate moistening rate at model layers          | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%upd_mf               |                                                                         | instantaneous convective updraft mass flux                      |               |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%dwn_mf               |                                                                         | instantaneous convective downdraft mass flux                    |               |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%det_mf               |                                                                         | instantaneous convective detrainment mass flux                  |               |    2 | real        | kind_phys | none   | F        |
!! | physics%Diag(i)%cldcov               |                                                                         | instantaneous 3D cloud fraction                                 |               |    2 | real        | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%slmsk     | sea_land_ice_mask_real                                                 | landmask: sea/land/ice=0/1/2                         | flag          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%tsfc      | surface_skin_temperature                                               | ocean surface skin temperature                       | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%tisfc     | sea_ice_temperature                                                    | sea uce surface skin temperature                     | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%snowd     | surface_snow_thickness_water_equivalent                                | water equivalent snow depth over land                | mm            |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%zorl      | surface_roughness_length                                               | surface roughness length                             | cm            |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%fice      | sea_ice_concentration                                                  | ice fraction over open water                         | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%hprim     |                                                                        | topographic standard deviation                       | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%hprime    | statistical_measures_of_subgrid_orography                              | orographic metrics                                   | various       |    2 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%sncovr    | surface_snow_area_fraction_for_diagnostics                             | surface snow area fraction                           | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%snoalb    | upper_bound_on_max_albedo_over_deep_snow                               | maximum snow albedo                                  | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%alvsf     |                                                                        | mean vis albedo with strong cosz dependency          | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%alnsf     |                                                                        | mean nir albedo with strong cosz dependency          | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%alvwf     |                                                                        | mean vis albedo with weak cosz dependency            | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%alnwf     |                                                                        | mean nir albedo with weak cosz dependency            | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%facsf     |                                                                        | fractional coverage with strong cosz dependency      | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%facwf     |                                                                        | fractional coverage with weak cosz dependency        | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%slope     |                                                                        | sfc slope type for lsm                               |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%shdmin    | minimum_vegetation_area_fraction                                       | min fractional coverage of green vegetation          | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%shdmax    | maximum_vegetation_area_fraction                                       | max fractional coverage of green vegetation          | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%tg3       | deep_soil_temperature                                                  | deep soil temperature                                | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%vfrac     |                                                                        | vegetation fraction for lsm                          | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%vtype     |                                                                        | vegetation type for lsm                              | index         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%stype     |                                                                        | soil type                                            | index         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%uustar    | surface_friction_velocity                                              | boundary layer parameter                             | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%oro       | orography                                                              | orography                                            | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%oro_uf    | orography_unfiltered                                                   | unfiltered orography                                 | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%hice      | sea_ice_thickness                                                      | sea ice thickness                                    | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%weasd     | water_equivalent_accumulated_snow_depth                                | water equiv of acc snow depth over land and sea ice  | mm            |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%canopy    | canopy_water_amount                                                    | canopy water amount                                  | kg m-2        |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%ffmm      | Monin-Obukhov_similarity_function_for_momentum                         | Monin-Obukhov similarity function for momentum       | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%ffhh      | Monin-Obukhov_similarity_function_for_heat                             | Monin-Obukhov similarity function for heat           | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%f10m      | ratio_of_wind_at_lowest_model_layer_and_wind_at_10m                    | ratio of sigma level 1 wind and 10m wind             | ratio         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%tprcp     | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep | total precipitation amount in each time step         | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%srflag    | flag_for_precipitation_type                                            | snow/rain flag for precipitation                     | flag          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%slc       | volume_fraction_of_unfrozen_soil_moisture                              | liquid soil moisture                                 | frac          |    2 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%smc       | volume_fraction_of_soil_moisture                                       | total soil moisture                                  | frac          |    2 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%stc       | soil_temperature                                                       | soil temperature                                     | K             |    2 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%t2m       | temperature_at_2m                                                      | 2 meter temperature                                  | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%q2m       | specific_humidity_at_2m                                                | 2 meter specific humidity                            | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%tref      | sea_surface_reference_temperature                                      | sea surface reference temperature                    | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%z_c       | sub-layer_cooling_thickness                                            | sub-layer cooling thickness                          | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%c_0       | coefficient_c_0                                                        | coefficient 1 to calculate d(Tz)/d(Ts)               | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%c_d       | coefficient_c_d                                                        | coefficient 2 to calculate d(Tz)/d(Ts)               | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%w_0       | coefficient_w_0                                                        | coefficient 3 to calculate d(Tz)/d(Ts)               | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%w_d       | coefficient_w_d                                                        | coefficient 4 to calculate d(Tz)/d(Ts)               | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xt        | diurnal_thermocline_layer_heat_content                                 | heat content in diurnal thermocline layer            | K m           |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xs        | sea_water_salinity                                                     | salinity  content in diurnal thermocline layer       | ppt m         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xu        | diurnal_thermocline_layer_x_current                                    | u-current content in diurnal thermocline layer       | m2 s-1        |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xv        | diurnal_thermocline_layer_y_current                                    | v-current content in diurnal thermocline layer       | m2 s-1        |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xz        | diurnal_thermocline_layer_thickness                                    | diurnal thermocline layer thickness                  | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%zm        | ocean_mixed_layer_thickness                                            | mixed layer thickness                                | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xtts      | sensitivity_of_dtl_heat_content_to_surface_temperature                 | d(xt)/d(ts)                                          | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%xzts      | sensitivity_of_dtl_thickness_to_surface_temperature                    | d(xz)/d(ts)                                          | m K-1         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%d_conv    | free_convection_layer_thickness                                        | thickness of free convection layer (FCL)             | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%ifd       | index_of_dtlm_start                                                    | index to start dtlm run or not                       | index         |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%dt_cool   | sub-layer_cooling_amount                                               | sub-layer cooling amount                             | K             |    1 | real    | kind_phys | none   | F        |
!! | physics%Sfcprop(i)%qrain     | sensible_heat_flux_due_to_rainfall                                     | sensible heat flux due to rainfall                   | W             |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%xlon         | longitude                                              | grid longitude in radians                          | radians       |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%xlat         | latitude                                               | grid latitude in radians                           | radians       |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%xlat_d       |                                                        | grid latitude in degrees                           | degrees       |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%sinlat       | sine_of_latitude                                       | sine of the grid latitude                          | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%coslat       | cosine_of_latitude                                     | cosine of the grid latitude                        | none          |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%area         | cell_area                                              | area of the grid cell                              | m2            |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%dx           | cell_size                                              | size of the grid cell                              | m             |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%ddy_o3       |                                                        | interpolation weight for ozone                     | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%jindx1_o3    |                                                        | interpolation low index for ozone                  | index         |    1 | integer |           | none   | F        |
!! | physics%Grid(i)%jindx2_o3    |                                                        | interpolation high index for ozone                 | index         |    1 | integer |           | none   | F        |
!! | physics%Grid(i)%ddy_h        |                                                        | interpolation weight for h2o                       | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Grid(i)%jindx1_h     |                                                        | interpolation low index for h2o                    | index         |    1 | integer |           | none   | F        |
!! | physics%Grid(i)%jindx2_h     |                                                        | interpolation high index for h2o                   | index         |    1 | integer |           | none   | F        |
!! | physics%Coupling(i)%nirbmdi        | surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step            | sfc nir beam sw downward flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nirdfdi        | surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step           | sfc nir diff sw downward flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%visbmdi        | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step  | sfc uv+vis beam sw downward flux                     | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%visdfdi        | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step | sfc uv+vis diff sw downward flux                     | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nirbmui        | surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step              | sfc nir beam sw upward flux                          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nirdfui        | surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step             | sfc nir diff sw upward flux                          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%visbmui        | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step    | sfc uv+vis beam sw upward flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%visdfui        | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step   | sfc uv+vis diff sw upward flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%sfcdsw         | surface_downwelling_shortwave_flux_on_radiation_time_step                                 | total sky sfc downward sw flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%sfcnsw         | surface_net_downwelling_shortwave_flux_on_radiation_time_step                             | total sky sfc netsw flx into ground                  | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%sfcdlw         | surface_downwelling_longwave_flux_on_radiation_time_step                                  | total sky sfc downward lw flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dusfcin_cpl    |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvsfcin_cpl    |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dtsfcin_cpl    |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dqsfcin_cpl    |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%ulwsfcin_cpl   |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%slimskin_cpl   |                                                                                           |                                                      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%rain_cpl       |                                                                                           | total rain precipitation                             |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%snow_cpl       |                                                                                           | total snow precipitation                             |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dusfc_cpl      |                                                                                           | sfc u momentum flux                                  |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvsfc_cpl      |                                                                                           | sfc v momentum flux                                  |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dtsfc_cpl      |                                                                                           | sfc sensible heat flux                               |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dqsfc_cpl      |                                                                                           | sfc latent heat flux                                 |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dlwsfc_cpl     |                                                                                           | sfc downward lw flux                                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dswsfc_cpl     |                                                                                           | sfc downward sw flux                                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dnirbm_cpl     |                                                                                           | sfc nir beam downward sw flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dnirdf_cpl     |                                                                                           | sfc nir diff downward sw flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvisbm_cpl     |                                                                                           | sfc uv+vis beam dnwd sw flux                         | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvisdf_cpl     |                                                                                           | sfc uv+vis diff dnwd sw flux                         | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nlwsfc_cpl     |                                                                                           | net downward lw flux                                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nswsfc_cpl     |                                                                                           | net downward sw flux                                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nnirbm_cpl     |                                                                                           | net nir beam downward sw flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nnirdf_cpl     |                                                                                           | net nir diff downward sw flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nvisbm_cpl     |                                                                                           | net uv+vis beam downward sw rad flux                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nvisdf_cpl     |                                                                                           | net uv+vis diff downward sw rad flux                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dusfci_cpl     |                                                                                           | instantaneous sfc u momentum flux                    |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvsfci_cpl     |                                                                                           | instantaneous sfc v momentum flux                    |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dtsfci_cpl     |                                                                                           | instantaneous sfc sensible heat flux                 |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dqsfci_cpl     |                                                                                           | instantaneous sfc latent heat flux                   |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dlwsfci_cpl    |                                                                                           | instantaneous sfc downward lw flux                   |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dswsfci_cpl    |                                                                                           | instantaneous sfc downward sw flux                   |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dnirbmi_cpl    |                                                                                           | instantaneous sfc nir beam downward sw flux          |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dnirdfi_cpl    |                                                                                           | instantaneous sfc nir diff downward sw flux          |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvisbmi_cpl    |                                                                                           | instantaneous sfc uv+vis beam downward sw flux       |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dvisdfi_cpl    |                                                                                           | instantaneous sfc uv+vis diff downward sw flux       |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nlwsfci_cpl    |                                                                                           | instantaneous net sfc downward lw flux               |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nswsfci_cpl    |                                                                                           | instantaneous net sfc downward sw flux               |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nnirbmi_cpl    |                                                                                           | instantaneous net nir beam sfc downward sw flux      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nnirdfi_cpl    |                                                                                           | instantaneous net nir diff sfc downward sw flux      |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nvisbmi_cpl    |                                                                                           | instantaneous net uv+vis beam downward sw flux       |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%nvisdfi_cpl    |                                                                                           | instantaneous net uv+vis diff downward sw flux       |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%t2mi_cpl       |                                                                                           | instantaneous T2m                                    |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%q2mi_cpl       |                                                                                           | instantaneous Q2m                                    |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%u10mi_cpl      |                                                                                           | instantaneous U10m                                   |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%v10mi_cpl      |                                                                                           | instantaneous V10m                                   |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%tsfci_cpl      |                                                                                           | instantaneous sfc temperature                        |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%psurfi_cpl     |                                                                                           | instantaneous sfc pressure                           |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%oro_cpl        |                                                                                           | orography                                            |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%slmsk_cpl      |                                                                                           | land/sea/ice mask                                    |               |    1 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%shum_wts       |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%sppt_wts       |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%skebu_wts      |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%skebv_wts      |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%vcu_wts        |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%vcv_wts        |                                                                                           |                                                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dqdti          |                                                                                           | instantaneous total moisture tendency                | kg kg-1 s-1   |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%cnvqci         |                                                                                           | instantaneous total convective conensate             | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%upd_mfi        |                                                                                           | instantaneous convective updraft mass flux           |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%dwn_mfi        |                                                                                           | instantaneous convective downdraft mass flux         |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%det_mfi        |                                                                                           | instantaneous convective detrainment mass flux       |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Coupling(i)%cldcovi        |                                                                                           | instantaneous 3D cloud fraction                      |               |    2 | real    | kind_phys | none   | F        |
!! | physics%Cldprop(i)%cv                   | fraction_of_convective_cloud                            | fraction of convective cloud                            | frac          |    1 | real    | kind_phys | none   | F        |
!! | physics%Cldprop(i)%cvt                  | pressure_at_top_of_convective_cloud                     | convective cloud top pressure                           | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%Cldprop(i)%cvb                  | pressure_at_bottom_of_convective_cloud                  | convective cloud bottom pressure                        | Pa            |    1 | real    | kind_phys | none   | F        |
!! | physics%n_ozone_layers       | vertical_dimension_of_ozone_forcing_data_from_host          | number of vertical layers in ozone forcing data coming from host        | count         |    0 | integer                       |           | none     | F        |
!! | physics%n_ozone_lats         | number_of_latitutde_points_in_ozone_forcing_data_from_host  | number of latitude points in ozone forcing data coming from host        | count         |    0 | integer                       |           | none     | F        |
!! | physics%n_ozone_times        | number_of_time_levels_in_ozone_forcing_data_from_host       | number of time levels in ozone forcing data coming from host            | count         |    0 | integer                       |           | none     | F        |
!! | physics%n_ozone_coefficients | number_of_coefficients_in_ozone_forcing_data_from_host      | number of coeffcients in ozone forcing data coming from host            | count         |    0 | integer                       |           | none     | F        |
!! | physics%ozone_lat           | latitude_of_ozone_forcing_data_from_host                    | latitude value of the ozone forcing data coming from host               | degree        |    1 | real                          | kind_phys | none     | F        |
!! | physics%ozone_pres           | natural_log_of_ozone_forcing_data_pressure_levels_from_host | natural logarithm of the pressure levels of the ozone forcing data      | Pa            |    1 | real                          | kind_phys | none     | F        |
!! | physics%ozone_time           | time_levels_in_ozone_forcing_data_from_host                 | time values of the ozone forcing data coming from host                  | day           |    1 | real                          | kind_phys | none     | F        |
!! | physics%ozone_forcing_in     | ozone_forcing_from_host                                     | ozone forcing data from host                                            | various       |    4 | real                          | kind_phys | none     | F        |

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
    type(GFS_sfccycle_type), allocatable     :: Sfccycle(:)
    type(GFS_interstitial_type), allocatable :: Interstitial(:)
    type(GFS_init_type), allocatable         :: Init_parm(:)
    integer                                  :: LTP

    integer :: n_ozone_coefficients
    integer :: n_ozone_layers
    integer :: n_ozone_times
    integer :: n_ozone_lats

    real(kind=kind_phys), allocatable :: ozone_lat(:), ozone_pres(:), ozone_time(:)
    real(kind=kind_phys), allocatable :: ozone_forcing_in(:,:,:,:)

    contains
      procedure :: create => physics_create
      procedure :: associate => physics_associate
  end type physics_type

  contains

  subroutine scm_state_create(scm_state, n_columns, n_levels, n_time_levels)
    class(scm_state_type)             :: scm_state
    integer, intent(in)               :: n_columns, n_levels, n_time_levels

    scm_state%experiment_name = clear_char
    scm_state%model_name = clear_char
    scm_state%output_dir = clear_char
    scm_state%physics_suite_dir = clear_char
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
    scm_state%n_tracers = 3
    allocate(scm_state%tracer_names(scm_state%n_tracers))
    scm_state%water_vapor_index = 1
    scm_state%ozone_index = 2
    scm_state%cloud_water_index = 3
    scm_state%tracer_names(1) = 'vap_wat'
    scm_state%tracer_names(2) = 'o3mr'
    scm_state%tracer_names(3) = 'liq_wat'
    scm_state%n_itt_swrad = int_zero
    scm_state%n_itt_lwrad = int_zero
    scm_state%n_itt_out = int_zero
    scm_state%n_levels_smooth = 5
    allocate(scm_state%blksz(n_columns))
    scm_state%blksz = int_one

    scm_state%sfc_flux_spec = .false.
    scm_state%sfc_type = int_zero
    scm_state%sfc_type_real = real_zero
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
    scm_state%temp_tracer = real_zero
    scm_state%temp_T = real_zero
    scm_state%temp_u = real_zero
    scm_state%temp_v = real_zero

    allocate(scm_state%w_ls(n_columns, n_levels), scm_state%omega(n_columns, 1, n_levels), scm_state%u_g(n_columns, n_levels), &
      scm_state%v_g(n_columns, n_levels), scm_state%dT_dt_rad(n_columns, n_levels), scm_state%h_advec_thil(n_columns, n_levels), &
      scm_state%h_advec_qt(n_columns, n_levels), scm_state%v_advec_thil(n_columns, n_levels), &
      scm_state%v_advec_qt(n_columns, n_levels), scm_state%pres_surf(n_columns,1), scm_state%T_surf(n_columns,1), &
      scm_state%u_nudge(n_columns, n_levels), scm_state%v_nudge(n_columns, n_levels), &
      scm_state%T_nudge(n_columns, n_levels), scm_state%thil_nudge(n_columns, n_levels), &
      scm_state%qt_nudge(n_columns, n_levels), scm_state%sh_flux(n_columns), scm_state%lh_flux(n_columns), &
      scm_state%u_force_tend(n_columns,n_levels), scm_state%v_force_tend(n_columns,n_levels), &
      scm_state%T_force_tend(n_columns,n_levels), scm_state%qv_force_tend(n_columns,n_levels))
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

  end subroutine scm_state_create

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

  subroutine physics_create(physics, n_columns, n_levels, lats, pres)
    class(physics_type) :: physics
    integer, intent(in) :: n_columns, n_levels
    real(kind=kind_phys), intent(in) :: lats(:), pres(:)

    real(kind=kind_phys) :: kind_phys_zero

    integer :: i
    integer, dimension(8) :: zeroes_8

    zeroes_8(:) = int_zero
    kind_phys_zero = real_zero

    allocate(physics%Model(n_columns), physics%Statein(n_columns), physics%Stateout(n_columns), physics%Sfcprop(n_columns), &
      physics%Coupling(n_columns), physics%Grid(n_columns), physics%Tbd(n_columns), physics%Cldprop(n_columns), &
      physics%Radtend(n_columns), physics%Diag(n_columns), physics%Sfccycle(n_columns), physics%Interstitial(n_columns), &
      physics%Init_parm(n_columns))

    do i=1, n_columns
      physics%Init_parm(i)%me = int_one
      physics%Init_parm(i)%master = int_one
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
      physics%Init_parm(i)%logunit= int_one
      physics%Init_parm(i)%bdat(:) = zeroes_8(:)
      physics%Init_parm(i)%cdat(:) = zeroes_8(:)
      physics%Init_parm(i)%dt_dycore = kind_phys_zero
      physics%Init_parm(i)%dt_phys = kind_phys_zero
      physics%Init_parm(i)%ak => null()
      physics%Init_parm(i)%bk => null()
      physics%Init_parm(i)%xlon => null()
      physics%Init_parm(i)%xlat => null()
      physics%Init_parm(i)%area => null()
      physics%Init_parm(i)%tracer_names => null()
      physics%Init_parm(i)%blksz => null()
    end do

    !set ozone forcing array dimensions
    physics%n_ozone_coefficients = 4
    physics%n_ozone_layers = n_levels
    physics%n_ozone_lats = n_columns
    physics%n_ozone_times = 2

    allocate(physics%ozone_lat(physics%n_ozone_lats), physics%ozone_pres(physics%n_ozone_layers), &
      physics%ozone_time(physics%n_ozone_times+1), &
      physics%ozone_forcing_in(physics%n_ozone_lats, physics%n_ozone_layers, physics%n_ozone_coefficients, physics%n_ozone_times))
    physics%ozone_lat = lats
    physics%ozone_pres = log(pres)
    physics%ozone_time = (/12.0, 13.0, 14.0/)
    physics%ozone_forcing_in = real_zero

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


  end subroutine physics_associate

end module gmtb_scm_type_defs
