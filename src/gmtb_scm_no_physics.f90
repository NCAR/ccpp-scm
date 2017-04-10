module gmtb_scm_main

contains

subroutine gmtb_scm_main_sub()

  use gmtb_scm_kinds, only: sp, dp, qp
  use gmtb_scm_input
  use gmtb_scm_utils
  use gmtb_scm_vgrid
  use gmtb_scm_setup
  use gmtb_scm_forcing
  use gmtb_scm_time_integration
  use gmtb_scm_output

  use            :: ccpp_types,                         &
                    only: STR_LEN, ccpp_t, ccpp_suite_t
  use            :: ccpp,                               &
                    only: ccpp_init
  use            :: ccpp_ipd,                           &
                    only: ccpp_ipd_run
  use            :: ccpp_fields,                        &
                    only: ccpp_fields_add

  implicit none

  character(len=80)                 :: experiment_name !> name of model configuration file
  character(len=80)                 :: model_name !< name of "host" model (must be "GFS" for prototype)
  character(len=80), allocatable    :: physics_suite_name(:) !< name of physics suite (must be "GFS_operational" for prototype)
  character(len=80)                 :: output_dir !< name of output directory to place netCDF file
  character(len=80)                 :: physics_suite_dir !< location of the physics suite XML files for the IPD (relative to the executable path)
  character(len=80)                 :: output_file !< name of output file (without the file extension)
  character(len=80)                 :: case_name !< name of case initialization and forcing to use (different than experiment name, which names the model run (as a control, experiment_1, etc.))
  character(len=5)                  :: i_string

  integer                           :: i, k, ioerror, allocate_status, grid_error !< dummy indices and error statuses
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
  integer, parameter                :: n_levels_smooth = 5 !< the number of levels over which the input profiles are smoothed into the reference profiles
  integer, parameter                :: n_tracers = 3 !< number of tracers
  integer, parameter                :: ozone_index = 2 !< index for ozone in the tracer array
  integer, parameter                :: cloud_water_index = 3 !< index for cloud water in the tracer array

  logical                           :: use_IPD !< flag for using the GMTB IPD
  logical                           :: use_leapfrog !< flag for using the leapfrog as the time-integration scheme for the forcing
  logical                           :: use_forward !< flag for using the forward Euler scheme for time-stepping the forcing
  logical                           :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
  integer                           :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer                           :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer                           :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere

  real(kind=dp)                           :: model_time !< elapsed model time (s)
  real(kind=dp)                           :: dt !< physics time step (s)
  real(kind=dp)                           :: dt_now !< time step currently being used (if it changes due to time-stepping scheme)
  real(kind=dp)                           :: runtime !< total runtime (s)
  real(kind=dp)                           :: output_frequency !< how often output is written (s)
  real(kind=dp)                           :: swrad_frequency !< how often SW radiation is called
  real(kind=dp)                           :: lwrad_frequency !< how often LW radiation is called
  real(kind=dp)                           :: relax_time !< time scale for hor. wind nudging (s)
  real(kind=dp)                           :: deg_to_rad_const !< conversion constant from degrees to radians
  real(kind=dp), parameter                :: c_filter = 0.15 !< parameter that controls the amount of damping in the leapfrog filter

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

  integer                                 :: init_year, init_month, init_day, init_hour

  !> - Define the reference profile variables.
  integer                           :: ref_nlev !< number of levels in the reference profile
  real(kind=dp), allocatable              :: ref_pres(:) !< pressure (Pa) of the reference profile levels
  real(kind=dp), allocatable              :: ref_T(:) !< absolute T (K) of the reference profile levels
  real(kind=dp), allocatable              :: ref_qv(:) !< water vapor specific humidity (kg/kg) of the reference profile levels
  real(kind=dp), allocatable              :: ref_ozone(:) !< ozone mass mixing ratio (kg/kg) of the reference profile levels

  !> - Define the SCM state variables; variables with appended "i" are interface; variables with appended "l" are layer-centered.
  !!  - index order for grid is (horizontal, vertical);
  !!  - index order for state variables is (horizontal, vertical, timesteps);
  !!  - index order for tracer is (horizontal, vertical, tracer_index, timesteps)
  real(kind=dp), allocatable              :: pres_i(:,:), pres_l(:,:) !< pressure on grid interfaces, centers (Pa)
  real(kind=dp), allocatable              :: si(:,:), sl(:,:) !< sigma on grid interfaces, centers
  real(kind=dp), allocatable              :: exner_i(:,:), exner_l(:,:) !< exner function on grid interfaces, centers
  real(kind=dp), allocatable              :: geopotential_i(:,:), geopotential_l(:,:) !< geopotential on grid interfaces, centers
  real(kind=dp), allocatable, target              :: state_T(:,:,:) !< model state absolute temperature at grid centers (K)
  real(kind=dp), allocatable, target              :: state_u(:,:,:), state_v(:,:,:) !< model state horizontal winds at grid centers (m/s)
  real(kind=dp), allocatable, target              :: state_tracer(:,:,:,:) !< model state tracer at grid centers
  real(kind=dp), allocatable              :: temp_T(:,:,:), temp_u(:,:,:), temp_v(:,:,:), temp_tracer(:,:,:,:) !< used for time-filtering
  real(kind=dp), allocatable              :: lat(:), lon(:) !< latitude and longitude (radians)

  !> - Define forcing-related variables (indexing is (horizontal, vertical)).
  real(kind=dp), allocatable              :: u_force_tend(:,:), v_force_tend(:,:), T_force_tend(:,:), qv_force_tend(:,:) !< total u, v, T, q forcing (units/s) (horizontal, vertical)
  real(kind=dp), allocatable              :: w_ls(:,:), omega(:,:), u_g(:,:), v_g(:,:), dT_dt_rad(:,:), h_advec_thil(:,:), &
    h_advec_qt(:,:), v_advec_thil(:,:), v_advec_qt(:,:), u_nudge(:,:), v_nudge(:,:), T_nudge(:,:), thil_nudge(:,:), qt_nudge(:,:) !< forcing terms interpolated to the model time and grid
  real(kind=dp), allocatable              :: T_surf(:), pres_surf(:) !< surface temperature and pressure interpolated to the model time
  real(kind=dp), allocatable              :: sh_flux(:), lh_flux(:) !< surface sensible and latent heat fluxes interpolated to the model time

  real(kind=dp), allocatable              :: a_k(:), b_k(:) !< used to determine grid sigma and pressure levels

  type(ccpp_t), allocatable, target                      :: cdata(:)
  type(ccpp_suite_t), allocatable                        :: suite(:)
  integer                                                :: ipd_loop

  real, target  :: gravity

  integer, allocatable                                   :: n_phy_fields(:)
  integer                                                :: cdata_time_index
  integer                                                :: ierr !< Integer error flag

  use_IPD = .true.


  call get_config_nml(experiment_name, model_name, n_cols, case_name, dt, time_scheme, runtime, output_frequency, &
    swrad_frequency, lwrad_frequency, n_levels, output_dir, output_file, thermo_forcing_type, mom_forcing_type, relax_time, &
    sfc_flux_spec, reference_profile_choice, init_year, init_month, init_day, init_hour, physics_suite_name, physics_suite_dir, &
    n_phy_fields)

  call get_case_init(case_name, sfc_flux_spec, input_nlev, input_ntimes, input_pres, input_time, input_height, input_thetail, &
    input_qt, input_ql, input_qi, input_u, input_v, input_tke, input_ozone, input_lat, input_lon, input_pres_surf, input_T_surf, &
    input_sh_flux_sfc, input_lh_flux_sfc, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge, &
    input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
    input_v_advec_thetail, input_v_advec_qt)

  call get_reference_profile(ref_nlev, ref_pres, ref_T, ref_qv, ref_ozone, reference_profile_choice)

  select case(trim(adjustl(model_name)))
    case("GFS")
      !>  - Call get_GFS_grid in \ref vgrid to read in the necessary coefficients and calculate the pressure-related variables on the grid.
      call get_GFS_vgrid(input_pres_surf(1), n_levels, n_cols, pres_i, pres_l, si, sl, exner_l, exner_i, a_k, b_k, grid_error)
      !>  - Exit if an unsupported number of levels is specified or the file with grid coefficients cannot be opened.
      if (grid_error == 1) then
        write(*,*) 'When using the GFS host model, only 28, 42, 60, 64, and 91 levels are currently supported. Exiting...'
        stop
      end if
      if (grid_error == 2) then
        write(*,*) 'The grid coefficient file could not be opened. Exiting...'
        stop
      end if
    case default
      write(*,*) 'Only the GFS model is currently supported. Exiting...'
      stop
  end select

  select case(time_scheme)
    case(1)
      n_time_levels = 1
    case(2)
      n_time_levels = 2
    case default
      n_time_levels = 2
  end select

  allocate(lat(n_cols), lon(n_cols))

  allocate(cdata(n_cols), suite(n_cols))

  allocate(state_T(n_cols, n_levels, n_time_levels), &
    state_u(n_cols, n_levels, n_time_levels), state_v(n_cols, n_levels, n_time_levels), &
    state_tracer(n_cols, n_levels, n_tracers, n_time_levels))

  allocate(temp_tracer(n_cols, n_levels, n_tracers, n_time_levels), temp_T(n_cols, n_levels, n_time_levels), &
    temp_u(n_cols, n_levels, n_time_levels), temp_v(n_cols, n_levels, n_time_levels))

  allocate(w_ls(n_cols, n_levels), omega(n_cols, n_levels), u_g(n_cols, n_levels), v_g(n_cols, n_levels), &
    dT_dt_rad(n_cols, n_levels), h_advec_thil(n_cols, n_levels), h_advec_qt(n_cols, n_levels), v_advec_thil(n_cols, n_levels), &
    v_advec_qt(n_cols, n_levels), pres_surf(n_cols), T_surf(n_cols), u_nudge(n_cols, n_levels), v_nudge(n_cols, n_levels), &
    T_nudge(n_cols, n_levels), thil_nudge(n_cols, n_levels), qt_nudge(n_cols, n_levels), sh_flux(n_cols), lh_flux(n_cols))

  allocate(u_force_tend(n_cols,n_levels), v_force_tend(n_cols,n_levels), T_force_tend(n_cols,n_levels), &
    qv_force_tend(n_cols,n_levels))
  u_force_tend = 0.0
  v_force_tend = 0.0
  T_force_tend = 0.0
  qv_force_tend = 0.0

  !read in physics suite namelist here

  call set_state(input_nlev, input_pres, input_qt, input_thetail, input_ql, input_qi, input_u, input_v, input_ozone, &
    n_levels, n_cols, ozone_index, cloud_water_index, pres_l, n_levels_smooth, ref_nlev, ref_pres, ref_qv, ref_T, ref_ozone, &
    state_tracer(:,:,:,1), state_T(:,:,1), state_u(:,:,1), state_v(:,:,1))

  call calc_GFS_geopotential(n_levels, n_cols, state_T(:,:,1), state_tracer(:,:,1,1), exner_i, exner_l, geopotential_i, &
    geopotential_l)

  model_time = 0.0
  itt = 1

  call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge, &
    input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
    input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
    input_sh_flux_sfc, input_lh_flux_sfc, n_levels, n_cols, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, &
    thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, &
    lh_flux)

  call output_init(output_dir, output_file, n_cols, n_levels, init_year, init_month, init_day, init_hour)

  itt_out = 1

  call output_append(output_dir, output_file, itt_out, model_time, pres_l, pres_i, sl, si, state_tracer(:,:,1,1), &
    state_T(:,:,1), state_u(:,:,1), state_v(:,:,1), state_tracer(:,:,3,1), u_force_tend, v_force_tend, T_force_tend, &
    qv_force_tend, w_ls, u_g, v_g, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, T_surf, pres_surf)

  !physics initialization section

  !set the array index of the time level of the state variables that the cdata
  !points to (this is the time level that will be updated during ipd_run;
  !if suite returns tendencies, SCM must apply them to this time level)
  select case(time_scheme)
    case(1)
      cdata_time_index = 1
    case(2)
      cdata_time_index = 2
    case default
      cdata_time_index = 2
  end select

  do i = 1, n_cols
    !set up each column's physics suite (which may be different)

    call ccpp_init( &
         trim(adjustl(physics_suite_dir))//trim(adjustl(physics_suite_name(i)))//'.xml', &
         cdata(i), ierr)

    select case(physics_suite_name(i))
      case ('suite_DUMMY_scm')
        call ccpp_fields_add(cdata(i), 'temperature', 'K', &
                             state_T(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'eastward_wind', 'm s-1', &
                             state_u(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'northward_wind', 'm s-1', &
                             state_v(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'water_vapor_specific_humidity', &
                             'kg kg-1', state_tracer(i,:,1,cdata_time_index), &
                             ierr)
      case ('suite_DUMMY_scm2')
        call ccpp_fields_add(cdata(i), 'temperature', 'K', &
                             state_T(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'eastward_wind', 'm s-1', &
                             state_u(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'northward_wind', 'm s-1', &
                             state_v(i,:,cdata_time_index), ierr)
        call ccpp_fields_add(cdata(i), 'water_vapor_specific_humidity', &
                             'kg kg-1', state_tracer(i,:,1,cdata_time_index), &
                             ierr)
      case default
        write(i_string,'(I5)') i
        write(*,*) 'The physics suite '//trim(physics_suite_name(i))//' specified for column #'//i_string&
          //' is not set up yet to be used in this model. Stopping...'
        STOP
    end select
  end do

  !first time step (call once)

  if (time_scheme == 1) then
    dt_now = dt
    model_time = dt_now

    call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge,&
      input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
      input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
      input_sh_flux_sfc, input_lh_flux_sfc, n_levels, n_cols, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
      thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, &
      lh_flux)

    call calc_GFS_pres_exner_geopotential(pres_surf, n_levels, n_cols, a_k, b_k, state_T(:,:,1), state_tracer(:,:,1,1), pres_i, &
      pres_l, si, sl, exner_l, exner_i, geopotential_l, geopotential_i)

    !pass in state variables to be modified by forcing and physics
    call do_time_step(n_levels, n_cols, time_scheme, state_tracer, state_T, state_u, state_v, cdata, w_ls, omega, u_g, v_g, &
      u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend,&
      v_force_tend, T_force_tend, qv_force_tend, exner_l, geopotential_i, pres_l, lat, dt_now, thermo_forcing_type, &
      mom_forcing_type, relax_time)

  else if (time_scheme == 2) then
    !if using the leapfrog scheme, we initialize by taking one half forward time step and one half (unfiltered) leapfrog time step to get to the end of the first time step
    dt_now = 0.5*dt
    model_time = dt_now

    !save initial state
    temp_tracer(:,:,:,1) = state_tracer(:,:,:,1)
    temp_T(:,:,1) = state_T(:,:,1)
    temp_u(:,:,1) = state_u(:,:,1)
    temp_v(:,:,1) = state_v(:,:,1)

    call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge,&
      input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
      input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
      input_sh_flux_sfc, input_lh_flux_sfc, n_levels, n_cols, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
      thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, &
      lh_flux)

    call calc_GFS_pres_exner_geopotential(pres_surf, n_levels, n_cols, a_k, b_k, state_T(:,:,1), state_tracer(:,:,1,1), pres_i, &
      pres_l, si, sl, exner_l, exner_i, geopotential_l, geopotential_i)

    call apply_forcing_forward_Euler(n_levels, n_cols, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, &
      u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, &
      geopotential_i, pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend,&
      qv_force_tend)

    !apply_forcing_forward_Euler updates state variables time level 1, so must copy this data to time_level 2 (where cdata points)
    state_T(:,:,2) = state_T(:,:,1)
    state_tracer(:,:,:,2) = state_tracer(:,:,:,1)
    state_u(:,:,2) = state_u(:,:,1)
    state_v(:,:,2) = state_v(:,:,1)

    do i=1, n_cols
      do ipd_loop = 1 , cdata(i)%suite%ipds_max
          cdata(i)%suite%ipd_n = ipd_loop
          call ccpp_ipd_run(cdata(i))
      end do
    end do

    !the filter routine (called after the following leapfrog time step) expects time level 2 in temp_tracer to be the updated, unfiltered state after the previous time step
    temp_tracer(:,:,:,2) = state_tracer(:,:,:,2)
    temp_T(:,:,2) = state_T(:,:,2)
    temp_u(:,:,2) = state_u(:,:,2)
    temp_v(:,:,2) = state_v(:,:,2)

    !do half a leapfrog time step to get to the end of one full time step
    model_time = dt
    call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge,&
      input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
      input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
      input_sh_flux_sfc, input_lh_flux_sfc, n_levels, n_cols, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
      thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, &
      lh_flux)

    call calc_GFS_pres_exner_geopotential(pres_surf, n_levels, n_cols, a_k, b_k, state_T(:,:,1), state_tracer(:,:,1,1), pres_i, &
      pres_l, si, sl, exner_l, exner_i, geopotential_l, geopotential_i)

    !calling do_time_step with the leapfrog scheme active expects state variables in time level 1 to have values from 2 time steps ago, so set them equal to the initial values
    state_T(:,:,1) = temp_T(:,:,1)
    state_u(:,:,1) = temp_u(:,:,1)
    state_v(:,:,1) = temp_v(:,:,1)
    state_tracer(:,:,:,1) = temp_tracer(:,:,:,1)

    !go forward one leapfrog time step
    call do_time_step(n_levels, n_cols, time_scheme, state_tracer, state_T, state_u, state_v, cdata, w_ls, omega, u_g, v_g, &
      u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend,&
      v_force_tend, T_force_tend, qv_force_tend, exner_l, geopotential_i, pres_l, lat, dt_now, thermo_forcing_type, &
      mom_forcing_type, relax_time)

    !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
    call filter(c_filter, temp_tracer, temp_T, temp_u, temp_v, state_tracer, state_T, state_u, state_v)

    !> \todo tracers besides water vapor do not need to be filtered (is this right?)
    state_tracer(:,:,2,1) = state_tracer(:,:,2,2)
    state_tracer(:,:,3,1) = state_tracer(:,:,3,2)
  end if

  !prepare for time loop
  n_timesteps = ceiling(runtime/dt)
  n_itt_swrad = floor(swrad_frequency/dt)
  n_itt_lwrad = floor(lwrad_frequency/dt)
  n_itt_out = floor(output_frequency/dt)

  dt_now = dt

  do itt = 2, n_timesteps
    !>  - Calculate the elapsed model time.
    model_time = itt*dt

    !>  - Save previously unfiltered state as temporary for use in the time filter.
    if(time_scheme == 2) then
      temp_tracer = state_tracer
      temp_T = state_T
      temp_u = state_u
      temp_v = state_v
    end if

    call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge,&
      input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
      input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
      input_sh_flux_sfc, input_lh_flux_sfc, n_levels, n_cols, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
      thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, &
      lh_flux)

    call calc_GFS_pres_exner_geopotential(pres_surf, n_levels, n_cols, a_k, b_k, state_T(:,:,1), state_tracer(:,:,1,1), pres_i, &
      pres_l, si, sl, exner_l, exner_i, geopotential_l, geopotential_i)

    !pass in state variables to be modified by forcing and physics
    call do_time_step(n_levels, n_cols, time_scheme, state_tracer, state_T, state_u, state_v, cdata, w_ls, omega, u_g, v_g, &
      u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend,&
      v_force_tend, T_force_tend, qv_force_tend, exner_l, geopotential_i, pres_l, lat, dt_now, thermo_forcing_type, &
      mom_forcing_type, relax_time)

    select case(time_scheme)
      case (1)
        !for forward Euler scheme, no filtering is done; simply transfer output state variables from slot 2 to slot 1
        ! state_T(:,:,1) = state_T(:,:,2)
        ! state_u(:,:,1) = state_u(:,:,2)
        ! state_v(:,:,1) = state_v(:,:,2)
        ! state_tracer(:,:,:,1) = state_tracer(:,:,:,2)
      case (2)
        !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
        call filter(c_filter, temp_tracer, temp_T, temp_u, temp_v, state_tracer, state_T, state_u, state_v)

        !> \todo tracers besides water vapor do not need to be filtered (is this right?)
        state_tracer(:,:,2,1) = state_tracer(:,:,2,2)
        state_tracer(:,:,3,1) = state_tracer(:,:,3,2)
    end select

    if(mod(itt, n_itt_out)==0) then
      itt_out = itt_out+1
      write(*,*) "itt = ",itt
      write(*,*) "model time (s) = ",model_time
      write(*,*) "calling output routine..."

      call output_append(output_dir, output_file, itt_out, model_time, pres_l, pres_i, sl, si, state_tracer(:,:,1,1), &
        state_T(:,:,1), state_u(:,:,1), state_v(:,:,1), state_tracer(:,:,3,1), u_force_tend, v_force_tend, T_force_tend, &
        qv_force_tend, w_ls, u_g, v_g, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, T_surf, pres_surf)

    end if
  end do

end subroutine gmtb_scm_main_sub

end module gmtb_scm_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref gmtb_scm_main_sub above.
program gmtb_scm
  use gmtb_scm_main
  call gmtb_scm_main_sub()
end program gmtb_scm
