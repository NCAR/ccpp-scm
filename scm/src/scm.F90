module scm_main

implicit none

contains

subroutine scm_main_sub()

  use scm_kinds, only: sp, dp, qp
  use scm_input
  use scm_utils
  use scm_vgrid
  use scm_setup
  use scm_forcing
  use scm_time_integration
  use scm_output
  use scm_type_defs
       
  use :: ccpp_static_api,                      &
         only: ccpp_physics_init,              &
               ccpp_physics_timestep_init,     &
               ccpp_physics_run,               &
               ccpp_physics_timestep_finalize, &
               ccpp_physics_finalize


  implicit none

  type(scm_state_type), target :: scm_state
  type(scm_input_type), target :: scm_input_instance
  type(scm_reference_type), target :: scm_reference

  integer      :: i, j, kdt_rad, idtend, itrac
  real(kind=8) :: rinc(5) !(DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
  integer      :: jdat(1:8)

  integer                           :: cdata_time_index
  integer                           :: ierr
  character(len=16) :: logfile_name
  logical                           :: in_spinup

  call get_config_nml(scm_state)
  
  select case(scm_state%input_type)
    case(0)
      call get_case_init(scm_state, scm_input_instance)
    case(1)
      call get_case_init_DEPHY(scm_state, scm_input_instance)
    case default
      write(*,*) 'An unrecognized specification of the input_type namelist variable is being used. Exiting...'
      stop
  end select
  
  call get_reference_profile(scm_state, scm_reference)

  call get_FV3_vgrid(scm_input_instance, scm_state)

  !allocate(cdata_cols(scm_state%n_cols))

  call set_state(scm_input_instance, scm_reference, scm_state)

  call calc_geopotential(1, scm_state)

  scm_state%model_time = 0.0
  scm_state%itt = 1

  in_spinup = (scm_state%do_spinup .and. scm_state%itt <= scm_state%spinup_timesteps)
  
  if (in_spinup) then
    call set_spinup_nudging(scm_state)
  end if
  
  call interpolate_forcing(scm_input_instance, scm_state, in_spinup)

  call physics%create(scm_state%n_cols)
  
  !physics initialization section

  !set the array index of the time level of the state variables that the cdata
  !points to (this is the time level that will be updated during ipd_run;
  !if suite returns tendencies, SCM must apply them to this time level)
  select case(scm_state%time_scheme)
    case(1)
      cdata_time_index = 1
    case default
      cdata_time_index = 2
  end select
  
  !open a logfile
  if (physics%Init_parm%me == physics%Init_parm%master .and. physics%Init_parm%logunit>=0) then
    write (logfile_name, '(A7,I0.5,A4)') 'logfile.out'
    open(unit=physics%Init_parm%logunit, file=trim(scm_state%output_dir)//'/'//logfile_name, action='write', status='replace')
  end if
  
  physics%Init_parm%levs = scm_state%n_levels
  physics%Init_parm%bdat(1) = scm_state%init_year
  physics%Init_parm%bdat(2) = scm_state%init_month
  physics%Init_parm%bdat(3) = scm_state%init_day
  physics%Init_parm%bdat(5) = scm_state%init_hour
  physics%Init_parm%cdat(:) = physics%Init_parm%bdat
  physics%Init_parm%dt_dycore = scm_state%dt
  physics%Init_parm%dt_phys = scm_state%dt
  physics%Init_parm%ak => scm_state%a_k
  physics%Init_parm%bk => scm_state%b_k
  !physics%Init_parm%xlon => scm_state%lon  !rank mismatch -> why does Init_parm%xlon have 2 dimensions?
  !physics%Init_parm%xlat => scm_state%lat  !rank mismatch -> why does Init_parm%xlat have 2 dimensions?
  !physics%Init_parm%area => scm_state%area !rank mismatch -> why does Init_parm%area have 2 dimensions?
  physics%Init_parm%tracer_names => scm_state%tracer_names
  physics%Init_parm%tracer_types => scm_state%tracer_types
  physics%Init_parm%fn_nml = scm_state%physics_nml
  physics%Init_parm%blksz => scm_state%blksz
  physics%Init_parm%tile_num = 1
  physics%Init_parm%hydrostatic = .true.
  physics%Init_parm%restart = .false.
  physics%Init_parm%nwat = scm_state%nwat
  
  ! Allocate and initialize DDTs
  call GFS_suite_setup(physics%Model, physics%Statein, physics%Stateout,           &
                       physics%Sfcprop, physics%Coupling, physics%Grid,            &
                       physics%Tbd, physics%Cldprop, physics%Radtend,              &
                       physics%Diag, physics%Interstitial, 0, 1, 1,                &
                       physics%Init_parm, scm_state%n_cols, scm_state%lon,         &
                       scm_state%lat, scm_state%area)
  
  !override radiation frequency
  if (scm_state%force_rad_T > 0) then
    !turn off radiation since it is already accounted for in the forcing
    !set radiation calling periods to -1 in order to prevent lsswr/lslwr from being be set to true in GFS_time_vary_pre
    physics%Model%nsswr = -1
    physics%Model%nslwr = -1
  end if
  
  !override fhzero in physics namelist if n_itt_diag is set in the case namelist
  if (scm_state%n_itt_diag >= 1) then
    physics%Model%nszero = scm_state%n_itt_diag
    physics%Model%fhzero = scm_state%n_itt_diag*scm_state%dt/3600.0
  end if
    
  !check for problematic diagnostic and radiation periods
  if (mod(physics%Model%nszero,scm_state%n_itt_out) /= 0) then
    write(*,*) "***ERROR***: The diagnostic output period must be a multiple of the output period."
    write(*,*) "From ", adjustl(trim(scm_state%physics_nml)), ", fhzero = ",physics%Model%fhzero
    write(*,*) "implying a diagnostic output period of ", physics%Model%nszero*scm_state%dt, "seconds."
    write(*,*) "The given output period in the case configuration namelist is ", scm_state%output_period,"seconds."
    STOP
  end if
  
  if (mod(physics%Model%nsswr,scm_state%n_itt_out) /= 0) then
    write(*,*) "***WARNING***: The shortwave radiation calling period is different than the output period."
    write(*,*) "From ", adjustl(trim(scm_state%physics_nml)), ", fhswr = ",physics%Model%fhswr
    write(*,*) "The given output period in the case configuration namelist is ", scm_state%output_period
    write(*,*) "This will cause the effective output period of variables that are only given values during shortwave calls to be ",&
      lcm(scm_state%n_itt_out,physics%Model%nsswr)*scm_state%dt," seconds."
  end if
  
  if (mod(physics%Model%nslwr,scm_state%n_itt_out) /= 0) then
    write(*,*) "***WARNING***: The longwave radiation calling period is different than the output period."
    write(*,*) "From ", adjustl(trim(scm_state%physics_nml)), ", fhlwr = ",physics%Model%fhlwr
    write(*,*) "The given output period in the case configuration namelist is ", scm_state%output_period
    write(*,*) "This will cause the effective output period of variables that are only given values during longwave calls to be ",&
      lcm(scm_state%n_itt_out,physics%Model%nslwr)*scm_state%dt," seconds."
  end if
  
  cdata%blk_no = 1
  cdata%thrd_no = 1
  
  call physics%associate(scm_state)
  call physics%set(scm_input_instance, scm_state)
  
  ! When asked to calculate 3-dim. tendencies, set Stateout variables to
  ! Statein variables here in order to capture the first call to dycore
  if (physics%Model%ldiag3d) then
    physics%Stateout%gu0 = physics%Statein%ugrs
    physics%Stateout%gv0 = physics%Statein%vgrs
    physics%Stateout%gt0 = physics%Statein%tgrs
    physics%Stateout%gq0 = physics%Statein%qgrs
  endif
  
  !initialize the column's physics

  write(0,'(a,i0,a)') "Calling ccpp_physics_init with suite '" // trim(trim(adjustl(scm_state%physics_suite_name))) // "'"
  call ccpp_physics_init(cdata, suite_name=trim(trim(adjustl(scm_state%physics_suite_name))), ierr=ierr)
  write(0,'(a,i0,a,i0)') "Called ccpp_physics_init with suite '" // trim(trim(adjustl(scm_state%physics_suite_name))) // "', ierr=", ierr
  if (ierr/=0) then
      write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_init: ' // trim(cdata%errmsg) // '. Exiting...'
      stop 1
  end if
  
  physics%Model%first_time_step = .true.

  call output_init(scm_state, physics)
  call output_append(scm_state, physics, force=.true.)

  !first time step (call once)
  if (.not. in_spinup) then
     scm_state%dt_now = scm_state%dt
     scm_state%model_time = scm_state%dt_now
  end if
  
  call interpolate_forcing(scm_input_instance, scm_state, in_spinup)
  
  if (.not. scm_state%model_ics) call calc_pres_exner_geopotential(1, scm_state)
  
  !pass in state variables to be modified by forcing and physics
  call do_time_step(scm_state, physics, cdata, in_spinup)

  if (.not. in_spinup) then
    call output_append(scm_state, physics)
  end if

  !prepare for time loop
  scm_state%n_timesteps = ceiling(scm_state%runtime/scm_state%dt) + scm_state%spinup_timesteps
  
  scm_state%dt_now = scm_state%dt
  
  !if (.not. in_spinup) then
    physics%Model%first_time_step = .false.
  !end if

  do i = 2, scm_state%n_timesteps
    !are we still in spinup mode?
    if (scm_state%do_spinup .and. i <= scm_state%spinup_timesteps) then
      in_spinup = .true.
    else
      in_spinup = .false.
      scm_state%itt = i - scm_state%spinup_timesteps
    end if
        
    !>  - Calculate the elapsed model time.
    if (.not. in_spinup) then
      scm_state%model_time = scm_state%itt*scm_state%dt

      rinc = 0
      rinc(4) = (scm_state%itt-1)*scm_state%dt
      !w3movdat is a GFS routine to calculate the current date (jdat) from an elapsed time and an initial date (rinc is single prec.)
      call w3movdat(rinc, physics%Model%idat, jdat)
      physics%Model%jdat = jdat
    end if

    call interpolate_forcing(scm_input_instance, scm_state, in_spinup)
    
    call calc_pres_exner_geopotential(1, scm_state)

    !zero out diagnostics output on EVERY time step - breaks diagnostics averaged over many timesteps
    !call physics%Diag%rad_zero(physics%Model)
    !call physics%Diag%phys_zero(physics%Model)

    !pass in state variables to be modified by forcing and physics
    call do_time_step(scm_state, physics, cdata, in_spinup)
    
    write(*,*) "itt = ",scm_state%itt
    write(*,*) "model time (s) = ",scm_state%model_time
    if (scm_state%lsm_ics .or. scm_state%model_ics) then
      write(*,*) "Bowen ratio: ",physics%Interstitial%dtsfc1(1)/physics%Interstitial%dqsfc1(1)
      write(*,*) "sensible heat flux (W m-2): ",physics%Interstitial%dtsfc1(1)
      write(*,*) "latent heat flux (W m-2): ",physics%Interstitial%dqsfc1(1)
    end if
    
    if (.not. in_spinup) then
      call output_append(scm_state, physics)
    end if
  end do

  call ccpp_physics_finalize(cdata, suite_name=trim(trim(adjustl(scm_state%physics_suite_name))), ierr=ierr)

  if (ierr/=0) then
      write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_finalize: ' // trim(cdata%errmsg) // '. Exiting...'
      stop 1
  end if

end subroutine scm_main_sub

end module scm_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref scm_main_sub above.
program scm
  use scm_main
  call scm_main_sub()
end program scm
