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
  use mpi_f08

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

  type(MPI_Comm)           :: fcst_mpi_comm

  integer      :: i, j, kdt_rad, idtend, itrac
  real(kind=8) :: rinc(5) !(DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
  integer      :: jdat(1:8)

  integer                           :: cdata_time_index
  integer                           :: ierr
  character(len=16) :: logfile_name
  logical                           :: in_spinup

  call MPI_INIT(ierr)
  if (ierr/=0) then
      write(*,*) 'An error occurred in MPI_INIT: ', ierr
      error stop
  end if
  fcst_mpi_comm = MPI_COMM_WORLD

  call get_config_nml(scm_state)

  select case(scm_state%input_type)
    case(0)
      call get_case_init(scm_state, scm_input_instance)
    case(1)
      call get_case_init_DEPHY(scm_state, scm_input_instance)
    case default
      write(*,*) 'An unrecognized specification of the input_type namelist variable is being used. Exiting...'
      error stop
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
    case(2)
      cdata_time_index = 2
    case default
      cdata_time_index = 2
  end select

  !open a logfile
  if (physics%Init_parm%me == physics%Init_parm%master .and. physics%Init_parm%logunit>=0) then
    write (logfile_name, '(A7,I0.5,A4)') 'logfile.out'
    open(unit=physics%Init_parm%logunit, file=trim(scm_state%output_dir)//'/'//logfile_name, action='write', status='replace')
  end if

  physics%Init_parm%fcst_mpi_comm   =  fcst_mpi_comm
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
                       physics%Diag, physics%Interstitial, 1, 1,                   &
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
    error stop
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
  cdata%thrd_cnt = 1

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
      error stop
  end if

  physics%Model%first_time_step = .true.

  call output_init(scm_state, physics)
  call output_append(scm_state, physics, force=.true.)

  !first time step (call once)

  if (scm_state%time_scheme == 1) then
     if (.not. in_spinup) then
       scm_state%dt_now = scm_state%dt
       scm_state%model_time = scm_state%dt_now
     end if

     call interpolate_forcing(scm_input_instance, scm_state, in_spinup)

     if (.not. scm_state%model_ics) call calc_pres_exner_geopotential(1, scm_state)

     !pass in state variables to be modified by forcing and physics
     call do_time_step(scm_state, physics, cdata, in_spinup)

  else if (scm_state%time_scheme == 2) then
  !   !if using the leapfrog scheme, we initialize by taking one half forward time step and one half (unfiltered) leapfrog time step to get to the end of the first time step
    if (.not. in_spinup) then
      scm_state%dt_now = 0.5*scm_state%dt
      scm_state%model_time = scm_state%dt_now
    end if

    !save initial state
    scm_state%temp_tracer(:,:,:,1) = scm_state%state_tracer(:,:,:,1)
    scm_state%temp_T(:,:,1) = scm_state%state_T(:,:,1)
    scm_state%temp_u(:,:,1) = scm_state%state_u(:,:,1)
    scm_state%temp_v(:,:,1) = scm_state%state_v(:,:,1)

    call interpolate_forcing(scm_input_instance, scm_state, in_spinup)

    call calc_pres_exner_geopotential(1, scm_state)

    if (scm_state%input_type == 0) then
      call apply_forcing_forward_Euler(scm_state, in_spinup)
    else
      call apply_forcing_DEPHY(scm_state, in_spinup)
    end if

    !apply_forcing_forward_Euler updates state variables time level 1, so must copy this data to time_level 2 (where cdata points)
    scm_state%state_T(:,:,2) = scm_state%state_T(:,:,1)
    scm_state%state_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,1)
    scm_state%state_u(:,:,2) = scm_state%state_u(:,:,1)
    scm_state%state_v(:,:,2) = scm_state%state_v(:,:,1)

    ! Calculate total non-physics tendencies by substracting old Stateout
    ! variables from new/updated Statein variables (gives the tendencies
    ! due to anything else than physics)
    do i=1, scm_state%n_cols
      if (physics%Model%ldiag3d) then
        idtend = physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_non_physics)
        if(idtend>=1) then
          physics%Diag%dtend(i,:,idtend) = physics%Diag%dtend(i,:,idtend) &
                 + (physics%Statein%ugrs(i,:) - physics%Stateout%gu0(i,:))
        endif

        idtend = physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_non_physics)
        if(idtend>=1) then
          physics%Diag%dtend(i,:,idtend) = physics%Diag%dtend(i,:,idtend) &
                 + (physics%Statein%vgrs(i,:) - physics%Stateout%gv0(i,:))
        endif

        idtend = physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_non_physics)
        if(idtend>=1) then
          physics%Diag%dtend(i,:,idtend) = physics%Diag%dtend(i,:,idtend) &
                 + (physics%Statein%tgrs(i,:) - physics%Stateout%gt0(i,:))
        endif

        if (physics%Model%qdiag3d) then
          do itrac=1,physics%Model%ntrac
            idtend = physics%Model%dtidx(itrac+100,physics%Model%index_of_process_non_physics)
            if(idtend>=1) then
              physics%Diag%dtend(i,:,idtend) = physics%Diag%dtend(i,:,idtend) &
                     + (physics%Statein%qgrs(i,:,itrac) - physics%Stateout%gq0(i,:,itrac))
            endif
          enddo
        endif
      endif
    end do

    call ccpp_physics_timestep_init(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
    if (ierr/=0) then
        write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_init: ' // trim(cdata%errmsg) // '. Exiting...'
        error stop
    end if

    !--- determine if radiation diagnostics buckets need to be cleared
    if (nint(physics%Model%fhzero*3600) >= nint(max(physics%Model%fhswr,physics%Model%fhlwr))) then
      if (mod(physics%Model%kdt,physics%Model%nszero) == 1 .or. physics%Model%nszero == 1) then
        call physics%Diag%rad_zero  (physics%Model)
      endif
    else
      kdt_rad = nint(min(physics%Model%fhswr,physics%Model%fhlwr)/physics%Model%dtp)
      if (mod(physics%Model%kdt,kdt_rad) == 1) then
        call physics%Diag%rad_zero  (physics%Model)
      endif
    endif

    !--- determine if physics diagnostics buckets need to be cleared
    if (mod(physics%Model%kdt,physics%Model%nszero) == 1 .or. physics%Model%nszero == 1) then
      call physics%Diag%phys_zero (physics%Model)
    endif

    call ccpp_physics_run(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
    if (ierr/=0) then
        write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_run: ' // trim(cdata%errmsg) // '. Exiting...'
        error stop
    end if

    call ccpp_physics_timestep_finalize(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
    if (ierr/=0) then
        write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_finalize: ' // trim(cdata%errmsg) // '. Exiting...'
        error stop
    end if

    !the filter routine (called after the following leapfrog time step) expects time level 2 in temp_tracer to be the updated, unfiltered state after the previous time step
    scm_state%temp_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,2)
    scm_state%temp_T(:,:,2) = scm_state%state_T(:,:,2)
    scm_state%temp_u(:,:,2) = scm_state%state_u(:,:,2)
    scm_state%temp_v(:,:,2) = scm_state%state_v(:,:,2)

    !do half a leapfrog time step to get to the end of one full time step
    if (.not. in_spinup) then
      scm_state%model_time = scm_state%dt
    end if
    call interpolate_forcing(scm_input_instance, scm_state, in_spinup)

    call calc_pres_exner_geopotential(1, scm_state)

    !calling do_time_step with the leapfrog scheme active expects state variables in time level 1 to have values from 2 time steps ago, so set them equal to the initial values
    scm_state%state_T(:,:,1) = scm_state%temp_T(:,:,1)
    scm_state%state_u(:,:,1) = scm_state%temp_u(:,:,1)
    scm_state%state_v(:,:,1) = scm_state%temp_v(:,:,1)
    scm_state%state_tracer(:,:,:,1) = scm_state%temp_tracer(:,:,:,1)

    !go forward one leapfrog time step
    call do_time_step(scm_state, physics, cdata, in_spinup)

    !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
    call filter(scm_state)

    !> \todo tracers besides water vapor do not need to be filtered (is this right?)
    scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,scm_state%cloud_water_index,2)
    scm_state%state_tracer(:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,scm_state%ozone_index,2)
  end if

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

    !>  - Save previously unfiltered state as temporary for use in the time filter.
    if(scm_state%time_scheme == 2) then
      scm_state%temp_tracer = scm_state%state_tracer
      scm_state%temp_T = scm_state%state_T
      scm_state%temp_u = scm_state%state_u
      scm_state%temp_v = scm_state%state_v
    end if

    call interpolate_forcing(scm_input_instance, scm_state, in_spinup)

    call calc_pres_exner_geopotential(1, scm_state)

    !zero out diagnostics output on EVERY time step - breaks diagnostics averaged over many timesteps
    !call physics%Diag%rad_zero(physics%Model)
    !call physics%Diag%phys_zero(physics%Model)

    !pass in state variables to be modified by forcing and physics
    call do_time_step(scm_state, physics, cdata, in_spinup)

    if (scm_state%time_scheme == 2) then
      !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
      call filter(scm_state)

      !> \todo tracers besides water vapor do not need to be filtered (is this right?)
      scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,scm_state%cloud_water_index,2)
      scm_state%state_tracer(:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,scm_state%ozone_index,2)
    end if

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
      error stop
  end if

  call MPI_FINALIZE(ierr)
  if (ierr/=0) then
      write(*,*) 'An error occurred in MPI_FINALIZE: ', ierr
      error stop
  end if

end subroutine scm_main_sub

end module scm_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref scm_main_sub above.
program scm
  use scm_main
  use mpi_f08

  call scm_main_sub()
end program scm
