module gmtb_scm_main

implicit none

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
  use gmtb_scm_type_defs
       
  use :: ccpp_static_api,                    &
         only: ccpp_physics_init,            &
               ccpp_physics_run,             &
               ccpp_physics_finalize


  implicit none

  type(scm_state_type), target :: scm_state
  type(scm_input_type), target :: scm_input
  type(scm_reference_type), target :: scm_reference

  integer      :: i, j, grid_error
  real(kind=8) :: rinc(5) !(DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)
  integer      :: jdat(1:8)

  integer                           :: cdata_time_index
  integer                           :: ierr
  character(len=16) :: logfile_name

  call get_config_nml(scm_state)

  call get_case_init(scm_state, scm_input)

  call get_reference_profile(scm_state, scm_reference)

  select case(trim(adjustl(scm_state%model_name)))
    case("GFS")
      !>  - Call get_GFS_grid in \ref vgrid to read in the necessary coefficients and calculate the pressure-related variables on the grid.
      call get_GFS_vgrid(scm_input, scm_state, grid_error)
      !>  - Exit if an unsupported number of levels is specified or the file with grid coefficients cannot be opened.
      if (grid_error == 1) then
        write(*,*) 'When using the GFS host model, only 28, 42, 60, 64, and 91 levels are currently supported. Exiting...'
        stop
      end if
      if (grid_error == 2) then
        write(*,*) 'The grid coefficient file could not be opened. Exiting...'
        stop
      end if
    case("FV3")
      call get_FV3_vgrid(scm_input, scm_state)
    case default
      write(*,*) 'Only the GFS (GSM) and FV3 vertical coordinates are currently supported. Exiting...'
      stop
  end select

  allocate(cdata_cols(scm_state%n_cols))

  call set_state(scm_input, scm_reference, scm_state)

  call calc_geopotential(1, scm_state)

  scm_state%model_time = 0.0
  scm_state%itt = 1

  call interpolate_forcing(scm_input, scm_state)

  call output_init(scm_state)

  scm_state%itt_out = 1

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

  do i = 1, scm_state%n_cols
     !set up each column's physics suite (which may be different)
      ! call ccpp_init(trim(adjustl(scm_state%physics_suite_name(i))), cdata_cols(i), ierr)
      ! if (ierr/=0) then
      !     write(*,'(a,i0,a)') 'An error occurred in ccpp_init for column ', i, '. Exiting...'
      !     stop
      ! end if

     !open a logfile for each column
      if (physics%Init_parm(i)%me == physics%Init_parm(i)%master .and. physics%Init_parm(i)%logunit>=0) then
          write (logfile_name, '(A7,I0.5,A4)') 'logfile', i, '.out'
          open(unit=physics%Init_parm(i)%logunit, file=trim(scm_state%output_dir)//'/'//logfile_name, action='write', status='replace')
      end if

      cdata_cols(i)%blk_no = i
      cdata_cols(i)%thrd_no = 1

      physics%Init_parm(i)%levs = scm_state%n_levels
      physics%Init_parm(i)%bdat(1) = scm_state%init_year
      physics%Init_parm(i)%bdat(2) = scm_state%init_month
      physics%Init_parm(i)%bdat(3) = scm_state%init_day
      physics%Init_parm(i)%bdat(5) = scm_state%init_hour
      physics%Init_parm(i)%cdat(:) = physics%Init_parm(i)%bdat
      physics%Init_parm(i)%dt_dycore = scm_state%dt
      physics%Init_parm(i)%dt_phys = scm_state%dt
      physics%Init_parm(i)%ak => scm_state%a_k(1,:)
      physics%Init_parm(i)%bk => scm_state%b_k(1,:)
      physics%Init_parm(i)%xlon => scm_state%lon
      physics%Init_parm(i)%xlat => scm_state%lat
      physics%Init_parm(i)%area => scm_state%area
      physics%Init_parm(i)%tracer_names => scm_state%tracer_names
      physics%Init_parm(i)%fn_nml = scm_state%physics_nml(1)
      physics%Init_parm(i)%blksz => scm_state%blksz
      physics%Init_parm(i)%tile_num = 1
      physics%Init_parm(i)%hydrostatic = .true.
      physics%Init_parm(i)%restart = .false.
      
      ! Allocate and initialize DDTs
      call GFS_suite_setup(physics%Model(i), physics%Statein(i), physics%Stateout(i),           &
                           physics%Sfcprop(i), physics%Coupling(i), physics%Grid(i),            &
                           physics%Tbd(i), physics%Cldprop(i), physics%Radtend(i),              &
                           physics%Diag(i), physics%Interstitial(i), 0, 1, 1,                   &
                           physics%Init_parm(i))
      
      call physics%associate(scm_state, i)

      !initialize each column's physics

      write(0,'(a,i0,a)') "Calling ccpp_physics_init for column ", i, " with suite '" // trim(trim(adjustl(scm_state%physics_suite_name(i)))) // "'"
      call ccpp_physics_init(cdata_cols(i), suite_name=trim(trim(adjustl(scm_state%physics_suite_name(i)))), ierr=ierr)
      write(0,'(a,i0,a,i0)') "Called ccpp_physics_init for column ", i, " with suite '" // trim(trim(adjustl(scm_state%physics_suite_name(i)))) // "', ierr=", ierr
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_init for column ', i, ': ' // trim(cdata_cols(i)%errmsg) // '. Exiting...'
          stop
      end if
      
      physics%Model(i)%first_time_step = .true.
  end do

  call output_append(scm_state, physics)

  !first time step (call once)

  if (scm_state%time_scheme == 1) then
     scm_state%dt_now = scm_state%dt
     scm_state%model_time = scm_state%dt_now

     call interpolate_forcing(scm_input, scm_state)

     if (.not. scm_state%model_ics) call calc_pres_exner_geopotential(1, scm_state)

     !pass in state variables to be modified by forcing and physics
     call do_time_step(scm_state, cdata_cols)

  else if (scm_state%time_scheme == 2) then
  !   !if using the leapfrog scheme, we initialize by taking one half forward time step and one half (unfiltered) leapfrog time step to get to the end of the first time step
    scm_state%dt_now = 0.5*scm_state%dt
    scm_state%model_time = scm_state%dt_now

    !save initial state
    scm_state%temp_tracer(:,:,:,:,1) = scm_state%state_tracer(:,:,:,:,1)
    scm_state%temp_T(:,:,:,1) = scm_state%state_T(:,:,:,1)
    scm_state%temp_u(:,:,:,1) = scm_state%state_u(:,:,:,1)
    scm_state%temp_v(:,:,:,1) = scm_state%state_v(:,:,:,1)

    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    call apply_forcing_forward_Euler(scm_state)

    !apply_forcing_forward_Euler updates state variables time level 1, so must copy this data to time_level 2 (where cdata points)
    scm_state%state_T(:,:,:,2) = scm_state%state_T(:,:,:,1)
    scm_state%state_tracer(:,:,:,:,2) = scm_state%state_tracer(:,:,:,:,1)
    scm_state%state_u(:,:,:,2) = scm_state%state_u(:,:,:,1)
    scm_state%state_v(:,:,:,2) = scm_state%state_v(:,:,:,1)

    do i=1, scm_state%n_cols
      call ccpp_physics_run(cdata_cols(i), suite_name=trim(trim(adjustl(scm_state%physics_suite_name(i)))), ierr=ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_run for column ', i, ': ' // trim(cdata_cols(i)%errmsg) // '. Exiting...'
          stop
      end if
    end do

    !the filter routine (called after the following leapfrog time step) expects time level 2 in temp_tracer to be the updated, unfiltered state after the previous time step
    scm_state%temp_tracer(:,:,:,:,2) = scm_state%state_tracer(:,:,:,:,2)
    scm_state%temp_T(:,:,:,2) = scm_state%state_T(:,:,:,2)
    scm_state%temp_u(:,:,:,2) = scm_state%state_u(:,:,:,2)
    scm_state%temp_v(:,:,:,2) = scm_state%state_v(:,:,:,2)

    !do half a leapfrog time step to get to the end of one full time step
    scm_state%model_time = scm_state%dt
    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    !calling do_time_step with the leapfrog scheme active expects state variables in time level 1 to have values from 2 time steps ago, so set them equal to the initial values
    scm_state%state_T(:,:,:,1) = scm_state%temp_T(:,:,:,1)
    scm_state%state_u(:,:,:,1) = scm_state%temp_u(:,:,:,1)
    scm_state%state_v(:,:,:,1) = scm_state%temp_v(:,:,:,1)
    scm_state%state_tracer(:,:,:,:,1) = scm_state%temp_tracer(:,:,:,:,1)

    !go forward one leapfrog time step
    call do_time_step(scm_state, cdata_cols)

    !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
    call filter(scm_state)

    !> \todo tracers besides water vapor do not need to be filtered (is this right?)
    scm_state%state_tracer(:,:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,:,scm_state%cloud_water_index,2)
    scm_state%state_tracer(:,:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,:,scm_state%ozone_index,2)
  end if

  scm_state%itt_out = scm_state%itt_out + 1
  call output_append(scm_state, physics)

  !prepare for time loop
  scm_state%n_timesteps = ceiling(scm_state%runtime/scm_state%dt)
  scm_state%n_itt_out = floor(scm_state%output_frequency/scm_state%dt)

  scm_state%dt_now = scm_state%dt

  do i=1, scm_state%n_cols
    physics%Model(i)%first_time_step = .false.
  end do

  do i = 2, scm_state%n_timesteps
    scm_state%itt = i
    !>  - Calculate the elapsed model time.
    scm_state%model_time = scm_state%itt*scm_state%dt

    rinc = 0
    rinc(4) = (scm_state%itt-1)*scm_state%dt
    !w3movdat is a GFS routine to calculate the current date (jdat) from an elapsed time and an initial date (rinc is single prec.)
    call w3movdat(rinc, physics%Model(1)%idat, jdat)
    do j=1, scm_state%n_cols
      physics%Model(j)%jdat = jdat
    end do

    !>  - Save previously unfiltered state as temporary for use in the time filter.
    if(scm_state%time_scheme == 2) then
      scm_state%temp_tracer = scm_state%state_tracer
      scm_state%temp_T = scm_state%state_T
      scm_state%temp_u = scm_state%state_u
      scm_state%temp_v = scm_state%state_v
    end if

    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    !zero out diagnostics output on EVERY time step - breaks diagnostics averaged over many timesteps
    do j=1, scm_state%n_cols
      call physics%Diag(j)%rad_zero(physics%Model(j))
      call physics%Diag(j)%phys_zero(physics%Model(j))
    end do

    !pass in state variables to be modified by forcing and physics
    call do_time_step(scm_state, cdata_cols)

    if (scm_state%time_scheme == 2) then
      !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
      call filter(scm_state)

      !> \todo tracers besides water vapor do not need to be filtered (is this right?)
      scm_state%state_tracer(:,:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,:,scm_state%cloud_water_index,2)
      scm_state%state_tracer(:,:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,:,scm_state%ozone_index,2)
    end if

    if(mod(scm_state%itt, scm_state%n_itt_out)==0) then
      scm_state%itt_out = scm_state%itt_out+1
      write(*,*) "itt = ",scm_state%itt
      write(*,*) "model time (s) = ",scm_state%model_time
      write(*,*) "calling output routine..."

      call output_append(scm_state, physics)

    end if
  end do

  do i=1, scm_state%n_cols
      call ccpp_physics_finalize(cdata_cols(i), suite_name=trim(trim(adjustl(scm_state%physics_suite_name(i)))), ierr=ierr)

      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_finalize for column ', i, ': ' // trim(cdata_cols(i)%errmsg) // '. Exiting...'
          stop
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
