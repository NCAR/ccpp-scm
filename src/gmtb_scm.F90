module gmtb_scm_main

implicit none

contains

subroutine gmtb_scm_main_sub()

  use gmtb_scm_kinds, only: sp, dp, qp
  use gmtb_scm_type_defs, only: scm_state_type, scm_input_type, scm_reference_type
  use gmtb_scm_input
  use gmtb_scm_utils
  use gmtb_scm_vgrid
  use gmtb_scm_setup
  use gmtb_scm_forcing
  use gmtb_scm_time_integration
  use gmtb_scm_output

  use            :: ccpp_types,                         &
                    only: ccpp_t
  use            :: ccpp,                               &
                    only: ccpp_init
  use            :: ccpp_fcall,                           &
                    only: ccpp_run
  use            :: ccpp_fields,                        &
                    only: ccpp_fields_add
  use ccpp_errors, only: ccpp_error

  implicit none

  type(scm_state_type), target :: scm_state
  type(scm_input_type), target :: scm_input
  type(scm_reference_type), target :: scm_reference

  integer                           :: i, grid_error !< dummy indices and error statuses

  type(ccpp_t), allocatable, target                      :: cdata(:)
  integer                                                :: ipd_index, subcycle_index, scheme_index

  integer                                                :: cdata_time_index
  integer                                                :: ierr !< Integer error flag

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

  allocate(cdata(scm_state%n_cols))

  call set_state(scm_input, scm_reference, scm_state)

  call calc_geopotential(1, scm_state)

  scm_state%model_time = 0.0
  scm_state%itt = 1

  call interpolate_forcing(scm_input, scm_state)

  call output_init(scm_state)

  scm_state%itt_out = 1

  call output_append(scm_state)

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
     call ccpp_init( &
          trim(adjustl(scm_state%physics_suite_dir))//trim(adjustl(scm_state%physics_suite_name(i)))//'.xml', &
          cdata(i), ierr)

    !use ccpp_fields.inc to call ccpp_fields_add for all variables to be exposed to CCPP (this is auto-generated from /src/ccpp/scripts/ccpp_prebuild.py - the script parses tables in gmtb_scm_type_defs.f90)
#include "ccpp_fields.inc"
  end do

  !first time step (call once)

  if (scm_state%time_scheme == 1) then
     scm_state%dt_now = scm_state%dt
     scm_state%model_time = scm_state%dt_now

     call interpolate_forcing(scm_input, scm_state)

     call calc_pres_exner_geopotential(1, scm_state)

     !pass in state variables to be modified by forcing and physics
     call do_time_step(scm_state, cdata)

  else if (scm_state%time_scheme == 2) then
  !   !if using the leapfrog scheme, we initialize by taking one half forward time step and one half (unfiltered) leapfrog time step to get to the end of the first time step
    scm_state%dt_now = 0.5*scm_state%dt
    scm_state%model_time = scm_state%dt_now

    !save initial state
    scm_state%temp_tracer(:,:,:,1) = scm_state%state_tracer(:,:,:,1)
    scm_state%temp_T(:,:,1) = scm_state%state_T(:,:,1)
    scm_state%temp_u(:,:,1) = scm_state%state_u(:,:,1)
    scm_state%temp_v(:,:,1) = scm_state%state_v(:,:,1)

    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    call apply_forcing_forward_Euler(scm_state)

    !apply_forcing_forward_Euler updates state variables time level 1, so must copy this data to time_level 2 (where cdata points)
    scm_state%state_T(:,:,2) = scm_state%state_T(:,:,1)
    scm_state%state_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,1)
    scm_state%state_u(:,:,2) = scm_state%state_u(:,:,1)
    scm_state%state_v(:,:,2) = scm_state%state_v(:,:,1)

    do i=1, scm_state%n_cols
      do ipd_index = 1 , cdata(i)%suite%ipds_max
        do subcycle_index = 1, cdata(i)%suite%ipds(ipd_index)%subcycles_max
          do scheme_index = 1, cdata(i)%suite%ipds(ipd_index)%subcycles(subcycle_index)%schemes_max
            call ccpp_run(cdata(i)%suite%ipds(ipd_index)%subcycles(subcycle_index)%schemes(scheme_index), cdata(i), ierr)
          end do !ipd parts
        end do !subcycles
      end do !schemes
    end do !columns

    !the filter routine (called after the following leapfrog time step) expects time level 2 in temp_tracer to be the updated, unfiltered state after the previous time step
    scm_state%temp_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,2)
    scm_state%temp_T(:,:,2) = scm_state%state_T(:,:,2)
    scm_state%temp_u(:,:,2) = scm_state%state_u(:,:,2)
    scm_state%temp_v(:,:,2) = scm_state%state_v(:,:,2)

    !do half a leapfrog time step to get to the end of one full time step
    scm_state%model_time = scm_state%dt
    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    !calling do_time_step with the leapfrog scheme active expects state variables in time level 1 to have values from 2 time steps ago, so set them equal to the initial values
    scm_state%state_T(:,:,1) = scm_state%temp_T(:,:,1)
    scm_state%state_u(:,:,1) = scm_state%temp_u(:,:,1)
    scm_state%state_v(:,:,1) = scm_state%temp_v(:,:,1)
    scm_state%state_tracer(:,:,:,1) = scm_state%temp_tracer(:,:,:,1)

    !go forward one leapfrog time step
    call do_time_step(scm_state, cdata)

    !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
    call filter(scm_state)

    !> \todo tracers besides water vapor do not need to be filtered (is this right?)
    scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,scm_state%cloud_water_index,2)
    scm_state%state_tracer(:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,scm_state%ozone_index,2)
  end if

  !prepare for time loop
  scm_state%n_timesteps = ceiling(scm_state%runtime/scm_state%dt)
  ! n_itt_swrad = floor(swrad_frequency/dt)
  ! n_itt_lwrad = floor(lwrad_frequency/dt)
  scm_state%n_itt_out = floor(scm_state%output_frequency/scm_state%dt)

  scm_state%dt_now = scm_state%dt

  do i = 2, scm_state%n_timesteps
    scm_state%itt = i
    !>  - Calculate the elapsed model time.
    scm_state%model_time = scm_state%itt*scm_state%dt

    !>  - Save previously unfiltered state as temporary for use in the time filter.
    if(scm_state%time_scheme == 2) then
      scm_state%temp_tracer = scm_state%state_tracer
      scm_state%temp_T = scm_state%state_T
      scm_state%temp_u = scm_state%state_u
      scm_state%temp_v = scm_state%state_v
    end if

    call interpolate_forcing(scm_input, scm_state)

    call calc_pres_exner_geopotential(1, scm_state)

    !pass in state variables to be modified by forcing and physics
    call do_time_step(scm_state, cdata)

    select case(scm_state%time_scheme)
      case (1)
        !for forward Euler scheme, no filtering is done; simply transfer output state variables from slot 2 to slot 1
        ! scm_state%state_T(:,:,1) = scm_state%state_T(:,:,2)
        ! scm_state%state_u(:,:,1) = scm_state%state_u(:,:,2)
        ! scm_state%state_v(:,:,1) = scm_state%state_v(:,:,2)
        ! scm_state%state_tracer(:,:,:,1) = scm_state%state_tracer(:,:,:,2)
      case (2)
        !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
        call filter(scm_state)

        !> \todo tracers besides water vapor do not need to be filtered (is this right?)
        scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) = scm_state%state_tracer(:,:,scm_state%cloud_water_index,2)
        scm_state%state_tracer(:,:,scm_state%ozone_index,1) = scm_state%state_tracer(:,:,scm_state%ozone_index,2)
    end select

    if(mod(scm_state%itt, scm_state%n_itt_out)==0) then
      scm_state%itt_out = scm_state%itt_out+1
      write(*,*) "itt = ",scm_state%itt
      write(*,*) "model time (s) = ",scm_state%model_time
      write(*,*) "calling output routine..."

      call output_append(scm_state)

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
