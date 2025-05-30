!> \file scm_time_integration.f90
!!  Contains subroutines to handle the SCM time stepping

module scm_time_integration

use iso_fortran_env, only: error_unit
use scm_kinds, only: sp, dp, qp
use scm_forcing

use ccpp_types,        only: ccpp_t
use :: ccpp_static_api,                      &
       only: ccpp_physics_timestep_init,     &
             ccpp_physics_run,               &
             ccpp_physics_timestep_finalize


implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup time_integration scm_time_integration
!! @{
!! Contains subroutines to handle the SCM time stepping.

!> This subroutine performs the Robert-Asselin time filtering of the state variables.
subroutine filter(scm_state)
  use scm_type_defs, only: scm_state_type

  type(scm_state_type), intent(inout)          :: scm_state

  !> \section filter_alg Algorithm
  !! The filtered state variables are calculated using
  !! \f[
  !! \overline{x^\tau}=(1-c)x^\tau + 0.5c\left(x^{\tau +1} + \overline{x^{\tau - 1}}\right)
  !! \f]
  !! where \f$\overline{x^\tau}\f$ is the filtered value of variable \f$x\f$ at the current iteration, \f$x^\tau\f$ is the unfiltered value of the previous time step, \f$x^{\tau +1}\f$ is the unfiltered
  !! value that was just updated by the forcing and physics, and \f$\overline{x^{\tau - 1}}\f$ is the filtered value of the variable from the previous iteration, and \f$c\f$ is the filtering constant.
  scm_state%state_tracer(:,:,scm_state%water_vapor_index,1) = &
    (1.0 - scm_state%c_filter)*scm_state%temp_tracer(:,:,scm_state%water_vapor_index,2) + &
    0.5*scm_state%c_filter*(scm_state%state_tracer(:,:,scm_state%water_vapor_index,2) + &
    scm_state%temp_tracer(:,:,scm_state%water_vapor_index,1))
  scm_state%state_T(:,:,1) = (1.0 - scm_state%c_filter)*scm_state%temp_T(:,:,2) + &
    0.5*scm_state%c_filter*(scm_state%state_T(:,:,2) + scm_state%temp_T(:,:,1))
  scm_state%state_u(:,:,1) = (1.0 - scm_state%c_filter)*scm_state%temp_u(:,:,2) + &
    0.5*scm_state%c_filter*(scm_state%state_u(:,:,2) + scm_state%temp_u(:,:,1))
  scm_state%state_v(:,:,1) = (1.0 - scm_state%c_filter)*scm_state%temp_v(:,:,2) + &
    0.5*scm_state%c_filter*(scm_state%state_v(:,:,2) + scm_state%temp_v(:,:,1))

end subroutine

!> This subroutine calls nuopc_rad_update and nuopc_rad_run in nuopc_physics.F90 (if necessary) and apply_forcing_leapfrog from \ref forcing and nuopc_phys_run, also from nuopc_physics.F90.
!! The subroutine nuopc_rad_update calculates the time-dependent parameters required to run radiation, and nuopc_rad_run calculates the radiative heating rate (but does not apply it). The
!! subroutine apply_forcing_leapfrog advances the state variables forward using the leapfrog method and nuopc_phys_run further changes the state variables using the forward method. By the end of
!! this subroutine, the unfiltered state variables will have been stepped forward in time.
subroutine do_time_step(scm_state, physics, cdata, in_spinup)
  use scm_type_defs, only: scm_state_type, physics_type

  type(scm_state_type), intent(inout)          :: scm_state
  type(physics_type), intent(inout)            :: physics
  type(ccpp_t), intent(inout)                  :: cdata
  logical, intent(in)                          :: in_spinup

  integer :: i, ierr, kdt_rad, idtend, itrac

  !> \section do_time_step_alg Algorithm
  !! @{

  !> - Call apply_forcing_* from \ref forcing. This routine updates the "input" state variables for the physics call (updates filtered values from previous timestep, if leapfrog scheme). It effectively replaces the change of the state variables due to dynamics.
  select case(scm_state%time_scheme)
    case(1)
      if (scm_state%input_type == 0) then
        call apply_forcing_forward_Euler(scm_state, in_spinup)
      else
        call apply_forcing_DEPHY(scm_state, in_spinup)
      end if
    case(2)
      if (scm_state%input_type == 0) then
        call apply_forcing_leapfrog(scm_state)
      else
        error stop 'The application of forcing terms from the DEPHY file format has not been implemented for the leapfrog time scheme.'
      end if

    case default
      if (scm_state%input_type == 0) then
        call apply_forcing_forward_Euler(scm_state, in_spinup)
      else
        call apply_forcing_DEPHY(scm_state, in_spinup)
      end if
  end select

  if (scm_state%time_scheme == 2) then
    !IPD cdata points to time level 2 for updating state variables; update time level 2 state variables with those where the forcing has been applied this time step
    scm_state%state_T(:,:,2) = scm_state%state_T(:,:,1)
    scm_state%state_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,1)
    scm_state%state_u(:,:,2) = scm_state%state_u(:,:,1)
    scm_state%state_v(:,:,2) = scm_state%state_v(:,:,1)
  end if

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
  enddo

  call ccpp_physics_timestep_init(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
  if (ierr/=0) then
      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_init: ' // trim(cdata%errmsg) // '. Exiting...'
      error stop trim(cdata%errmsg)
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

  !CCPP run phase
  ! time_vary group doesn't have any run phase (omitted)
  ! radiation group
  call physics%Interstitial(1)%rad_reset(physics%Model)
  call ccpp_physics_run(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), group_name="radiation", ierr=ierr)
  if (ierr/=0) then
      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group radiation: ' // trim(cdata%errmsg) // '. Exiting...'
      error stop
  end if
  ! process-split physics
  call physics%Interstitial(1)%phys_reset(physics%Model)
  call ccpp_physics_run(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), group_name="phys_ps", ierr=ierr)
  if (ierr/=0) then
      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ps: ' // trim(cdata%errmsg) // '. Exiting...'
      error stop
  end if
  ! time-split physics
  call ccpp_physics_run(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), group_name="phys_ts", ierr=ierr)
  if (ierr/=0) then
      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_run for group phys_ts: ' // trim(cdata%errmsg) // '. Exiting...'
      error stop
  end if

  call ccpp_physics_timestep_finalize(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
  if (ierr/=0) then
      write(error_unit,'(a,i0,a)') 'An error occurred in ccpp_physics_timestep_finalize: ' // trim(cdata%errmsg) // '. Exiting...'
      error stop trim(cdata%errmsg)
  end if

  !if no physics call, need to transfer state_variables(:,:,1) to state_variables (:,:,2)
  ! scm_state%state_T(:,:,2) = scm_state%state_T(:,:,1)
  ! scm_state%state_tracer(:,:,:,2) = scm_state%state_tracer(:,:,:,1)
  ! scm_state%state_u(:,:,2) = scm_state%state_u(:,:,1)
  ! scm_state%state_v(:,:,2) = scm_state%state_v(:,:,1)

  !> @}
end subroutine

!> @}
!> @}
end module scm_time_integration
