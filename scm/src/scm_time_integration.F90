!> \file scm_time_integration.f90
!!  Contains subroutines to handle the SCM time stepping

module scm_time_integration

use scm_kinds, only: sp, dp, qp
use scm_forcing

use ccpp_api,        only: ccpp_t
use ccpp_static_api, only: ccpp_physics_run

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
subroutine do_time_step(scm_state, physics, cdata)
  use scm_type_defs, only: scm_state_type, physics_type

  type(scm_state_type), intent(inout)          :: scm_state
  type(physics_type), intent(inout)            :: physics
  type(ccpp_t), intent(inout)                  :: cdata

  integer :: i, ierr

  !> \section do_time_step_alg Algorithm
  !! @{

  !> - Call apply_forcing_* from \ref forcing. This routine updates the "input" state variables for the physics call (updates filtered values from previous timestep, if leapfrog scheme). It effectively replaces the change of the state variables due to dynamics.
  select case(scm_state%time_scheme)
    case(1)
      call apply_forcing_forward_Euler(scm_state)
    case(2)
      call apply_forcing_leapfrog(scm_state)
    case default
      call apply_forcing_forward_Euler(scm_state)
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
      physics%Diag%du3dt(i,:,8)  = physics%Diag%du3dt(i,:,8)  &
                                      + (physics%Statein%ugrs(i,:) - physics%Stateout%gu0(i,:))
      physics%Diag%dv3dt(i,:,8)  =   physics%Diag%dv3dt(i,:,8)  &
                                      + (  physics%Statein%vgrs(i,:) - physics%Stateout%gv0(i,:))
      physics%Diag%dt3dt(i,:,11) = physics%Diag%dt3dt(i,:,11) &
                                      + (physics%Statein%tgrs(i,:) - physics%Stateout%gt0(i,:))
      if (physics%Model%qdiag3d) then
        physics%Diag%dq3dt(i,:,12) = physics%Diag%dq3dt(i,:,12) &
              + (physics%Statein%qgrs(i,:,physics%Model%ntqv) - physics%Stateout%gq0(i,:,physics%Model%ntqv))
        physics%Diag%dq3dt(i,:,13) = physics%Diag%dq3dt(i,:,13) &
              + (physics%Statein%qgrs(i,:,physics%Model%ntoz) - physics%Stateout%gq0(i,:,physics%Model%ntoz))
      endif
    endif
  end do

  call ccpp_physics_run(cdata, suite_name=trim(adjustl(scm_state%physics_suite_name)), ierr=ierr)
  if (ierr/=0) then
      write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_run: ' // trim(cdata%errmsg) // '. Exiting...'
      stop
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
