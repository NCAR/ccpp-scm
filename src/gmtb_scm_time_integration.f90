!> \file gmtb_scm_time_integration.f90
!!  Contains subroutines to handle the SCM time stepping

module gmtb_scm_time_integration

use gmtb_scm_kinds, only: sp, dp, qp
!use nuopc_physics
use gmtb_scm_forcing

use            :: ccpp_types, only: ccpp_t
use            :: ccpp_ipd, only: ccpp_ipd_run

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup time_integration gmtb_scm_time_integration
!! @{
!! Contains subroutines to handle the SCM time stepping.

!> This subroutine performs the Robert-Asselin time filtering of the state variables.
subroutine filter(c_filter, temp_tracer, temp_T, temp_u, temp_v, state_tracer, state_T, state_u, state_v)
  real(kind=dp), intent(in)     :: c_filter !< filter constant (equal to 1 - the constant used in the GFS)
  real(kind=dp), intent(in)     :: temp_tracer(:,:,:,:) !< temporary copy of the state tracer variables before the latest forcing/physics update (horizontal, vertical, tracer_index, iteration)
  real(kind=dp), intent(in)     :: temp_T(:,:,:) !< temporary copy of the state temperature before the latest forcing/physics update (horizontal, vertical, iteration)
  real(kind=dp), intent(in)     :: temp_u(:,:,:) !< temporary copy of the state zonal wind before the latest forcing/physics update (horizontal, vertical, iteration)
  real(kind=dp), intent(in)     :: temp_v(:,:,:) !< temporary copy of the state meridional wind before the latest forcing/physics update (horizontal, vertical, iteration)
  real(kind=dp), intent(inout)  :: state_tracer(:,:,:,:) !< state tracer variables containing the section of the array to be updated and the section that contains the latest unfiltered values (horizontal, vertical, tracer_index, iteration)
  real(kind=dp), intent(inout)  :: state_T(:,:,:) !< state temperature containing the section of the array to be updated and the section that contains the latest unfiltered values (horizontal, vertical, iteration)
  real(kind=dp), intent(inout)  :: state_u(:,:,:) !< state zonal wind containing the section of the array to be updated and the section that contains the latest unfiltered values (horizontal, vertical, iteration)
  real(kind=dp), intent(inout)  :: state_v(:,:,:) !< state merdional wind containing the section of the array to be updated and the section that contains the latest unfiltered values (horizontal, vertical, iteration)

  !> \section filter_alg Algorithm
  !! The filtered state variables are calculated using
  !! \f[
  !! \overline{x^\tau}=(1-c)x^\tau + 0.5c\left(x^{\tau +1} + \overline{x^{\tau - 1}}\right)
  !! \f]
  !! where \f$\overline{x^\tau}\f$ is the filtered value of variable \f$x\f$ at the current iteration, \f$x^\tau\f$ is the unfiltered value of the previous time step, \f$x^{\tau +1}\f$ is the unfiltered
  !! value that was just updated by the forcing and physics, and \f$\overline{x^{\tau - 1}}\f$ is the filtered value of the variable from the previous iteration, and \f$c\f$ is the filtering constant.
  state_tracer(:,:,1,1) = (1.0 - c_filter)*temp_tracer(:,:,1,2) + 0.5*c_filter*(state_tracer(:,:,1,2) + temp_tracer(:,:,1,1))
  state_T(:,:,1) = (1.0 - c_filter)*temp_T(:,:,2) + 0.5*c_filter*(state_T(:,:,2) + temp_T(:,:,1))
  state_u(:,:,1) = (1.0 - c_filter)*temp_u(:,:,2) + 0.5*c_filter*(state_u(:,:,2) + temp_u(:,:,1))
  state_v(:,:,1) = (1.0 - c_filter)*temp_v(:,:,2) + 0.5*c_filter*(state_v(:,:,2) + temp_v(:,:,1))

end subroutine

!> This subroutine calls nuopc_rad_update and nuopc_rad_run in nuopc_physics.F90 (if necessary) and apply_forcing_leapfrog from \ref forcing and nuopc_phys_run, also from nuopc_physics.F90.
!! The subroutine nuopc_rad_update calculates the time-dependent parameters required to run radiation, and nuopc_rad_run calculates the radiative heating rate (but does not apply it). The
!! subroutine apply_forcing_leapfrog advances the state variables forward using the leapfrog method and nuopc_phys_run further changes the state variables using the forward method. By the end of
!! this subroutine, the unfiltered state variables will have been stepped forward in time.
! subroutine do_time_step(n_levels, n_columns, time_scheme, state_fldin, state_fldout, sfc_prop, diags, intrfc_fld, cld_prop, &
!   rad_tend, mdl_parm, tbddata, dyn_parm, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
!   thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend, v_force_tend, T_force_tend, &
!   qv_force_tend, exner_l, phii, pres_l, lat, dt, lsswr, lslwr, thermo_forcing_type, mom_forcing_type, relax_time)
subroutine do_time_step(n_levels, n_columns, time_scheme, &
  state_tracer, state_T, state_u, state_v, cdata, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge,&
  thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend, v_force_tend, T_force_tend, &
  qv_force_tend, exner_l, phii, pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time)

  integer, intent(in)                    :: n_levels !< number of model levels
  integer, intent(in)                    :: n_columns !< number of model columns
  integer, intent(in)                    :: time_scheme !< time-stepping scheme choice (1=forward, 2=filtered leapfrog)
  ! type(model_parameters), intent(in)     :: mdl_parm !< NUOPC DDT containing model configuration parameters
  ! type(dynamic_parameters), intent(inout)   :: dyn_parm !< NUOPC DDT containing parameters that change throughout the integration often
  ! type(state_fields_in), intent(in)      :: state_fldin !< NUOPC DDT containing the model state variables before the time step
  ! type(sfc_properties), intent(in)       :: sfc_prop !< NUOPC DDT containing surface variables
  ! type(diagnostics), intent(in)          :: diags !< NUOPC DDT containing diagnostic variables
  ! type(cloud_properties), intent(in)     :: cld_prop !< NUOPC DDT containing cloud property variables
  ! type(radiation_tendencies), intent(in) :: rad_tend !< NUOPC DDT containing radiaitve tendencies
  ! type(interface_fields), intent(in)     :: intrfc_fld !< NUOPC DDT containing interfacial variables
  ! type(state_fields_out), intent(inout)  :: state_fldout !< NUOPC DDT that are updated during the time step
  ! type(tbd_ddt), intent(in)              :: tbddata !< NUOPC DDT containing unclassified variables

  real(kind=dp), intent(inout)                 :: state_tracer(:,:,:,:) !< model state tracers
  real(kind=dp), intent(inout)                 :: state_T(:,:,:) !< model state temperature
  real(kind=dp), intent(inout)                 :: state_u(:,:,:) !< model state zonal wind
  real(kind=dp), intent(inout)                 :: state_v(:,:,:) !< model state meridional wind
  type(ccpp_t), intent(inout)                  :: cdata(:)
  real(kind=dp), intent(in)                    :: w_ls(:,:) !< large scale w (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: omega(:,:) !< large scale pressure velocity (Pa/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: u_g(:,:) !< geostrophic zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: v_g(:,:) !< geostrophic meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: u_nudge(:,:) !< nudging zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: v_nudge(:,:) !< nudging meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: T_nudge(:,:) !< nudging absolute temperature (K) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: thil_nudge(:,:) !< nudging liq. pot. temperatuer (K) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: qt_nudge(:,:) !< nudging specific humidity (K) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: dT_dt_rad(:,:) !< radiative heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: h_advec_thil(:,:) !< change in theta_il due to horizontal advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: h_advec_qt(:,:) !< change in q_t due to horizontal advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: v_advec_thil(:,:) !< change in theta_il due to vertical advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: v_advec_qt(:,:) !< change in q_t due to vertical advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in)                    :: lat(:) !< latitude in radians
  real(kind=dp), intent(in)                    :: dt !< time step (s)
  real(kind=dp), intent(in)                    :: exner_l(:,:) !< exner function at model level centers (horizontal, vertical)
  real(kind=dp), intent(in)                    :: phii(:,:) !< geopotential at model level interfaces (horizontal, vertical)
  real(kind=dp), intent(in)                    :: pres_l(:,:) !< pressure at model level centers (horizontal, vertical)
  real(kind=dp), intent(in)                    :: relax_time !< time scale for relaxation to observations (s)
  ! logical, intent(in)                    :: lsswr !< logical flag for calling shortwave radiation
  ! logical, intent(in)                    :: lslwr !< logical flag for calling longwave radiation
  integer, intent(in)                    :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer, intent(in)                    :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"

  real(kind=dp), intent(out)                   :: u_force_tend(:,:) !< sum of all forcing for zonal momentum (horizontal, vertical)
  real(kind=dp), intent(out)                   :: v_force_tend(:,:) !< sum of all forcing for meridional momentum (horizontal, vertical)
  real(kind=dp), intent(out)                   :: T_force_tend(:,:) !< sum of all forcing for temperature (horizontal, vertical)
  real(kind=dp), intent(out)                   :: qv_force_tend(:,:) !< sum of all forcing for water vapor (horizontal, vertical)

  integer :: i, ipd_index, subcycle_index, scheme_index, ierr

  !> \section do_time_step_alg Algorithm
  !! @{
  !! - Call nuopc_rad_update and nuopr_rad_run as necessary. State variables are not affected by these calls. Only radiative heating rates (and other rad. diagnostics) are returned.
  ! if(lsswr .or. lslwr) then
    ! call nuopc_rad_update (mdl_parm, dyn_parm)
    !
    ! call nuopc_rad_run (state_fldin, sfc_prop, diags,intrfc_fld,&
    !         cld_prop, rad_tend, mdl_parm, dyn_parm)
  ! end if

  !> - Call apply_forcing_* from \ref forcing. This routine updates the "input" state variables for the physics call (updates filtered values from previous timestep, if leapfrog scheme). It effectively replaces the change of the state variables due to dynamics.
  select case(time_scheme)
    case(1)
      call apply_forcing_forward_Euler(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, &
        u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, &
        phii, pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, &
        qv_force_tend)
    case(2)
      ! call apply_forcing_forward_Euler(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, &
      !   u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, &
      !   phii, pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, &
      !   qv_force_tend)
      call apply_forcing_leapfrog(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, u_nudge, &
        v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, phii, &
        pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, qv_force_tend)
    case default
      call apply_forcing_forward_Euler(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, &
        u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, &
        phii, pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, &
        qv_force_tend)
  end select

  if (time_scheme == 2) then
    !IPD cdata points to time level 2 for updating state variables; update time level 2 state variables with those where the forcing has been applied this time step
    state_T(:,:,2) = state_T(:,:,1)
    state_tracer(:,:,:,2) = state_tracer(:,:,:,1)
    state_u(:,:,2) = state_u(:,:,1)
    state_v(:,:,2) = state_v(:,:,1)
  end if

  do i=1, n_columns
    do ipd_index = 1 , cdata(i)%suite%ipds_max
      do subcycle_index = 1, cdata(i)%suite%ipds(ipd_index)%subcycles_max
        do scheme_index = 1, cdata(i)%suite%ipds(ipd_index)%subcycles(subcycle_index)%schemes_max
          call ccpp_ipd_run(cdata(i)%suite%ipds(ipd_index)%subcycles(subcycle_index)%schemes(scheme_index), cdata(i), ierr)
        end do !ipd parts
      end do !subcycles
    end do !schemes
  end do !columns


  !> - Call nuopc_phys_run to complete the forward time step. The unfiltered, updated state variables are in the NUOPC DDT state_fldout, which points to state_var(:,:,2).
  ! call nuopc_phys_run(state_fldin, state_fldout, sfc_prop,&
  !           diags, intrfc_fld, cld_prop, rad_tend,&
  !           mdl_parm, tbddata, dyn_parm )

  !if no physics call, need to transfer state_variables(:,:,1) to state_variables (:,:,2)
  ! state_T(:,:,2) = state_T(:,:,1)
  ! state_tracer(:,:,:,2) = state_tracer(:,:,:,1)
  ! state_u(:,:,2) = state_u(:,:,1)
  ! state_v(:,:,2) = state_v(:,:,1)

  !> @}
end subroutine

!> @}
!> @}
end module gmtb_scm_time_integration
