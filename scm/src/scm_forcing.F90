!> \file scm_forcing.f90
!!  Contains subroutines to handle the SCM forcing -- interpolating in space and time, etc.

module scm_forcing

use scm_kinds, only: sp, dp, qp
use scm_utils, only: interpolate_to_grid_centers

use scm_physical_constants, only: con_pi, con_omega, con_g, con_cp, con_rd

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup forcing scm_forcing
!! @{
!! Contains subroutines to handle the SCM forcing -- interpolating in space and time, etc.

!> This subroutine interpolates the model forcing, column position, and surface properties to the current model time and to the model grid.
!! \note The input forcing file contains forcing for one column, yet the SCM code can accommodate more than one column. Right now, all columns are assumed to have the same forcing.
subroutine interpolate_forcing(scm_input, scm_state)
  use scm_type_defs, only: scm_input_type, scm_state_type

  type(scm_input_type), intent(in) :: scm_input
  type(scm_state_type), intent(inout) :: scm_state

  integer :: i, n
  integer :: low_t_index, top_index !< index of the time in the input file immediately before the current model time, index of the last calculated level
  real(kind=dp) :: lifrac

  real(kind=dp) :: w_ls_bracket(2,scm_state%n_levels), omega_bracket(2,scm_state%n_levels), u_g_bracket(2,scm_state%n_levels), &
    v_g_bracket(2,scm_state%n_levels), u_nudge_bracket(2,scm_state%n_levels), v_nudge_bracket(2,scm_state%n_levels), &
    T_nudge_bracket(2,scm_state%n_levels), thil_nudge_bracket(2,scm_state%n_levels), qt_nudge_bracket(2,scm_state%n_levels), &
    dT_dt_rad_bracket(2,scm_state%n_levels), h_advec_thil_bracket(2,scm_state%n_levels),h_advec_qt_bracket(2,scm_state%n_levels), &
    v_advec_thil_bracket(2,scm_state%n_levels), v_advec_qt_bracket(2,scm_state%n_levels) !< forcing terms that "bracket" around the model time

  !> \section interpolate_forcing_alg Algorithm
  !! @{

  !> - Check for the case where the elapsed model time extends beyond the supplied forcing.
  if(scm_state%model_time >= scm_input%input_time(scm_input%input_ntimes)) then
    !>  - If so, hold the forcing terms constant at the last supplied values. The forcing still needs to be interpolated to the grid.
    write(*,*) "The model_time has exceeded the specifed period of forcing. Forcing will now be held constant at the last &
      specified values."

      !>  - For all forcing terms, call interpolate_to_grid_centers from \ref utils for each variable. This subroutine returns the last vertical index calculated in case forcing terms above the case input needs to be specified.
      do i=1, scm_state%n_cols
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_w_ls(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          w_ls_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_omega(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          omega_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_u_g(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, u_g_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_v_g(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, v_g_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_u_nudge(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          u_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_v_nudge(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          v_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_T_nudge(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          T_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_thil_nudge(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          thil_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_qt_nudge(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          qt_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_dT_dt_rad(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          dT_dt_rad_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_h_advec_thetail(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          h_advec_thil_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_h_advec_qt(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          h_advec_qt_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_v_advec_thetail(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          v_advec_thil_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
          scm_input%input_v_advec_qt(scm_input%input_ntimes,:), scm_state%pres_l(i,:), scm_state%n_levels, &
          v_advec_qt_bracket(1,:), top_index, 3)

        !>  - If the input forcing file does not reach to the model domain top, fill in values above the input forcing file domain with those from the top level.
        if (top_index < scm_state%n_levels .AND. top_index.GT.0) then
          w_ls_bracket(1,top_index+1:scm_state%n_levels) = 0.0!w_ls_bracket(1,top_index)
          omega_bracket(1,top_index+1:scm_state%n_levels) = 0.0!omega_bracket(1,top_index)
          u_g_bracket(1,top_index+1:scm_state%n_levels) = u_g_bracket(1,top_index)
          v_g_bracket(1,top_index+1:scm_state%n_levels) = v_g_bracket(1,top_index)
          u_nudge_bracket(1,top_index+1:scm_state%n_levels) = u_nudge_bracket(1,top_index)
          v_nudge_bracket(1,top_index+1:scm_state%n_levels) = v_nudge_bracket(1,top_index)
          T_nudge_bracket(1,top_index+1:scm_state%n_levels) = T_nudge_bracket(1,top_index)
          thil_nudge_bracket(1,top_index+1:scm_state%n_levels) = thil_nudge_bracket(1,top_index)
          qt_nudge_bracket(1,top_index+1:scm_state%n_levels) = qt_nudge_bracket(1,top_index)
          dT_dt_rad_bracket(1,top_index+1:scm_state%n_levels) = dT_dt_rad_bracket(1,top_index)
          h_advec_thil_bracket(1,top_index+1:scm_state%n_levels) = h_advec_thil_bracket(1,top_index)
          h_advec_qt_bracket(1,top_index+1:scm_state%n_levels) = h_advec_qt_bracket(1,top_index)
          v_advec_thil_bracket(1,top_index+1:scm_state%n_levels) = v_advec_thil_bracket(1,top_index)
          v_advec_qt_bracket(1,top_index+1:scm_state%n_levels) = v_advec_qt_bracket(1,top_index)
        end if

        !>  - For this case, no time interpolation is necessary; just set the forcing terms to the vertically-interpolated values.
        scm_state%w_ls(i,:) = w_ls_bracket(1,:)
        scm_state%omega(i,:) = omega_bracket(1,:)
        scm_state%u_g(i,:) = u_g_bracket(1,:)
        scm_state%v_g(i,:) = v_g_bracket(1,:)
        scm_state%u_nudge(i,:) = u_nudge_bracket(1,:)
        scm_state%v_nudge(i,:) = v_nudge_bracket(1,:)
        scm_state%T_nudge(i,:) = T_nudge_bracket(1,:)
        scm_state%thil_nudge(i,:) = thil_nudge_bracket(1,:)
        scm_state%qt_nudge(i,:) = qt_nudge_bracket(1,:)
        scm_state%dT_dt_rad(i,:) = dT_dt_rad_bracket(1,:)
        scm_state%h_advec_thil(i,:) = h_advec_thil_bracket(1,:)
        scm_state%h_advec_qt(i,:) = h_advec_qt_bracket(1,:)
        scm_state%v_advec_thil(i,:) = v_advec_thil_bracket(1,:)
        scm_state%v_advec_qt(i,:) = v_advec_qt_bracket(1,:)

        !>  - Set the surface parameters to the last available data.
        scm_state%pres_surf(i) = scm_input%input_pres_surf(scm_input%input_ntimes)
        scm_state%T_surf(i) = scm_input%input_T_surf(scm_input%input_ntimes)
        scm_state%sh_flux(i) = scm_input%input_sh_flux_sfc(scm_input%input_ntimes)
        scm_state%lh_flux(i) = scm_input%input_lh_flux_sfc(scm_input%input_ntimes)
      end do
  else
  !> - When the model elapsed time is within the time-frame specified by the input forcing, the forcing must be interpolated in time and space.
    !>  - Determine the index in the input file for the time immediately preceeding the model time and determine the linear interpolation value.
    low_t_index = 0
    do n=1, scm_input%input_ntimes
      if (scm_input%input_time(n) > scm_state%model_time) then
        low_t_index = n-1
        lifrac = (scm_state%model_time - scm_input%input_time(low_t_index))/&
          (scm_input%input_time(low_t_index+1) - scm_input%input_time(low_t_index))
        exit
      end if
    end do

    do i=1, scm_state%n_cols
      !>  - For all forcing terms, call interpolate_to_grid_centers from \ref utils for each variable for each time level that "bracket" around
      !>    the current model time. This subroutine returns the last vertical index calculated in case forcing terms above the case input needs
      !>    to be specified.
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_w_ls(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, w_ls_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_w_ls(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, w_ls_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_omega(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, omega_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_omega(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, omega_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_u_g(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, u_g_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_u_g(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, u_g_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_g(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_g_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_g(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_g_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_u_nudge(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, u_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_u_nudge(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, u_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_nudge(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_nudge(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_T_nudge(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, T_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_T_nudge(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, T_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_thil_nudge(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, thil_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_thil_nudge(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, thil_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_qt_nudge(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, qt_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_qt_nudge(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, qt_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_dT_dt_rad(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, dT_dt_rad_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_dT_dt_rad(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, dT_dt_rad_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_h_advec_thetail(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, h_advec_thil_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
        scm_input%input_h_advec_thetail(low_t_index+1,:), scm_state%pres_l(i,:), scm_state%n_levels, &
        h_advec_thil_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_h_advec_qt(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, h_advec_qt_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_h_advec_qt(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, h_advec_qt_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_advec_thetail(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_advec_thil_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, &
        scm_input%input_v_advec_thetail(low_t_index+1,:), scm_state%pres_l(i,:), scm_state%n_levels, &
        v_advec_thil_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_advec_qt(low_t_index,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_advec_qt_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v_advec_qt(low_t_index+1,:), &
        scm_state%pres_l(i,:), scm_state%n_levels, v_advec_qt_bracket(2,:), top_index, 3)

      !>  - If the input forcing file does not reach to the model domain top, fill in values above the input forcing file domain with those from the top level.
      if (top_index < scm_state%n_levels) then
        w_ls_bracket(1,top_index+1:scm_state%n_levels) = 0.0!w_ls_bracket(1,top_index)
        w_ls_bracket(2,top_index+1:scm_state%n_levels) = 0.0!w_ls_bracket(2,top_index)
        omega_bracket(1,top_index+1:scm_state%n_levels) = 0.0!omega_bracket(1,top_index)
        omega_bracket(2,top_index+1:scm_state%n_levels) = 0.0!omega_bracket(2,top_index)
        u_g_bracket(1,top_index+1:scm_state%n_levels) = u_g_bracket(1,top_index)
        u_g_bracket(2,top_index+1:scm_state%n_levels) = u_g_bracket(2,top_index)
        v_g_bracket(1,top_index+1:scm_state%n_levels) = v_g_bracket(1,top_index)
        v_g_bracket(2,top_index+1:scm_state%n_levels) = v_g_bracket(2,top_index)
        u_nudge_bracket(1,top_index+1:scm_state%n_levels) = u_nudge_bracket(1,top_index)
        u_nudge_bracket(2,top_index+1:scm_state%n_levels) = u_nudge_bracket(2,top_index)
        v_nudge_bracket(1,top_index+1:scm_state%n_levels) = v_nudge_bracket(1,top_index)
        v_nudge_bracket(2,top_index+1:scm_state%n_levels) = v_nudge_bracket(2,top_index)
        T_nudge_bracket(1,top_index+1:scm_state%n_levels) = T_nudge_bracket(1,top_index)
        T_nudge_bracket(2,top_index+1:scm_state%n_levels) = T_nudge_bracket(2,top_index)
        thil_nudge_bracket(1,top_index+1:scm_state%n_levels) = thil_nudge_bracket(1,top_index)
        thil_nudge_bracket(2,top_index+1:scm_state%n_levels) = thil_nudge_bracket(2,top_index)
        qt_nudge_bracket(1,top_index+1:scm_state%n_levels) = qt_nudge_bracket(1,top_index)
        qt_nudge_bracket(2,top_index+1:scm_state%n_levels) = qt_nudge_bracket(2,top_index)
        dT_dt_rad_bracket(1,top_index+1:scm_state%n_levels) = dT_dt_rad_bracket(1,top_index)
        dT_dt_rad_bracket(2,top_index+1:scm_state%n_levels) = dT_dt_rad_bracket(2,top_index)
        h_advec_thil_bracket(1,top_index+1:scm_state%n_levels) = h_advec_thil_bracket(1,top_index)
        h_advec_thil_bracket(2,top_index+1:scm_state%n_levels) = h_advec_thil_bracket(2,top_index)
        h_advec_qt_bracket(1,top_index+1:scm_state%n_levels) = h_advec_qt_bracket(1,top_index)
        h_advec_qt_bracket(2,top_index+1:scm_state%n_levels) = h_advec_qt_bracket(2,top_index)
        v_advec_thil_bracket(1,top_index+1:scm_state%n_levels) = v_advec_thil_bracket(1,top_index)
        v_advec_thil_bracket(2,top_index+1:scm_state%n_levels) = v_advec_thil_bracket(2,top_index)
        v_advec_qt_bracket(1,top_index+1:scm_state%n_levels) = v_advec_qt_bracket(1,top_index)
        v_advec_qt_bracket(2,top_index+1:scm_state%n_levels) = v_advec_qt_bracket(2,top_index)
      end if

      !>  - Interpolate the forcing terms in time.
      scm_state%w_ls(i,:) = (1.0 - lifrac)*w_ls_bracket(1,:) + lifrac*w_ls_bracket(2,:)
      scm_state%omega(i,:) = (1.0 - lifrac)*omega_bracket(1,:) + lifrac*omega_bracket(2,:)
      scm_state%u_g(i,:) = (1.0 - lifrac)*u_g_bracket(1,:) + lifrac*u_g_bracket(2,:)
      scm_state%v_g(i,:) = (1.0 - lifrac)*v_g_bracket(1,:) + lifrac*v_g_bracket(2,:)
      scm_state%u_nudge(i,:) = (1.0 - lifrac)*u_nudge_bracket(1,:) + lifrac*u_nudge_bracket(2,:)
      scm_state%v_nudge(i,:) = (1.0 - lifrac)*v_nudge_bracket(1,:) + lifrac*v_nudge_bracket(2,:)
      scm_state%T_nudge(i,:) = (1.0 - lifrac)*T_nudge_bracket(1,:) + lifrac*T_nudge_bracket(2,:)
      scm_state%thil_nudge(i,:) = (1.0 - lifrac)*thil_nudge_bracket(1,:) + lifrac*thil_nudge_bracket(2,:)
      scm_state%qt_nudge(i,:) = (1.0 - lifrac)*qt_nudge_bracket(1,:) + lifrac*qt_nudge_bracket(2,:)
      scm_state%dT_dt_rad(i,:) = (1.0 - lifrac)*dT_dt_rad_bracket(1,:) + lifrac*dT_dt_rad_bracket(2,:)
      scm_state%h_advec_thil(i,:) = (1.0 - lifrac)*h_advec_thil_bracket(1,:) + lifrac*h_advec_thil_bracket(2,:)
      scm_state%h_advec_qt(i,:) = (1.0 - lifrac)*h_advec_qt_bracket(1,:) + lifrac*h_advec_qt_bracket(2,:)
      scm_state%v_advec_thil(i,:) = (1.0 - lifrac)*v_advec_thil_bracket(1,:) + lifrac*v_advec_thil_bracket(2,:)
      scm_state%v_advec_qt(i,:) = (1.0 - lifrac)*v_advec_qt_bracket(1,:) + lifrac*v_advec_qt_bracket(2,:)

      !>  - Interpolate the surface parameters in time.
      scm_state%pres_surf(i) = (1.0 - lifrac)*scm_input%input_pres_surf(low_t_index) + &
        lifrac*scm_input%input_pres_surf(low_t_index+1)
      scm_state%T_surf(i) = (1.0 - lifrac)*scm_input%input_T_surf(low_t_index) + lifrac*scm_input%input_T_surf(low_t_index+1)
      scm_state%sh_flux(i) = (1.0 - lifrac)*scm_input%input_sh_flux_sfc(low_t_index) + &
        lifrac*scm_input%input_sh_flux_sfc(low_t_index+1)
      scm_state%lh_flux(i) = (1.0 - lifrac)*scm_input%input_lh_flux_sfc(low_t_index) + &
        lifrac*scm_input%input_lh_flux_sfc(low_t_index+1)
    end do
  end if

  !> @}
end subroutine interpolate_forcing

!> This subroutine updates the state variables due to the model forcing using the leapfrog time integration scheme.
!! It overwrites the state variable arrays where the filtered values from the previous time step are kept. These parts of the state
!! variable arrays are pointed to by the state_fields_in DDT and are later updated by calling nuopc_phys_run.
subroutine apply_forcing_leapfrog(scm_state)
  use scm_type_defs, only: scm_state_type

  type(scm_state_type), intent(inout) :: scm_state

  real(kind=dp) :: old_u(scm_state%n_cols, scm_state%n_levels), old_v(scm_state%n_cols, scm_state%n_levels), &
    old_T(scm_state%n_cols, scm_state%n_levels), old_qv(scm_state%n_cols, scm_state%n_levels)
  real(kind=dp) :: w_ls_i(scm_state%n_cols, scm_state%n_levels+1), zi(scm_state%n_cols, scm_state%n_levels+1)
  real(kind=dp) :: theta(scm_state%n_cols, scm_state%n_levels)

  integer :: i,k
  real(kind=dp) :: f_coriolis, grav_inv, g_over_cp, omega_plus, omega_minus, dth_dp_plus, dth_dp_minus, &
    dqv_dp_plus, dqv_dp_minus

  !> \section apply_leapfrog_forcing_alg Algorithm
  !! @{

  grav_inv = 1.0/con_g
  g_over_cp = con_g/con_cp

  !> - Save old state variables (filtered from previous time step)
  old_u = scm_state%state_u(:,:,1)
  old_v = scm_state%state_v(:,:,1)
  old_T = scm_state%state_T(:,:,1)
  old_qv = scm_state%state_tracer(:,:,scm_state%water_vapor_index,1)

  theta = old_T/scm_state%exner_l(:,:)

  !> - Initialize forcing sums to zero.
  scm_state%u_force_tend = 0.0
  scm_state%v_force_tend = 0.0
  scm_state%T_force_tend = 0.0
  scm_state%qv_force_tend = 0.0

  !if(.not. nudge_wind .or. .not. nudge_thermo) then
    !>  - Calculate w_ls and z (height) at model layer interfaces.
    do i=1, scm_state%n_cols
      w_ls_i(i,1) = 0.0
      zi(i,1) = scm_state%geopotential_i(i,1)*grav_inv
      do k=2, scm_state%n_levels
        w_ls_i(i,k) = 0.5*(scm_state%w_ls(i,k-1) + scm_state%w_ls(i,k))
        zi(i,k) = scm_state%geopotential_i(i,k)*grav_inv
      end do
      w_ls_i(i,scm_state%n_levels+1) = w_ls_i(i,scm_state%n_levels)
      zi(i,scm_state%n_levels+1) = scm_state%geopotential_i(i,scm_state%n_levels+1)*grav_inv
    end do
  !end if

  select case(scm_state%mom_forcing_type)
    case (1)
      write(*,*) 'momentum forcing type = 1 is not implemented. Pick 2 or 3. Stopping...'
      stop
    case (2)
      !> - Calculate change in state momentum variables due to vertical advection (subsidence).

      !>  - Calculate tendencies due to vertical advection using same discretization as in previous GFS SCM implmentation (staggered central difference)
      !!    \f[
      !!    \frac{\partial x}{\partial t}|_{vert. advection} = \frac{w_{k+1}\left(x_{k+1} - x_{k}\right) + w_k\left(x_k - x_{k-1}\right)}{-2\left(z_{k+1}-z_{k}\right)}
      !!    \f]
      !!    where \f$w\f$ is the vertical velocity on interface levels, \f$x\f$ is some state variable on grid centers, \f$z\f$ is the height of interface levels. An model layer shares the same index as the interface below.
      do i=1, scm_state%n_cols
        do k=2, scm_state%n_levels-1
          scm_state%u_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_u(i,k+1) - old_u(i,k)) + w_ls_i(i,k)*(old_u(i,k) - old_u(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
          scm_state%v_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_v(i,k+1) - old_v(i,k)) + w_ls_i(i,k)*(old_v(i,k) - old_v(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
        end do
        !>  - Handle the top and bottom levels separately using special discretizations.
        scm_state%u_force_tend(i,1) = -w_ls_i(i,2)*(old_u(i,2) - old_u(i,1))/(zi(i,2)-zi(i,1))
        scm_state%u_force_tend(i,scm_state%n_levels) = -w_ls_i(i,scm_state%n_levels)*&
          (old_u(i,scm_state%n_levels) - old_u(i,scm_state%n_levels-1))/(zi(i,scm_state%n_levels+1)-zi(i,scm_state%n_levels))
        scm_state%v_force_tend(i,1) = -w_ls_i(i,2)*(old_v(i,2) - old_v(i,1))/(zi(i,2)-zi(i,1))
        scm_state%v_force_tend(i,scm_state%n_levels) = -w_ls_i(i,scm_state%n_levels)*&
          (old_v(i,scm_state%n_levels) - old_v(i,scm_state%n_levels-1))/(zi(i,scm_state%n_levels+1)-zi(i,scm_state%n_levels))

        !> - Add forcing due to geostrophic wind
        !>  - Calculate Coriolis parameter.
        f_coriolis = 2.0*con_omega*sin(scm_state%lat(i))
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%u_force_tend(i,k) = scm_state%u_force_tend(i,k) +  f_coriolis*(old_v(i,k) - scm_state%v_g(i,k))
          scm_state%v_force_tend(i,k) = scm_state%v_force_tend(i,k) -  f_coriolis*(old_u(i,k) - scm_state%u_g(i,k))
        end do
      end do
    case (3)
      !> - Calculate change in state momentum variables due to nudging.
      do i=1, scm_state%n_cols
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%u_force_tend(i,k) = (scm_state%u_nudge(i,k) - old_u(i,k))/scm_state%relax_time
          scm_state%v_force_tend(i,k) = (scm_state%v_nudge(i,k) - old_v(i,k))/scm_state%relax_time
        end do
      end do
    case default
      scm_state%u_force_tend = 0.0
      scm_state%v_force_tend = 0.0
  end select

  select case (scm_state%thermo_forcing_type)
    case (1)
      do i=1, scm_state%n_cols
        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, scm_state%n_levels
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + scm_state%dT_dt_rad(i,k) + &
            scm_state%exner_l(i,k)*(scm_state%h_advec_thil(i,k) +scm_state%v_advec_thil(i,k))
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) + scm_state%h_advec_qt(i,k) + scm_state%v_advec_qt(i,k)
        end do
      end do
    case (2)
      do i=1, scm_state%n_cols
        do k=2, scm_state%n_levels-1
          !upstream scheme (for boundaries, assume vertical derivatives are 0 => no vertical advection)
          omega_plus = MAX(scm_state%omega(i,k), 0.0)
          omega_minus = MIN(scm_state%omega(i,k), 0.0)
          dth_dp_plus = (theta(i,k) - theta(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dth_dp_minus = (theta(i,k+1) - theta(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          dqv_dp_plus = (old_qv(i,k)-old_qv(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dqv_dp_minus = (old_qv(i,k+1)-old_qv(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          scm_state%qv_force_tend(i,k) = -omega_plus*dqv_dp_minus - omega_minus*dqv_dp_plus
          scm_state%T_force_tend(i,k) = scm_state%exner_l(i,k)*(-omega_plus*dth_dp_minus - omega_minus*dth_dp_plus)

        end do

        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, scm_state%n_levels
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + scm_state%dT_dt_rad(i,k) + &
            scm_state%exner_l(i,k)*scm_state%h_advec_thil(i,k)
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) + scm_state%h_advec_qt(i,k)
        end do
      end do
    case (3)
      !> - Calculate change in state temperature/moisture variables due to nudging.
      do i=1, scm_state%n_cols
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%T_force_tend(i,k) = (scm_state%T_nudge(i,k) - old_T(i,k))/scm_state%relax_time
          scm_state%qv_force_tend(i,k) = (scm_state%qt_nudge(i,k) - old_qv(i,k))/scm_state%relax_time
        end do

        do k=2, scm_state%n_levels-1
          !upstream scheme (for boundaries, assume vertical derivatives are 0 => no vertical advection)
          omega_plus = MAX(scm_state%omega(i,k), 0.0)
          omega_minus = MIN(scm_state%omega(i,k), 0.0)
          dth_dp_plus = (theta(i,k) - theta(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dth_dp_minus = (theta(i,k+1) - theta(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          dqv_dp_plus = (old_qv(i,k)-old_qv(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dqv_dp_minus = (old_qv(i,k+1)-old_qv(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) -omega_plus*dqv_dp_minus - omega_minus*dqv_dp_plus
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + &
            scm_state%exner_l(i,k)*(-omega_plus*dth_dp_minus - omega_minus*dth_dp_plus)
        end do


      end do
    case default
      scm_state%T_force_tend = 0.0
      scm_state%qv_force_tend = 0.0
  end select

  do i=1, scm_state%n_cols
    do k=1, scm_state%n_levels
      !> - Update the state variables using the leapfrog scheme:
      !!   \f[
      !!   x^{\tau + 1} = \overline{x^{\tau - 1}} + 2\Delta t\frac{\partial x}{\partial t}|^\tau_{forcing}
      !!   \f]
      !!   \f$\overline{x^{\tau - 1}}\f$ is the filtered value at the previous time step and \f$\frac{\partial x}{\partial t}|^\tau_{forcing}\f$ is the sum of forcing terms calculated in this time step.
      scm_state%state_u(i,k,1) = old_u(i,k) + 2.0*scm_state%dt*scm_state%u_force_tend(i,k)
      scm_state%state_v(i,k,1) = old_v(i,k) + 2.0*scm_state%dt*scm_state%v_force_tend(i,k)
      scm_state%state_T(i,k,1) = scm_state%state_T(i,k,1) + 2.0*scm_state%dt*(scm_state%T_force_tend(i,k))
      scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) = scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) + &
        2.0*scm_state%dt*(scm_state%qv_force_tend(i,k))
      ! scm_state%state_u(i,k,1) = old_u(i,k) + scm_state%dt*scm_state%u_force_tend(i,k)
      ! scm_state%state_v(i,k,1) = old_v(i,k) + scm_state%dt*scm_state%v_force_tend(i,k)
      ! scm_state%state_T(i,k,1) = scm_state%state_T(i,k,1) + scm_state%dt*(scm_state%T_force_tend(i,k))
      ! scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) = scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) + scm_state%dt*(scm_state%qv_force_tend(i,k))
    end do
  end do
  !> @}

end subroutine apply_forcing_leapfrog

subroutine apply_forcing_forward_Euler(scm_state)
  use scm_type_defs, only: scm_state_type

  type(scm_state_type), intent(inout) :: scm_state

  real(kind=dp) :: old_u(scm_state%n_cols, scm_state%n_levels), old_v(scm_state%n_cols, scm_state%n_levels), &
    old_T(scm_state%n_cols, scm_state%n_levels), old_qv(scm_state%n_cols, scm_state%n_levels)
  real(kind=dp) :: w_ls_i(scm_state%n_cols, scm_state%n_levels+1), zi(scm_state%n_cols, scm_state%n_levels+1)
  real(kind=dp) :: theta(scm_state%n_cols, scm_state%n_levels)

  integer :: i,k
  real(kind=dp) :: f_coriolis, grav_inv, g_over_cp, omega_plus, omega_minus, dth_dp_plus, dth_dp_minus, &
    dqv_dp_plus, dqv_dp_minus

  !> \section apply_leapfrog_forcing_alg Algorithm
  !! @{

  grav_inv = 1.0/con_g
  g_over_cp = con_g/con_cp

  !> - Save old state variables (filtered from previous time step)
  old_u = scm_state%state_u(:,:,1)
  old_v = scm_state%state_v(:,:,1)
  old_T = scm_state%state_T(:,:,1)
  old_qv = scm_state%state_tracer(:,:,scm_state%water_vapor_index,1)

  theta = old_T/scm_state%exner_l(:,:)

  !> - Initialize forcing sums to zero.
  scm_state%u_force_tend = 0.0
  scm_state%v_force_tend = 0.0
  scm_state%T_force_tend = 0.0
  scm_state%qv_force_tend = 0.0

  !if(.not. nudge_wind .or. .not. nudge_thermo) then
    !>  - Calculate w_ls and z (height) at model layer interfaces.
    do i=1, scm_state%n_cols
      w_ls_i(i,1) = 0.0
      zi(i,1) = scm_state%geopotential_i(i,1)*grav_inv
      do k=2, scm_state%n_levels
        w_ls_i(i,k) = 0.5*(scm_state%w_ls(i,k-1) + scm_state%w_ls(i,k))
        zi(i,k) = scm_state%geopotential_i(i,k)*grav_inv
      end do
      w_ls_i(i,scm_state%n_levels+1) = w_ls_i(i,scm_state%n_levels)
      zi(i,scm_state%n_levels+1) = scm_state%geopotential_i(i,scm_state%n_levels+1)*grav_inv
    end do
  !end if

  select case(scm_state%mom_forcing_type)
    case (1)
      write(*,*) 'momentum forcing type = 1 is not implemented. Pick 2 or 3. Stopping...'
      stop
    case (2)
      !> - Calculate change in state momentum variables due to vertical advection (subsidence).

      !>  - Calculate tendencies due to vertical advection using same discretization as in previous GFS SCM implmentation (staggered central difference)
      !!    \f[
      !!    \frac{\partial x}{\partial t}|_{vert. advection} = \frac{w_{k+1}\left(x_{k+1} - x_{k}\right) + w_k\left(x_k - x_{k-1}\right)}{-2\left(z_{k+1}-z_{k}\right)}
      !!    \f]
      !!    where \f$w\f$ is the vertical velocity on interface levels, \f$x\f$ is some state variable on grid centers, \f$z\f$ is the height of interface levels. An model layer shares the same index as the interface below.
      do i=1, scm_state%n_cols
        do k=2, scm_state%n_levels-1
          scm_state%u_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_u(i,k+1) - old_u(i,k)) + w_ls_i(i,k)*(old_u(i,k) - old_u(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
          scm_state%v_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_v(i,k+1) - old_v(i,k)) + w_ls_i(i,k)*(old_v(i,k) - old_v(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
        end do
        !>  - Handle the top and bottom levels separately using special discretizations.
        scm_state%u_force_tend(i,1) = -w_ls_i(i,2)*(old_u(i,2) - old_u(i,1))/(zi(i,2)-zi(i,1))
        scm_state%u_force_tend(i,scm_state%n_levels) = -w_ls_i(i,scm_state%n_levels)*&
          (old_u(i,scm_state%n_levels) - old_u(i,scm_state%n_levels-1))/(zi(i,scm_state%n_levels+1)-zi(i,scm_state%n_levels))
        scm_state%v_force_tend(i,1) = -w_ls_i(i,2)*(old_v(i,2) - old_v(i,1))/(zi(i,2)-zi(i,1))
        scm_state%v_force_tend(i,scm_state%n_levels) = -w_ls_i(i,scm_state%n_levels)*&
          (old_v(i,scm_state%n_levels) - old_v(i,scm_state%n_levels-1))/(zi(i,scm_state%n_levels+1)-zi(i,scm_state%n_levels))

        !> - Add forcing due to geostrophic wind
        !>  - Calculate Coriolis parameter.
        f_coriolis = 2.0*con_omega*sin(scm_state%lat(i))
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%u_force_tend(i,k) = scm_state%u_force_tend(i,k) +  f_coriolis*(old_v(i,k) - scm_state%v_g(i,k))
          scm_state%v_force_tend(i,k) = scm_state%v_force_tend(i,k) -  f_coriolis*(old_u(i,k) - scm_state%u_g(i,k))
        end do
      end do
    case (3)
      !> - Calculate change in state momentum variables due to nudging.
      do i=1, scm_state%n_cols
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%u_force_tend(i,k) = (scm_state%u_nudge(i,k) - old_u(i,k))/scm_state%relax_time
          scm_state%v_force_tend(i,k) = (scm_state%v_nudge(i,k) - old_v(i,k))/scm_state%relax_time
        end do
      end do
    case default
      scm_state%u_force_tend = 0.0
      scm_state%v_force_tend = 0.0
  end select

  select case (scm_state%thermo_forcing_type)
    case (1)
      do i=1, scm_state%n_cols
        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, scm_state%n_levels
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + scm_state%dT_dt_rad(i,k) + &
            scm_state%exner_l(i,k)*(scm_state%h_advec_thil(i,k) +scm_state%v_advec_thil(i,k))
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) + scm_state%h_advec_qt(i,k) + scm_state%v_advec_qt(i,k)
        end do
      end do
    case (2)
      do i=1, scm_state%n_cols
        do k=2, scm_state%n_levels-1
          !upstream scheme (for boundaries, assume vertical derivatives are 0 => no vertical advection)
          omega_plus = MAX(scm_state%omega(i,k), 0.0)
          omega_minus = MIN(scm_state%omega(i,k), 0.0)
          dth_dp_plus = (theta(i,k) - theta(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dth_dp_minus = (theta(i,k+1) - theta(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          dqv_dp_plus = (old_qv(i,k)-old_qv(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dqv_dp_minus = (old_qv(i,k+1)-old_qv(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          scm_state%qv_force_tend(i,k) = -omega_plus*dqv_dp_minus - omega_minus*dqv_dp_plus
          scm_state%T_force_tend(i,k) = scm_state%exner_l(i,k)*(-omega_plus*dth_dp_minus - omega_minus*dth_dp_plus)
        end do

        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, scm_state%n_levels
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + scm_state%dT_dt_rad(i,k) + &
            scm_state%exner_l(i,k)*scm_state%h_advec_thil(i,k)
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) + scm_state%h_advec_qt(i,k)
        end do
      end do
    case (3)
      !> - Calculate change in state temperature/moisture variables due to nudging.
      do i=1, scm_state%n_cols
        do k=1, scm_state%n_levels
          !accumulate forcing tendencies
          scm_state%T_force_tend(i,k) = (scm_state%T_nudge(i,k) - old_T(i,k))/scm_state%relax_time
          scm_state%qv_force_tend(i,k) = (scm_state%qt_nudge(i,k) - old_qv(i,k))/scm_state%relax_time
        end do

        do k=2, scm_state%n_levels-1
          !upstream scheme (for boundaries, assume vertical derivatives are 0 => no vertical advection)
          omega_plus = MAX(scm_state%omega(i,k), 0.0)
          omega_minus = MIN(scm_state%omega(i,k), 0.0)
          dth_dp_plus = (theta(i,k) - theta(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dth_dp_minus = (theta(i,k+1) - theta(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          dqv_dp_plus = (old_qv(i,k)-old_qv(i,k-1))/(scm_state%pres_l(i,k)-scm_state%pres_l(i,k-1))
          dqv_dp_minus = (old_qv(i,k+1)-old_qv(i,k))/(scm_state%pres_l(i,k+1)-scm_state%pres_l(i,k))
          scm_state%qv_force_tend(i,k) = scm_state%qv_force_tend(i,k) -omega_plus*dqv_dp_minus - omega_minus*dqv_dp_plus
          scm_state%T_force_tend(i,k) = scm_state%T_force_tend(i,k) + &
            scm_state%exner_l(i,k)*(-omega_plus*dth_dp_minus - omega_minus*dth_dp_plus)
        end do
      end do
    case default
      scm_state%T_force_tend = 0.0
      scm_state%qv_force_tend = 0.0
  end select

  do i=1, scm_state%n_cols
    do k=1, scm_state%n_levels
      !> - Update the state variables using the forward Euler scheme:
      !!   \f[
      !!   x^{\tau + 1} = x^{\tau} + \Delta t\frac{\partial x}{\partial t}|^\tau_{forcing}
      !!   \f]
      !!   \f$x^{\tau}\f$ is the value at the previous time step and \f$\frac{\partial x}{\partial t}|^\tau_{forcing}\f$ is the sum of forcing terms calculated in this time step.
      scm_state%state_u(i,k,1) = old_u(i,k) + scm_state%dt*scm_state%u_force_tend(i,k)
      scm_state%state_v(i,k,1) = old_v(i,k) + scm_state%dt*scm_state%v_force_tend(i,k)
      scm_state%state_T(i,k,1) = scm_state%state_T(i,k,1) + scm_state%dt*(scm_state%T_force_tend(i,k))
      scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) = scm_state%state_tracer(i,k,scm_state%water_vapor_index,1) + &
        scm_state%dt*(scm_state%qv_force_tend(i,k))
    end do
  end do
  !> @}

end subroutine apply_forcing_forward_Euler

!> @}
!> @}
end module scm_forcing
