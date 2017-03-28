!> \file gmtb_scm_forcing.f90
!!  Contains subroutines to handle the SCM forcing -- interpolating in space and time, etc.

module gmtb_scm_forcing

use gmtb_scm_kinds, only: sp, dp, qp
use gmtb_scm_utils, only: interpolate_to_grid_centers

use gmtb_scm_physical_constants, only: con_pi, con_omega, con_g, con_cp, con_rd

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup forcing gmtb_scm_forcing
!! @{
!! Contains subroutines to handle the SCM forcing -- interpolating in space and time, etc.

!> This subroutine interpolates the model forcing, column position, and surface properties to the current model time and to the model grid.
!! \note The input forcing file contains forcing for one column, yet the SCM code can accommodate more than one column. Right now, all columns are assumed to have the same forcing.
subroutine interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, &
    input_v_nudge, input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
    input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
    input_sh_flux_sfc, input_lh_flux_sfc, n_levels,&
    n_columns, grid_pres, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, &
    h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, lat, lon, pres_surf, T_surf, sh_flux, lh_flux)

  integer, intent(in) :: input_ntimes !< number of forcing profiles in the input file corresponding to changing forcing in time
  integer, intent(in) :: input_nlev !< number of levels in each input forcing profile
  real(kind=dp), intent(in) :: input_w_ls(:,:) !< input large-scale w (m/s) (time, levels)
  real(kind=dp), intent(in) :: input_omega(:,:) !< input large-scale omega (Pa/s) (time, levels)
  real(kind=dp), intent(in) :: input_u_g(:,:) !< input geostrophic zonal wind (m/s) (time, levels)
  real(kind=dp), intent(in) :: input_v_g(:,:) !< input geostrophic meridional wind (m/s) (time, levels)
  real(kind=dp), intent(in) :: input_u_nudge(:,:) !< input nudging zonal wind (m/s) (time, levels)
  real(kind=dp), intent(in) :: input_v_nudge(:,:) !< input nudging meridional wind (m/s) (time, levels)
  real(kind=dp), intent(in) :: input_T_nudge(:,:) !< input nudging abs. temperature (K)) (time, levels)
  real(kind=dp), intent(in) :: input_thil_nudge(:,:) !< input nudging liq. pot. temperature (K) (time, levels)
  real(kind=dp), intent(in) :: input_qt_nudge(:,:) !< input nudging specific humidity (kg/kg) (time, levels)
  real(kind=dp), intent(in) :: input_dT_dt_rad(:,:) !< input radiative heating rate (K/s) (time, levels)
  real(kind=dp), intent(in) :: input_h_advec_thetail(:,:) !< input change in theta_il due to horizontal advection (K/s) (time, levels)
  real(kind=dp), intent(in) :: input_h_advec_qt(:,:) !< input change in q_t due to horizontal advection (kg/kg /s) (time, levels)
  real(kind=dp), intent(in) :: input_v_advec_thetail(:,:) !< input change in theta_il due to vertical advection (K/s) (time, levels)
  real(kind=dp), intent(in) :: input_v_advec_qt(:,:) !< input change in q_t due to vertical advection (kg/kg /s) (time, levels)
  real(kind=dp), intent(in) :: input_pres(:) !< input pressure levels (Pa) (levels)
  real(kind=dp), intent(in) :: input_time(:) !< input times (s) (time)
  real(kind=dp), intent(in) :: input_lat(:) !< input latitude (degrees) (time)
  real(kind=dp), intent(in) :: input_lon(:) !< input longitude (degrees) (time)
  real(kind=dp), intent(in) :: input_pres_surf(:) !< input surface pressure (Pa) (time)
  real(kind=dp), intent(in) :: input_T_surf(:) !< input surface temperature (K) (time)
  real(kind=dp), intent(in) :: input_sh_flux_sfc(:) !< input surface sensible heat flux (K m s^-1) (time)
  real(kind=dp), intent(in) :: input_lh_flux_sfc(:) !< input surface latent heat flux (kg kg^-1 m s^-1) (time)
  integer, intent(in) :: n_levels !< number of model levels
  integer, intent(in) :: n_columns !< number of model columns
  real(kind=dp), intent(in) :: grid_pres(:,:) !< pressure at grid centers (Pa) (horizontal, vertical)
  real(kind=dp), intent(in) :: model_time !< current model elapsed time (s)
  real(kind=dp), intent(out) :: w_ls(:,:) !< output large scale w (m/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: omega(:,:) !< output large scale pressure velocity (Pa/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: u_g(:,:) !< output geostrophic zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: v_g(:,:) !< output geostrophic meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: u_nudge(:,:) !< output nudging zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: v_nudge(:,:) !< output nudging meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: T_nudge(:,:) !< output nudging abs. temperature (K) (horizontal, vertical)
  real(kind=dp), intent(out) :: thil_nudge(:,:) !< output nudging liq. pot. temperature (K) (horizontal, vertical)
  real(kind=dp), intent(out) :: qt_nudge(:,:) !< output nudging specific humidity (kg/kg) (horizontal, vertical)
  real(kind=dp), intent(out) :: dT_dt_rad(:,:) !< output radiative heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: h_advec_thil(:,:) !< output change in theta_il due to horizontal advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: h_advec_qt(:,:) !< output change in q_t due to horizontal advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(out) :: v_advec_thil(:,:) !< output change in theta_il due to vertical advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(out) :: v_advec_qt(:,:) !< output change in q_t due to vertical advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(out) :: lat(:) !< output latitude (horizontal)
  real(kind=dp), intent(out) :: lon(:) !< output longitude (horizontal)
  real(kind=dp), intent(out) :: pres_surf(:) !< output surface pressure (Pa) (horizontal)
  real(kind=dp), intent(out) :: T_surf(:) !< output surface temperature (K) (horizontal)
  real(kind=dp), intent(out) :: sh_flux(:) !< output surface sensible heat flux (K m s^-1) (horizontal)
  real(kind=dp), intent(out) :: lh_flux(:) !< output surface latent heat flux (kg kg^-1 m s^-1) (K) (horizontal)

  integer :: i, n
  integer :: low_t_index, top_index !< index of the time in the input file immediately before the current model time, index of the last calculated level
  real(kind=dp) :: lifrac
  real(kind=dp) :: deg_to_rad_const

  real(kind=dp) :: w_ls_bracket(2,n_levels), omega_bracket(2,n_levels), u_g_bracket(2,n_levels), v_g_bracket(2,n_levels), &
    u_nudge_bracket(2,n_levels), v_nudge_bracket(2,n_levels), T_nudge_bracket(2,n_levels), thil_nudge_bracket(2,n_levels), &
    qt_nudge_bracket(2,n_levels), dT_dt_rad_bracket(2,n_levels), h_advec_thil_bracket(2,n_levels),h_advec_qt_bracket(2,n_levels), &
    v_advec_thil_bracket(2,n_levels), v_advec_qt_bracket(2,n_levels) !< forcing terms that "bracket" around the model time

  !> \section interpolate_forcing_alg Algorithm
  !! @{

  deg_to_rad_const = con_pi/180.0

  !> - Check for the case where the elapsed model time extends beyond the supplied forcing.
  if(model_time >= input_time(input_ntimes)) then
    !>  - If so, hold the forcing terms constant at the last supplied values. The forcing still needs to be interpolated to the grid.
    write(*,*) "The model_time has exceeded the specifed period of forcing. Forcing will now be held constant at the last &
      specified values."

      !>  - For all forcing terms, call interpolate_to_grid_centers from \ref utils for each variable. This subroutine returns the last vertical index calculated in case forcing terms above the case input needs to be specified.
      do i=1, n_columns
        call interpolate_to_grid_centers(input_nlev, input_pres, input_w_ls(input_ntimes,:), grid_pres(i,:), n_levels, &
          w_ls_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_omega(input_ntimes,:), grid_pres(i,:), n_levels, &
          omega_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_u_g(input_ntimes,:), grid_pres(i,:), n_levels, &
          u_g_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_v_g(input_ntimes,:), grid_pres(i,:), n_levels, &
          v_g_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_u_nudge(input_ntimes,:), grid_pres(i,:), n_levels, &
          u_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_v_nudge(input_ntimes,:), grid_pres(i,:), n_levels, &
          v_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_T_nudge(input_ntimes,:), grid_pres(i,:), n_levels, &
          T_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_thil_nudge(input_ntimes,:), grid_pres(i,:), n_levels, &
          thil_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_qt_nudge(input_ntimes,:), grid_pres(i,:), n_levels, &
          qt_nudge_bracket(1,:), top_index, 1)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_dT_dt_rad(input_ntimes,:), grid_pres(i,:), n_levels, &
          dT_dt_rad_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_thetail(input_ntimes,:), grid_pres(i,:), n_levels, &
          h_advec_thil_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_qt(input_ntimes,:), grid_pres(i,:), n_levels, &
          h_advec_qt_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_thetail(input_ntimes,:), grid_pres(i,:), n_levels, &
          v_advec_thil_bracket(1,:), top_index, 3)
        call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_qt(input_ntimes,:), grid_pres(i,:), n_levels, &
          v_advec_qt_bracket(1,:), top_index, 3)

        !>  - If the input forcing file does not reach to the model domain top, fill in values above the input forcing file domain with those from the top level.
        if (top_index < n_levels) then
          w_ls_bracket(1,top_index+1:n_levels) = w_ls_bracket(1,top_index)
          omega_bracket(1,top_index+1:n_levels) = omega_bracket(1,top_index)
          u_g_bracket(1,top_index+1:n_levels) = u_g_bracket(1,top_index)
          v_g_bracket(1,top_index+1:n_levels) = v_g_bracket(1,top_index)
          u_nudge_bracket(1,top_index+1:n_levels) = u_nudge_bracket(1,top_index)
          v_nudge_bracket(1,top_index+1:n_levels) = v_nudge_bracket(1,top_index)
          T_nudge_bracket(1,top_index+1:n_levels) = T_nudge_bracket(1,top_index)
          thil_nudge_bracket(1,top_index+1:n_levels) = thil_nudge_bracket(1,top_index)
          qt_nudge_bracket(1,top_index+1:n_levels) = qt_nudge_bracket(1,top_index)
          dT_dt_rad_bracket(1,top_index+1:n_levels) = dT_dt_rad_bracket(1,top_index)
          h_advec_thil_bracket(1,top_index+1:n_levels) = h_advec_thil_bracket(1,top_index)
          h_advec_qt_bracket(1,top_index+1:n_levels) = h_advec_qt_bracket(1,top_index)
          v_advec_thil_bracket(1,top_index+1:n_levels) = v_advec_thil_bracket(1,top_index)
          v_advec_qt_bracket(1,top_index+1:n_levels) = v_advec_qt_bracket(1,top_index)
        end if

        !>  - For this case, no time interpolation is necessary; just set the forcing terms to the vertically-interpolated values.
        w_ls(i,:) = w_ls_bracket(1,:)
        omega(i,:) = omega_bracket(1,:)
        u_g(i,:) = u_g_bracket(1,:)
        v_g(i,:) = v_g_bracket(1,:)
        u_nudge(i,:) = u_nudge_bracket(1,:)
        v_nudge(i,:) = v_nudge_bracket(1,:)
        T_nudge(i,:) = T_nudge_bracket(1,:)
        thil_nudge(i,:) = thil_nudge_bracket(1,:)
        qt_nudge(i,:) = qt_nudge_bracket(1,:)
        dT_dt_rad(i,:) = dT_dt_rad_bracket(1,:)
        h_advec_thil(i,:) = h_advec_thil_bracket(1,:)
        h_advec_qt(i,:) = h_advec_qt_bracket(1,:)
        v_advec_thil(i,:) = v_advec_thil_bracket(1,:)
        v_advec_qt(i,:) = v_advec_qt_bracket(1,:)

        !>  - Set the surface parameters to the last available data.
        lat(i) = input_lat(input_ntimes)
        lon(i) = input_lon(input_ntimes)
        pres_surf(i) = input_pres_surf(input_ntimes)
        T_surf(i) = input_T_surf(input_ntimes)
        sh_flux(i) = input_sh_flux_sfc(input_ntimes)
        lh_flux(i) = input_lh_flux_sfc(input_ntimes)
      end do
      !>  - Convert lat, lon to radians.
      lat = lat*deg_to_rad_const
      lon = lon*deg_to_rad_const
  else
  !> - When the model elapsed time is within the time-frame specified by the input forcing, the forcing must be interpolated in time and space.
    !>  - Determine the index in the input file for the time immediately preceeding the model time and determine the linear interpolation value.
    low_t_index = 0
    do n=1, input_ntimes
      if (input_time(n) > model_time) then
        low_t_index = n-1
        lifrac = (model_time - input_time(low_t_index))/(input_time(low_t_index+1) - input_time(low_t_index))
        exit
      end if
    end do

    do i=1, n_columns
      !>  - For all forcing terms, call interpolate_to_grid_centers from \ref utils for each variable for each time level that "bracket" around the current model time. This subroutine returns the last vertical index calculated in case forcing terms above the case input needs to be specified.
      call interpolate_to_grid_centers(input_nlev, input_pres, input_w_ls(low_t_index,:), grid_pres(i,:), n_levels, &
        w_ls_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_w_ls(low_t_index+1,:), grid_pres(i,:), n_levels, &
        w_ls_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_omega(low_t_index,:), grid_pres(i,:), n_levels, &
        omega_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_omega(low_t_index+1,:), grid_pres(i,:), n_levels, &
        omega_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_u_g(low_t_index,:), grid_pres(i,:), n_levels, &
        u_g_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_u_g(low_t_index+1,:), grid_pres(i,:), n_levels, &
        u_g_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_g(low_t_index,:), grid_pres(i,:), n_levels, &
        v_g_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_g(low_t_index+1,:), grid_pres(i,:), n_levels, &
        v_g_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_u_nudge(low_t_index,:), grid_pres(i,:), n_levels, &
        u_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_u_nudge(low_t_index+1,:), grid_pres(i,:), n_levels, &
        u_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_nudge(low_t_index,:), grid_pres(i,:), n_levels, &
        v_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_nudge(low_t_index+1,:), grid_pres(i,:), n_levels, &
        v_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_T_nudge(low_t_index,:), grid_pres(i,:), n_levels, &
        T_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_T_nudge(low_t_index+1,:), grid_pres(i,:), n_levels, &
        T_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_thil_nudge(low_t_index,:), grid_pres(i,:), n_levels, &
        thil_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_thil_nudge(low_t_index+1,:), grid_pres(i,:), n_levels, &
        thil_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_qt_nudge(low_t_index,:), grid_pres(i,:), n_levels, &
        qt_nudge_bracket(1,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_qt_nudge(low_t_index+1,:), grid_pres(i,:), n_levels, &
        qt_nudge_bracket(2,:), top_index, 1)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_dT_dt_rad(low_t_index,:), grid_pres(i,:), n_levels, &
        dT_dt_rad_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_dT_dt_rad(low_t_index+1,:), grid_pres(i,:), n_levels, &
        dT_dt_rad_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_thetail(low_t_index,:), grid_pres(i,:), n_levels, &
        h_advec_thil_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_thetail(low_t_index+1,:), grid_pres(i,:), n_levels, &
        h_advec_thil_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_qt(low_t_index,:), grid_pres(i,:), n_levels, &
        h_advec_qt_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_h_advec_qt(low_t_index+1,:), grid_pres(i,:), n_levels, &
        h_advec_qt_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_thetail(low_t_index,:), grid_pres(i,:), n_levels, &
        v_advec_thil_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_thetail(low_t_index+1,:), grid_pres(i,:), n_levels, &
        v_advec_thil_bracket(2,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_qt(low_t_index,:), grid_pres(i,:), n_levels, &
        v_advec_qt_bracket(1,:), top_index, 3)
      call interpolate_to_grid_centers(input_nlev, input_pres, input_v_advec_qt(low_t_index+1,:), grid_pres(i,:), n_levels, &
        v_advec_qt_bracket(2,:), top_index, 3)

      !>  - If the input forcing file does not reach to the model domain top, fill in values above the input forcing file domain with those from the top level.
      if (top_index < n_levels) then
        w_ls_bracket(1,top_index+1:n_levels) = w_ls_bracket(1,top_index)
        w_ls_bracket(2,top_index+1:n_levels) = w_ls_bracket(2,top_index)
        omega_bracket(1,top_index+1:n_levels) = omega_bracket(1,top_index)
        omega_bracket(2,top_index+1:n_levels) = omega_bracket(2,top_index)
        u_g_bracket(1,top_index+1:n_levels) = u_g_bracket(1,top_index)
        u_g_bracket(2,top_index+1:n_levels) = u_g_bracket(2,top_index)
        v_g_bracket(1,top_index+1:n_levels) = v_g_bracket(1,top_index)
        v_g_bracket(2,top_index+1:n_levels) = v_g_bracket(2,top_index)
        u_nudge_bracket(1,top_index+1:n_levels) = u_nudge_bracket(1,top_index)
        u_nudge_bracket(2,top_index+1:n_levels) = u_nudge_bracket(2,top_index)
        v_nudge_bracket(1,top_index+1:n_levels) = v_nudge_bracket(1,top_index)
        v_nudge_bracket(2,top_index+1:n_levels) = v_nudge_bracket(2,top_index)
        T_nudge_bracket(1,top_index+1:n_levels) = T_nudge_bracket(1,top_index)
        T_nudge_bracket(2,top_index+1:n_levels) = T_nudge_bracket(2,top_index)
        thil_nudge_bracket(1,top_index+1:n_levels) = thil_nudge_bracket(1,top_index)
        thil_nudge_bracket(2,top_index+1:n_levels) = thil_nudge_bracket(2,top_index)
        qt_nudge_bracket(1,top_index+1:n_levels) = qt_nudge_bracket(1,top_index)
        qt_nudge_bracket(2,top_index+1:n_levels) = qt_nudge_bracket(2,top_index)
        dT_dt_rad_bracket(1,top_index+1:n_levels) = dT_dt_rad_bracket(1,top_index)
        dT_dt_rad_bracket(2,top_index+1:n_levels) = dT_dt_rad_bracket(2,top_index)
        h_advec_thil_bracket(1,top_index+1:n_levels) = h_advec_thil_bracket(1,top_index)
        h_advec_thil_bracket(2,top_index+1:n_levels) = h_advec_thil_bracket(2,top_index)
        h_advec_qt_bracket(1,top_index+1:n_levels) = h_advec_qt_bracket(1,top_index)
        h_advec_qt_bracket(2,top_index+1:n_levels) = h_advec_qt_bracket(2,top_index)
        v_advec_thil_bracket(1,top_index+1:n_levels) = v_advec_thil_bracket(1,top_index)
        v_advec_thil_bracket(2,top_index+1:n_levels) = v_advec_thil_bracket(2,top_index)
        v_advec_qt_bracket(1,top_index+1:n_levels) = v_advec_qt_bracket(1,top_index)
        v_advec_qt_bracket(2,top_index+1:n_levels) = v_advec_qt_bracket(2,top_index)
      end if

      !>  - Interpolate the forcing terms in time.
      w_ls(i,:) = (1.0 - lifrac)*w_ls_bracket(1,:) + lifrac*w_ls_bracket(2,:)
      omega(i,:) = (1.0 - lifrac)*omega_bracket(1,:) + lifrac*omega_bracket(2,:)
      u_g(i,:) = (1.0 - lifrac)*u_g_bracket(1,:) + lifrac*u_g_bracket(2,:)
      v_g(i,:) = (1.0 - lifrac)*v_g_bracket(1,:) + lifrac*v_g_bracket(2,:)
      u_nudge(i,:) = (1.0 - lifrac)*u_nudge_bracket(1,:) + lifrac*u_nudge_bracket(2,:)
      v_nudge(i,:) = (1.0 - lifrac)*v_nudge_bracket(1,:) + lifrac*v_nudge_bracket(2,:)
      T_nudge(i,:) = (1.0 - lifrac)*T_nudge_bracket(1,:) + lifrac*T_nudge_bracket(2,:)
      thil_nudge(i,:) = (1.0 - lifrac)*thil_nudge_bracket(1,:) + lifrac*thil_nudge_bracket(2,:)
      qt_nudge(i,:) = (1.0 - lifrac)*qt_nudge_bracket(1,:) + lifrac*qt_nudge_bracket(2,:)
      dT_dt_rad(i,:) = (1.0 - lifrac)*dT_dt_rad_bracket(1,:) + lifrac*dT_dt_rad_bracket(2,:)
      h_advec_thil(i,:) = (1.0 - lifrac)*h_advec_thil_bracket(1,:) + lifrac*h_advec_thil_bracket(2,:)
      h_advec_qt(i,:) = (1.0 - lifrac)*h_advec_qt_bracket(1,:) + lifrac*h_advec_qt_bracket(2,:)
      v_advec_thil(i,:) = (1.0 - lifrac)*v_advec_thil_bracket(1,:) + lifrac*v_advec_thil_bracket(2,:)
      v_advec_qt(i,:) = (1.0 - lifrac)*v_advec_qt_bracket(1,:) + lifrac*v_advec_qt_bracket(2,:)

      !>  - Interpolate the surface parameters in time.
      lat(i) = (1.0 - lifrac)*input_lat(low_t_index) + lifrac*input_lat(low_t_index+1)
      lon(i) = (1.0 - lifrac)*input_lon(low_t_index) + lifrac*input_lon(low_t_index+1)
      pres_surf(i) = (1.0 - lifrac)*input_pres_surf(low_t_index) + lifrac*input_pres_surf(low_t_index+1)
      T_surf(i) = (1.0 - lifrac)*input_T_surf(low_t_index) + lifrac*input_T_surf(low_t_index+1)
      sh_flux(i) = (1.0 - lifrac)*input_sh_flux_sfc(low_t_index) + lifrac*input_sh_flux_sfc(low_t_index+1)
      lh_flux(i) = (1.0 - lifrac)*input_lh_flux_sfc(low_t_index) + lifrac*input_lh_flux_sfc(low_t_index+1)
    end do
    !>  - Convert lat, lon to radians.
    lat = lat*deg_to_rad_const
    lon = lon*deg_to_rad_const
  end if

  !> @}
end subroutine interpolate_forcing

!> This subroutine updates the state variables due to the model forcing using the leapfrog time integration scheme.
!! It overwrites the state variable arrays where the filtered values from the previous time step are kept. These parts of the state
!! variable arrays are pointed to by the state_fields_in DDT and are later updated by calling nuopc_phys_run.
subroutine apply_forcing_leapfrog(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, u_nudge, &
  v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, phii, pres_l, &
  lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, qv_force_tend)
  integer, intent(in) :: n_columns !< number of model columns
  integer, intent(in) :: n_levels !< number of model levels
  real(kind=dp), intent(inout) :: state_tracer(:,:,:,:) !< model state tracers (horizontal, vertical, tracer, iteration)
  real(kind=dp), intent(inout) :: state_T(:,:,:) !< model state temperature (horizontal, vertical, iteration)
  real(kind=dp), intent(inout) :: state_u(:,:,:) !< model state zonal wind (horizontal, vertical, iteration)
  real(kind=dp), intent(inout) :: state_v(:,:,:) !< model state meridional wind (horizontal, vertical, iteration)
  real(kind=dp), intent(in) :: w_ls(:,:) !< large scale w (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: omega(:,:) !< large scale pressure velocity (Pa/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: u_g(:,:) !< geostrophic zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_g(:,:) !< geostrophic meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: u_nudge(:,:) !< nudging zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_nudge(:,:) !< nudging meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: T_nudge(:,:) !< nudging abs. temperatuere (K) (horizontal, vertical)
  real(kind=dp), intent(in) :: thil_nudge(:,:) !< nudging liq. pot. temperature (K) (horizontal, vertical)
  real(kind=dp), intent(in) :: qt_nudge(:,:) !< nudging specific humidity (kg/kg) (horizontal, vertical)
  real(kind=dp), intent(in) :: dT_dt_rad(:,:) !< radiative heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: h_advec_thil(:,:) !< change in theta_il due to horizontal advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: h_advec_qt(:,:) !< change in q_t due to horizontal advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_advec_thil(:,:) !< change in theta_il due to vertical advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_advec_qt(:,:) !< change in q_t due to vertical advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in) :: dt !< time step (s)
  real(kind=dp), intent(in) :: lat(:) !< latitude in radians
  real(kind=dp), intent(in) :: exner_l(:,:) !< exner function at model level centers (horizontal, vertical)
  real(kind=dp), intent(in) :: phii(:,:) !< geopotential at model level interfaces (horizontal, vertical)
  real(kind=dp), intent(in) :: pres_l(:,:) !< pressure at model level centers (horizontal, vertical)
  integer, intent(in) :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer, intent(in) :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  real(kind=dp), intent(in) :: relax_time !< time scale for hor. wind nudging
  real(kind=dp), intent(out) :: u_force_tend(:,:) !< sum of all forcing for zonal momentum (horizontal, vertical)
  real(kind=dp), intent(out) :: v_force_tend(:,:) !< sum of all forcing for meridional momentum (horizontal, vertical)
  real(kind=dp), intent(out) :: T_force_tend(:,:) !< sum of all forcing for temperature (horizontal, vertical)
  real(kind=dp), intent(out) :: qv_force_tend(:,:) !< sum of all forcing for water vapor (horizontal, vertical)

  real(kind=dp) :: old_u(n_columns, n_levels), old_v(n_columns, n_levels), old_T(n_columns, n_levels), old_qv(n_columns, n_levels)
  real(kind=dp) :: w_ls_i(n_columns, n_levels+1), zi(n_columns, n_levels+1)
  real(kind=dp) :: theta(n_columns, n_levels)

  integer :: i,k
  real(kind=dp) :: f_coriolis, grav_inv, g_over_cp

  !> \section apply_leapfrog_forcing_alg Algorithm
  !! @{

  grav_inv = 1.0/con_g
  g_over_cp = con_g/con_cp

  !> - Save old state variables (filtered from previous time step)
  old_u = state_u(:,:,1)
  old_v = state_v(:,:,1)
  old_T = state_T(:,:,1)
  old_qv = state_tracer(:,:,1,1)

  theta = old_T/exner_l

  !> - Initialize forcing sums to zero.
  u_force_tend = 0.0
  v_force_tend = 0.0
  T_force_tend = 0.0
  qv_force_tend = 0.0

  !if(.not. nudge_wind .or. .not. nudge_thermo) then
    !>  - Calculate w_ls and z (height) at model layer interfaces.
    do i=1, n_columns
      w_ls_i(i,1) = 0.0
      zi(i,1) = phii(i,1)*grav_inv
      do k=2, n_levels
        w_ls_i(i,k) = 0.5*(w_ls(i,k-1) + w_ls(i,k))
        zi(i,k) = phii(i,k)*grav_inv
      end do
      w_ls_i(i,n_levels+1) = w_ls_i(i,n_levels)
      zi(i,n_levels+1) = phii(i,n_levels+1)*grav_inv
    end do
  !end if

  select case(mom_forcing_type)
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
      do i=1, n_columns
        do k=2, n_levels-1
          u_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_u(i,k+1) - old_u(i,k)) + w_ls_i(i,k)*(old_u(i,k) - old_u(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
          v_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_v(i,k+1) - old_v(i,k)) + w_ls_i(i,k)*(old_v(i,k) - old_v(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
        end do
        !>  - Handle the top and bottom levels separately using special discretizations.
        u_force_tend(i,1) = -w_ls_i(i,2)*(old_u(i,2) - old_u(i,1))/(zi(i,2)-zi(i,1))
        u_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_u(i,n_levels) - old_u(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))
        v_force_tend(i,1) = -w_ls_i(i,2)*(old_v(i,2) - old_v(i,1))/(zi(i,2)-zi(i,1))
        v_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_v(i,n_levels) - old_v(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))

        !> - Add forcing due to geostrophic wind
        !>  - Calculate Coriolis parameter.
        f_coriolis = 2.0*con_omega*sin(lat(i))
        do k=1, n_levels
          !accumulate forcing tendencies
          u_force_tend(i,k) = u_force_tend(i,k) +  f_coriolis*(old_v(i,k) - v_g(i,k))
          v_force_tend(i,k) = v_force_tend(i,k) -  f_coriolis*(old_u(i,k) - u_g(i,k))
        end do
      end do
    case (3)
      !> - Calculate change in state momentum variables due to nudging.
      do i=1, n_columns
        do k=1, n_levels
          !accumulate forcing tendencies
          u_force_tend(i,k) = (u_nudge(i,k) - old_u(i,k))/relax_time
          v_force_tend(i,k) = (v_nudge(i,k) - old_v(i,k))/relax_time
        end do
      end do
    case default
      u_force_tend = 0.0
      v_force_tend = 0.0
  end select

  select case (thermo_forcing_type)
    case (1)
      do i=1, n_columns
        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, n_levels
          T_force_tend(i,k) = T_force_tend(i,k) + dT_dt_rad(i,k) + exner_l(i,k)*(h_advec_thil(i,k) +v_advec_thil(i,k))
          qv_force_tend(i,k) = qv_force_tend(i,k) + h_advec_qt(i,k) + v_advec_qt(i,k)
        end do
      end do
    case (2)
      do i=1, n_columns
        ! do k=2, n_levels-1
        !
        !   qv_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_qv(i,k+1) - old_qv(i,k)) + w_ls_i(i,k)*(old_qv(i,k) - old_qv(i,k-1)))/&
        !     (zi(i,k+1)-zi(i,k))
        !   !T_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_T(i,k+1) - old_T(i,k)) + w_ls_i(i,k)*(old_T(i,k) - old_T(i,k-1)))/(zi(i,k+1)-zi(i,k)) - &
        !   !  0.5*(w_ls_i(i,k+1) + w_ls_i(i,k))*g_over_cp
        !   T_force_tend(i,k) = exner_l(i,k)*(-0.5*(w_ls_i(i,k+1)*(theta(i,k+1) - theta(i,k)) + &
        !     w_ls_i(i,k)*(theta(i,k) - theta(i,k-1))))/(zi(i,k+1)-zi(i,k))
        ! end do
        ! !>  - Handle the top and bottom levels separately using special discretizations.
        ! qv_force_tend(i,1) = -w_ls_i(i,2)*(old_qv(i,2) - old_qv(i,1))/(zi(i,2)-zi(i,1))
        ! qv_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_qv(i,n_levels) - old_qv(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))
        ! T_force_tend(i,1) = -w_ls_i(i,2)*(old_T(i,2) - old_T(i,1))/(zi(i,2)-zi(i,1)) !- 0.5*(w_ls_i(2) + w_ls_i(1))*g_over_cp
        ! T_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_T(i,n_levels) - old_T(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels)) - &
        !   0.5*(w_ls_i(i,n_levels+1) + w_ls_i(i,n_levels))*g_over_cp

        do k=2, n_levels-1
          !upstream scheme
          if(omega(i,k) < 0.0) then
            qv_force_tend(i,k) = -omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
            T_force_tend(i,k) =  exner_l(i,k)*(-omega(i,k)*(theta(i,k) - theta(i,k-1))/(pres_l(i,k)-pres_l(i,k-1)))

          else
            qv_force_tend(i,k) = -omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
            T_force_tend(i,k) = exner_l(i,k)*(-omega(i,k)*(theta(i,k+1) - theta(i,k))/(pres_l(i,k+1)-pres_l(i,k)))

          end if
        end do
        T_force_tend(i,1) = exner_l(i,1)*(-0.5*(omega(i,2)+omega(i,1))*(theta(i,2) - theta(i,1))/(pres_l(i,2)-pres_l(i,1)))
        T_force_tend(i,n_levels) = exner_l(i,n_levels)*(-0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (theta(i,n_levels) - theta(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1)))
        qv_force_tend(i,1) = -0.5*(omega(i,2)+omega(i,1))*(old_qv(i,2) - old_qv(i,1))/(pres_l(i,2)-pres_l(i,1))
        qv_force_tend(i,n_levels) = -0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (old_qv(i,n_levels) - old_qv(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1))

        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, n_levels
          T_force_tend(i,k) = T_force_tend(i,k) + dT_dt_rad(i,k) + exner_l(i,k)*h_advec_thil(i,k)
          qv_force_tend(i,k) = qv_force_tend(i,k) + h_advec_qt(i,k)
        end do
      end do
    case (3)
      !> - Calculate change in state temperature/moisture variables due to nudging.
      do i=1, n_columns
        do k=1, n_levels
          !accumulate forcing tendencies
          T_force_tend(i,k) = (T_nudge(i,k) - old_T(i,k))/relax_time
          qv_force_tend(i,k) = (qt_nudge(i,k) - old_qv(i,k))/relax_time
        end do

        do k=2, n_levels-1
          !upstream scheme
          if(omega(i,k) < 0.0) then
            qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
            T_force_tend(i,k) =  T_force_tend(i,k) + &
              exner_l(i,k)*(-omega(i,k)*(theta(i,k) - theta(i,k-1))/(pres_l(i,k)-pres_l(i,k-1)))

          else
            qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
            T_force_tend(i,k) = T_force_tend(i,k) + &
              exner_l(i,k)*(-omega(i,k)*(theta(i,k+1) - theta(i,k))/(pres_l(i,k+1)-pres_l(i,k)))

          end if
        end do
        T_force_tend(i,1) = T_force_tend(i,1) + &
          exner_l(i,1)*(-0.5*(omega(i,2)+omega(i,1))*(theta(i,2) - theta(i,1))/(pres_l(i,2)-pres_l(i,1)))
        T_force_tend(i,n_levels) = T_force_tend(i,n_levels) + exner_l(i,n_levels)*(-0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (theta(i,n_levels) - theta(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1)))
        qv_force_tend(i,1) = qv_force_tend(i,1) - 0.5*(omega(i,2)+omega(i,1))*(old_qv(i,2) - old_qv(i,1))/(pres_l(i,2)-pres_l(i,1))
        qv_force_tend(i,n_levels) = qv_force_tend(i,n_levels) - 0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (old_qv(i,n_levels) - old_qv(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1))

      end do

      ! do i=1, n_columns
      !   do k=2, n_levels-1
      !     if(omega(i,k) < 0.0) then !need to bring in pres_l
      !       qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
      !       T_force_tend(i,k) = T_force_tend(i,k) - omega(i,k)*(old_T(i,k) - old_T(i,k-1))/(pres_l(i,k)-pres_l(i,k-1) + &
      !         old_T(i,k)*con_rd/(pres_l(i,k)*con_cp))
      !     else
      !       qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
      !       T_force_tend(i,k) = T_force_tend(i,k) - omega(i,k)*(old_T(i,k+1) - old_T(i,k))/(pres_l(i,k+1)-pres_l(i,k) + &
      !         old_T(i,k)*con_rd/(pres_l(i,k)*con_cp))
      !     end if
      !     ! qv_force_tend(i,k) = qv_force_tend(i,k) - 0.5*(w_ls_i(i,k+1)*(old_qv(i,k+1) - old_qv(i,k)) + &
      !     !   w_ls_i(i,k)*(old_qv(i,k) - old_qv(i,k-1)))/(zi(i,k+1)-zi(i,k))
      !     ! !T_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_T(i,k+1) - old_T(i,k)) + w_ls_i(i,k)*(old_T(i,k) - old_T(i,k-1)))/(zi(i,k+1)-zi(i,k)) - &
      !     ! !  0.5*(w_ls_i(i,k+1) + w_ls_i(i,k))*g_over_cp
      !     ! T_force_tend(i,k) = T_force_tend(i,k) + exner_l(i,k)*(-0.5*(w_ls_i(i,k+1)*(theta(i,k+1) - theta(i,k)) + &
      !     !   w_ls_i(i,k)*(theta(i,k) - theta(i,k-1))))/(zi(i,k+1)-zi(i,k))
      !   end do
      !   !>  - Handle the top and bottom levels separately using special discretizations.
      !   ! qv_force_tend(i,1) = qv_force_tend(i,1) - w_ls_i(i,2)*(old_qv(i,2) - old_qv(i,1))/(zi(i,2)-zi(i,1))
      !   ! qv_force_tend(i,n_levels) = qv_force_tend(i,n_levels) - w_ls_i(i,n_levels)*(old_qv(i,n_levels) - old_qv(i,n_levels-1))/&
      !   !   (zi(i,n_levels+1)-zi(i,n_levels))
      !   ! T_force_tend(i,1) = T_force_tend(i,1) - w_ls_i(i,2)*(old_T(i,2) - old_T(i,1))/(zi(i,2)-zi(i,1)) !- 0.5*(w_ls_i(2) + w_ls_i(1))*g_over_cp
      !   ! T_force_tend(i,n_levels) = T_force_tend(i,n_levels) - w_ls_i(i,n_levels)*(old_T(i,n_levels) - old_T(i,n_levels-1))/&
      !   !   (zi(i,n_levels+1)-zi(i,n_levels)) - 0.5*(w_ls_i(i,n_levels+1) + w_ls_i(i,n_levels))*g_over_cp
      !
      ! end do
    case default
      T_force_tend = 0.0
      qv_force_tend = 0.0
  end select

  do i=1, n_columns
    do k=1, n_levels
      !> - Update the state variables using the leapfrog scheme:
      !!   \f[
      !!   x^{\tau + 1} = \overline{x^{\tau - 1}} + 2\Delta t\frac{\partial x}{\partial t}|^\tau_{forcing}
      !!   \f]
      !!   \f$\overline{x^{\tau - 1}}\f$ is the filtered value at the previous time step and \f$\frac{\partial x}{\partial t}|^\tau_{forcing}\f$ is the sum of forcing terms calculated in this time step.
      state_u(i,k,1) = old_u(i,k) + 2.0*dt*u_force_tend(i,k)
      state_v(i,k,1) = old_v(i,k) + 2.0*dt*v_force_tend(i,k)
      state_T(i,k,1) = state_T(i,k,1) + 2.0*dt*(T_force_tend(i,k))
      state_tracer(i,k,1,1) = state_tracer(i,k,1,1) + 2.0*dt*(qv_force_tend(i,k))
      ! state_u(i,k,1) = old_u(i,k) + dt*u_force_tend(i,k)
      ! state_v(i,k,1) = old_v(i,k) + dt*v_force_tend(i,k)
      ! state_T(i,k,1) = state_T(i,k,1) + dt*(T_force_tend(i,k))
      ! state_tracer(i,k,1,1) = state_tracer(i,k,1,1) + dt*(qv_force_tend(i,k))
    end do
  end do
  !> @}

end subroutine apply_forcing_leapfrog

subroutine apply_forcing_forward_Euler(n_levels, n_columns, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, &
  u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, exner_l, phii, &
  pres_l, lat, dt, thermo_forcing_type, mom_forcing_type, relax_time, u_force_tend, v_force_tend, T_force_tend, qv_force_tend)

  integer, intent(in) :: n_columns !< number of model columns
  integer, intent(in) :: n_levels !< number of model levels
  real(kind=dp), intent(inout) :: state_tracer(:,:,:,:) !< model state tracers (horizontal, vertical, tracer, iteration)
  real(kind=dp), intent(inout) :: state_T(:,:,:) !< model state temperature (horizontal, vertical, iteration)
  real(kind=dp), intent(inout) :: state_u(:,:,:) !< model state zonal wind (horizontal, vertical, iteration)
  real(kind=dp), intent(inout) :: state_v(:,:,:) !< model state meridional wind (horizontal, vertical, iteration)
  real(kind=dp), intent(in) :: w_ls(:,:) !< large scale w (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: omega(:,:) !< large scale pressure velocity (Pa/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: u_g(:,:) !< geostrophic zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_g(:,:) !< geostrophic meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: u_nudge(:,:) !< nudging zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_nudge(:,:) !< nudging meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: T_nudge(:,:) !< nudging abs. temperatuere (K) (horizontal, vertical)
  real(kind=dp), intent(in) :: thil_nudge(:,:) !< nudging liq. pot. temperature (K) (horizontal, vertical)
  real(kind=dp), intent(in) :: qt_nudge(:,:) !< nudging specific humidity (kg/kg) (horizontal, vertical)
  real(kind=dp), intent(in) :: dT_dt_rad(:,:) !< radiative heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: h_advec_thil(:,:) !< change in theta_il due to horizontal advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: h_advec_qt(:,:) !< change in q_t due to horizontal advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_advec_thil(:,:) !< change in theta_il due to vertical advection (K/s) (horizontal, vertical)
  real(kind=dp), intent(in) :: v_advec_qt(:,:) !< change in q_t due to vertical advection (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in) :: dt !< time step (s)
  real(kind=dp), intent(in) :: lat(:) !< latitude in radians
  real(kind=dp), intent(in) :: exner_l(:,:) !< exner function at model level centers (horizontal, vertical)
  real(kind=dp), intent(in) :: phii(:,:) !< geopotential at model level interfaces (horizontal, vertical)
  real(kind=dp), intent(in) :: pres_l(:,:) !< pressure at model level centers (horizontal, vertical)
  integer, intent(in) :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer, intent(in) :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  real(kind=dp), intent(in) :: relax_time !< time scale for hor. wind nudging
  real(kind=dp), intent(out) :: u_force_tend(:,:) !< sum of all forcing for zonal momentum (horizontal, vertical)
  real(kind=dp), intent(out) :: v_force_tend(:,:) !< sum of all forcing for meridional momentum (horizontal, vertical)
  real(kind=dp), intent(out) :: T_force_tend(:,:) !< sum of all forcing for temperature (horizontal, vertical)
  real(kind=dp), intent(out) :: qv_force_tend(:,:) !< sum of all forcing for water vapor (horizontal, vertical)

  real(kind=dp) :: old_u(n_columns, n_levels), old_v(n_columns, n_levels), old_T(n_columns, n_levels), old_qv(n_columns, n_levels)
  real(kind=dp) :: w_ls_i(n_columns, n_levels+1), zi(n_columns, n_levels+1)
  real(kind=dp) :: theta(n_columns, n_levels)

  integer :: i,k
  real(kind=dp) :: f_coriolis, grav_inv, g_over_cp

  !> \section apply_leapfrog_forcing_alg Algorithm
  !! @{

  grav_inv = 1.0/con_g
  g_over_cp = con_g/con_cp

  !> - Save old state variables (filtered from previous time step)
  old_u = state_u(:,:,1)
  old_v = state_v(:,:,1)
  old_T = state_T(:,:,1)
  old_qv = state_tracer(:,:,1,1)

  theta = old_T/exner_l

  !> - Initialize forcing sums to zero.
  u_force_tend = 0.0
  v_force_tend = 0.0
  T_force_tend = 0.0
  qv_force_tend = 0.0

  !if(.not. nudge_wind .or. .not. nudge_thermo) then
    !>  - Calculate w_ls and z (height) at model layer interfaces.
    do i=1, n_columns
      w_ls_i(i,1) = 0.0
      zi(i,1) = phii(i,1)*grav_inv
      do k=2, n_levels
        w_ls_i(i,k) = 0.5*(w_ls(i,k-1) + w_ls(i,k))
        zi(i,k) = phii(i,k)*grav_inv
      end do
      w_ls_i(i,n_levels+1) = w_ls_i(i,n_levels)
      zi(i,n_levels+1) = phii(i,n_levels+1)*grav_inv
    end do
  !end if

  select case(mom_forcing_type)
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
      do i=1, n_columns
        do k=2, n_levels-1
          u_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_u(i,k+1) - old_u(i,k)) + w_ls_i(i,k)*(old_u(i,k) - old_u(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
          v_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_v(i,k+1) - old_v(i,k)) + w_ls_i(i,k)*(old_v(i,k) - old_v(i,k-1)))/&
            (zi(i,k+1)-zi(i,k))
        end do
        !>  - Handle the top and bottom levels separately using special discretizations.
        u_force_tend(i,1) = -w_ls_i(i,2)*(old_u(i,2) - old_u(i,1))/(zi(i,2)-zi(i,1))
        u_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_u(i,n_levels) - old_u(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))
        v_force_tend(i,1) = -w_ls_i(i,2)*(old_v(i,2) - old_v(i,1))/(zi(i,2)-zi(i,1))
        v_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_v(i,n_levels) - old_v(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))

        !> - Add forcing due to geostrophic wind
        !>  - Calculate Coriolis parameter.
        f_coriolis = 2.0*con_omega*sin(lat(i))
        do k=1, n_levels
          !accumulate forcing tendencies
          u_force_tend(i,k) = u_force_tend(i,k) +  f_coriolis*(old_v(i,k) - v_g(i,k))
          v_force_tend(i,k) = v_force_tend(i,k) -  f_coriolis*(old_u(i,k) - u_g(i,k))
        end do
      end do
    case (3)
      !> - Calculate change in state momentum variables due to nudging.
      do i=1, n_columns
        do k=1, n_levels
          !accumulate forcing tendencies
          u_force_tend(i,k) = (u_nudge(i,k) - old_u(i,k))/relax_time
          v_force_tend(i,k) = (v_nudge(i,k) - old_v(i,k))/relax_time
        end do
      end do
    case default
      u_force_tend = 0.0
      v_force_tend = 0.0
  end select

  select case (thermo_forcing_type)
    case (1)
      do i=1, n_columns
        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, n_levels
          T_force_tend(i,k) = T_force_tend(i,k) + dT_dt_rad(i,k) + exner_l(i,k)*(h_advec_thil(i,k) +v_advec_thil(i,k))
          qv_force_tend(i,k) = qv_force_tend(i,k) + h_advec_qt(i,k) + v_advec_qt(i,k)
        end do
      end do
    case (2)
      do i=1, n_columns
        ! do k=2, n_levels-1
        !
        !   qv_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_qv(i,k+1) - old_qv(i,k)) + w_ls_i(i,k)*(old_qv(i,k) - old_qv(i,k-1)))/&
        !     (zi(i,k+1)-zi(i,k))
        !   !T_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_T(i,k+1) - old_T(i,k)) + w_ls_i(i,k)*(old_T(i,k) - old_T(i,k-1)))/(zi(i,k+1)-zi(i,k)) - &
        !   !  0.5*(w_ls_i(i,k+1) + w_ls_i(i,k))*g_over_cp
        !   T_force_tend(i,k) = exner_l(i,k)*(-0.5*(w_ls_i(i,k+1)*(theta(i,k+1) - theta(i,k)) + &
        !     w_ls_i(i,k)*(theta(i,k) - theta(i,k-1))))/(zi(i,k+1)-zi(i,k))
        ! end do
        ! !>  - Handle the top and bottom levels separately using special discretizations.
        ! qv_force_tend(i,1) = -w_ls_i(i,2)*(old_qv(i,2) - old_qv(i,1))/(zi(i,2)-zi(i,1))
        ! qv_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_qv(i,n_levels) - old_qv(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels))
        ! T_force_tend(i,1) = -w_ls_i(i,2)*(old_T(i,2) - old_T(i,1))/(zi(i,2)-zi(i,1)) !- 0.5*(w_ls_i(2) + w_ls_i(1))*g_over_cp
        ! T_force_tend(i,n_levels) = -w_ls_i(i,n_levels)*(old_T(i,n_levels) - old_T(i,n_levels-1))/(zi(i,n_levels+1)-zi(i,n_levels)) - &
        !   0.5*(w_ls_i(i,n_levels+1) + w_ls_i(i,n_levels))*g_over_cp

        do k=2, n_levels-1
          !upstream scheme
          if(omega(i,k) < 0.0) then
            qv_force_tend(i,k) = -omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
            T_force_tend(i,k) =  exner_l(i,k)*(-omega(i,k)*(theta(i,k) - theta(i,k-1))/(pres_l(i,k)-pres_l(i,k-1)))

          else
            qv_force_tend(i,k) = -omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
            T_force_tend(i,k) = exner_l(i,k)*(-omega(i,k)*(theta(i,k+1) - theta(i,k))/(pres_l(i,k+1)-pres_l(i,k)))

          end if
        end do
        T_force_tend(i,1) = exner_l(i,1)*(-0.5*(omega(i,2)+omega(i,1))*(theta(i,2) - theta(i,1))/(pres_l(i,2)-pres_l(i,1)))
        T_force_tend(i,n_levels) = exner_l(i,n_levels)*(-0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (theta(i,n_levels) - theta(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1)))
        qv_force_tend(i,1) = -0.5*(omega(i,2)+omega(i,1))*(old_qv(i,2) - old_qv(i,1))/(pres_l(i,2)-pres_l(i,1))
        qv_force_tend(i,n_levels) = -0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (old_qv(i,n_levels) - old_qv(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1))

        !> - Add forcing due to prescribed radiation and horizontal advection
        do k=1, n_levels
          T_force_tend(i,k) = T_force_tend(i,k) + dT_dt_rad(i,k) + exner_l(i,k)*h_advec_thil(i,k)
          qv_force_tend(i,k) = qv_force_tend(i,k) + h_advec_qt(i,k)
        end do
      end do
    case (3)
      !> - Calculate change in state temperature/moisture variables due to nudging.
      do i=1, n_columns
        do k=1, n_levels
          !accumulate forcing tendencies
          T_force_tend(i,k) = (T_nudge(i,k) - old_T(i,k))/relax_time
          qv_force_tend(i,k) = (qt_nudge(i,k) - old_qv(i,k))/relax_time
        end do

        do k=2, n_levels-1
          !upstream scheme
          if(omega(i,k) < 0.0) then
            qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
            T_force_tend(i,k) =  T_force_tend(i,k) + &
              exner_l(i,k)*(-omega(i,k)*(theta(i,k) - theta(i,k-1))/(pres_l(i,k)-pres_l(i,k-1)))

          else
            qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
            T_force_tend(i,k) = T_force_tend(i,k) + &
              exner_l(i,k)*(-omega(i,k)*(theta(i,k+1) - theta(i,k))/(pres_l(i,k+1)-pres_l(i,k)))

          end if
        end do
        T_force_tend(i,1) = T_force_tend(i,1) + &
          exner_l(i,1)*(-0.5*(omega(i,2)+omega(i,1))*(theta(i,2) - theta(i,1))/(pres_l(i,2)-pres_l(i,1)))
        T_force_tend(i,n_levels) = T_force_tend(i,n_levels) + exner_l(i,n_levels)*(-0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (theta(i,n_levels) - theta(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1)))
        qv_force_tend(i,1) = qv_force_tend(i,1) - 0.5*(omega(i,2)+omega(i,1))*(old_qv(i,2) - old_qv(i,1))/(pres_l(i,2)-pres_l(i,1))
        qv_force_tend(i,n_levels) = qv_force_tend(i,n_levels) - 0.5*(omega(i,n_levels-1)+omega(i,n_levels))*&
          (old_qv(i,n_levels) - old_qv(i,n_levels-1))/(pres_l(i,n_levels)-pres_l(i,n_levels-1))

      end do

      ! do i=1, n_columns
      !   do k=2, n_levels-1
      !     if(omega(i,k) < 0.0) then !need to bring in pres_l
      !       qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k)-old_qv(i,k-1))/(pres_l(i,k)-pres_l(i,k-1))
      !       T_force_tend(i,k) = T_force_tend(i,k) - omega(i,k)*(old_T(i,k) - old_T(i,k-1))/(pres_l(i,k)-pres_l(i,k-1) + &
      !         old_T(i,k)*con_rd/(pres_l(i,k)*con_cp))
      !     else
      !       qv_force_tend(i,k) = qv_force_tend(i,k) - omega(i,k)*(old_qv(i,k+1)-old_qv(i,k))/(pres_l(i,k+1)-pres_l(i,k))
      !       T_force_tend(i,k) = T_force_tend(i,k) - omega(i,k)*(old_T(i,k+1) - old_T(i,k))/(pres_l(i,k+1)-pres_l(i,k) + &
      !         old_T(i,k)*con_rd/(pres_l(i,k)*con_cp))
      !     end if
      !     ! qv_force_tend(i,k) = qv_force_tend(i,k) - 0.5*(w_ls_i(i,k+1)*(old_qv(i,k+1) - old_qv(i,k)) + &
      !     !   w_ls_i(i,k)*(old_qv(i,k) - old_qv(i,k-1)))/(zi(i,k+1)-zi(i,k))
      !     ! !T_force_tend(i,k) = -0.5*(w_ls_i(i,k+1)*(old_T(i,k+1) - old_T(i,k)) + w_ls_i(i,k)*(old_T(i,k) - old_T(i,k-1)))/(zi(i,k+1)-zi(i,k)) - &
      !     ! !  0.5*(w_ls_i(i,k+1) + w_ls_i(i,k))*g_over_cp
      !     ! T_force_tend(i,k) = T_force_tend(i,k) + exner_l(i,k)*(-0.5*(w_ls_i(i,k+1)*(theta(i,k+1) - theta(i,k)) + &
      !     !   w_ls_i(i,k)*(theta(i,k) - theta(i,k-1))))/(zi(i,k+1)-zi(i,k))
      !   end do
      !   !>  - Handle the top and bottom levels separately using special discretizations.
      !   ! qv_force_tend(i,1) = qv_force_tend(i,1) - w_ls_i(i,2)*(old_qv(i,2) - old_qv(i,1))/(zi(i,2)-zi(i,1))
      !   ! qv_force_tend(i,n_levels) = qv_force_tend(i,n_levels) - w_ls_i(i,n_levels)*(old_qv(i,n_levels) - old_qv(i,n_levels-1))/&
      !   !   (zi(i,n_levels+1)-zi(i,n_levels))
      !   ! T_force_tend(i,1) = T_force_tend(i,1) - w_ls_i(i,2)*(old_T(i,2) - old_T(i,1))/(zi(i,2)-zi(i,1)) !- 0.5*(w_ls_i(2) + w_ls_i(1))*g_over_cp
      !   ! T_force_tend(i,n_levels) = T_force_tend(i,n_levels) - w_ls_i(i,n_levels)*(old_T(i,n_levels) - old_T(i,n_levels-1))/&
      !   !   (zi(i,n_levels+1)-zi(i,n_levels)) - 0.5*(w_ls_i(i,n_levels+1) + w_ls_i(i,n_levels))*g_over_cp
      !
      ! end do
    case default
      T_force_tend = 0.0
      qv_force_tend = 0.0
  end select

  do i=1, n_columns
    do k=1, n_levels
      !> - Update the state variables using the forward Euler scheme:
      !!   \f[
      !!   x^{\tau + 1} = x^{\tau} + \Delta t\frac{\partial x}{\partial t}|^\tau_{forcing}
      !!   \f]
      !!   \f$x^{\tau}\f$ is the value at the previous time step and \f$\frac{\partial x}{\partial t}|^\tau_{forcing}\f$ is the sum of forcing terms calculated in this time step.
      state_u(i,k,1) = old_u(i,k) + dt*u_force_tend(i,k)
      state_v(i,k,1) = old_v(i,k) + dt*v_force_tend(i,k)
      state_T(i,k,1) = state_T(i,k,1) + dt*(T_force_tend(i,k))
      state_tracer(i,k,1,1) = state_tracer(i,k,1,1) + dt*(qv_force_tend(i,k))
    end do
  end do
  !> @}

end subroutine apply_forcing_forward_Euler

!> @}
!> @}
end module gmtb_scm_forcing
