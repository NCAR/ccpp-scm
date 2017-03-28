!> \file gmtb_scm_setup.f90
!!  Contains subroutines to initialize the SCM, including setting the atmospheric state, interpolating initial conditions to the
!!  model grid, and patching in a reference sounding above the provided initial conditions.

module gmtb_scm_setup

use gmtb_scm_kinds, only: sp, dp, qp
use gmtb_scm_physical_constants, only: con_hvap, con_hfus, con_cp, con_rocp
use gmtb_scm_utils, only: interpolate_to_grid_centers

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup setup gmtb_scm_setup
!! @{
!! Contains subroutines to initialize the SCM, including setting the atmospheric state, interpolating initial conditions to the
!!  model grid, and patching in a reference sounding above the provided initial conditions.

!> Subroutine to interpolate the initial conditions to the model grid and set the state variables.
subroutine set_state(n_input_levels, input_pres, input_qt, input_thetail, input_ql, input_qi, input_u, input_v, input_ozone, &
  n_model_levels, n_columns, ozone_index, cloud_water_index, pres_l, n_smooth_levels, n_ref_levels, ref_pres, ref_qv, ref_T, &
  ref_ozone, state_tracer, state_T, state_u, state_v)

  integer, intent(in) :: n_input_levels !< number of levels in the case input data
  integer, intent(in) :: n_model_levels !< number of model levels
  integer, intent(in) :: n_columns !< number of columns (should be 1 for SCM)
  integer, intent(in) :: n_smooth_levels !< number of model levels over which to smoothly transition to the reference profile
  integer, intent(in) :: n_ref_levels !< number of levels in the reference profile
  integer, intent(in) :: ozone_index !< tracer array index for ozone
  integer, intent(in) :: cloud_water_index !< tracer array index for cloud water

  !input profile variables
  real(kind=dp), intent(in) :: input_pres(:) !< pressure levels for the case input (Pa)
  real(kind=dp), intent(in) :: input_qt(:) !< total water specific humidity for the case input (kg kg^-1)
  real(kind=dp), intent(in) :: input_thetail(:) !< liquid-ice water potential temperature for the case input (K)
  real(kind=dp), intent(in) :: input_ql(:) !< liquid water specific humidity for the case input (kg kg^-1)
  real(kind=dp), intent(in) :: input_qi(:) !< ice water specific humidity for the case input (kg kg^-1)
  real(kind=dp), intent(in) :: input_u(:) !< east-west wind for the case input (m s^-1)
  real(kind=dp), intent(in) :: input_v(:) !< north-south wind for the case input (m s^-1)
  real(kind=dp), intent(in) :: input_ozone(:) !< ozone concentration for the case input (kg kg^-1)


  real(kind=dp), intent(in) :: pres_l(:,:) !< model pressure levels (Pa) (n_columns, n_model_levels)

  !reference profile variables
  real(kind=dp), intent(in) :: ref_pres(:) !< pressure levels for the reference profile (Pa)
  real(kind=dp), intent(in) :: ref_qv(:) !< water vapor specific humidity for the reference profile (kg kg^-1)
  real(kind=dp), intent(in) :: ref_T(:) !< temperature for the reference profile (K)
  real(kind=dp), intent(in) :: ref_ozone(:) !< ozone concentratin for the reference profile (kg kg^-1)

  !state variables to be calculated (n_columns, n_model_levels)
  real(kind=dp), intent(out) :: state_tracer(:,:,:) !< water vapor specific humidity for the model state, ozone, cloud water [all (kg kg^-1)]
  real(kind=dp), intent(out) :: state_T(:,:) !< temperature for the model state (K)
  real(kind=dp), intent(out) :: state_u(:,:)  !< east-west wind for the model state (m s^-1)
  real(kind=dp), intent(out) :: state_v(:,:) !< north-south wind for the model state (m s^-1)

  integer :: i,j, last_index_init, grid_error
  real(kind=dp), dimension(n_input_levels) :: input_qv, input_T
  real(kind=dp), parameter :: p0 = 1.0E5

  !> \section set_state_alg Algorithm
  !! @{

  !> - Calculate water vapor from total water, suspended liquid water, and suspended ice.
  input_qv = input_qt - input_ql - input_qi

  !> - \todo When patching in a reference sounding, need to handle the case when the reference sounding is too short; patch_in_ref
  !! checks for the case, but as of now, it just extrapolates where it needs to and returns an error code; error should be handled
  !! here or passed up to the main program.

  !> - For each column, interpolate the water vapor to the model grid.
  do i=1, n_columns
    call interpolate_to_grid_centers(n_input_levels, input_pres, input_qv, pres_l(i,:), n_model_levels, state_tracer(i,:,1), &
      last_index_init, 1)
    !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
    if(last_index_init < n_model_levels) THEN
      call patch_in_ref(last_index_init, n_smooth_levels, n_ref_levels, ref_pres, ref_qv, pres_l(i,:), n_model_levels, &
        state_tracer(i,:,1), grid_error)
    end if
  end do

  !> - Calculate the input absolute temperature from input pressure, theta_il, ql, and qi.
  input_T = (input_pres/p0)**con_rocp*(input_thetail + (con_hvap/con_cp)*input_ql + (con_hfus/con_cp)*input_qi)

  !> - For each column, interpolate the temperature to the model grid.
  do i=1, n_columns
    call interpolate_to_grid_centers(n_input_levels, input_pres, input_T, pres_l(i,:), n_model_levels, state_T(i,:), &
      last_index_init, 1)
    !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
    if(last_index_init < n_model_levels) THEN
      call patch_in_ref(last_index_init, n_smooth_levels, n_ref_levels, ref_pres, ref_T, pres_l(i,:), n_model_levels, state_T(i,:),&
        grid_error)
    end if
  end do

  !> - For each column, interpolate the u-wind to the model grid.
  do i=1, n_columns
    call interpolate_to_grid_centers(n_input_levels, input_pres, input_u, pres_l(i,:), n_model_levels, state_u(i,:), &
      last_index_init, 1)
    if(last_index_init < n_model_levels) THEN
      do j=last_index_init + 1, n_model_levels
        !>  - The standard atmosphere doesn't have wind data; assume zero-gradient above the input data.
        state_u(i,j) = state_u(i,last_index_init)
      end do
    end if
  end do

  !> - For each column, interpolate the v-wind to the model grid.
  do i=1, n_columns
    call interpolate_to_grid_centers(n_input_levels, input_pres, input_v, pres_l(i,:), n_model_levels, state_v(i,:), &
      last_index_init, 1)
    if(last_index_init < n_model_levels) THEN
      do j=last_index_init + 1, n_model_levels
        !>  - The standard atmosphere doesn't have wind data; assume zero-gradient above the input data.
        state_v(i,j) = state_v(i,last_index_init)
      end do
    end if
  end do

  !> - For each column, interpolate the ozone to the model grid.
  if(ozone_index > 0) then
    do i=1, n_columns
      call interpolate_to_grid_centers(n_input_levels, input_pres, input_ozone, pres_l(i,:), n_model_levels, &
        state_tracer(i,:,ozone_index), last_index_init, 1)
      !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
      if(last_index_init < n_model_levels) THEN
        call patch_in_ref(last_index_init, n_smooth_levels, n_ref_levels, ref_pres, ref_ozone, pres_l(i,:), n_model_levels, &
          state_tracer(i,:,ozone_index), grid_error)
      end if
    end do
  end if

  state_tracer(:,:,cloud_water_index) = 0.0
  !> @}
end subroutine set_state

!> Subroutine to patch in a reference profile (smoothly) above a given model level.
subroutine patch_in_ref(last_index_init, n_levels_smooth, n_ref_levels, ref_pres, ref_field, pres_l, n_model_levels, model_field, &
  err)
  integer, intent(inout) :: last_index_init !< last vertical model index that was initialized
  integer, intent(in) :: n_levels_smooth !< number of levels below the last initialized model level to start smoothing in the reference profile
  integer, intent(in) :: n_ref_levels !< number of reference profile levels
  integer, intent(in) :: n_model_levels !< number of model levels
  real(kind=dp), intent(in) :: ref_pres(:) !< reference profile pressure levels (Pa)
  real(kind=dp), intent(in) :: ref_field(:) !< reference profile data
  real(kind=dp), intent(in) :: pres_l(:) !< model pressure levels (Pa)
  real(kind=dp), intent(inout) :: model_field(:) !< model profile data
  integer, intent(out) :: err !< error code to return to calling procedure (0: no error, 1: reference profile shorter than input, 2: reference profile shorter than model domain)

  integer i, j, last_initialized
  real(kind=dp) :: smooth_frac, lifrac, ref_val_model, gradient
  logical :: found

  !> \section patch_in_ref_alg Algorithm
  !! @{

  !> - Check for the reference profile being shorter than the input (error code 1), or shorter than the model domain (error code 2).
  err = 0
  if (ref_pres(n_ref_levels) > pres_l(last_index_init)) then
    !reference profile shorter than input
    err = 1
  elseif (ref_pres(n_ref_levels) > pres_l(n_model_levels)) then
    !reference profile shorter than model domain
    err = 2
  end if

  !> - Recalculate the model data profile n_levels_smooth below the last index calculated in order to provided a smoother transition to the reference profile.
  if (err /= 1) then
    do i=last_index_init - (n_levels_smooth-1), last_index_init
      !>  - smooth_frac is the fraction of the reference value to use (remainder is model field already calculated)
      smooth_frac = 1 - (last_index_init - i + 1)/REAL(n_levels_smooth+1)

      !>  - First, interpolate reference profile to model levels.
      found = .false.
      do j=1, n_ref_levels
        if (ref_pres(j) <= pres_l(i)) then
          found = .true.
          lifrac = (pres_l(i) - ref_pres(j-1))/(ref_pres(j)-ref_pres(j-1))
          ref_val_model = ref_field(j-1) + lifrac*(ref_field(j)-ref_field(j-1))
          exit
        end if
      end do

      ! !if the reference pressure levels do not extend to the top of the input profile, return an error code (user should ensure that the
      ! ! reference profile extends past the input profile at a bare minimum.)
      ! if (.not. found) then
      !   !reference sounding does not specify values up to highest model level
      !   write(*,*) 'The reference sounding patched in to the initial profile does not specify values above model level',i-1
      !   write(*,*) 'The SCM can not be initialized since the combination of the case initial profile and the reference profile do not&
      !     extend to the top of the model domain.'
      !   err = 1
      !   return
      ! end if

      !>  - Then, recalculate field at model level using smooth_frac
      model_field(i) = smooth_frac*ref_val_model + (1.0 - smooth_frac)*model_field(i)
    end do

    !> - Above the highest model level specified in the input, just interpolate to the reference profile.
    last_initialized = n_model_levels
    do i=last_index_init + 1, n_model_levels
      found = .false.
      do j=1, n_ref_levels
        if (ref_pres(j) <= pres_l(i)) then
          found = .true.
          lifrac = (pres_l(i) - ref_pres(j-1))/(ref_pres(j)-ref_pres(j-1))
          model_field(i) = ref_field(j-1) + lifrac*(ref_field(j)-ref_field(j-1))
          exit
        end if
      end do
      !>  - If the reference profile does not proceed to the top of the model domain, keep track of the highest level initialized.
      if (.not. found) then
        write(*,*) 'The reference sounding patched in to the initial profile does not specify values above model level',i-1
        write(*,*) 'lowest reference pressure level =',ref_pres(n_ref_levels)
        write(*,*) 'lowest model pressure level =',pres_l(n_model_levels)
        last_initialized = i-1

        exit
      end if
    end do
  end if

  last_index_init = last_initialized
  !>  - If the reference profile was shorter than the model domain, extrapolate to fill in the missing values (temporary solution).
  if(err /= 0) then
    gradient = (model_field(last_initialized) - model_field(last_initialized-1))/ &
      (pres_l(last_initialized) - pres_l(last_initialized-1))
    write(*,*) 'temporarily extrapolating profile above level',last_initialized
    do i=last_initialized+1, n_model_levels
      !determine gradient over the last two initialized levels
      model_field(i) = model_field(i-1) + gradient*(pres_l(i)-pres_l(i-1))
    end do
  end if
  !> @}
end subroutine patch_in_ref
!> @}
!> @}
end module gmtb_scm_setup
