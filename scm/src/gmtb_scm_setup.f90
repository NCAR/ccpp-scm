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
subroutine set_state(scm_input, scm_reference, scm_state)
  use gmtb_scm_type_defs, only : scm_input_type, scm_reference_type, scm_state_type

  type(scm_input_type), intent(in) :: scm_input
  type(scm_reference_type), intent(in) :: scm_reference
  type(scm_state_type), intent(inout) :: scm_state

  integer :: i,j, last_index_init, grid_error
  real(kind=dp), dimension(scm_input%input_nlev) :: input_qv, input_T
  real(kind=dp), parameter :: p0 = 1.0E5

  !> \section set_state_alg Algorithm
  !! @{

  !> - Calculate water vapor from total water, suspended liquid water, and suspended ice.
  input_qv = scm_input%input_qt - scm_input%input_ql - scm_input%input_qi

  !> - \todo When patching in a reference sounding, need to handle the case when the reference sounding is too short; patch_in_ref
  !! checks for the case, but as of now, it just extrapolates where it needs to and returns an error code; error should be handled
  !! here or passed up to the main program.

  !> - For each column, interpolate the water vapor to the model grid.
  do i=1, scm_state%n_cols
    call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, input_qv, scm_state%pres_l(i,1,:), &
      scm_state%n_levels, scm_state%state_tracer(i,1,:,scm_state%water_vapor_index,1), last_index_init, 1)
    !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
    if(last_index_init < scm_state%n_levels) THEN
      call patch_in_ref(last_index_init, scm_state%n_levels_smooth, scm_reference%ref_nlev, scm_reference%ref_pres, &
        scm_reference%ref_qv, scm_state%pres_l(i,1,:), scm_state%n_levels, &
        scm_state%state_tracer(i,1,:,scm_state%water_vapor_index,1), grid_error)
    end if
  end do

  !> - Calculate the input absolute temperature from input pressure, theta_il, ql, and qi.
  input_T = (scm_input%input_pres/p0)**con_rocp*(scm_input%input_thetail + (con_hvap/con_cp)*scm_input%input_ql + &
    (con_hfus/con_cp)*scm_input%input_qi)

  !> - For each column, interpolate the temperature to the model grid.
  do i=1, scm_state%n_cols
    call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, input_T, scm_state%pres_l(i,1,:), &
      scm_state%n_levels, scm_state%state_T(i,1,:,1), last_index_init, 1)
    !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
    if(last_index_init < scm_state%n_levels) THEN
      call patch_in_ref(last_index_init, scm_state%n_levels_smooth, scm_reference%ref_nlev, scm_reference%ref_pres, &
        scm_reference%ref_T, scm_state%pres_l(i,1,:), scm_state%n_levels, scm_state%state_T(i,1,:,1), grid_error)
    end if
  end do

  !> - For each column, interpolate the u-wind to the model grid.
  do i=1, scm_state%n_cols
    call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_u, scm_state%pres_l(i,1,:), &
      scm_state%n_levels, scm_state%state_u(i,1,:,1), last_index_init, 1)
    if(last_index_init < scm_state%n_levels) THEN
      do j=last_index_init + 1, scm_state%n_levels
        !>  - The standard atmosphere doesn't have wind data; assume zero-gradient above the input data.
        scm_state%state_u(i,1,j,1) = scm_state%state_u(i,1,last_index_init,1)
      end do
    end if
  end do

  !> - For each column, interpolate the v-wind to the model grid.
  do i=1, scm_state%n_cols
    call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_v, scm_state%pres_l(i,1,:), &
      scm_state%n_levels, scm_state%state_v(i,1,:,1), last_index_init, 1)
    if(last_index_init < scm_state%n_levels) THEN
      do j=last_index_init + 1, scm_state%n_levels
        !>  - The standard atmosphere doesn't have wind data; assume zero-gradient above the input data.
        scm_state%state_v(i,1,j,1) = scm_state%state_v(i,1,last_index_init,1)
      end do
    end if
  end do

  !> - For each column, interpolate the ozone to the model grid.
  if(scm_state%ozone_index > 0) then
    do i=1, scm_state%n_cols
      call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_ozone, scm_state%pres_l(i,1,:), &
        scm_state%n_levels, scm_state%state_tracer(i,1,:,scm_state%ozone_index,1), last_index_init, 1)
      !>  - If the input domain does not span the model domain, patch in McClatchey tropical standard atmosphere (smoothly over a number of levels) above.
      if(last_index_init < scm_state%n_levels) THEN
        call patch_in_ref(last_index_init, scm_state%n_levels_smooth, scm_reference%ref_nlev, scm_reference%ref_pres, &
          scm_reference%ref_ozone, scm_state%pres_l(i,1,:), scm_state%n_levels, &
          scm_state%state_tracer(i,1,:,scm_state%ozone_index,1), grid_error)
      end if
    end do
  end if

  scm_state%state_tracer(:,1,:,scm_state%cloud_water_index,1) = 0.0
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

!--------------
! GFS initialze
!--------------
subroutine GFS_suite_setup (Model, Statein, Stateout, Sfcprop,                   &
                            Coupling, Grid, Tbd, Cldprop, Radtend, Diag,         &
                            Interstitial, Init_parm, n_ozone_lats,               &
                            n_ozone_layers, n_ozone_times, n_ozone_coefficients, &
                            ozone_lat, ozone_pres, ozone_time, ozone_forcing_in)

  use machine,             only: kind_phys
  use GFS_typedefs,        only: GFS_init_type,                          &
                                 GFS_statein_type,  GFS_stateout_type,   &
                                 GFS_sfcprop_type,  GFS_coupling_type,   &
                                 GFS_control_type,  GFS_grid_type,       &
                                 GFS_tbd_type,      GFS_cldprop_type,    &
                                 GFS_radtend_type,  GFS_diag_type,       &
                                 GFS_interstitial_type
  use funcphys,            only: gfuncphys
  use module_microphysics, only: gsmconst
  use cldwat2m_micro,      only: ini_micro
  use aer_cloud,           only: aer_cloud_init
  use module_ras,          only: ras_init
  use ozne_def,            only: latsozp, levozp, timeoz, oz_coeff, oz_lat, oz_pres, oz_time, ozplin
  use GFS_rrtmg_setup,     only: GFS_rrtmg_setup_init

  !--- interface variables
  type(GFS_control_type),      intent(inout) :: Model
  type(GFS_statein_type),      intent(inout) :: Statein
  type(GFS_stateout_type),     intent(inout) :: Stateout
  type(GFS_sfcprop_type),      intent(inout) :: Sfcprop
  type(GFS_coupling_type),     intent(inout) :: Coupling
  type(GFS_grid_type),         intent(inout) :: Grid
  type(GFS_tbd_type),          intent(inout) :: Tbd
  type(GFS_cldprop_type),      intent(inout) :: Cldprop
  type(GFS_radtend_type),      intent(inout) :: Radtend
  type(GFS_diag_type),         intent(inout) :: Diag
  type(GFS_interstitial_type), intent(inout) :: Interstitial
  type(GFS_init_type),         intent(in)    :: Init_parm

  integer, intent(in) :: n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times
  real(kind=kind_phys), intent(in) :: ozone_lat(:), ozone_pres(:), ozone_time(:), ozone_forcing_in(:,:,:,:)

  !--- local variables
  real(kind=kind_phys), parameter   :: p_ref = 101325.0d0

  !--- set control properties (including namelist read)
  call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                   Init_parm%me, Init_parm%master,               &
                   Init_parm%logunit, Init_parm%isc,             &
                   Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                   Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                   Init_parm%gnx, Init_parm%gny,                 &
                   Init_parm%dt_dycore, Init_parm%dt_phys,       &
                   Init_parm%bdat, Init_parm%cdat,               &
                   Init_parm%tracer_names, Init_parm%blksz)

  !--- allocate memory for the variables stored in ozne_def and set them
  allocate(oz_lat(n_ozone_lats), oz_pres(n_ozone_layers), oz_time(n_ozone_times+1))
  allocate(ozplin(n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times))
  latsozp = n_ozone_lats
  levozp = n_ozone_layers
  timeoz = n_ozone_times
  oz_coeff = n_ozone_coefficients
  oz_lat = ozone_lat
  oz_pres = ozone_pres
  oz_time = ozone_time
  ozplin = ozone_forcing_in

  !--- initialize DDTs
  call Statein%create(1, Model)
  call Stateout%create(1, Model)
  call Sfcprop%create(1, Model)
  call Coupling%create(1, Model)
  call Grid%create(1, Model)
  call Tbd%create(1, 1, Model)
  call Cldprop%create(1, Model)
  call Radtend%create(1, Model)
  !--- internal representation of diagnostics
  call Diag%create(1, Model)
  !--- internal representation of interstitials for CCPP physics
  call Interstitial%create(1, Model)

  !--- populate the grid components
  call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

  !--- read in and initialize ozone and water
  if (Model%ntoz > 0) then
    call setindxoz (Init_parm%blksz, Grid%xlat_d, Grid%jindx1_o3, &
                      Grid%jindx2_o3, Grid%ddy_o3)
  endif

  if (Model%h2o_phys) then
    call setindxh2o (Init_parm%blksz, Grid%xlat_d, Grid%jindx1_h, &
                       Grid%jindx2_h, Grid%ddy_h)
  endif

  !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
  call gfuncphys ()

  call gsmconst (Model%dtp, Model%me, .TRUE.)

  !--- define sigma level for radiation initialization
  !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
  !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
  !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
  Model%si = (Init_parm%ak + Init_parm%bk * p_ref - Init_parm%ak(Model%levr+1)) &
           / (p_ref - Init_parm%ak(Model%levr+1))

  !--- initialize Morrison-Gettleman microphysics
  if (Model%ncld == 2) then
    call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice)
    call aer_cloud_init ()
  endif

  !--- initialize ras
  if (Model%ras) call ras_init (Model%levs, Model%me)

  !--- initialize soil vegetation
  call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)

  !--- lsidea initialization
  if (Model%lsidea) then
    print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
    stop
    !--- NEED TO get the logic from the old phys/gloopb.f initialization area
  endif

  !--- sncovr may not exist in ICs from chgres.
  !--- FV3GFS handles this as part of the IC ingest
  !--- this not is placed here to alert users to the need to study
  !--- the FV3GFS_io.F90 module

end subroutine GFS_suite_setup

!------------------
! GFS_grid_populate
!------------------
subroutine GFS_grid_populate (Grid, xlon, xlat, area)
  use machine,             only: kind_phys
  use physcons,            only: pi => con_pi
  use GFS_typedefs,        only: GFS_grid_type

  implicit none

  type(GFS_grid_type)              :: Grid
  real(kind=kind_phys), intent(in) :: xlon(:,:)
  real(kind=kind_phys), intent(in) :: xlat(:,:)
  real(kind=kind_phys), intent(in) :: area(:,:)

  !--- local variables
  integer :: n_columns, i

  n_columns = size(Grid%xlon)

  do i=1, n_columns
   Grid%xlon = xlon(i,1)
   Grid%xlat = xlat(i,1)
   Grid%xlat_d(i) = xlat(i,1) * 180.0_kind_phys/pi
   Grid%sinlat(i) = sin(Grid%xlat(i))
   Grid%coslat(i) = sqrt(1.0_kind_phys - Grid%sinlat(i)*Grid%sinlat(i))
   Grid%area(i)   = area(i,1)
   Grid%dx(i)     = sqrt(area(i,1))
  end do

end subroutine GFS_grid_populate

end module gmtb_scm_setup
