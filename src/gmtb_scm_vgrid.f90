!> \file gmtb_scm_vgrid.f90
!!  Contains the vertical grid setup routines.

module gmtb_scm_vgrid

use gmtb_scm_kinds, only: sp, dp, qp
use gmtb_scm_physical_constants, only : con_cp, con_rocp, con_fvirt

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup vgrid gmtb_scm_vgrid
!! @{
!! Contains the vertical grid setup routines.

!> Subroutine for setting up the GFS hybrid coordinate vertical grid. Files with precalculated coefficients (A_k and B_k) from the
!! "fix" directory in Patrick Tripp's V2 code are used; their filenames start with "global_hyblev". There are only files for 28, 42,
!! 60, 64, and 91 levels. If a different number of levels are specified, an error is returned. Using the A_k and B_k coefficients, the
!! model level pressures, sigma values, and exner function (at interfaces and layer centers) are calculated and returned to the calling
!! procedure.
subroutine get_GFS_vgrid(pres_sfc, n_levels, n_columns, pres_i, pres_l, sigi, sigl, exner_l, exner_i, a_k, b_k, error)
  !create GFS hybrid coordinate vertical grid
  integer, intent(in) :: n_levels !< number of model levels
  integer, intent(in) :: n_columns !< number of columns
  real(kind=dp), intent(in) :: pres_sfc !< surface pressure (Pa)
  real(kind=dp), allocatable, intent(out) :: pres_l(:,:)  !< pressure at model level centers
  real(kind=dp), allocatable, intent(out) :: pres_i(:,:) !< pressure at model interfaces
  real(kind=dp), allocatable, intent(out) :: sigi(:,:) !< sigma at model interfaces
  real(kind=dp), allocatable, intent(out) :: sigl(:,:) !< sigma at model layers
  real(kind=dp), allocatable, intent(out) :: exner_l(:,:) !< exner function at model level centers
  real(kind=dp), allocatable, intent(out) :: exner_i(:,:) !< exner function at model interfaces
  real(kind=dp), allocatable, intent(out) :: a_k(:) !< GFS hybrid pressure coordinate coefficients
  real(kind=dp), allocatable, intent(out) :: b_k(:) !< GFS hybrid pressure coordinate coefficients
  integer, intent(out) :: error !< error code

  integer               :: i, ierr

  real(kind=dp)               :: pres_sfc_inv, p0
  character(len=37)     :: filename
  character(len=80)     :: line
  character(len=16)     :: file_format

  !>  \section get_GFS_vgrid_alg Algorithm
  !!  @{

  error = 0

  !file format for all but the 91 level file is the same
  file_format = '(1F12.3, 1F12.8)'


  !> - Check to see if the desired number of grid levels is available. If not, return an error code.
  select case (n_levels)
    case(28)
      filename = "../model_config/global_hyblev.l28.txt"
    case(42)
      filename = "../model_config/global_hyblev.l42.txt"
    case(60)
      filename = "../model_config/global_hyblev.l60.txt"
    case(64)
      filename = "../model_config/global_hyblev.l64.txt"
    case(91)
      filename = "../model_config/global_hyblev.l91.txt"
      !file format for the 91-level file is different
      file_format = '(1F14.6, 1F10.6)'
    case default
      error = 1
      return
  end select

  allocate(pres_l(n_columns, n_levels), pres_i(n_columns, n_levels+1), exner_l(n_columns, n_levels), &
    exner_i(n_columns, n_levels+1))
  allocate(a_k(n_levels+1), b_k(n_levels+1), sigl(n_columns, n_levels), sigi(n_columns, n_levels+1))

  !> - Open the appropriate file.
  open(unit=1, file=filename, status='old', action='read', iostat=ierr)
  if(ierr /= 0) then
    write(*,*) 'There was an error opening the file ', filename, ' in the model_config directory. &
      Error code = ',ierr
    error = 2
  endif

  !> - The first line contains the number of coefficients and number of levels (these should already by known; discard this info)
  read(1,'(a)',iostat=ierr) line
  !> - Read in the coefficient data.
  do i=1, n_levels+1
    read(1,file_format)  a_k(i), b_k(i)
  end do
  close(1)

  !> - Calculate interface pressures, sigma, and exner function.

  p0 = pres_sfc
  pres_sfc_inv = 1.0/p0
  do i=1, n_levels+1
    pres_i(:,i) = a_k(i) + b_k(i)*p0
    sigi(:,i) = a_k(i)*pres_sfc_inv + b_k(i)
    exner_i(:,i) = (pres_i(:,i)/1.0E5)**con_rocp
  end do

  !> - Calculate layer center pressures, sigma, and exner function.
  do i=1, n_levels
    pres_l(:,i) = ((1.0/(con_rocp+1.0))*(pres_i(:,i)**(con_rocp+1.0) - pres_i(:,i+1)**(con_rocp+1.0))/ &
      (pres_i(:,i) - pres_i(:,i+1)))**(1.0/con_rocp)
    sigl(:,i) = 0.5*(sigi(:,i) + sigi(:,i+1))

    exner_l(:,i) = (pres_l(:,i)/1.0E5)**con_rocp

  end do
  !> @}
end subroutine get_GFS_vgrid

!> This subroutine calculates the pressure and exner function at grid centers and interface levels given a surface pressure and interface-level GFS grid coefficients.
!! This subroutine should be called to update the pressures of the model levels as the surface pressure of the column changes.
subroutine calc_GFS_pres_exner_geopotential(pres_sfc, n_levels, n_columns, a_k, b_k, T, qv, pres_i, pres_l, sigi, sigl, exner_l, &
  exner_i, geopotential_l, geopotential_i)

  integer, intent(in) :: n_levels !< number of model levels
  integer, intent(in) :: n_columns !< number of columns
  real(kind=dp), intent(in) :: pres_sfc(:) !< surface pressure (Pa)
  real(kind=dp), intent(in) :: a_k(:), b_k(:) !< GFS hybrid pressure coordinate coefficients
  real(kind=dp), intent(in) :: T(:,:), qv(:,:)
  real(kind=dp), intent(inout) :: pres_l(:,:)  !< pressure at model level centers
  real(kind=dp), intent(inout) :: pres_i(:,:) !< pressure at model interfaces
  real(kind=dp), intent(inout) :: sigi(:,:) !< sigma at model interfaces
  real(kind=dp), intent(inout) :: sigl(:,:) !< sigma at model layers
  real(kind=dp), intent(inout) :: exner_l(:,:) !< exner function at model level centers
  real(kind=dp), intent(inout) :: exner_i(:,:) !< exner function at model interfaces
  real(kind=dp), intent(inout) :: geopotential_l(:,:) !< geopotential at model level centers
  real(kind=dp), intent(inout) :: geopotential_i(:,:) !< geopotential function at model interfaces

  real(kind=dp)               :: pres_sfc_inv, tem, dgeopotential_lower_half, dgeopotential_upper_half
  integer               :: i,k

  !> - Calculate interface pressures, sigma, and exner function.
  do i=1, n_columns
    pres_sfc_inv = 1.0/pres_sfc(i)
    do k=1, n_levels+1
      pres_i(i,k) = a_k(k) + b_k(k)*pres_sfc(i)
      sigi(i,k) = a_k(k)*pres_sfc_inv + b_k(k)
      exner_i(i,k) = (pres_i(i,k)*1.0E-5)**con_rocp
    end do
  end do

  !> - Calculate layer center pressures, sigma, and exner function.
  do i=1, n_columns
    do k=1, n_levels
      pres_l(i,k) = ((1.0/(con_rocp+1.0))*(pres_i(i,k)**(con_rocp+1.0) - pres_i(i,k+1)**(con_rocp+1.0))/ &
        (pres_i(i,k) - pres_i(i,k+1)))**(1.0/con_rocp)
      sigl(i,k) = 0.5*(sigi(i,k) + sigi(i,k+1))
      exner_l(i,k) = (pres_l(i,k)*1.0E-5)**con_rocp
    end do
  end do

  do i=1, n_columns
    geopotential_i(i,1) = 0.0
  end do

  do i=1, n_columns
    do k=1, n_levels
      tem = con_cp*T(i,k)*(1.0 + con_fvirt*max(qv(i,k), 0.0))/exner_l(i,k)
      dgeopotential_lower_half = (exner_i(i,k) - exner_l(i,k))*tem
      dgeopotential_upper_half = (exner_l(i,k) - exner_i(i,k+1))*tem
      geopotential_l(i,k) = geopotential_i(i,k) + dgeopotential_lower_half
      geopotential_i(i,k+1) = geopotential_l(i,k) + dgeopotential_upper_half
    end do
  end do

end subroutine calc_GFS_pres_exner_geopotential

subroutine calc_GFS_geopotential(n_levels, n_columns, T, qv, exner_i, exner_l, geopotential_i, geopotential_l)
  integer, intent(in) :: n_levels, n_columns
  real(kind=dp), intent(in) :: T(:,:), qv(:,:)
  real(kind=dp), intent(in) :: exner_i(:,:), exner_l(:,:)
  real(kind=dp), allocatable, intent(out) :: geopotential_i(:,:), geopotential_l(:,:)

  integer i,k
  real(kind=dp) :: tem, dgeopotential_lower_half, dgeopotential_upper_half

  allocate(geopotential_l(n_columns, n_levels), geopotential_i(n_columns, n_levels+1))

  do i=1, n_columns
    geopotential_i(i,1) = 0.0
  end do

  do i=1, n_columns
    do k=1, n_levels
      tem = con_cp*T(i,k)*(1.0 + con_fvirt*max(qv(i,k), 0.0))/exner_l(i,k)
      dgeopotential_lower_half = (exner_i(i,k) - exner_l(i,k))*tem
      dgeopotential_upper_half = (exner_l(i,k) - exner_i(i,k+1))*tem
      geopotential_l(i,k) = geopotential_i(i,k) + dgeopotential_lower_half
      geopotential_i(i,k+1) = geopotential_l(i,k) + dgeopotential_upper_half
    end do
  end do

end subroutine calc_GFS_geopotential
!> @}
!> @}
end module gmtb_scm_vgrid
