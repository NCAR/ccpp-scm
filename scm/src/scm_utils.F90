!> \file scm_utils.f90
!!  Contains miscellaneous helper subroutines.

module scm_utils

use scm_kinds, only: sp, dp, qp
use scm_physical_constants, only: con_rd, con_g

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup utils scm_utils
!! @{
!! Contains miscellaneous helper subroutines.

!> Subroutine to interpolate a generic input field to the model vertical grid. The index of the last model level initialized is returned
!! so that the calling procedure can patch in a reference profile above.
subroutine interpolate_to_grid_centers(n_input_levels, input_pres, input_field, pres_l, n_model_levels, model_field, &
  last_initialized, below_input_level_method)
  integer, intent(in) :: n_input_levels !< number of levels in the input profile
  integer, intent(in) :: n_model_levels !< number of model levels
  real(kind=dp), intent(in) :: input_pres(:)  !< pressure levels for the input profile (Pa)
  real(kind=dp), intent(in) :: input_field(:) !< data from the input profile
  real(kind=dp), intent(in) :: pres_l(:) !< model pressure levels (Pa)
  real(kind=dp), intent(out):: model_field(:) !< interpolated model data profile
  integer, intent(out):: last_initialized !< last vertical model index that was properly initialized
  integer, intent(in) :: below_input_level_method !< specified method for handling model levels below the input levels (1= zero gradient, 2=constant gradient, 3=set to zero)

  integer i, j
  real(kind=dp) :: lifrac, gradient
  logical :: found

  !> \todo Given the exponential decrease in pressure, this routine should probably be converted to using a different method of interpolation (rather tha linear).

  !> \section interpolate_to_grid_centers_alg Algorithm
  !! @{
  do i=1, n_model_levels
    !> - Find level where input_pres < pres_l(i)
    found = .false.
    do j=1, n_input_levels
      if (input_pres(j) <= pres_l(i)) then
        found = .true.
        !> - If there are model levels below the input levels, set them according to the method specified
        if (j == 1) then
          select case (below_input_level_method)
            case (1)
              !zero gradient method
              model_field(i) = input_field(1)
            case (2)
              !constant gradient method
              gradient = (input_field(2) - input_field(1))/(input_pres(2) - input_pres(1))
              model_field(i) = input_field(1) - gradient*(input_pres(1) - pres_l(i))
            case (3)
              !zero method
              model_field(i) = 0.0
          end select
        else
          !> - Once found, calculate the linear interpolation fraction.
          lifrac = (pres_l(i) - input_pres(j-1))/(input_pres(j)-input_pres(j-1))

          !> - Calculate the model field at level i using linear interpolation.
          model_field(i) = input_field(j-1) + lifrac*(input_field(j)-input_field(j-1))
        end if
        last_initialized = i
        exit !the inner loop, since model level i is calculated
      end if
    end do

    !> - If no level is found, the input sounding does not reach to top of model domain; return the index of the last level initialized.
    if (.not. found) then
      last_initialized = i-1
      exit !the outer (model layers) loop, since there are no layers above this one that can be calculated
    end if
  end do
  !> @}
end subroutine

subroutine w_to_omega(n_col, n_lev, w, p, T, omega)
  integer, intent(in)    :: n_col
  integer, intent(in)    :: n_lev
  real(kind=dp), intent(in) :: w(:,:)
  real(kind=dp), intent(in) :: p(:,:)
  real(kind=dp), intent(in) :: T(:,:)
  real(kind=dp), intent(out) :: omega(:,:)

  real(kind=dp) :: rho(n_col,n_lev)

  rho = p/(con_rd*T)
  omega = -1*w*rho*con_g

end subroutine w_to_omega
!> @}
!> @}
end module scm_utils

module NetCDF_read
  use scm_kinds, only : sp, dp, qp
  use netcdf
  
  implicit none
  
  real(kind=dp) :: missing_value = -9999.0
  
  interface NetCDF_read_var
    module procedure NetCDF_read_var_0d_int
    module procedure NetCDF_read_var_0d
    module procedure NetCDF_read_var_1d
    module procedure NetCDF_read_var_2d
  end interface
  
  contains
  
  subroutine NetCDF_read_var_0d_int(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    integer, intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID))
      call check(NF90_GET_VAR(ncid,varID,var_data))
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data))
      end if
    end if
    
  end subroutine NetCDF_read_var_0d_int
  
  subroutine NetCDF_read_var_0d(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID))
      call check(NF90_GET_VAR(ncid,varID,var_data))
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data))
      end if
    end if
    
  end subroutine NetCDF_read_var_0d
  
  subroutine NetCDF_read_var_1d(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID))
      call check(NF90_GET_VAR(ncid,varID,var_data))
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data))
      end if
    end if
    
  end subroutine NetCDF_read_var_1d
  
  subroutine NetCDF_read_var_2d(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID))
      call check(NF90_GET_VAR(ncid,varID,var_data))
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data))
      end if
    end if
    
  end subroutine NetCDF_read_var_2d
  
  !Generic subroutine to check for netCDF I/O errors
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= NF90_NOERR) then
      print *, trim(nf90_strerror(status))
      stop "stopped"
    end if
  end subroutine check
  
end module NetCDF_read  
  
module data_qc
  use scm_kinds, only : sp, dp, qp
  use NetCDF_read, only: missing_value
  
  implicit none
  
  interface check_missing
    module procedure check_missing_0d
    module procedure check_missing_1d
  end interface
  
  interface conditionally_set_var
    module procedure conditionally_set_var_0d
    module procedure conditionally_set_var_1d
  end interface
  
  contains
    
  subroutine check_missing_0d(var, missing)
    real(kind=dp), intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if (var == missing_value) missing = .true.
  end subroutine check_missing_0d
  
  subroutine check_missing_1d(var, missing)
    real(kind=dp), dimension(:), intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if ( ANY(var == missing_value)) missing = .true.
  end subroutine check_missing_1d
  
  subroutine conditionally_set_var_0d(input, set_var, input_name, req, missing)
    real(kind=dp), intent(in) :: input
    real(kind=dp), intent(inout) :: set_var
    character (*), intent(in) :: input_name
    logical, intent(in) :: req
    logical, intent(out) :: missing
    
    call check_missing(input, missing)
    
    if (.not. missing) then
      set_var = input
    else
      if (req) then
        write(0,'(a,i0,a)') "The variable '" // input_name // "' in the case data file had missing data, but it is required for the given physics configuration. Stopping..."
        STOP
      end if
    end if
    
  end subroutine conditionally_set_var_0d
  
  subroutine conditionally_set_var_1d(input, set_var, input_name, req, missing)
    real(kind=dp), dimension(:), intent(in) :: input
    real(kind=dp), dimension(:), intent(inout) :: set_var
    character (*), intent(in) :: input_name
    logical, intent(in) :: req
    logical, intent(out) :: missing
    
    call check_missing(input, missing)
    
    if (.not. missing) then
      set_var = input
    else
      if (req) then
        write(0,'(a,i0,a)') "The variable '" // input_name // "' in the case data file had missing data, but it is required for the given physics configuration. Stopping..."
        STOP
      end if
    end if
    
  end subroutine conditionally_set_var_1d
end module data_qc
