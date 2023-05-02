!> \file scm_utils.f90
!!  Contains miscellaneous helper subroutines.

module scm_utils

use scm_kinds, only: sp, dp, qp
use scm_physical_constants, only: con_rd, con_g

implicit none

interface find_vertical_index_pressure
  module procedure find_vertical_index_pressure_sp
  module procedure find_vertical_index_pressure_dp
end interface

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

subroutine find_vertical_index_pressure_sp(p_thresh, pres, k_out)
  real(kind=sp), intent(in) :: p_thresh
  real(kind=sp), intent(in) :: pres(:)
  integer, intent(out) :: k_out
  
  integer :: k
  
  k_out = -999
  do k=1, size(pres)
    if (pres(k) <= p_thresh) then
      k_out = k
      exit
    end if
  end do
  
end subroutine find_vertical_index_pressure_sp

subroutine find_vertical_index_pressure_dp(p_thresh, pres, k_out)
  real(kind=dp), intent(in) :: p_thresh
  real(kind=dp), intent(in) :: pres(:)
  integer, intent(out) :: k_out
  
  integer :: k
  
  k_out = -999
  do k=1, size(pres)
    if (pres(k) <= p_thresh) then
      k_out = k
      exit
    end if
  end do
  
end subroutine find_vertical_index_pressure_dp

subroutine find_vertical_index_height(z_thresh, height, k_out)
  real(kind=sp), intent(in) :: z_thresh
  real(kind=sp), intent(in) :: height(:)
  integer, intent(out) :: k_out
  
  integer k
  
  k_out = -999
  do k=1, size(height)
    if (height(k) >= z_thresh) then
      k_out = k
      exit
    end if
  end do
  
end subroutine find_vertical_index_height

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

integer function lcm(a,b)
    integer, intent(in) :: a,b
        lcm = a*b / gcd(a,b)
end function lcm
 
integer function gcd(a,b)
    integer, intent(in) :: a,b
    integer :: c,d,t
    c = a
    d = b
    do while (d/=0)
      t = d
      d = mod(c,d)
      c = t
    end do
    gcd = abs(c)
end function gcd

!> @}
!> @}
end module scm_utils

module NetCDF_read
  use scm_kinds, only : sp, dp, qp
  use netcdf
  
  implicit none
  
  real(kind=sp) :: missing_value = -9999.0
  integer       :: missing_value_int = -9999
  character (3) :: missing_value_char = "mis"
  
  interface NetCDF_read_var
    module procedure NetCDF_read_var_0d_int
    module procedure NetCDF_read_var_1d_int
    module procedure NetCDF_read_var_2d_int
    module procedure NetCDF_read_var_3d_int
    module procedure NetCDF_read_var_0d_sp
    module procedure NetCDF_read_var_1d_sp
    module procedure NetCDF_read_var_2d_sp
    module procedure NetCDF_read_var_3d_sp
    module procedure NetCDF_read_var_4d_sp
    module procedure NetCDF_read_var_0d_dp
    module procedure NetCDF_read_var_1d_dp
    module procedure NetCDF_read_var_2d_dp
    module procedure NetCDF_read_var_3d_dp
    module procedure NetCDF_read_var_4d_dp
  end interface
  
  interface NetCDF_conditionally_read_var
    module procedure NetCDF_conditionally_read_var_1d_sp
    module procedure NetCDF_conditionally_read_var_2d_sp
    module procedure NetCDF_conditionally_read_var_3d_sp
    module procedure NetCDF_conditionally_read_var_4d_sp
  end interface NetCDF_conditionally_read_var
  
  interface NetCDF_read_att
    module procedure NetCDF_read_att_char
    module procedure NetCDF_read_att_int
    module procedure NetCDF_read_att_char_or_int
    module procedure NetCDF_read_att_sp
  end interface
  
  contains
  
  subroutine NetCDF_read_var_0d_int(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    integer, intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_0d_int
  
  subroutine NetCDF_read_var_1d_int(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    integer, dimension(:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_1d_int
  
  subroutine NetCDF_read_var_2d_int(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    integer, dimension(:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_2d_int
  
  subroutine NetCDF_read_var_3d_int(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    integer, dimension(:,:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_3d_int
  
  subroutine NetCDF_read_var_0d_sp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=sp), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_0d_sp
  
  subroutine NetCDF_read_var_1d_sp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=sp), dimension(:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_1d_sp
  
  subroutine NetCDF_read_var_2d_sp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=sp), dimension(:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_2d_sp
  
  subroutine NetCDF_read_var_3d_sp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=sp), dimension(:,:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_3d_sp
  
  subroutine NetCDF_read_var_4d_sp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=sp), dimension(:,:,:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_4d_sp
  
  subroutine NetCDF_read_var_0d_dp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_0d_dp
  
  subroutine NetCDF_read_var_1d_dp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_1d_dp
  
  subroutine NetCDF_read_var_2d_dp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_2d_dp
  
  subroutine NetCDF_read_var_3d_dp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:,:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_3d_dp
  
  subroutine NetCDF_read_var_4d_dp(ncid, var_name, req, var_data)
    
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    logical, intent(in) :: req
    real(kind=dp), dimension(:,:,:,:), intent(out) :: var_data
    
    integer :: varID, ierr
    
    if (req) then
      call check(NF90_INQ_VARID(ncid,var_name,varID),var_name)
      call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
    else
      ierr = NF90_INQ_VARID(ncid,var_name,varID)
      if (ierr /= NF90_NOERR) then
        var_data = missing_value
      else
        call check(NF90_GET_VAR(ncid,varID,var_data),var_name)
      end if
    end if
    
  end subroutine NetCDF_read_var_4d_dp
  
  !Generic subroutine to check for netCDF I/O errors
  subroutine check(status, var_name)
    integer, intent ( in) :: status
    character(*), intent(in) :: var_name

    if(status /= NF90_NOERR) then
      print *, trim(nf90_strerror(status)),' ',var_name
      stop "stopped"
    end if
  end subroutine check
  
  subroutine NetCDF_read_att_int(ncid, var_id, att_name, req, att_data)
  
    integer, intent(in) :: ncid, var_id
    character (*), intent(in) :: att_name
    logical, intent(in) :: req
    integer, intent(out) :: att_data
    
    integer :: varID, ierr
    
    if (req) then
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        write(*,*) 'There was an error reading the required '//adjustl(trim(att_name))//' attribute. Stopping...'
        stop
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    else
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        att_data = missing_value_int
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    end if
    
  end subroutine NetCDF_read_att_int
  
  subroutine NetCDF_read_att_sp(ncid, var_id, att_name, req, att_data)
  
    integer, intent(in) :: ncid, var_id
    character (*), intent(in) :: att_name
    logical, intent(in) :: req
    real(kind=sp), intent(out) :: att_data
    
    integer :: varID, ierr
    
    if (req) then
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        write(*,*) 'There was an error reading the required '//adjustl(trim(att_name))//' attribute. Stopping...'
        stop
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    else
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        att_data = missing_value
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    end if
    
  end subroutine NetCDF_read_att_sp
  
  subroutine NetCDF_read_att_char(ncid, var_id, att_name, req, att_data)
  
    integer, intent(in) :: ncid, var_id
    character (*), intent(in) :: att_name
    logical, intent(in) :: req
    character (*), intent(out) :: att_data
    
    integer :: varID, ierr
    
    if (req) then
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        write(*,*) 'There was an error reading the required '//adjustl(trim(att_name))//' attribute. Stopping...'
        stop
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    else
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name)
      if (ierr /= NF90_NOERR) then
        att_data = missing_value_char
      else
        call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
      end if
    end if
    
  end subroutine NetCDF_read_att_char
  
  subroutine NetCDF_read_att_char_or_int(ncid, var_id, att_name, req, att_data, char_att_data)
  
    integer, intent(in) :: ncid, var_id
    character (*), intent(in) :: att_name
    logical, intent(in) :: req
    character (*), intent(out) :: char_att_data
    integer, intent(out) :: att_data
    
    integer :: varID, ierr, type
    
    if (req) then
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name, xtype = type)
      if (ierr /= NF90_NOERR) then
        write(*,*) 'There was an error reading the required '//adjustl(trim(att_name))//' attribute. Stopping...'
        stop
      else
        if (type == NF90_CHAR) then
          call check(NF90_GET_ATT(ncid, var_id, att_name, char_att_data),att_name)
          att_data = missing_value_int
        else if (type == NF90_INT) then
          call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
          char_att_data = missing_value_char
        else
          att_data = missing_value_int
          char_att_data = missing_value_char
        end if
      end if
    else
      ierr = NF90_INQUIRE_ATTRIBUTE(ncid, var_id, att_name, xtype = type)
      if (ierr /= NF90_NOERR) then
        char_att_data = missing_value_char
        att_data = missing_value_int
      else
        if (type == NF90_CHAR) then
          call check(NF90_GET_ATT(ncid, var_id, att_name, char_att_data),att_name)
          att_data = missing_value_int
        else if (type == NF90_INT) then
          call check(NF90_GET_ATT(ncid, var_id, att_name, att_data),att_name)
          char_att_data = missing_value_char
        else
          att_data = missing_value_int
          char_att_data = missing_value_char
        end if
      end if
    end if
    
  end subroutine NetCDF_read_att_char_or_int

  subroutine NetCDF_conditionally_read_var_1d_sp(var_ctl, var_att, var_name, filename, ncid, var_data)
    integer, intent(in) :: var_ctl, ncid
    character (*), intent(in) :: var_att, var_name, filename
    real(kind=sp), dimension(:), intent(out) :: var_data
    real(kind=sp) :: missing_value_eps

    missing_value_eps = missing_value + 0.01

    if (var_ctl > 0) then
      call NetCDF_read_var(ncid, var_name, .False., var_data)
      if (maxval(var_data) < missing_value_eps) then
        write(*,*) 'The global attribute '//var_att//' in '//filename//' indicates that the variable '//var_name//' should be present, but it is missing. Stopping ...'
        stop
      end if
    else
      var_data = missing_value
    end if
  end subroutine NetCDF_conditionally_read_var_1d_sp

  subroutine NetCDF_conditionally_read_var_2d_sp(var_ctl, var_att, var_name, filename, ncid, var_data)
    integer, intent(in) :: var_ctl, ncid
    character (*), intent(in) :: var_att, var_name, filename
    real(kind=sp), dimension(:,:), intent(out) :: var_data
    real(kind=sp) :: missing_value_eps

    missing_value_eps = missing_value + 0.01

    if (var_ctl > 0) then
      call NetCDF_read_var(ncid, var_name, .False., var_data)
      if (maxval(var_data) < missing_value_eps) then
        write(*,*) 'The global attribute '//var_att//' in '//filename//' indicates that the variable '//var_name//' should be present, but it is missing. Stopping ...'
        stop
      end if
    else
      var_data = missing_value
    end if
  end subroutine NetCDF_conditionally_read_var_2d_sp  

  subroutine NetCDF_conditionally_read_var_3d_sp(var_ctl, var_att, var_name, filename, ncid, var_data)
    integer, intent(in) :: var_ctl, ncid
    character (*), intent(in) :: var_att, var_name, filename
    real(kind=sp), dimension(:,:,:), intent(out) :: var_data
    real(kind=sp) :: missing_value_eps
    
    missing_value_eps = missing_value + 0.01
    
    if (var_ctl > 0) then
      call NetCDF_read_var(ncid, var_name, .False., var_data)
      if (maxval(var_data) < missing_value_eps) then
        write(*,*) 'The global attribute '//var_att//' in '//filename//' indicates that the variable '//var_name//' should be present, but it is missing. Stopping ...'
        stop
      end if
    else
      var_data = missing_value
    end if
  end subroutine NetCDF_conditionally_read_var_3d_sp
  
  subroutine NetCDF_conditionally_read_var_4d_sp(var_ctl, var_att, var_name, filename, ncid, var_data)
    integer, intent(in) :: var_ctl, ncid
    character (*), intent(in) :: var_att, var_name, filename
    real(kind=sp), dimension(:,:,:,:), intent(out) :: var_data
    real(kind=sp) :: missing_value_eps
    
    missing_value_eps = missing_value + 0.01
    
    if (var_ctl > 0) then
      call NetCDF_read_var(ncid, var_name, .False., var_data)
      if (maxval(var_data) < missing_value_eps) then
        write(*,*) 'The global attribute '//var_att//' in '//filename//' indicates that the variable '//var_name//' should be present, but it is missing. Stopping ...'
        stop
      end if
    else
      var_data = missing_value
    end if
  end subroutine NetCDF_conditionally_read_var_4d_sp
  
end module NetCDF_read  

module NetCDF_def
  use NetCDF_read, only : check
  use netcdf
  
  implicit none
  
  contains
  
  subroutine NetCDF_def_var(ncid, var_name, var_type, desc, unit, varid, dims)
    use NetCDF_read, only: missing_value, missing_value_int
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    integer, intent(in) :: var_type
    integer, optional, dimension(:), intent(in) :: dims
    character (*), intent(in) :: desc
    character (*), intent(in) :: unit
    integer, intent(out) :: varid
    
    if (present(dims)) then
      call CHECK(NF90_DEF_VAR(NCID=ncid,NAME=var_name,XTYPE=var_type,DIMIDS=dims,VARID=varid),var_name)
    else
      call CHECK(NF90_DEF_VAR(NCID=ncid,NAME=var_name,XTYPE=var_type,VARID=varid),var_name)
    end if
    call CHECK(NF90_PUT_ATT(NCID=ncid,VARID=varid,NAME="description",VALUES=desc),var_name)
    call CHECK(NF90_PUT_ATT(NCID=ncid,VARID=varid,NAME="units",VALUES=unit),var_name)
    
    if (var_type == NF90_FLOAT) then
      call CHECK(NF90_PUT_ATT(NCID=ncid,VARID=varid,NAME="_FillValue",VALUES=missing_value),var_name)
    elseif (var_type == NF90_INT) then
      call CHECK(NF90_PUT_ATT(NCID=ncid,VARID=varid,NAME="_FillValue",VALUES=missing_value_int),var_name)
    else
      write(0,'(a,i0,a)') "The variable '" // var_name // "' is defined as a type other than NF90_FLOAT or NF90_INT. Stopping..."
      STOP
    end if
  
  end subroutine NetCDF_def_var
end module NetCDF_def

module NetCDF_put
  use NetCDF_read, only : check
  use netcdf
  use scm_kinds, only : sp, dp, qp
  
  implicit none
  
  interface NetCDF_put_var
    module procedure NetCDF_put_var_int_0d
    module procedure NetCDF_put_var_1d
    module procedure NetCDF_put_var_2d
  end interface
  
  contains
  
  subroutine NetCDF_put_var_int_0d(ncid, var_name, var, known_varid)
    integer, intent(in) :: ncid
    character (*), intent(in) :: var_name
    integer, intent(in) :: var
    integer, intent(in), optional :: known_varid
    
    integer :: var_id
    
    !write(*,*) 'Putting variable: ',var_name
    if (present(known_varid)) then
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=known_varid,VALUES=var),var_name)
    else
      CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME=var_name,VARID=var_id),var_name)
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=var),var_name)
    end if
    
  end subroutine NetCDF_put_var_int_0d
  
  subroutine NetCDF_put_var_1d(ncid, var_name, var, itt, mult_const)
    integer, intent(in) :: ncid, itt
    character (*), intent(in) :: var_name
    real(kind=dp), intent(in), dimension(:) :: var
    real(kind=dp), intent(in), optional :: mult_const
    
    integer :: var_id
    
    !write(*,*) 'Putting variable: ',var_name
    CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME=var_name,VARID=var_id),var_name)
    if (present(mult_const)) then
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=var*mult_const,START=(/1,itt /)),var_name)
    else
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=var,START=(/1,itt /)),var_name)
    end if
    
  end subroutine NetCDF_put_var_1d
  
  subroutine NetCDF_put_var_2d(ncid, var_name, var, itt, mult_const)
    integer, intent(in) :: ncid, itt
    character (*), intent(in) :: var_name
    real(kind=dp), intent(in), dimension(:,:) :: var
    real(kind=dp), intent(in), optional :: mult_const
    
    integer :: var_id
    
    !write(*,*) 'Putting variable: ',var_name
    CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME=var_name,VARID=var_id),var_name)
    if (present(mult_const)) then
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=var*mult_const,START=(/1,1,itt /)),var_name)
    else
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=var,START=(/1,1,itt /)),var_name)
    end if
    
  end subroutine NetCDF_put_var_2d
end module NetCDF_put

module data_qc
  use scm_kinds, only : sp, dp, qp
  use NetCDF_read, only: missing_value, missing_value_int
  
  implicit none
  
  interface check_missing
    module procedure check_missing_int_0d
    module procedure check_missing_int_1d
    module procedure check_missing_real_dp_0d
    module procedure check_missing_real_dp_1d
  end interface
  
  interface conditionally_set_var
    module procedure conditionally_set_int_var_0d
    module procedure conditionally_set_int_var_1d
    module procedure conditionally_set_real_dp_var_0d
    module procedure conditionally_set_real_dp_var_1d
  end interface
  
  contains
  
  subroutine check_missing_int_0d(var, missing)
    integer, intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if (var == missing_value_int) missing = .true.
  end subroutine check_missing_int_0d
  
  subroutine check_missing_int_1d(var, missing)
    integer, dimension(:), intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if ( ANY(var == missing_value_int)) missing = .true.
  end subroutine check_missing_int_1d
  
  subroutine check_missing_real_dp_0d(var, missing)
    real(kind=dp), intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if (var == missing_value) missing = .true.
  end subroutine check_missing_real_dp_0d
  
  subroutine check_missing_real_dp_1d(var, missing)
    real(kind=dp), dimension(:), intent(in) :: var
    logical, intent(out) :: missing
    
    missing = .false.
    if ( ANY(var == missing_value)) missing = .true.
  end subroutine check_missing_real_dp_1d
  
  subroutine conditionally_set_int_var_0d(input, set_var, input_name, req, missing)
    integer, intent(in) :: input
    integer, intent(inout) :: set_var
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
    
  end subroutine conditionally_set_int_var_0d
  
  subroutine conditionally_set_int_var_1d(input, set_var, input_name, req, missing)
    integer, dimension(:), intent(in) :: input
    integer, dimension(:), intent(inout) :: set_var
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
    
  end subroutine conditionally_set_int_var_1d
  
  subroutine conditionally_set_real_dp_var_0d(input, set_var, input_name, req, missing)
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
    
  end subroutine conditionally_set_real_dp_var_0d
  
  subroutine conditionally_set_real_dp_var_1d(input, set_var, input_name, req, missing)
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
    
  end subroutine conditionally_set_real_dp_var_1d
end module data_qc
