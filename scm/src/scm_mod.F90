module scm_mod
!> \section arg_table_scm_mod
!! \htmlinclude scm_mod.html
!!
  use ccpp_config,   only: ty_ccpp_config
  use scm_type_defs, only: physics_type
  implicit none
  
  integer, parameter :: character_length = 80
  integer, parameter :: int_zero = 0
  integer, parameter :: int_one = 1
  integer, parameter :: int_neg_one = -1
  real(kind=dp), parameter :: real_zero = 0.0
  real(kind=dp), parameter :: real_one = 1.0

  character(len = 80) :: clear_char = ''

  type(physics_type),   target :: physics
  type(ty_ccpp_config), target :: ccpp_cfg

end module scm_mod
