!> \file scm_mod.F90
!!  Contains variable declarations for CCPP SCM.
!!

module scm_mod
!> \section arg_table_scm_mod
!! \htmlinclude scm_mod.html
!!
  use ccpp_config,   only: ty_ccpp_config
  use scm_type_defs, only: physics_type
  implicit none

  type(physics_type),   target :: physics
  type(ty_ccpp_config), target :: ccpp_cfg

end module scm_mod
