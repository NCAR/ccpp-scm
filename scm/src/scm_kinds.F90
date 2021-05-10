!> \file scm_kinds.f90
!!  Contains definition of kinds used in the CCPP SCM.

module scm_kinds

!! \section arg_table_scm_kinds
!! \htmlinclude scm_kinds.html
!!

  integer, parameter :: sp  = selected_real_kind(P= 6,R=37)
  integer, parameter :: dp  = selected_real_kind(P=13,R=300)
  integer, parameter :: qp = selected_real_kind(P=27,R=2400)

end module scm_kinds
