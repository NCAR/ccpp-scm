!> \file scm_kinds.f90
!!  Contains definition of kinds used in the CCPP SCM.

module scm_kinds

!! \section arg_table_scm_kinds
!! \htmlinclude scm_kinds.html
!!
  integer, parameter :: sp = selected_real_kind(P= 6,R=37)
#ifndef SINGLE_PREC
  integer, parameter :: dp = selected_real_kind(P=13,R=300)
  integer, parameter :: qp = selected_real_kind(P=27,R=2400)
#else
  integer, parameter :: dp = sp
  integer, parameter :: qp = sp
#endif

  ! these types exists to allow generic interface compilation in scm_utils.F90
  integer, parameter :: kind_scm_sp = selected_real_kind(P= 6,R=37)
  integer, parameter :: kind_scm_dp = selected_real_kind(P=13,R=300)
  integer, parameter :: kind_scm_qp = selected_real_kind(P=27,R=2400)
end module scm_kinds
