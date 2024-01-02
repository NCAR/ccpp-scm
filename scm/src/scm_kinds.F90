!> \file scm_kinds.f90
!!  Contains definition of kinds used in the CCPP SCM.

module scm_kinds
  use, intrinsic :: iso_fortran_env

!! \section arg_table_scm_kinds
!! \htmlinclude scm_kinds.html
!!
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
end module scm_kinds
