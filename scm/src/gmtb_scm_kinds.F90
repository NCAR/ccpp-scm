!> \file gmtb_scm_kinds.f90
!!  Contains definition of kinds used in the GMTB SCM.

module gmtb_scm_kinds

  integer, parameter :: sp  = selected_real_kind(P= 6,R=37)
	integer, parameter :: dp  = selected_real_kind(P=13,R=300)
	integer, parameter :: qp = selected_real_kind(P=27,R=2400)

end module gmtb_scm_kinds
