module gmtb_scm_physical_constants

use gmtb_scm_kinds, only: dp

implicit none

public
  real(kind=dp),parameter:: con_pi     =3.1415926535897931

  real(kind=dp),parameter:: con_g      =9.80665e+0
  real(kind=dp),parameter:: con_omega  =7.2921e-5

  real(kind=dp),parameter:: con_rd     =2.8705e+2
  real(kind=dp),parameter:: con_rv     =4.6150e+2
  real(kind=dp),parameter:: con_cp     =1.0046e+3
  real(kind=dp),parameter:: con_hvap   =2.5000e+6
  real(kind=dp),parameter:: con_hfus   =3.3358e+5

  real(kind=dp),parameter:: con_rocp   =con_rd/con_cp
  real(kind=dp),parameter:: con_fvirt  =con_rv/con_rd - 1

end module gmtb_scm_physical_constants
