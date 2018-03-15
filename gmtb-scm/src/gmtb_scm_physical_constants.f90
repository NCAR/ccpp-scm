module gmtb_scm_physical_constants

use gmtb_scm_kinds, only: dp

implicit none

public

!> \section arg_table_gmtb_scm_physical_constants
!! | local_name        | standard_name                                               | long_name                                             | units         | rank | type              |    kind   | intent | optional |
!! |-----------------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | con_cp                | specific_heat_of_dry_air_at_constant_pressure          | specific heat of dry air at constant pressure           | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_fvirt             | ratio_of_vapor_to_dry_air_gas_constants_minus_one      | rv/rd - 1 (rv = ideal gas constant for water vapor)     | none          |    0 | real              | kind_phys | none   | F        |
!! | con_g                 | gravitational_acceleration                             | gravitational acceleration                              | m s-2         |    0 | real              | kind_phys | none   | F        |
!! | con_pi                | pi                                                     | ratio of a circle's circumference to its diameter       | radians       |    0 | real              | kind_phys | none   | F        |
!! | con_rd                | gas_constant_dry_air                                   | ideal gas constant for dry air                          | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_rv                | gas_constant_water_vapor                               | ideal gas constant for water vapor                      | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!!
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
