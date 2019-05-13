module gmtb_scm_physical_constants

use gmtb_scm_kinds, only: dp

implicit none

public

!> \section arg_table_gmtb_scm_physical_constants
!! | local_name             | standard_name                                            | long_name                                               | units         | rank | type              |    kind   | intent | optional |
!! |------------------------|----------------------------------------------------------|---------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | con_cliq               | specific_heat_of_liquid_water_at_constant_pressure       | specific heat of liquid water at constant pressure      | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_cp                 | specific_heat_of_dry_air_at_constant_pressure            | specific heat of dry air at constant pressure           | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_cvap               | specific_heat_of_water_vapor_at_constant_pressure        | specific heat of water vapor at constant pressure       | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_eps                | ratio_of_dry_air_to_water_vapor_gas_constants            | rd/rv                                                   | none          |    0 | real              | kind_phys | none   | F        |
!! | con_epsm1              | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one  | (rd/rv) - 1                                             | none          |    0 | real              | kind_phys | none   | F        |
!! | con_fvirt              | ratio_of_vapor_to_dry_air_gas_constants_minus_one        | (rv/rd) - 1 (rv = ideal gas constant for water vapor)   | none          |    0 | real              | kind_phys | none   | F        |
!! | con_g                  | gravitational_acceleration                               | gravitational acceleration                              | m s-2         |    0 | real              | kind_phys | none   | F        |
!! | con_hfus               | latent_heat_of_fusion_of_water_at_0C                     | latent heat of fusion                                   | J kg-1        |    0 | real              | kind_phys | none   | F        |
!! | con_hvap               | latent_heat_of_vaporization_of_water_at_0C               | latent heat of evaporation/sublimation                  | J kg-1        |    0 | real              | kind_phys | none   | F        |
!! | con_pi                 | pi                                                       | ratio of a circle's circumference to its diameter       | radians       |    0 | real              | kind_phys | none   | F        |
!! | con_rd                 | gas_constant_dry_air                                     | ideal gas constant for dry air                          | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_rv                 | gas_constant_water_vapor                                 | ideal gas constant for water vapor                      | J kg-1 K-1    |    0 | real              | kind_phys | none   | F        |
!! | con_t0c                | temperature_at_zero_celsius                              | temperature at 0 degrees Celsius                        | K             |    0 | real              | kind_phys | none   | F        |
!! | con_ttp                | triple_point_temperature_of_water                        | triple point temperature of water                       | K             |    0 | real              | kind_phys | none   | F        |
!! | con_vonKarman          | vonKarman_constant                                       | vonKarman constant                                      | none          |    0 | real              | kind_phys | none   | F        |
!! | cimin                  | minimum_sea_ice_concentration                            | minimum sea ice concentration                           | frac          |    0 | real              | kind_phys | none   | F        |
!! | rlapse                 | air_temperature_lapse_rate_constant                      | environmental air temperature lapse rate constant       | K m-1         |    0 | real              | kind_phys | none   | F        |
!! | con_jcal               | joules_per_calorie_constant                              | joules per calorie constant                             | J cal-1       |    0 | real              | kind_phys | none   | F        |
!! | con_rhw0               | sea_water_reference_density                              | sea water reference density                             | kg m-3        |    0 | real              | kind_phys | none   | F        |
!! | con_sbc                | steffan_boltzmann_constant                               | Steffan-Boltzmann constant                              | W m-2 K-4     |    0 | real              | kind_phys | none   | F        |
!! | con_tice               | freezing_point_temperature_of_seawater                   | freezing point temperature of seawater                  | K             |    0 | real              | kind_phys | none   | F        |
!!
  real(kind=dp),parameter:: con_pi     =3.1415926535897931

  real(kind=dp),parameter:: con_g      =9.80665e+0
  real(kind=dp),parameter:: con_omega  =7.2921e-5

  real(kind=dp),parameter:: con_rd     =2.8705e+2
  real(kind=dp),parameter:: con_rv     =4.6150e+2
  real(kind=dp),parameter:: con_cp     =1.0046e+3
  real(kind=dp),parameter:: con_cliq   =4.1855e+3
  real(kind=dp),parameter:: con_cvap   =1.8460e+3
  real(kind=dp),parameter:: con_hvap   =2.5000e+6
  real(kind=dp),parameter:: con_hfus   =3.3358e+5
  real(kind=dp),parameter:: con_t0c    =2.7315e+2
  real(kind=dp),parameter:: con_ttp    =2.7316e+2


  real(kind=dp),parameter:: con_rocp   =con_rd/con_cp
  real(kind=dp),parameter:: con_fvirt  =con_rv/con_rd - 1
  real(kind=dp),parameter:: con_eps    =con_rd/con_rv
  real(kind=dp),parameter:: con_epsm1  =con_rd/con_rv-1.

  real(kind=dp),parameter:: con_vonKarman = 0.4
  
  real(kind=kind_phys),parameter:: cimin      =0.15
  real(kind=kind_phys), parameter:: rlapse  = 0.65e-2
  real(kind=kind_phys),parameter:: con_jcal   =4.1855E+0
  real(kind=kind_phys),parameter:: con_rhw0   =1022.0
  real(kind=kind_phys),parameter:: con_sbc    =5.670400e-8
  real(kind=kind_phys),parameter:: con_tice   =2.7120e+2

end module gmtb_scm_physical_constants
