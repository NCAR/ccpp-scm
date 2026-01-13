module scm_physical_constants

use scm_kinds, only: dp

implicit none

public

!> \section arg_table_scm_physical_constants
!! \htmlinclude scm_physical_constants.html
!!
  integer      ,parameter:: con_zero   =0
  real(kind=dp),parameter:: con_pi     =3.1415926535897931

  real(kind=dp),parameter:: con_rerth  =6.3712e+6
  real(kind=dp),parameter:: con_g      =9.80665e+0
  real(kind=dp),parameter:: con_omega  =7.2921e-5
  real(kind=dp),parameter:: con_p0     =1.01325e+5

  real(kind=dp),parameter:: con_rgas   =8.314472_dp
  real(kind=dp),parameter:: con_rd     =2.8705e+2
  real(kind=dp),parameter:: con_rv     =4.6150e+2
  real(kind=dp),parameter:: con_cp     =1.0046e+3
  real(kind=dp),parameter:: con_cliq   =4.1855e+3
  real(kind=dp),parameter:: con_csol   =2.1060e+3
  real(kind=dp),parameter:: con_cvap   =1.8460e+3
  real(kind=dp),parameter:: con_hvap   =2.5000e+6
  real(kind=dp),parameter:: con_hfus   =3.3358e+5
  real(kind=dp),parameter:: con_psat   =6.1078e+2_dp                 !< pres at H2O 3pt (\f$Pa\f$)
  real(kind=dp),parameter:: con_t0c    =2.7315e+2
  real(kind=dp),parameter:: con_ttp    =2.7316e+2
  real(kind=dp),parameter:: con_epsq   =1.0E-12_dp
  real(kind=dp),parameter:: con_epsqs  =1.0E-10_dp

  real(kind=dp),parameter:: con_rocp   =con_rd/con_cp
  real(kind=dp),parameter:: con_rog    =con_rd/con_g
  real(kind=dp),parameter:: con_fvirt  =con_rv/con_rd - 1._dp
  real(kind=dp),parameter:: con_eps    =con_rd/con_rv
  real(kind=dp),parameter:: con_epsm1  =con_rd/con_rv - 1._dp
  real(kind=dp),parameter:: con_1ovg   =1._dp/con_g

  real(kind=dp),parameter:: con_avgd   =6.0221415e23_dp
  real(kind=dp),parameter:: con_amd    =28.9644_dp                   !< molecular wght of dry air (\f$g/mol\f$)
  real(kind=dp),parameter:: con_amw    =18.0154_dp
  real(kind=dp),parameter:: con_amo3   =47.9982_dp
  real(kind=dp),parameter:: karman     =0.4_dp

  real(kind=dp),parameter:: cimin      =0.15
  !> minimum rain amount
  real(kind=dp),parameter:: rainmin    =1.e-13_dp
  real(kind=dp),parameter:: rlapse     =0.65e-2
  real(kind=dp),parameter:: con_jcal   =4.1855E+0
  real(kind=dp),parameter:: con_rhw0   =1022.0
  real(kind=dp),parameter:: con_sbc    =5.670400e-8
  real(kind=dp),parameter:: con_tice   =2.7120e+2

  real(kind=dp),parameter:: rhowater   =1000._dp
  real(kind=dp),parameter:: rholakeice = 0.917e3_dp          !< density of ice on lake (kg/m^3)


  real(kind=dp),parameter:: con_c         = 2.99792458e+8_dp !< speed of light (\f$m/s\f$)
  real(kind=dp),parameter:: con_plnk      = 6.6260693e-34_dp !< planck constant (\f$J/s\f$)
  real(kind=dp),parameter:: con_boltz     = 1.3806505e-23_dp !< boltzmann constant (\f$J/K\f$)
  real(kind=dp),parameter:: con_solr_2002 = 1.3660e+3_dp     !< solar constant (\f$W/m^{2}\f$)-Liu(2002)
  real(kind=dp),parameter:: con_solr_2008 = 1.3608e+3_dp     !< solar constant (\f$W/m^{2}\f$)-nasa-sorce Tim(2008)
  real(kind=dp),parameter:: con_thgni     = -38.15_dp        !< temperature the H.G.Nuc. ice starts

  ! for gfdlmp v3
  real(kind=dp), parameter:: con_rhoair_IFS = 1.0  ! reference air density (kg/m^3), ref: IFS
  real(kind=dp), parameter:: con_rhosnow = 100.0   ! density of snow (kg/m^3)
  real(kind=dp), parameter :: con_visd  = 1.717e-5 ! dynamics viscosity of air at 0 deg C and 1000 hPa (Mason, 1971) (kg/m/s)
  real(kind=dp), parameter :: con_visk  = 1.35e-5  ! kinematic viscosity of air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
  real(kind=dp), parameter :: con_vdifu = 2.25e-5  ! diffusivity of water vapor in air at 0 deg C  and 1000 hPa (Mason, 1971) (m^2/s)
  real(kind=dp), parameter :: con_tcond = 2.40e-2  ! thermal conductivity of air at 0 deg C  and 1000 hPa (Mason, 1971) (J/m/s/K)
  real(kind=dp), parameter :: con_cdg   = 3.15121  ! drag coefficient of graupel (Locatelli and Hobbs, 1974)
  real(kind=dp), parameter :: con_cdh   = 0.5      ! drag coefficient of hail (Heymsfield and Wright, 2014)
  real(kind=dp), parameter :: con_rhocw = 1.0e3    ! density of cloud water (kg/m^3)
  real(kind=dp), parameter :: con_rhoci = 9.17e2   ! density of cloud ice (kg/m^3)
  real(kind=dp), parameter :: con_rhocr = 1.0e3    ! density of rain (Lin et al. 1983) (kg/m^3)
  real(kind=dp), parameter :: con_rhocg = 4.0e2    ! density of graupel (Rutledge and Hobbs 1984) (kg/m^3)
  real(kind=dp), parameter :: con_rhoch = 9.17e2   ! density of hail (Lin et al. 1983) (kg/m^3)
  real(kind=dp), parameter :: con_qcmin = 1.0e-15  ! min value for cloud condensates (kg/kg)
  real(kind=dp), parameter :: con_qfmin = 1.0e-8   ! min value for sedimentation (kg/kg)
  real(kind=dp), parameter :: con_one      = 1_dp
  real(kind=dp), parameter :: con_p001     = 0.001_dp
  real(kind=dp), parameter :: con_secinday = 86400._dp

  ! --- constants from physcons.F90 ---
  real(kind=dp),parameter:: decorr_con = 2.50_dp      !< Decorrelation length constant (km) for iovr = 4 or 5 and idcor = 0
  real(kind=dp),parameter:: qamin = 1.e-16_dp !< Minimum aerosol concentration
end module scm_physical_constants
