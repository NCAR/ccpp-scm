!> \file gmtb_scm.f90
!!  Contains the main program for the GMTB SCM.

!> \defgroup SCM GMTB Single Column Model
!! @{

module gmtb_scm_main

contains

  !> \brief Main SCM program
  !! \author Grant J. Firl
  !! \date March-June 2016
  !!
  !! The prototype GMTB SCM can be thought of as an update to the existing GFS SCM to 1) use the NGGPS Interoperable Physics Driver (IPD)
  !! (and to keep up with its development) and 2) have the capability to run GCSS-style idealized cases. Its initialization draws heavily
  !! from the standalone IPD developed by Patrick Tripp at NOAA EMC. Many of the algorithms are heavily influenced by the original GFS SCM,
  !! and clues for implementation were also gleaned from how the IPD is implemented (called) in the GSM. The standalone driver demonstrated how to
  !! use the nuopc_physics module to calculate tendencies due to radiation and the rest of the physics suite once, including
  !! initialization. It uses 8 columns and it initializes GFS state variables from unformatted binary output files that are (I believe)
  !! generated from a GFS model run.
  !!
  !! The SCM prototype differs from the standalone driver in several important ways. First, although the state arrays will maintain
  !! a "horizontal" dimension (be capable of representing more than one column), the code is designed to use a single column. Second,
  !! in order to avoid becoming a GFS-only SCM, input is achieved through host-model agnostic netCDF files. In addition, variables needed
  !! to initialize or configure the physics suite are obtained through external namelist files to avoid using GFS output files. Third,
  !! the model will have the ability to use different vertical grids and coordinates as well as different time-stepping schemes (initial
  !! capability follows the GFS, however). Fourth, the SCM includes functionality to advance through time, calling physics repeatedly,
  !! and to replace dynamics/advection with specified large-scale forcing.
subroutine gmtb_scm_main_sub()
  !> \section alg Main Algorithm
  !! @{

  !> \subsection modules Load the necessary modules for the GTMB SCM.
  !! @{

  !> - Load the GMTB SCM modules called from the main program.
  use gmtb_scm_kinds, only: sp, dp, qp
  use gmtb_scm_input
  use gmtb_scm_utils
  use gmtb_scm_vgrid
  use gmtb_scm_setup
  use gmtb_scm_forcing
  use gmtb_scm_time_integration
  use gmtb_scm_output

  !> - Load modules that are needed for using the interoperable physics driver (IPD); includes derived data types (DDTs) and subroutines that
  !! are called to drive the physics suite.
  use nuopc_physics, only: state_fields_in, state_fields_out, sfc_properties, diagnostics, interface_fields, cloud_properties, &
        radiation_tendencies, model_parameters, nuopc_phys_init, nuopc_rad_update, nuopc_rad_run, nuopc_phys_run, tbd_ddt, &
        use_nuopc, dynamic_parameters

  !> - Load the modules that are needed for the GFS physics suite
  use module_radsw_parameters,    only : topfsw_type, sfcfsw_type
  use module_radlw_parameters,    only : topflw_type, sfcflw_type
  use funcphys
  use namelist_soilveg
  use physcons,                   only : dxmax, dxmin, dxinv, con_pi, max_lat, max_lon, min_lat, min_lon
  use physparam,                  only : ipsd0
  use machine
  use mersenne_twister,           only: random_setseed, random_index, random_stat
  use module_ras,                 only : ras_init
  use ozne_def,                   only : latsozp, levozp, timeoz, pl_coeff, pl_lat, pl_time, pl_pres
  use h2o_def,                    only : latsh2o, levh2o, timeh2o,  h2o_coeff, h2o_lat, h2o_pres, h2o_time
  use module_radiation_driver,    only : radupdate

  !> @}

  implicit none

  !> \subsection var_dec Variable Declarations
  !! @{

  !> - Define the variables necessary for the SCM infrastructure.

  character(len=80)                 :: experiment_name !> name of model configuration file
  character(len=80)                 :: model_name !< name of "host" model (must be "GFS" for prototype)
  character(len=80)                 :: physics_suite !< name of physics suite (must be "GFS_operational" for prototype)
  character(len=80)                 :: output_dir !< name of output directory to place netCDF file
  character(len=80)                 :: output_file !< name of output file (without the file extension)
  character(len=80)                 :: case_name !< name of case initialization and forcing to use (different than experiment name, which names the model run (as a control, experiment_1, etc.))

  integer                           :: i, k, ioerror, allocate_status, grid_error !< dummy indices and error statuses
  integer                           :: n_levels !< number of model levels (must be 64 for prototype)
  integer                           :: itt !< current model iteration
  integer                           :: itt_out  !< output iteration counter
  integer                           :: time_scheme !< 1=> forward Euler, 2=> filtered leapfrog
  integer                           :: n_timesteps !< number of timesteps needed to integrate over runtime
  integer                           :: n_time_levels !< number of time levels to keep track of for time-integration scheme (2 for leapfrog)
  integer                           :: n_itt_swrad !< number of iterations between calls to SW rad
  integer                           :: n_itt_lwrad !< number of iterations between calls to LW rad
  integer                           :: n_itt_out !< number of iterations between calls to write the output
  integer, parameter                :: n_levels_smooth = 5 !< the number of levels over which the input profiles are smoothed into the reference profiles
  integer, parameter                :: ngptc = 1 !< number of columns

  logical                           :: use_leapfrog !< flag for using the leapfrog as the time-integration scheme for the forcing
  logical                           :: use_forward !< flag for using the forward Euler scheme for time-stepping the forcing
  logical                           :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
  integer                           :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer                           :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer                           :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere

  real(kind=dp)                           :: model_time !< elapsed model time (s)
  real(kind=dp)                           :: model_time_m1 !< elapsed model time at the end of the previous timestep (s)
  real(kind=dp)                           :: dt !< physics time step (s)
  real(kind=dp)                           :: runtime !< total runtime (s)
  real(kind=dp)                           :: output_frequency !< how often output is written (s)
  real(kind=dp)                           :: swrad_frequency !< how often SW radiation is called
  real(kind=dp)                           :: lwrad_frequency !< how often LW radiation is called
  real(kind=dp)                           :: relax_time !< time scale for hor. wind nudging (s)
  real(kind=dp)                           :: deg_to_rad_const
  real(kind=dp), parameter                :: c_filter = 0.15 !< parameter that controls the amount of damping in the leapfrog filter
  real(kind=dp), parameter                :: rhoh2o = 1000.0


  !> - Define the case-specific initialization and forcing variables.
  integer                           :: input_nlev !< number of levels in the input file
  integer                           :: input_ntimes !< number of times in the input file where forcing is available
  real(kind=dp), allocatable              :: input_pres(:) !< pressure (Pa) of input levels
  real(kind=dp), allocatable              :: input_time(:) !< time (s) since beginning of forcing file
  real(kind=dp), allocatable              :: input_height(:) !< height of input levels (m) (initial)
  real(kind=dp), allocatable              :: input_thetail(:) !< ice-liquid water potential temperature(K) (initial)
  real(kind=dp), allocatable              :: input_qt(:) !< total water specific humidity (kg/kg) (initial)
  real(kind=dp), allocatable              :: input_ql(:) !< suspended liquid specific humidity (kg/kg) (initial)
  real(kind=dp), allocatable              :: input_qi(:) !< suspended ice specific humidity (kg/kg) (initial)
  real(kind=dp), allocatable              :: input_u(:) !< zonal wind (m/s) (initial)
  real(kind=dp), allocatable              :: input_v(:) !< meridional wind (m/s) (initial)
  real(kind=dp), allocatable              :: input_tke(:) !< turbulence kinetic energy (m^2/s^2) (initial)
  real(kind=dp), allocatable              :: input_ozone(:) !< ozone mass mixing ratio (kg/kg) (initial)
  real(kind=dp), allocatable              :: input_lat(:) !< time-series of latitude of column center
  real(kind=dp), allocatable              :: input_lon(:) !< time-series of longitude of column center
  real(kind=dp), allocatable              :: input_pres_surf(:) !< time-series of surface pressure (Pa)
  real(kind=dp), allocatable              :: input_T_surf(:) !< time-series of surface temperture
  real(kind=dp), allocatable              :: input_w_ls(:,:) !< large-scale vertical velocity (m/s) (time, levels)
  real(kind=dp), allocatable              :: input_omega(:,:) !< large-scale pressure vertical velocity (Pa/s) (time, levels)
  real(kind=dp), allocatable              :: input_u_g(:,:) !< geostrophic zonal wind (m/s) (time, levels)
  real(kind=dp), allocatable              :: input_v_g(:,:) !< geostrophic meridional wind (m/s) (time, levels)
  real(kind=dp), allocatable              :: input_u_nudge(:,:) !< E-W wind profile to nudge towards (m/s) (time, levels)
  real(kind=dp), allocatable              :: input_v_nudge(:,:) !< N-S wind profile to nudge towards (m/s) (time, levels)
  real(kind=dp), allocatable              :: input_T_nudge(:,:) !< absolute temperature profile to nudge towards (K) (time, levels)
  real(kind=dp), allocatable              :: input_thil_nudge(:,:) !< liquid potential temperature profile to nudge towards (K) (time, levels)
  real(kind=dp), allocatable              :: input_qt_nudge(:,:) !< specific humidity profile to nudge towards (kg/kg) (time, levels)
  real(kind=dp), allocatable              :: input_dT_dt_rad(:,:) !< large-scale T-tendency (K/s) (time, levels)
  real(kind=dp), allocatable              :: input_h_advec_thetail(:,:) !< large-scale tendency of ice-liquid potential temperature due to horizontal advection (K/s) (time, levels)
  real(kind=dp), allocatable              :: input_h_advec_qt(:,:) !< large-scale tendency of total water specific humidity due to horizontal advection (kg/kg /s) (time, levels)
  real(kind=dp), allocatable              :: input_v_advec_thetail(:,:) !< large-scale tendency of ice-liquid potential temperature due to vertical advection (K/s) (time, levels)
  real(kind=dp), allocatable              :: input_v_advec_qt(:,:) !< large-scale tendency of total water specific humidity due to vertical advection (kg/kg /s) (time, levels)
  real(kind=dp), allocatable              :: input_sh_flux_sfc(:) !< time-series of surface sensible heat flux (K m s^-1)
  real(kind=dp), allocatable              :: input_lh_flux_sfc(:) !< time-series of surface latent heat flux (kg kg^-1 m s^-1)

  !> - Define the reference profile variables.
  integer                           :: ref_nlev !< number of levels in the reference profile
  real(kind=dp), allocatable              :: ref_pres(:) !< pressure (Pa) of the reference profile levels
  real(kind=dp), allocatable              :: ref_T(:) !< absolute T (K) of the reference profile levels
  real(kind=dp), allocatable              :: ref_qv(:) !< water vapor specific humidity (kg/kg) of the reference profile levels
  real(kind=dp), allocatable              :: ref_ozone(:) !< ozone mass mixing ratio (kg/kg) of the reference profile levels

  !> - Define the SCM state variables; variables with appended "i" are interface; variables with appended "l" are layer-centered.
  !!  - index order for grid is (horizontal, vertical);
  !!  - index order for state variables is (horizontal, vertical, timesteps);
  !!  - index order for tracer is (horizontal, vertical, tracer_index, timesteps)
  real(kind=dp), allocatable              :: pres_i(:,:), pres_l(:,:) !< pressure on grid interfaces, centers (Pa)
  real(kind=dp), allocatable              :: si(:,:), sl(:,:) !< sigma on grid interfaces, centers
  real(kind=dp), allocatable              :: exner_i(:,:), exner_l(:,:) !< exner function on grid interfaces, centers
  real(kind=dp), allocatable              :: state_T(:,:,:) !< model state absolute temperature at grid centers (K)
  real(kind=dp), allocatable              :: state_u(:,:,:), state_v(:,:,:) !< model state horizontal winds at grid centers (m/s)
  real(kind=dp), allocatable              :: state_tracer(:,:,:,:) !< model state tracer at grid centers
  real(kind=dp), allocatable              :: temp_T(:,:,:), temp_u(:,:,:), temp_v(:,:,:), temp_tracer(:,:,:,:) !< used for time-filtering

  !> - Define forcing-related variables (indexing is (horizontal, vertical)).
  real(kind=dp), allocatable              :: u_force_tend(:,:), v_force_tend(:,:), T_force_tend(:,:), qv_force_tend(:,:) !< total u, v, T, q forcing (units/s) (horizontal, vertical)
  real(kind=dp), allocatable              :: w_ls(:,:), omega(:,:), u_g(:,:), v_g(:,:), dT_dt_rad(:,:), h_advec_thil(:,:), &
    h_advec_qt(:,:), v_advec_thil(:,:), v_advec_qt(:,:), u_nudge(:,:), v_nudge(:,:), T_nudge(:,:), thil_nudge(:,:), qt_nudge(:,:) !< forcing terms interpolated to the model time and grid
  real(kind=dp), allocatable              :: T_surf(:), pres_surf(:) !< surface temperature and pressure interpolated to the model time
  real(kind=dp), allocatable              :: sh_flux(:), lh_flux(:) !< surface sensible and latent heat fluxes interpolated to the model time

  !> - Define ozone forcing variables (the input forcing terms reside on their own (latitude, levels, time) grid and have 'pl_coeff' terms)
  real(kind=kind_phys), allocatable :: ozone_forcing_in(:,:,:,:) !<array of ozone forcing to be read from file generated by chgres script (in GFS); zeroed out for SCM; (latsozp, levozp, pl_coeff, timeoz); using climatological ozone uses (latsozc, levozc, timeozc)
  real(kind=kind_phys), allocatable :: ozone_forcing_out(:,:,:) !<ozone forcing array in format used by ozphys (after ozinterpol); 'tbd%prdout' in nuopc_physics, 'prdoz' in gbphys, 'prdout' in ozphys
  real(kind=kind_phys), allocatable :: o3out(:,:,:) !<temporary ozone forcing array output by ozinterpol (values are passed to ozone_forcing_out in nuopc_phys_run immediately after ozinterpol)
  real(kind=kind_phys), allocatable :: ozone_pres(:) !< pressure levels for input ozone forcing data
  real(kind=kind_phys), allocatable :: ozone_ddy(:) !< fractional distance between latitude bands for ozone forcing data?
  integer, allocatable              :: ozone_jindx1(:), ozone_jindx2(:) !< latitudnal indices for ozone forcing?

  !> - Define h2o forcing variables (the input forcing terms reside on their own (latitude, levels, time) grid and have 'h2o_coeff' terms)
  real(kind=kind_phys), allocatable :: h2o_forcing_in(:,:,:,:) !<array of h2o forcing to be read from file generated by chgres script (in GFS); zeroed out for SCM; (latsh2o, levh2o, h2o_coeff, timeh2o)
  real(kind=kind_phys), allocatable :: h2o_forcing_out(:,:,:) !<h2o forcing array in format used by h2ophys (after h2ointerpol); 'tbd%prdout' in nuopc_physics, 'prdoz' in gbphys, 'prdout' in ozphys
  real(kind=kind_phys), allocatable :: h2oout(:,:,:) !<temporaryh2o forcing array output by h2ointerpol (values are passed to h2o_forcing_out in nuopc_phys_run immediately after h2ointerpol)
  real(kind=kind_phys), allocatable :: h2o_ddy(:) !< fractional distance between latitude bands for h2o forcing data?
  integer, allocatable              :: h2o_jindx1(:), h2o_jindx2(:) !< latitudnal indices for h2o forcing?



  !> - Define GFS-related grid coefficients.
  real(kind=dp), allocatable              :: a_k(:), b_k(:) !< used to determine grid sigma and pressure levels

  !> - Define variables from the standalone driver.
  !>  - Define namelist variables (see namelists in model_config).
  integer :: ntcw
  integer :: ntiw
  integer :: ntlnc
  integer :: ntinc
  integer :: ncld
  integer :: ntoz
  integer :: NTRAC
  integer :: me
  integer :: lsoil
  integer :: lsm
  integer :: nmtvr
  integer :: nrcm
!  integer :: levozp !defined in ozne_def; currently read in as namelist variable
!  integer :: levh2o !defined in h2o_def
  integer :: lonr
  integer :: latr
  integer :: lats_node_r
  integer :: jcap
  integer :: num_p3d
  integer :: num_p2d
  integer :: npdf3d
  integer :: ncnvcld3d
! integer :: pl_coeff !defined in ozne_def; currently read in as namelist variable
!  integer :: h2o_coeff !defined in h2o_def
  integer :: ncw(2)
  integer :: nlunit      ! local namelist unit for set_soilveg

  real (kind=kind_phys) :: crtrh(3)
  real (kind=kind_phys) :: cdmbgwd(2)
  real (kind=kind_phys) :: ccwf(2)
  real (kind=kind_phys) :: dlqf(2)
  real (kind=kind_phys) :: ctei_rm(2)
  real (kind=kind_phys) :: cgwf(2)
  real (kind=kind_phys) :: prslrd0
  real (kind=kind_phys) :: ral_ts

  integer :: imfshalcnv, imfdeepcnv
  integer :: nctp ! number of cloud types in CS scheme
  integer :: ntke ! tke location in the tracer array
  integer :: ntot3d ! number of total 3d fields for phy_f3d
  integer :: ntot2d ! number of total 2d fields for phy_f2d

  logical :: ras
  logical :: pre_rad
  logical :: ldiag3d
  logical :: lgocart
  logical :: cplflx
  logical :: lssav_cpl
  logical :: flipv
  logical :: old_monin
  logical :: cnvgwd
  logical :: shal_cnv
  logical :: cal_pre
  logical :: mom4ice
  logical :: mstrat
  logical :: trans_trac
  integer :: nstf_name(5)
  logical :: moist_adj
  integer :: thermodyn_id
  integer :: sfcpress_id
  logical :: gen_coord_hybrid
  logical :: lsidea
  logical :: pdfcld
  logical :: shcnvcw
  logical :: redrag
  logical :: hybedmf
  logical :: dspheat

  logical :: shoc_cld           ! flag for SHOC in grrad
  logical :: cscnv              ! flag for Chikira-Sugiyama convection
  logical :: do_shoc            ! flag for SHOC
  logical :: shocaftcnv         ! flag for SHOC
  logical :: h2o_phys           ! flag for stratospheric h2o forcing

  ! Radiation option control parameters
  !real (kind=kind_phys), allocatable :: si_loc(:)
  integer :: ictm, isol, ico2, iaer, ialb, iems, isot, ivegsrc
  integer :: iovr_sw, iovr_lw, isubc_sw, isubc_lw
  integer :: idate(4), iflip
  logical :: sas_shal, crick_proof, ccnorm, norad_precip
  integer :: levr_loc

  !>  - Define DDT variables used in nuopc_physics.F90
  type(model_parameters)     :: mdl_parm
  type(dynamic_parameters)   :: dyn_parm
  type(state_fields_in)      :: state_fldin
  type(sfc_properties)       :: sfc_prop
  type(diagnostics)          :: diags
  type(cloud_properties)     :: cld_prop
  type(radiation_tendencies) :: rad_tend
  type(interface_fields)     :: intrfc_fld
  type(state_fields_out)     :: state_fldout
  type(tbd_ddt)              :: tbddata

  integer :: njeff !< defined as "number of used points" in gbphys; the number of used points can be less than the horizontal dimension
  integer :: lat !< defined as "latitude index used for debug prints" in nuopc_physics.F90

  !most of these variables are defined in the namelists in model_config
  real (kind=kind_phys), dimension(ngptc) :: xlon, xlat, sinlat, coslat
  real (kind=kind_phys)                   :: solhr, solcon, dtlw, dtsw
  integer, dimension(ngptc)               :: icsdsw, icsdlw
  integer, dimension(8)                   :: idat, jdat
  logical                                 :: lsswr     ! logical flags for sw radiation calls
  logical                                 :: lslwr     ! logical flags for lw radiation calls
  logical                                 :: lssav     ! logical flag for store 3-d cloud field
  integer                                 :: ipt       ! index for diagnostic printout point
  logical                                 :: lprnt     ! control flag for diagnostic print out
  logical                                 :: lmfshal              ! mass-flux shallow conv scheme flag
  logical                                 :: lmfdeep2             ! scale-aware mass-flux deep conv scheme flag

  real(kind=kind_phys)                    :: deltim, slag, sdec, cdec
  real(kind=4)                            :: rinc(5) !(DAYS, HOURS, MINUTES, SECONDS, MILLISECONDS)

  !>  - The remaining varialbes are defined in the DDT sections of nuopc_physics.F90
  !  for ddt statein setrad

  !real(kind=kind_phys), allocatable :: tracer(:,:,:)

  ! for ddt sfc_prop
  real (kind=kind_phys), dimension(ngptc) :: slmsk, tsfc, snowd, sncovr, snoalb, zorl,   &
                                             hprim,alvsf, alnsf, alvwf, alnwf, facsf, facwf,  &
                                             fice,tisfc,slope,shdmin,shdmax,tg3,vfrac,vtype,&
                                             stype,uustar,oro,oro_uf,hice,sheleg,canopy,&
                                             ffmm,ffhh,f10m,t2m,q2m
  ! nmtvr=14
  real (kind=kind_phys)                   :: hprim_v(ngptc,14)

  ! for ddt diags setrad
  integer :: NFXR           ! second dimension of input/output array fluxr
  real (kind=kind_phys) :: fluxr_v(ngptc,39) ! 39=number of fluxes
  type (topfsw_type), dimension(ngptc) :: topfsw
  type (sfcfsw_type), dimension(ngptc) :: sfcfsw
  type (topflw_type), dimension(ngptc) :: topflw
  type (sfcflw_type), dimension(ngptc) :: sfcflw
  real (kind=kind_phys), dimension(ngptc,4) :: dswcmp, uswcmp

  ! for ddt cld_prop setrad

  real (kind=kind_phys), dimension(ngptc)               :: cv, cvt, cvb,flgmin_v
  real (kind=kind_phys), dimension(2)                   :: flgmin
  real (kind=kind_phys), allocatable, dimension(:, :)   :: fcice, frain, rrime, deltaq, cnvw,cnvc, cldcov_v, cnvqc_v
  real (kind=kind_phys)                                 :: sup

  ! for ddt rad_tend setrad
  real (kind=kind_phys), allocatable, dimension(:,:)    :: swh, swhc, hlw, hlwc, dtdt,dtdtr
  real (kind=kind_phys), dimension (ngptc)              :: sfalb,  coszen_v, tsflw, sfcemis, coszdg,rqtk
  real (kind=kind_phys), allocatable                    :: hlwd(:,:,:)

  ! for setphys setrad
  real (kind=kind_phys) :: dtp,dtf,clstp
  integer               :: nnp, nlons(ngptc)

  !  for ddt statein setphys
  real(kind=kind_phys)              :: adjtrc(3)
  real(kind=kind_phys), allocatable :: phii(:,:), phil(:,:)

  ! for ddt diags setphys
  real (kind=kind_phys),dimension(ngptc) :: srunoff, evbsa, evcwa, snohfa, transa, sbsnoa, snowca, soilm, tmpmin, tmpmax, dusfc,&
                                            dvsfc, dtsfc, dqsfc, geshem, gflux, dlwsfc, ulwsfc, suntim, runoff, ep, cldwrk, dugwd, &
                                            dvgwd, psmean, bengsh, spfhmin, spfhmax,rain, rainc, u10m, v10m, zlvl, psurf, &
                                            hpbl, pwat, t1, q1, u1, v1, chh, cmm, dlwsfci, ulwsfci, dswsfci, uswsfci, dusfci, &
                                            dvsfci, dtsfci, dqsfci, gfluxi, epi, smcwlt2, smcref2, wet1, sr
  real(kind=kind_phys), allocatable :: forcet(:,:), forceq(:,:)

  real(kind=kind_phys), allocatable :: dt3dt_v(:,:,:), du3dt_v(:,:,:), dv3dt_v(:,:,:),dq3dt_v(:,:,:), dqdt_v(:,:)

  real (kind=kind_phys),dimension(ngptc)  :: sfcdsw,sfcnsw,sfcdlw
  real(kind=kind_phys), dimension(ngptc)  :: nirbmdi, nirdfdi, visbmdi, visdfdi, nirbmui, nirdfui,visbmui, visdfui, &
                                            aoi_du, aoi_dv, aoi_dt, aoi_dq, aoi_dlw, aoi_dsw, aoi_dnirbm,&
                                            aoi_dnirdf, aoi_dvisbm, aoi_dvisdf, aoi_rain, aoi_snow, aoi_nlw, aoi_nsw,&
                                            aoi_nnirbm, aoi_nnirdf, aoi_nvisbm, aoi_nvisdf, aoi_dusfci, aoi_dvsfci, aoi_dtsfci,&
                                            aoi_dqsfci,  aoi_dlwsfci, aoi_dswsfci, aoi_dnirbmi, aoi_dnirdfi, aoi_dvisbmi,&
                                            aoi_dvisdfi, aoi_nlwsfci, aoi_nswsfci, aoi_nnirbmi, aoi_nnirdfi, aoi_nvisbmi,&
                                            aoi_nvisdfi, aoi_t2mi,    aoi_q2mi, aoi_u10mi,   aoi_v10mi,   aoi_tseai, aoi_psurfi,&
                                            nst_xt,  nst_xs,   nst_xu,  nst_xv, nst_xz, nst_zm,  nst_xtts, nst_xzts, nst_d_conv,&
                                            nst_ifd, nst_dt_cool, nst_qrain, slimsk
  real(kind=kind_phys), dimension(ngptc)  :: aoi_slimskin,aoi_ulwsfcin,aoi_dusfcin, aoi_dvsfcin, aoi_dqsfcin, aoi_dtsfcin
  real(kind=kind_phys), allocatable :: s_sppt_wts(:,:), s_qphys(:,:), s_cplrain0(:), s_cplsnow0(:), s_totprcp0(:), s_cnvprcp0(:), &
                                       s_gu0(:,:), s_gv0(:,:), s_gt0(:,:), s_gr0(:,:,:)
  logical :: s_do_sppt, s_do_shum, s_do_skeb, s_do_vc


  !  for ddt tbddata
  real(kind=kind_phys), dimension(ngptc)    :: dpshc,acv,acvb,acvt,tprcp,srflag, nst_tref, nst_z_c, nst_c_0, nst_c_d, nst_w_0, &
                                                nst_w_d
  real(kind=kind_phys), allocatable         :: fscav(:), fswtr(:)
  real(kind=kind_phys), allocatable         :: phy_fctd(:,:)

  real(kind=kind_phys), allocatable         :: rannum_v(:,:)
  real(kind=kind_phys)                      :: bkgd_vdif_m, bkgd_vdif_h, bkgd_vdif_s, psautco(2), prautco(2), evpco, wminco(2)
  real(kind=kind_phys), dimension(ngptc,4)  :: smc_v,stc_v,slc_v

  real(kind=kind_phys), allocatable   :: upd_mfv(:,:), dwn_mfv(:,:), det_mfv(:,:)
  !  ngptc=8, levs=64, num_p3d, npdf3d, nblck,lats_node_r,lonr
  ! num_p3d=4,npdf3d=0,nblck=25,lats_node_r=3,num_p2d=3
  real(kind=kind_phys), allocatable   :: phy_f3dv(:,:,:)
  real(kind=kind_phys)                ::  phy_f2dv(ngptc,3)

  !  ---  variables used for random number generator (thread safe mode)
  type (random_stat) :: stat
  integer :: numrdm(ngptc*2)
  integer :: ipseed1, ipseed2
  integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds

  real(kind=kind_phys) :: work1, work2, rannum(ngptc*2)

  !> @}

  !> \subsection init SCM Initialization
  !! @{
  !! - Define the namelists that get read in to configure the GFS physics suite. There are separate namelists for GFS physics (general), GFS radiation, GFS surface scheme, GFS clouds, and miscellaneous GFS variables (GFS_tbd)
  NAMELIST /GFS_physics/ ntcw, ncld, ntoz, ntrac, me, lsoil, lsm, nmtvr, nrcm, levozp, lonr, latr, jcap, &
    num_p3d, num_p2d, npdf3d, ncnvcld3d, pl_coeff, ncw, nstf_name, thermodyn_id, sfcpress_id, crtrh, cdmbgwd, ccwf, dlqf, ctei_rm, &
    cgwf, prslrd0, ral_ts, ras, pre_rad, ldiag3d, lgocart, cplflx, flipv, old_monin, cnvgwd, shal_cnv, imfshalcnv, imfdeepcnv, &
    cal_pre, mom4ice, mstrat, trans_trac, moist_adj, gen_coord_hybrid, lsidea, pdfcld, shcnvcw, redrag, hybedmf, dspheat, nfxr, &
    dtp, dtf, clstp, nlons, cscnv, nctp, ntke, do_shoc, shocaftcnv, ntot3d, ntot2d, ntiw, ntlnc, ntinc, levh2o, h2o_phys,&
    h2o_coeff, lats_node_r
  NAMELIST /GFS_rad/ ictm, isol, ico2, iaer, ialb, iems, iovr_sw, iovr_lw, isubc_sw, isubc_lw, shoc_cld, sas_shal, crick_proof, &
    ccnorm, norad_precip, idate, iflip, nlunit, icsdsw, icsdlw, isot, ivegsrc
  NAMELIST /GFS_surface/ slmsk, tsfc, snowd, sncovr, snoalb, zorl, hprim, alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc, &
    slope, shdmin, shdmax, tg3, vfrac, vtype, stype, uustar, oro, oro_uf
  NAMELIST /GFS_clouds/ cv, cvb, cvt, flgmin, sup
  NAMELIST /GFS_tbd/ dpshc, bkgd_vdif_m, bkgd_vdif_h, bkgd_vdif_s, psautco, prautco, evpco, wminco

  !> - Call get_config in \ref input to fetch the basic model configuration variables (from file or command line)
  call get_config_nml(experiment_name, model_name, physics_suite, case_name, dt, time_scheme, runtime, output_frequency, &
    swrad_frequency, lwrad_frequency, n_levels, output_dir, output_file, thermo_forcing_type, mom_forcing_type, relax_time, &
    sfc_flux_spec, reference_profile_choice)

  !> - Call get_case_init in \ref input to read and return the case input data.
  call get_case_init(case_name, input_nlev, input_ntimes, input_pres, input_time, input_height, input_thetail, &
    input_qt, input_ql, input_qi, input_u, input_v, input_tke, input_ozone, input_lat, input_lon, input_pres_surf, input_T_surf, &
    input_sh_flux_sfc, input_lh_flux_sfc, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge, &
    input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
    input_v_advec_thetail, input_v_advec_qt)

  !> - Call get_reference_profile in \ref input to return a reference profile to use above the case profile.
  call get_reference_profile(ref_nlev, ref_pres, ref_T, ref_qv, ref_ozone, reference_profile_choice)

  !> - Set up the vertical grid and time-integration scheme based on the host model choice.
  select case(trim(adjustl(model_name)))
    case("GFS")
      !>  - Call get_GFS_grid in \ref vgrid to read in the necessary coefficients and calculate the pressure-related variables on the grid.
      call get_GFS_vgrid(input_pres_surf(1), n_levels, ngptc, pres_i, pres_l, si, sl, exner_l, exner_i, a_k, b_k, grid_error)
      !>  - Exit if an unsupported number of levels is specified or the file with grid coefficients cannot be opened.
      if (grid_error == 1) then
        write(*,*) 'When using the GFS host model, only 28, 42, 60, 64, and 91 levels are currently supported. Exiting...'
        stop
      end if
      if (grid_error == 2) then
        write(*,*) 'The grid coefficient file could not be opened. Exiting...'
        stop
      end if
    case default
      write(*,*) 'Only the GFS model is currently supported. Exiting...'
      stop
  end select

  !>  - Set the number of time levels needed to be stored in memory based on the choice of time-stepping scheme for the forcing
  select case(time_scheme)
    ! case(1)
    !   n_time_levels = 1
    ! case(2)
    !   n_time_levels = 2
    case default
      n_time_levels = 2
  end select


  !> - Once all dimensions are known, allocate the necessary arrays.

  allocate(forcet(ngptc, n_levels), forceq(ngptc, n_levels))
  forcet = 0.0
  forceq = 0.0

  allocate(state_T(ngptc, n_levels, n_time_levels), &
    state_u(ngptc, n_levels, n_time_levels), state_v(ngptc, n_levels, n_time_levels), &
    state_tracer(ngptc, n_levels, 3, n_time_levels))
  !> - \todo state_tracer and temp_tracer should use a variable for dimension 3.

  allocate(temp_tracer(ngptc, n_levels, 3, n_time_levels), temp_T(ngptc, n_levels, n_time_levels), &
    temp_u(ngptc, n_levels, n_time_levels), temp_v(ngptc, n_levels, n_time_levels))

  allocate(w_ls(ngptc, n_levels), omega(ngptc, n_levels), u_g(ngptc, n_levels), v_g(ngptc, n_levels), dT_dt_rad(ngptc, n_levels), &
    h_advec_thil(ngptc, n_levels), h_advec_qt(ngptc, n_levels), v_advec_thil(ngptc, n_levels), v_advec_qt(ngptc, n_levels), &
    pres_surf(ngptc), T_surf(ngptc), u_nudge(ngptc, n_levels), v_nudge(ngptc, n_levels), T_nudge(ngptc, n_levels), &
    thil_nudge(ngptc, n_levels), qt_nudge(ngptc, n_levels), sh_flux(ngptc), lh_flux(ngptc))

  allocate(fcice(ngptc, n_levels), frain(ngptc, n_levels), rrime(ngptc, n_levels), deltaq(ngptc, n_levels), cnvw(ngptc, n_levels), &
    cnvc(ngptc, n_levels), cldcov_v(ngptc, n_levels), cnvqc_v(ngptc, n_levels))
  fcice = 0.0
  frain = 0.0
  rrime = 1.0
  deltaq = 0.0
  cnvw = 0.0
  cnvc = 0.0
  cldcov_v = 0.0
  cnvqc_v = 0.0

  allocate(swh(ngptc, n_levels), swhc(ngptc, n_levels), hlw(ngptc, n_levels), hlwc(ngptc, n_levels), dtdt(ngptc, n_levels), &
    dtdtr(ngptc, n_levels), hlwd(ngptc,n_levels,6))
  swh = 0.0
  swhc = 0.0
  hlw = 0.0
  hlwc = 0.0
  dtdt = 0.0
  dtdtr = 0.0
  hlwd = 0.0

  allocate(phii(ngptc, n_levels+1), phil(ngptc, n_levels))
  phii = 0.0
  phil = 0.0

  allocate(dt3dt_v(ngptc,n_levels,6), du3dt_v(ngptc,n_levels,4),dv3dt_v(ngptc,n_levels,4),dq3dt_v(ngptc,n_levels,9),&
    dqdt_v(ngptc, n_levels))
  dt3dt_v = 0.0
  du3dt_v = 0.0
  dv3dt_v = 0.0
  dq3dt_v = 0.0
  dqdt_v = 0.0



  allocate(u_force_tend(ngptc,n_levels), v_force_tend(ngptc,n_levels), T_force_tend(ngptc,n_levels), qv_force_tend(ngptc,n_levels))
  u_force_tend = 0.0
  v_force_tend = 0.0
  T_force_tend = 0.0
  qv_force_tend = 0.0

  allocate(s_sppt_wts(ngptc,n_levels), s_qphys(ngptc,n_levels), s_cplrain0(ngptc), s_cplsnow0(ngptc), s_totprcp0(ngptc), &
    s_cnvprcp0(ngptc), s_gu0(ngptc,n_levels), s_gv0(ngptc,n_levels), s_gt0(ngptc,n_levels), s_gr0(ngptc,n_levels,3))
  s_do_sppt = .false.
  s_do_shum = .false.
  s_do_skeb = .false.
  s_do_vc = .false.

  !state_tracer(:,:,1,1) = state_qv(:,:,1)
  !state_tracer(:,:,3,1) = 0.0 !cloud water

  !> - Read in namelists for GFS model control variables rather than read a GFS output file as in the standalone driver.
  !! - \todo These namelists are GFS physics-specific. This should probably be placed in the host model switch above.
  write(*,*) '../model_config/'//trim(experiment_name)//'.nml'
  open(unit=1, file='../model_config/'//trim(experiment_name)//'.nml', status='old', action='read', iostat=ioerror)
  if(ioerror /= 0) then
    write(*,*) 'There was an error opening the file GFS_standalone in the model_config directory. &
      Error code = ',ioerror
    stop
  endif
  read(1, NML=GFS_physics)
  read(1, NML=GFS_rad)
  read(1, NML=GFS_surface)
  read(1, NML=GFS_clouds)
  read(1, NML=GFS_tbd)
  close(1)

  allocate(upd_mfv(ngptc,n_levels), dwn_mfv(ngptc,n_levels), det_mfv(ngptc,n_levels), phy_f3dv(ngptc,n_levels,ntot3d))
  upd_mfv = 0.0
  dwn_mfv = 0.0
  det_mfv = 0.0
  phy_f3dv = 0.0
  phy_f2dv = 0.0

  !for interactive ozone, we need to allocate and define ozone forcing arrays; for now, ozone forcing is always zero; eventually, this can be read from file created by chgres

  !horizontal extent of ozone forcing is same as number of columns
  latsozp = ngptc

  !there is only one ozone forcing time
  timeoz = 2

  !levozp set to same as model levels
  levozp = n_levels

  !allocate variables needed for ozone forcing
  allocate(ozone_forcing_in(latsozp,levozp,pl_coeff,timeoz), ozone_forcing_out(lats_node_r,levozp, pl_coeff), ozone_pres(levozp), &
    ozone_jindx1(lats_node_r), ozone_jindx2(lats_node_r), ozone_ddy(lats_node_r), o3out(levozp,lats_node_r,pl_coeff))

  !allocate variables that are defined in ozne_def
  allocate (pl_lat(latsozp), pl_pres(levozp),pl_time(timeoz+1))

  !pressure levels of ozone forcing data; for SCM are equal to model levels; set pl_pres in ozne_def so that ozone routines can access

  ozone_pres = pres_l(1,:)
  pl_pres = ozone_pres

  !set ozone forcing terms to 0
  ozone_forcing_in = 0.0

  !set ozone latitudnal indices to 1
  ozone_jindx1 = 1
  ozone_jindx2 = 1

  !set this fractional distance to 0
  ozone_ddy = 0.0

  !set ozone latitude array equal to input latitude for SCM
  pl_lat = input_lat(:)

  !set local time of ozone forcing data to 12
  pl_time = (/12.0, 13.0, 14.0/)

  !for interactive stratospheric h2o, we need to allocate and define h2o forcing arrays; for now, h2o forcing is always zero; eventually, this can be read from file created by chgres

  !horizontal extent of ozone forcing is same as number of columns
  latsh2o = ngptc

  !there is only one ozone forcing time
  timeh2o = 2

  !allocate variables needed for ozone forcing
  allocate(h2o_forcing_in(latsh2o,levh2o,h2o_coeff,timeh2o), h2o_forcing_out(lats_node_r,levh2o, h2o_coeff), h2o_pres(levh2o), &
    h2o_jindx1(lats_node_r), h2o_jindx2(lats_node_r), h2o_ddy(lats_node_r), h2oout(levh2o,lats_node_r,h2o_coeff))

  !allocate variables that are defined in ozne_def
  allocate (h2o_lat(latsh2o), h2o_time(timeh2o+1))

  !pressure levels of ozone forcing data; for SCM are equal to model levels; set pl_pres in ozne_def so that ozone routines can access
  h2o_pres = pres_l(1,:)

  !set ozone forcing terms to 0
  h2o_forcing_in = 0.0

  !set ozone latitudnal indices to 1
  h2o_jindx1 = 1
  h2o_jindx2 = 1

  !set this fractional distance to 0
  h2o_ddy = 0.0

  !set ozone latitude array equal to input latitude for SCM
  h2o_lat = input_lat(:)

  !set local time of ozone forcing data to 12
  h2o_time = (/12.0, 13.0, 14.0/)


  !calculate dxmax, dxmin, dxinv
  dxmax = log(1.0/(max_lat*max_lon))
  dxmin = log(1.0/(latr*lonr))
  dxinv = 1.0/(dxmax - dxmin)

  allocate(fscav(ntrac-ncld+2), fswtr(ntrac-ncld+2), phy_fctd(ngptc, nctp))
  fscav = 0.0
  fswtr = 0.0
  phy_fctd = 0.0

  allocate(rannum_v(ngptc, nrcm))
  rannum_v = 0.0

  write(*,*) 'hybedmf=',hybedmf
  !> - \todo The physics and dynamics timesteps are equal in the SCM. This may need to be changed in the future.
  dtp = dt
  dtf = dt

  !> - Since some GFS initialization routines look for datasets in the working directory, use a system call to copy runtime data into the working directory.
  call copy_data_to_working_dir('../standalone_data')


  !> - Call set_state in \ref setup to interpolate the case input data to the model grid and patch in the reference sounding above the case data as necessary.
  call set_state(input_nlev, input_pres, input_qt, input_thetail, input_ql, input_qi, input_u, input_v, input_ozone, &
    n_levels, ngptc, ntoz, ntcw, pres_l, n_levels_smooth, ref_nlev, ref_pres, ref_qv, ref_T, ref_ozone, state_tracer(:,:,:,1), &
    state_T(:,:,1), state_u(:,:,1), state_v(:,:,1))

  !> - Initialize the model elapsed time and model time step iteration counter.
  model_time = 0.0
  itt = 1

  !> - Call interpolate_forcing in \ref forcing to interpolate the input forcing to the model grid and model time.
  call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, input_v_nudge, &
    input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
    input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
    input_sh_flux_sfc, input_lh_flux_sfc, n_levels, ngptc, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, &
    thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, xlat, xlon, pres_surf, T_surf, sh_flux, &
    lh_flux)

  !> - \todo need to convert w_ls to omega and pass it into state_fldin (is pointer)
  !call w_to_omega(ngptc, n_levels, w_ls, pres_l, state_T(:,:,1), omega)

  !> - Put initial date (idate) into idat format for rad_update. (Is this needed? Done in nuopc_phys_init...)
  idat = 0
  idat(1) = idate(4)
  idat(2) = idate(2)
  idat(3) = idate(3)
  idat(5) = idate(1)

  lprnt = .false.
  suntim = 0.0
  tmpmax = 0.0
  spfhmax = 0.0
  tmpmin = 1.0e6
  spfhmin = 1.0e6


  !> - Set jdat (current date) to idat
  jdat = idat

  !> - Set timesteps for longwave and shortwave radiation.
  dtlw = lwrad_frequency
  dtsw = swrad_frequency

  !> - Set the lmfshal and lmfdeep2 flags as done in gloopr (with comment: "set up parameters for Xu & Randall's cloudiness computation")

  lmfshal = shal_cnv .and. (imfshalcnv > 0)
  lmfdeep2 = (imfdeepcnv == 2)

  !> - Initialize the model time step.
  deltim = dt

  !> - Initialize the longwave/shortwave radiation flags to True.
  lsswr = .true.
  lslwr = .true.
  lssav = .true.

  !> - Convert input latitude/longitude to radians and calculate sin and cos of latitude.
  deg_to_rad_const = con_pi/180.0
  xlon = input_lon(1)*deg_to_rad_const
  xlat = input_lat(1)*deg_to_rad_const
  sinlat = sin(xlat)
  coslat = cos(xlat)
  !> - \todo Check that solhr is initialized correctly.
  solhr = idate(1)
  lat = ngptc

  !> - Initialize physics substep to 1, since physics time step = dynamics time step and there are no substeps.
  nnp = 1
  ipt = ngptc

  !> - Initialize a random number array needed for RAS or old SAS
  if (imfdeepcnv <= 0 .or. cal_pre) then
    ipseed1 = idate(1) + idate(2) + idate(3) + idate(4)
    call random_setseed(ipseed1)
    call random_number(rannum)
    do k=1, nrcm
      do i=1, ngptc
        rannum_v(i,k) = rannum(i + ngptc*(k-1))
      end do
    end do
    ipseed1 =  ipseed1 + nint(rannum(1)*1000)
  end if

  !> - Call ras_init if RAS is the chosen convection scheme. \todo Not tested yet.
  if (ras) call ras_init(n_levels, me)


  !> - \todo hprim_v is an array of surface-related variables that are geographic location dependent that are used in the gravity wave drag parameterization.
  !! in Fletcher's SCM, these are read in; presumably generated with outside scripts/code based on an orography dataset and a given location. These have
  !! no influence on maritime cases, but will need to be addressed for land cases.
  hprim_v(1,:) = (/ 105.25647735595703, 1.8643060922622681, -0.34177213907241821, 5.2631579339504242E-002, -0.12643678486347198, &
    -0.18934911489486694, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 1.0000000000000000, 44.882732391357422, &
    0.29146522283554077, 2.4096453562378883E-003, 191.56666564941406 /) ! orographic std dev
  hprim_v = 0.0

  !> - \todo Need to understand what this (adjtrc) does...
  adjtrc = (/ 1.0, 1.0, 1.0 /)

  !> - Call output_init in \ref output to create and initialize the output file.
  call output_init(output_dir, output_file, ngptc, n_levels, idate(4), idate(2), idate(3), idate(1))

  !> - Write out the initial profiles and forcing.
  itt_out = 1
  cldcov_v(:,:) = 0.0
  dqsfci = 0.0
  dtsfci = 0.0
  dusfci = 0.0
  dvsfci = 0.0
  tprcp = 0.0
  rain = 0.0
  rainc = 0.0
  call output_append(output_dir, output_file, itt_out, model_time, pres_l, pres_i, sl, si, phil, phii, state_tracer(:,:,1,1), &
    state_T(:,:,1), state_u(:,:,1), state_v(:,:,1), state_tracer(:,:,3,1), u_force_tend, v_force_tend, T_force_tend, &
    qv_force_tend, w_ls, u_g, v_g, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, T_surf, pres_surf, dqsfci, &
    dtsfci, dusfci, dvsfci, cldcov_v, swh, hlw, rhoh2o*tprcp/dt, rain/dt, rainc/dt, pwat, dt3dt_v(:,:,1)/dt, dt3dt_v(:,:,2)/dt, &
    dt3dt_v(:,:,3)/dt, dt3dt_v(:,:,4)/dt, dt3dt_v(:,:,5)/dt, dt3dt_v(:,:,6)/dt, dq3dt_v(:,:,1)/dt, dq3dt_v(:,:,2)/dt, &
    dq3dt_v(:,:,3)/dt, dq3dt_v(:,:,4)/dt, upd_mfv/dt, dwn_mfv/dt, det_mfv/dt, hpbl, topfsw%upfxc, topfsw%dnfxc, topfsw%upfx0, &
    sfcfsw%upfxc, sfcfsw%dnfxc, sfcfsw%upfx0, sfcfsw%dnfx0, topflw%upfxc, topflw%upfx0, sfcflw%upfxc, sfcflw%upfx0, sfcflw%dnfxc, &
    sfcflw%dnfx0)

  !> \subsubsection standalone The following code is mostly copied from the NUOPC standalone driver.
  !! @{
  !! - The "*_readin" subroutines in the standalone driver are not used in the SCM since the variables are initialized in another way. Similarly, the "*_saveout" routines are not called since the SCM is not being used for regression testing.

  !standalone driver calls a routine to read in these variables from a GFS output file; they have been moved to an editable namelist

  ! call phys_init_readin ( n_levels, ntcw, ncld, ntoz, NTRAC, levs, me, lsoil, lsm, nmtvr, nrcm, levozp,  &
  !                          lonr, latr, jcap, num_p3d, num_p2d, npdf3d, pl_coeff, ncw, nst_fcst, &
  !                          thermodyn_id, sfcpress_id, crtrh, cdmbgwd,  &
  !                          ccwf, dlqf, ctei_rm, cgwf, prslrd0, ras, pre_rad, ldiag3d, lgocart,  &
  !                          cplflx, flipv, old_monin, cnvgwd, shal_cnv, sashal, newsas, cal_pre, mom4ice,  &
  !                          mstrat, trans_trac, moist_adj, &
  !                          gen_coord_hybrid, lsidea, pdfcld, shcnvcw, redrag, hybedmf, dspheat, &
  !                          dxmax, dxmin, dxinv,  &
  !                          ! For radiation
  !                          si_loc, ictm, isol, ico2, iaer, ialb, iems,                    &
  !                          iovr_sw,iovr_lw,isubc_sw,isubc_lw,   &
  !                          sas_shal,crick_proof,ccnorm,norad_precip,idate,iflip,nlunit )
  !
  ! print *, "after phys_init_readin",idate


   !call gfuncphys !> gfuncphys has been moved to nuopc_phys_init in Hang Lei's April 2016 commit

   !for interactive ozone and h2o forcing, you must call ozoneini and h2oini (in gfs_physics_init, this is done before call to nuopc_phys_init); this points variables in mdl_parm DDT to those defined here
   call ozoneini(mdl_parm, ozone_jindx1, ozone_jindx2, ozone_forcing_in, ozone_ddy)
   call h2oini(mdl_parm, h2o_jindx1, h2o_jindx2, h2o_forcing_in, h2o_ddy, h2o_pres)

   !> - Call nuopc_phys_init in nuopc_physics.F90.
   !>  - This routine initializes idat, fills in mdl_parm DDT, calls gfuncphys (funcphys.f), rad_initialize (rad_initialize.f), and set_soilveg (set_soilveg.f)
   call nuopc_phys_init (mdl_parm,ntcw, ncld, ntoz, NTRAC, n_levels, me, lsoil, ntiw, ntlnc, ntinc,            &
                        lsm, nmtvr, nrcm, levozp, levh2o, lonr, latr, jcap, num_p3d, num_p2d,  &
                        npdf3d, ncnvcld3d, pl_coeff, ncw, crtrh, cdmbgwd, ccwf, dlqf, ctei_rm,    &
                        cgwf, prslrd0, ral_ts, ras, pre_rad, ldiag3d, lgocart, cplflx, flipv,  &
                        old_monin,cnvgwd, shal_cnv, imfshalcnv, imfdeepcnv, cal_pre, mom4ice,  &
                        mstrat,trans_trac,nstf_name,moist_adj,thermodyn_id,sfcpress_id, &
                        gen_coord_hybrid, n_levels, lsidea, pdfcld, shcnvcw, redrag,       &
                        hybedmf, dspheat,                                              &
                        dxmax, dxmin, dxinv, h2o_phys, h2o_coeff,  &
                        ! NEW from nems_slg_shoc
                        cscnv, nctp, ntke, do_shoc, shocaftcnv, ntot3d, ntot2d,   &
                        ! For radiation
                        si, ictm, isol, ico2, iaer, ialb, iems, isot, ivegsrc,                   &
                        iovr_sw,iovr_lw,isubc_sw,isubc_lw, shoc_cld,                      &
                        crick_proof,ccnorm,norad_precip,idate,iflip,nlunit, lats_node_r)

   !print statements from standalone driver
   print*, "after nuopc_phys_init"

   print*,'***deltim=',deltim,sinlat, coslat, solhr, NGPTC, njeff
   njeff = ngptc


  !> - Set deltim = 0.5*dt for first forward time step.
  if (time_scheme == 2) then
    deltim = 0.5*dt
  end if

  call radupdate ( idat, jdat, dtsw, deltim, lsswr, me, slag, sdec, cdec, solcon )

  !> - Call the "setrad" methods of the NUOPC DDTs in preparation for nuopc_rad_update and nuopc_rad_run.
  !>  - The dynamic_parameters DDT is the only DDT that has mostly non-pointers. This setrad routine must be called before nuopc_rad_update and nuopc_rad_run if the underlying parameters have changed.
  call dyn_parm%setrad(xlon, xlat,&
                sinlat, coslat, solhr, NGPTC, njeff, itt,&
                jdat, solcon, icsdsw, icsdlw, dtlw, dtsw,&
                lsswr, lslwr, lssav, lmfshal, lmfdeep2, ipt, lprnt,deltim,&
                slag, sdec, cdec )

  !>  - For the first time step, the state_fields_in DDT points to the filtered values.
  call state_fldin%setrad(pres_i,pres_l,exner_l,state_T(:,:,1),state_tracer(:,:,1,1),state_tracer(:,:,2:3,1),omega,ntrac,shoc_cld)

  call sfc_prop%setrad(slmsk, T_surf,snowd,sncovr,snoalb, zorl,   &
    hprim, fice,tisfc,alvsf, alnsf, alvwf, alnwf, facsf, facwf)

  call diags%setrad(NFXR, fluxr_v, topfsw, topflw,&
    dswcmp=dswcmp, uswcmp=uswcmp )

  call cld_prop%setrad(cv, cvt, cvb, fcice, frain, rrime, flgmin_v, &
    cldcov_v, deltaq, sup, cnvw, cnvc)

  call rad_tend%setrad( swh, sfalb, coszen_v, hlw,&
     tsflw, sfcemis, coszdg, hlwc, swhc)

  call intrfc_fld%setrad( sfcfsw, sfcflw )

  !standalone driver calls a subroutine to read in radiation initialization from GFS output file; now, initialization data is
  !read in from the case file and other parameters are read in via a namelist

  ! call rad_run_readin(state_fldin, sfc_prop, diags, intrfc_fld,&
  !                     cld_prop, rad_tend, mdl_parm, dyn_parm)

  !tracer(:,:,ntcw-1) = 0.0
  !> - \todo vvl is initialized to zero, but should be omega converted from large-scale w forcing.
  !omega(:,:) = 0.0

  !> - Call nuopc_rad_update from nuopc_physics.F90. This is a wrapper to the radupdate found in grrad.f used to "update time-sensitive data used by radiation" code
  call nuopc_rad_update (mdl_parm, dyn_parm)

  !> - Call nuopc_rad_run from nuopc_physics.F90. This is a wrapper to the grrad subroutine in grrad.f that is used to "set up and invoke main radiation calls."
  print*, 'after nuopc_rad_update ', deltim, dxmax, dxmin, dxinv
  call nuopc_rad_run (state_fldin, sfc_prop, diags,intrfc_fld,&
          cld_prop, rad_tend, mdl_parm, dyn_parm)

  !no need to call the saveout routine

  !call rad_run_saveout(diags, intrfc_fld, cld_prop, rad_tend)

  !> - \todo Calling dyn_parm%setphys resets slag, sdec, cdec to 0; need to call nuopc_rad_update again? -no way to get these values from dyn_parm (private)
  call radupdate ( idat, jdat, dtsw, deltim, lsswr, me, slag, sdec, cdec, solcon )


  !> - Call the "setphys" methods of the NUOPC DDTs in preparation for the nuopc_phys_run call.
  call dyn_parm%setphys(xlon, xlat,&
          sinlat, coslat, solhr, ngptc, njeff, itt,&
          lssav, lat, dtp, dtf, clstp, nnp, nlons, model_time/3600.0,&
          slag, sdec, cdec )

  !>  - The state_fields_in DDT is filled with the time-filtered values of the state variables.
  call state_fldin%setphys(pres_i, pres_l, exner_l, state_T(:,:,1), state_tracer(:,:,:,1),&
          omega, pres_surf, state_u(:,:,1), state_v(:,:,1), exner_i, phii, phil, adjtrc)


  ! if(time_scheme == 2) then
    !>  - For the leapfrog scheme, the state_fields_out DDT is filled with the unfiltered values of the state variables. The nuopc_phys_run call updates the unfiltered state variables.
    call state_fldout%setphys(state_T(:,:,2), state_tracer(:,:,:,2), state_u(:,:,2), state_v(:,:,2))
  ! else if (time_scheme == 1) then
  !   !>  - For the foward Euler scheme, the state_fields_out DDT is filled with time level 1 values.
  !   call state_fldout%setphys(state_T(:,:,1), state_tracer(:,:,:,1), state_u(:,:,1), state_v(:,:,1))
  ! end if

  call diags%setphys ( srunoff, evbsa, evcwa, snohfa, &
               transa, sbsnoa, snowca, soilm, tmpmin, tmpmax, dusfc,&
               dvsfc, dtsfc, dqsfc, geshem, gflux, dlwsfc, ulwsfc, &
               suntim, runoff, ep, cldwrk, dugwd, dvgwd, psmean, bengsh,&
               spfhmin, spfhmax,rain, rainc, &
               dt3dt_v, dq3dt_v, du3dt_v, dv3dt_v, dqdt_v,&
               u10m, v10m, zlvl, psurf, &
               hpbl, pwat, t1, q1, u1, v1, chh, cmm, dlwsfci, ulwsfci, &
               dswsfci, uswsfci, dusfci, dvsfci, dtsfci, dqsfci, gfluxi, &
               epi, smcwlt2, smcref2, wet1, sr)!, forcet, forceq)


  call intrfc_fld%setphys(sfcdsw, sfcnsw, sfcdlw, nirbmui, nirdfui, visbmui,&
               visdfui, nirbmdi, nirdfdi, visbmdi, visdfdi, aoi_du, aoi_dv,&
               aoi_dt, aoi_dq, aoi_dlw, aoi_dsw, aoi_dnirbm, aoi_dnirdf,&
               aoi_dvisbm, aoi_dvisdf, aoi_rain, aoi_nlw, aoi_nsw, aoi_nnirbm,&
               aoi_nnirdf, aoi_nvisbm, aoi_nvisdf, aoi_slimskin,aoi_ulwsfcin,&
               aoi_dusfcin, aoi_dvsfcin, aoi_dqsfcin, aoi_dtsfcin, aoi_snow,&
               nst_xt, nst_xs, nst_xu,&
               nst_xv, nst_xz, nst_zm, nst_xtts, nst_xzts, nst_d_conv,&
               nst_ifd, nst_dt_cool, nst_Qrain, aoi_dusfci, aoi_dvsfci,&
               aoi_dtsfci, aoi_dqsfci, aoi_dlwsfci, aoi_dswsfci, aoi_dnirbmi,&
               aoi_dnirdfi, aoi_dvisbmi, aoi_dvisdfi, aoi_nlwsfci, aoi_nswsfci,&
               aoi_nnirbmi, aoi_nnirdfi, aoi_nvisbmi, aoi_nvisdfi, aoi_t2mi,&
              aoi_q2mi, aoi_u10mi, aoi_v10mi, aoi_tseai, aoi_psurfi, & !oro,slimsk, &
              s_sppt_wts, s_qphys, s_cplrain0, s_cplsnow0, s_totprcp0, s_cnvprcp0, s_gu0, s_gv0, s_gt0, s_gr0, &
              s_do_sppt, s_do_shum, s_do_skeb, s_do_vc )

  call rad_tend%setphys( swh, sfalb, coszen_v,hlw, tsflw, sfcemis,&
     rqtk=rqtk, hlwd=hlwd, dtdtr=dtdt, swhc=swhc, hlwc=hlwc)

  call sfc_prop%setphys(hprim_v, slope, shdmin, shdmax, snoalb, tg3,&
                slmsk, vfrac, vtype, stype, uustar, oro,&
                oro_uf, hice, fice, tisfc, T_surf, snowd,&
                sheleg,  sncovr, zorl, canopy, ffmm, ffhh,&
                f10m, t2m, q2m)

  call cld_prop%setphys (flgmin,cv, cvt, cvb, cnvqc_v, sup )

  call tbddata%set ( dpshc, o3out, lats_node_r, ozone_forcing_out, h2oout, h2o_forcing_out, ozone_pres, rannum_v,&
               bkgd_vdif_m, bkgd_vdif_h, bkgd_vdif_s,&
               psautco, prautco, evpco, wminco, &
               acv, acvb, acvt, slc_v, smc_v, stc_v,&
               upd_mfv, dwn_mfv, det_mfv, &
               phy_f3dv, phy_f2dv, tprcp, srflag,&
               nst_Tref, nst_z_c, nst_c_0, nst_c_d, nst_w_0, nst_w_d, fscav, fswtr, phy_fctd)

  !standalone driver reads in a GFS output file to initialize variables for the physics call; now, variables are initialized using
  !namelists and case data

   ! call phys_run_readin (state_fldin, sfc_prop,&
   !          diags, intrfc_fld, cld_prop, rad_tend,&
   !          mdl_parm, tbddata, dyn_parm)

  !> - Call nuopc_phys_run from nuopc_physics.F90. This is simply a wrapper to gbphys which is the GFS atmospheric physics driver routine.
  call nuopc_phys_run(state_fldin, state_fldout, sfc_prop,&
            diags, intrfc_fld, cld_prop, rad_tend,&
            mdl_parm, tbddata, dyn_parm )

  ! end standalone driver code
  !> @}
  !! @}

  !> \subsection time_loop Time Loop
  !! @{
  !! \subsubsection time_loop_description Description
  !! @{
  !! The filtered leapfrog scheme can be represented as follows:
  !! \f[
  !! \frac{x^{\tau + 1} - \overline{x^{\tau - 1}}}{2\Delta t}=F^\tau
  !! \f]
  !! where \f$x^{\tau + 1}\f$ is the value of variable \f$x\f$ at time \f$\tau + 1\f$, \f$\overline{x^{\tau - 1}}\f$ is its filtered
  !! value at time \f$\tau - 1\f$, \f$\Delta t\f$ is the time step, and \f$F^\tau\f$ is the collection of processes that change \f$x\f$ at time \f$\tau\f$.
  !! The Robert-Asselin filtered value of \f$x\f$ is given by
  !! \f[
  !! \overline{x^\tau}=(1-c)x^\tau + 0.5c\left(x^{\tau +1} + \overline{x^{\tau - 1}}\right)
  !! \f]
  !! where \f$c\f$ is the filtering constant.
  !! This scheme is implemented in the code as follows:
  !! - For the current time step (going to time \f$\tau + 1\f$), we have available at the beginning of the time step the filtered value \f$\overline{x^{\tau-1}}\f$ and the unfiltered value \f$x^\tau\f$. These are represented by state_x(:,:,1) and state_x(:,:,2) in the code, respectively.
  !! - At the beginning of the time step, these values are saved in a temporary array for use in the filtering step later.
  !! - In the do_time_step routine, nuopc_rad_run uses state_x(:,:,1) (or \f$\overline{x^{\tau -1}}\f$) as input to calculate the heating rate at time \f$\tau\f$.
  !! - In the routine apply_forcing_leapfrog (in \ref forcing), two things happen:
  !!  -# Tendencies of the state variables due to the forcing are calculated for time \f$\tau\f$ (\f$F^\tau\f$ above).
  !!  -# The state variables are stepped forward in time using
  !!  \f[
  !!  x^{\tau+1}=\overline{x^{\tau -1}} + 2\Delta tF^\tau
  !!  \f]
  !!  In the code, state_x(:,:,1) is simply updated. This variable is pointed to by the state_fields_in NUOPC DDT and is used as input for the physics call. At this point in the code, this variable no longer represents the filtered values from the previous time step, but the unfiltered state variables after their change due to the forcing and before their changes due to the physics.
  !! - During the call to nuopc_phys_run, state_x(:,:,1) is used as the starting point for the state variable and its unfiltered value is updated by the physics and stored in state_x(:,:,2) (which state_fields_out points to). This is considered \f$x^{\tau +1}\f$.
  !! - Finally, using \f$x^{\tau +1}\f$ just calculated, and the saved values of \f$\overline{x^{\tau-1}}\f$ and \f$x^\tau\f$, the filtered value \f$\overline{x^\tau}\f$ is calculated by the subroutine filter in \ref time_integration to be used in the next iteration.
  !! @}
  !! \subsubsection time_loop_algorithm Algorithm
  !! @{
  !! - Convert various calling frequencies to time step numbers.
  n_timesteps = ceiling(runtime/dt)
  n_itt_swrad = floor(swrad_frequency/dt)
  n_itt_lwrad = floor(lwrad_frequency/dt)
  n_itt_out = floor(output_frequency/dt)

  !> - Set all time steps equal to the "master" SCM time step
  deltim = dt
  dtp = dt
  dtf = dt

  !> - Update dyn_parm DDT with new deltim.
  call dyn_parm%setrad(xlon, xlat,&
      sinlat, coslat, solhr, NGPTC, njeff, itt,&
      jdat, solcon, icsdsw, icsdlw, dtlw, dtsw,&
      lsswr, lslwr, lssav, lmfshal, lmfdeep2, ipt, lprnt,deltim,&
      slag, sdec, cdec )

  !> - Based on the time-integration method chosen, start the model time loop. Currently, only filtered-leapfrog is implemented.
  !if (use_leapfrog) then

    !> - Start the main time loop.
    do itt = 1, n_timesteps
      !>  - Calculate the elapsed model time.
      model_time = itt*dt
      model_time_m1 = (itt-1)*dt

      !>  - Save previously unfiltered state as temporary for use in the time filter.
      if(time_scheme == 2) then
        temp_tracer = state_tracer
        temp_T = state_T
        temp_u = state_u
        temp_v = state_v
      end if

      dt3dt_v = 0.0
      du3dt_v = 0.0
      dv3dt_v = 0.0
      dq3dt_v = 0.0
      dqdt_v = 0.0

      !>  - Call interpolate_forcing from \ref forcing to interpolate the forcing terms to the model grid and time.
      call interpolate_forcing(input_ntimes, input_nlev, input_w_ls, input_omega, input_u_g, input_v_g, input_u_nudge, &
        input_v_nudge, input_T_nudge, input_thil_nudge, input_qt_nudge, input_dT_dt_rad, input_h_advec_thetail, input_h_advec_qt, &
        input_v_advec_thetail, input_v_advec_qt, input_pres, input_time, input_lat, input_lon, input_pres_surf, input_T_surf, &
        n_levels, ngptc, pres_l, model_time, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, dT_dt_rad, &
        h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, xlat, xlon, pres_surf, T_surf)

      !>  - Update the pressure (and related variables) on the grid levels since some cases have a time-dependent surface pressure.
      call calc_GFS_pres_and_exner(pres_surf, n_levels, ngptc, a_k, b_k, pres_i, pres_l, si, sl, exner_l, exner_i)

      !call w_to_omega(ngptc, n_levels, w_ls, pres_l, state_T(:,:,1), omega)

      !>  - Recalculate sin and cos of latitude since xlon and xlat are time-dependent for some cases. (Column is moving in geographic space)
      sinlat = sin(xlat)
      coslat = cos(xlat)

      !>  - Before the time step is executed, perform various operations that are found in GLOOPR and GLOOPB before the physics time steps are executed in the GSM.
      !>   - GLOOPR code:

      !>    - Determine if SW, LW radiation needs to be called and set logical flags.
      lsswr = mod(itt, n_itt_swrad)==0
      lslwr = mod(itt, n_itt_lwrad)==0

      !>    - Update jdat (yr, mon, day, t-zone, hr, min, sec, mil-sec) given the elapsed model time and the initial time.
      rinc = 0
      rinc(4) = model_time_m1
      !w3movdat is a GFS routine to calculate the current date (jdat) from an elapsed time and an initial date (rinc is single prec.)
      call w3movdat(rinc, idat, jdat)
      solhr = mod(model_time_m1/3600.0 + idate(1),24.0)

      !>    - Set up random numbers for subgrid cloudiness in radiation
      if (lsswr .or. lslwr) then
        if (isubc_lw == 2 .or. isubc_sw ==2) then
          ipseed2 = mod(nint(100.0*sqrt(model_time)), ipsdlim) + 1 + ipsd0
          call random_setseed(ipseed2, stat)
          call random_index(ipsdlim, numrdm, stat)
          do i=1, ngptc
            icsdsw(i) = numrdm(i)
            icsdlw(i) = numrdm(i+ngptc)
          end do
        end if
      end if

      !>    - Grab variables from "3D arrays saved for restart", calculate minimum large ice fraction ?
      !>    - \todo Go back through GLOOPR to determine what needs to be included here.
      if (NUM_P3D == 3) then
        do k = 1, n_levels !was levr
          do i = 1, ngptc
            fcice(i,k) = phy_f3dv(i,k,1)
            frain(i,k) = phy_f3dv(i,k,2)
            rrime(i,k) = phy_f3dv(i,k,3)
          enddo
        enddo

        work1 = (log(coslat(1)/(ngptc)) - dxmin) * dxinv
        work1 = max(0.0, min(1.0,work1))
        work2 = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
        do i=1,ngptc
          flgmin_v(i) = work2
        enddo
      else
        do i=1,ngptc
          flgmin_v(i) = 0.0
        enddo
      endif

      if(num_p3d == 4 .and. npdf3d == 3) then
        do k = 1, n_levels
          do i = 1, ngptc
            deltaq(i,k) = phy_f3dv(i,k,5)
            cnvw(i,k)   = phy_f3dv(i,k,6)
            cnvc(i,k)   = phy_f3dv(i,k,7)
          enddo
        enddo
      else
        do k = 1, n_levels
          do i = 1, ngptc
            deltaq(i,k) = 0.0
            cnvw(i,k)   = 0.0
            cnvc(i,k)   = 0.0
          enddo
        enddo
      endif

      !### end GLOOPR code

      !### start GLOOPB code
      !>   - GLOOPB code:
      !!    - Set up random number array needed for RAS and old SAS (should probably only do this if these schemes are active)
      if (imfdeepcnv <= 0 .or. cal_pre) then
        call random_setseed(ipseed1)
        call random_number(rannum)
        do k=1, nrcm
          do i=1, ngptc
            rannum_v(i,k) = rannum(i + ngptc*(k-1))
          end do
        end do
        ipseed1 =  ipseed1 + nint(rannum(1)*1000)
      end if

      !>    - Set geopotential and exner functions to zero to force gbphys to recalculate? Set maximum depth of shallow convection.
      do i=1,ngptc
        phil(i,n_levels) = 0.0 ! will force calculation of geopotential in gbphys.
        dpshc(i)     = 0.3 * pres_i(i,1)
        !exner_l(i,1) = 0.0 ! forces calculation of (p/p00)**kappa in gbphys
        !exner_i(i,1) = 0.0 ! forces calculation of (p/p00)**kappa in gbphys
      enddo

      !>    - Zero out convective scheme quantities before the call to gbphys.
      !if (lgocart) then
        do k=1,n_levels
          do i=1,ngptc
            dqdt_v(i,k) = 0.
            upd_mfv(i,k) = 0.
            dwn_mfv(i,k) = 0.
            det_mfv(i,k) = 0.
            cnvqc_v(i,k) = 0.
          enddo
        enddo
      !endif

      !>    - Zero out "temperature change due to radiative heating per time step (K)" before call to gbphys.
      do k=1,n_levels
        do i=1,ngptc
          dtdt(i,k) = 0.0
        enddo
      enddo

      !>    - Zero out "mass change due to moisture variation" before call to gbphys.
      do i=1,ngptc
        rqtk(i) = 0.0
      enddo

      call radupdate ( idat, jdat, dtsw, deltim, lsswr, me, slag, sdec, cdec, solcon )
      !>    - Update dynparm physics parameters before call to gbphys.
      call dyn_parm%setphys(xlon, xlat, &
          sinlat, coslat, solhr, ngptc, njeff, itt, &
          lssav, lat, dtp, dtf, clstp, nnp, nlons, model_time/3600.0, &
          slag, sdec, cdec )
      !### end GLOOPB code

      !>  - Update dynparm radiation parameters before call to grrad. (should this only be called on a radiation time step?)
      call dyn_parm%setrad(xlon, xlat,&
          sinlat, coslat, solhr, NGPTC, njeff, itt,&
          jdat, solcon, icsdsw, icsdlw, dtlw, dtsw,&
          lsswr, lslwr, lssav, lmfshal, lmfdeep2, ipt, lprnt,deltim,&
          slag, sdec, cdec )

      !>  - Call do_time_step in \ref time_integration. This routine calls nuopc_rad_update and nuopc_rad_run (if necessary), apply_forcing_leapfrog from \ref forcing, and nuopc_phys_run. Each updates var(:,:,2), or the unfiltered values of state variables.
      call do_time_step(n_levels, ngptc, time_scheme, state_fldin, state_fldout, sfc_prop, diags, intrfc_fld, cld_prop, rad_tend, &
        mdl_parm, tbddata, dyn_parm, state_tracer, state_T, state_u, state_v, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, &
        thil_nudge, qt_nudge, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, u_force_tend, v_force_tend, &
        T_force_tend, qv_force_tend, exner_l, phii, pres_l, xlat, dt, lsswr, lslwr, thermo_forcing_type, mom_forcing_type, &
        relax_time)

      !>  - Call filter in \ref time_integration to perform the Asselin time filtering of the state variables.
      select case(time_scheme)
        case (1)
          !for forward Euler scheme, no filtering is done; simply transfer output state variables from slot 2 to slot 1
          state_T(:,:,1) = state_T(:,:,2)
          state_u(:,:,1) = state_u(:,:,2)
          state_v(:,:,1) = state_v(:,:,2)
          state_tracer(:,:,:,1) = state_tracer(:,:,:,2)
        case (2)
          !for filtered-leapfrog scheme, call the filtering routine to calculate values of the state variables to save in slot 1 using slot 2 vars (updated, unfiltered) output from the physics
          call filter(c_filter, temp_tracer, temp_T, temp_u, temp_v, state_tracer, state_T, state_u, state_v)

          !> \todo tracers besides water vapor do not need to be filtered (is this right?)
          state_tracer(:,:,2,1) = state_tracer(:,:,2,2)
          state_tracer(:,:,3,1) = state_tracer(:,:,3,2)
      end select


      !>  - Execute code that is found in GLOOPR, GLOOPB after the time step is executed.
      !!   - The following code is from GLOOPR after the call to radiation. Since the calls to radiation and the rest of the physics in GSM are separate, this code might need to be executed before the call to nuopc_phys_run (gbphys), so I'm not sure if it is OK to put here...
      !!   - Set some variables from the radiation output variables.
      if (lsswr) then
        do i = 1, ngptc
          sfcdsw(i) = sfcfsw(i)%dnfxc
          sfcnsw(i) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
        enddo
        if (cplflx) then
          do i = 1, ngptc
            nirbmdi(i) = dswcmp(i,1)
            nirdfdi(i) = dswcmp(i,2)
            visbmdi(i) = dswcmp(i,3)
            visdfdi(i) = dswcmp(i,4)
            nirbmui(i) = uswcmp(i,1)
            nirdfui(i) = uswcmp(i,2)
            visbmui(i) = uswcmp(i,3)
            visdfui(i) = uswcmp(i,4)
          enddo
        endif
      endif
      if (lslwr) then
        do i = 1, ngptc
          sfcdlw(i) = sfcflw(i)%dnfxc
        enddo
      endif

      !>  - Write output (snapshot of filtered state variables) on the output frequency by calling output_append from \ref output.
      if(mod(itt, n_itt_out)==0) then
        itt_out = itt_out+1
        write(*,*) "itt = ",itt
        write(*,*) "model time (s) = ",model_time
        write(*,*) "calling output routine..."


        call output_append(output_dir, output_file, itt_out, model_time, pres_l, pres_i, sl, si, phil, phii, state_tracer(:,:,1,1),&
          state_T(:,:,1), state_u(:,:,1), state_v(:,:,1), state_tracer(:,:,3,1), u_force_tend, v_force_tend, T_force_tend, &
          qv_force_tend, w_ls, u_g, v_g, dT_dt_rad, h_advec_thil, h_advec_qt, v_advec_thil, v_advec_qt, T_surf, pres_surf, dqsfci,&
          dtsfci, dusfci, dvsfci, cldcov_v, swh, hlw, rhoh2o*tprcp/dt, rain/dt, rainc/dt, pwat, dt3dt_v(:,:,1)/dt, &
          dt3dt_v(:,:,2)/dt, dt3dt_v(:,:,3)/dt, dt3dt_v(:,:,4)/dt, dt3dt_v(:,:,5)/dt, dt3dt_v(:,:,6)/dt, dq3dt_v(:,:,1)/dt, &
          dq3dt_v(:,:,2)/dt, dq3dt_v(:,:,3)/dt, dq3dt_v(:,:,4)/dt, upd_mfv/dt, dwn_mfv/dt, det_mfv/dt, hpbl, topfsw%upfxc, &
          topfsw%dnfxc, topfsw%upfx0, sfcfsw%upfxc, sfcfsw%dnfxc, sfcfsw%upfx0, sfcfsw%dnfx0, topflw%upfxc, topflw%upfx0, &
          sfcflw%upfxc, sfcflw%upfx0, sfcflw%dnfxc, sfcflw%dnfx0)

      end if

    end do
  !end if

  !> @}
  !! @}

  !> \subsection clean_up Clean Up
  !! @{
  !! - Finally, clean up the data copied to the working directory.
  call remove_data_from_working_dir()
  !> @}
  !! @}
end subroutine gmtb_scm_main_sub
end module gmtb_scm_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref gmtb_scm_main_sub above.
program gmtb_scm
  use gmtb_scm_main
  call gmtb_scm_main_sub()
end program gmtb_scm
 !> @}
