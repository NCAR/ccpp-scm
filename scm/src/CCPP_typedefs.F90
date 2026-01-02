module CCPP_typedefs

!> \section arg_table_CCPP_typedefs Argument Table
!! \htmlinclude CCPP_typedefs.html
!!

    ! Physics kind defininitions needed for interstitial DDTs
    use machine,  only: kind_grid, kind_dyn, kind_phys

    ! Physics type defininitions needed for interstitial DDTs
    use module_radsw_parameters,  only: profsw_type, cmpfsw_type
    use module_radlw_parameters,  only: proflw_type
    use GFS_typedefs,             only: GFS_control_type

    implicit none

    private

    ! GFS_interstitial_type         !< fields required to replace interstitial code in GFS_{physics,radiation}_driver.F90 in CCPP
    public GFS_interstitial_type

!! \section arg_table_GFS_interstitial_type
!! \htmlinclude GFS_interstitial_type.html
!!
  type GFS_interstitial_type

    real (kind=kind_phys), pointer      :: adjsfculw_land(:)  => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_ice(:)   => null()  !<
    real (kind=kind_phys), pointer      :: adjsfculw_water(:) => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfu(:)       => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfd(:)       => null()  !<
    real (kind=kind_phys), pointer      :: aerodp(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: alb1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: alpha(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: bexp1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd(:)              => null()  !<
    real (kind=kind_phys), pointer      :: cd_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cd_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cd_water(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq(:)             => null()  !<
    real (kind=kind_phys), pointer      :: cdq_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cdq_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cdq_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cf_upi(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: chh_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: chh_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: clcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: cldf(:)            => null()  !<
    real (kind=kind_phys), pointer      :: cldsa(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: cldtaulw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cldtausw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cld1d(:)           => null()  !<
    real (kind=kind_phys), pointer      :: clouds(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: clw(:,:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: clx(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: cmm_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: cmm_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cmm_water(:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_dqldt(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_fice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnv_mfd(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: cnv_ndrop(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_nice(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: cnvc(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_r(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ctei_rml(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cumabs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dd_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: de_lgth(:)         => null()  !<
    real (kind=kind_phys), pointer      :: del(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: del_gz(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: delr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dlength(:)         => null()  !<
    real (kind=kind_phys), pointer      :: dqdt(:,:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dqsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: drain(:)           => null()  !<
    real (kind=kind_phys), pointer      :: dtdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dtzm(:)            => null()  !<
    real (kind=kind_phys), pointer      :: dt_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dudt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dusfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvdftra(:,:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: dvdt(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfcg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dvsfc1(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dzlyr(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: elvmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: ep1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ep1d_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evap_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: evap_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: evap_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: evbs(:)            => null()  !<
    real (kind=kind_phys), pointer      :: evcw(:)            => null()  !<
    real (kind=kind_phys), pointer      :: pah(:)             => null()  !<
    real (kind=kind_phys), pointer      :: ecan(:)            => null()  !<
    real (kind=kind_phys), pointer      :: etran(:)           => null()  !<
    real (kind=kind_phys), pointer      :: edir(:)            => null()  !<
    real (kind=kind_phys), pointer      :: faerlw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: faersw(:,:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffhh_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fh2(:)             => null()  !<
    real (kind=kind_phys), pointer      :: fh2_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: fh2_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fh2_water(:)       => null()  !<
    logical,               pointer      :: flag_cice(:)       => null()  !<
    logical,               pointer      :: flag_guess(:)      => null()  !<
    logical,               pointer      :: flag_iter(:)       => null()  !<
    logical,               pointer      :: flag_lakefreeze(:) => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ffmm_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: fm10(:)            => null()  !<
    real (kind=kind_phys), pointer      :: fm10_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fm10_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: fm10_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: frland(:)          => null()  !<
    real (kind=kind_phys), pointer      :: fscav(:)           => null()  !<
    real (kind=kind_phys), pointer      :: fswtr(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_ice(:)    => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_land(:)   => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw_water(:)  => null()  !<
    real (kind=kind_phys), pointer      :: gamma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: gamq(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gamt(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gasvmr(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: gflx(:)            => null()  !<
    real (kind=kind_phys), pointer      :: gflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: gflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gflx_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: graupelmp(:)       => null()  !<
    real (kind=kind_phys), pointer      :: gwdcu(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: gwdcv(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: zvfun(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hffac(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hflxq(:)           => null()  !<
    real (kind=kind_phys), pointer      :: hflx_ice(:)        => null()  !<
    real (kind=kind_phys), pointer      :: hflx_land(:)       => null()  !<
    real (kind=kind_phys), pointer      :: hflx_water(:)      => null()  !<
    !--- radiation variables that need to be carried over from radiation to physics
    real (kind=kind_phys), pointer      :: htlwc(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htlw0(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htswc(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: htsw0(:,:)         => null()  !<
    !
    real (kind=kind_phys), pointer      :: icemp(:)           => null()  !<
    logical,               pointer      :: dry(:)             => null()  !<
    integer,               pointer      :: idxday(:)          => null()  !<
    logical,               pointer      :: icy(:)             => null()  !<
    logical,               pointer      :: lake(:)            => null()  !<
    logical,               pointer      :: ocean(:)           => null()  !<
    integer,               pointer      :: islmsk(:)          => null()  !<
    integer,               pointer      :: islmsk_cice(:)     => null()  !<
    logical,               pointer      :: wet(:)             => null()  !<
    integer                             :: kb                            !<
    integer,               pointer      :: kbot(:)            => null()  !<
    integer,               pointer      :: kcnv(:)            => null()  !<
    integer                             :: kd                            !<
    integer,               pointer      :: kinver(:)          => null()  !<
    integer,               pointer      :: kpbl(:)            => null()  !<
    integer                             :: kt                            !<
    integer,               pointer      :: ktop(:)            => null()  !<
    integer,               pointer      :: mbota(:,:)         => null()  !<
    integer,               pointer      :: mtopa(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: ncgl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncpr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: ncps(:,:)          => null()  !<
    integer                             :: nday                          !<
    real (kind=kind_phys), pointer      :: oa4(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: oc(:)              => null()  !<
    real (kind=kind_phys), pointer      :: olyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: plvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: plyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: prcpmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: prnum(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: q2mp(:)            => null()  !<
    real (kind=kind_phys), pointer      :: qgl(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qicn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlcn(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qrn(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: qsnw(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: qss_ice(:)         => null()  !<
    real (kind=kind_phys), pointer      :: qss_land(:)        => null()  !<
    real (kind=kind_phys), pointer      :: qss_water(:)       => null()  !<
    logical                             :: fullradar_diag                !<
    real (kind=kind_phys)               :: raddt                         !<
    real (kind=kind_phys), pointer      :: rainmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincd(:)          => null()  !<
    real (kind=kind_phys), pointer      :: raincs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rainmcadj(:)       => null()  !<
    real (kind=kind_phys), pointer      :: rainp(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb(:)              => null()  !<
    real (kind=kind_phys), pointer      :: rb_ice(:)          => null()  !<
    real (kind=kind_phys), pointer      :: rb_land(:)         => null()  !<
    real (kind=kind_phys), pointer      :: rb_water(:)        => null()  !<
    logical                             :: max_hourly_reset              !<
    logical                             :: ext_diag_thompson_reset       !<
    real (kind=kind_phys), pointer      :: rhc(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: runoff(:)          => null()  !<
    real (kind=kind_phys), pointer      :: save_q(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_t(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_tcp(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_u(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_v(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sbsno(:)           => null()  !<
    type (cmpfsw_type),    pointer      :: scmpsw(:)          => null()  !<
    real (kind=kind_phys), pointer      :: sfcalb(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: sigma(:)           => null()  !<
    real (kind=kind_phys), pointer      :: sigmaf(:)          => null()  !<
    real (kind=kind_phys), pointer      :: sigmafrac(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: sigmatot(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: snowc(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snohf(:)           => null()  !<
    real (kind=kind_phys), pointer      :: snowmp(:)          => null()  !<
    real (kind=kind_phys), pointer      :: snowmt(:)          => null()  !<
    real (kind=kind_phys), pointer      :: stress(:)          => null()  !<
    real (kind=kind_phys), pointer      :: stress_ice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: stress_land(:)     => null()  !<
    real (kind=kind_phys), pointer      :: stress_water(:)    => null()  !<
    real (kind=kind_phys), pointer      :: t2mmp(:)           => null()  !<
    real (kind=kind_phys), pointer      :: theta(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tlvl(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tkeh(:,:)          => null()  !< vertical turbulent kinetic energy (m2/s2) at the model layer interfaces
    real (kind=kind_phys), pointer      :: tlyr(:,:)          => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tprcp_water(:)     => null()  !<
    integer                             :: tracers_start_index           !<
    integer                             :: tracers_water                 !<
    real (kind=kind_phys), pointer      :: trans(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tseal(:)           => null()  !<
    real (kind=kind_phys), pointer      :: tsfa(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsfc_water(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tsfg(:)            => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: tsurf_water(:)     => null()  !<
    real (kind=kind_phys), pointer      :: ud_mf(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: uustar_ice(:)      => null()  !<
    real (kind=kind_phys), pointer      :: uustar_land(:)     => null()  !<
    real (kind=kind_phys), pointer      :: uustar_water(:)    => null()  !<
    real (kind=kind_phys), pointer      :: vdftra(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: vegf1d(:)          => null()  !<
    real (kind=kind_phys)               :: lndp_vgf                      !<

    real (kind=kind_phys), pointer      :: w_upi(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: wcbmax(:)          => null()  !<
    real (kind=kind_phys), pointer      :: wind(:)            => null()  !<
    real (kind=kind_phys), pointer      :: work1(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work2(:)           => null()  !<
    real (kind=kind_phys), pointer      :: work3(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xcosz(:)           => null()  !<
    real (kind=kind_phys), pointer      :: xlai1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: xmu(:)             => null()  !<
    real (kind=kind_phys), pointer      :: z01d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: zt1d(:)            => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_ice(:)       => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_land(:)      => null()  !<
    real (kind=kind_phys), pointer      :: ztmax_water(:)     => null()  !<
!==================================================================================================
! UGWP - five mechnanisms of momentum deposition due to various types of GWs
! (oss, ofd, obl, ogw) + ngw = sum( sso + ngw)
!==================================================================================================
! nGWs
    real (kind=kind_phys), pointer      :: dudt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: dvdt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: dtdt_ngw(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: kdis_ngw(:,:)      => null()  !<

    real (kind=kind_phys), pointer      :: tau_oss(: )        => null()  !< instantaneous momentum flux due to OSS
    real (kind=kind_phys), pointer      :: tau_tofd(:)        => null()  !< instantaneous momentum flux due to TOFD
    real (kind=kind_phys), pointer      :: tau_mtb(:)         => null()  !< instantaneous momentum of mountain blocking drag
    real (kind=kind_phys), pointer      :: tau_ogw(:)         => null()  !< instantaneous momentum flux of OGWs
    real (kind=kind_phys), pointer      :: tau_ngw(:)         => null()  !< instantaneous momentum flux of NGWs

    real (kind=kind_phys), pointer      :: zngw(:)            => null()  !< launch levels of NGWs
    real (kind=kind_phys), pointer      :: zmtb(:)            => null()  !< mountain blocking height
    real (kind=kind_phys), pointer      :: zlwb(:)            => null()  !< low level wave breaking height
    real (kind=kind_phys), pointer      :: zogw(:)            => null()  !< height of OGW-launch

    real (kind=kind_phys), pointer      :: dudt_mtb(:,:)      => null()  !< daily aver u-wind tend due to mountain blocking
    real (kind=kind_phys), pointer      :: dudt_tms(:,:)      => null()  !< daily aver u-wind tend due to TMS

    ! RRTMGP
    real (kind=kind_phys), pointer      :: p_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: p_lev(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: t_lev(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: t_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: relhum(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: tv_lay(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: qs_lay(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: q_lay(:,:)                => null()  !<
    real (kind=kind_phys), pointer      :: deltaZ(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: deltaZc(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: deltaP(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: cloud_overlap_param(:,:)  => null()  !< Cloud overlap parameter
    real (kind=kind_phys), pointer      :: cnv_cloud_overlap_param(:,:) => null()  !< Convective cloud overlap parameter
    real (kind=kind_phys), pointer      :: precip_overlap_param(:,:) => null()  !< Precipitation overlap parameter
    real (kind=kind_phys), pointer      :: tracer(:,:,:)             => null()  !<
    real (kind=kind_phys), pointer      :: aerosolslw(:,:,:,:)       => null()  !< Aerosol radiative properties in each LW band.
    real (kind=kind_phys), pointer      :: aerosolssw(:,:,:,:)       => null()  !< Aerosol radiative properties in each SW band.
    real (kind=kind_phys), pointer      :: precip_frac(:,:)          => null()  !< Precipitation fraction
    real (kind=kind_phys), pointer      :: cld_cnv_frac(:,:)         => null()  !< SGS convective cloud fraction
    real (kind=kind_phys), pointer      :: cld_cnv_lwp(:,:)          => null()  !< SGS convective cloud liquid water path
    real (kind=kind_phys), pointer      :: cld_cnv_reliq(:,:)        => null()  !< SGS convective cloud liquid effective radius
    real (kind=kind_phys), pointer      :: cld_cnv_iwp(:,:)          => null()  !< SGS convective cloud ice water path
    real (kind=kind_phys), pointer      :: cld_cnv_reice(:,:)        => null()  !< SGS convective cloud ice effecive radius
    real (kind=kind_phys), pointer      :: cld_pbl_lwp(:,:)          => null()  !< SGS PBL        cloud liquid water path
    real (kind=kind_phys), pointer      :: cld_pbl_reliq(:,:)        => null()  !< SGS PBL        cloud liquid effective radius
    real (kind=kind_phys), pointer      :: cld_pbl_iwp(:,:)          => null()  !< SGS PBL        cloud ice water path
    real (kind=kind_phys), pointer      :: cld_pbl_reice(:,:)        => null()  !< SGS PBL        cloud ice effecive radius
    real (kind=kind_phys), pointer      :: fluxlwUP_allsky(:,:)      => null()  !< RRTMGP upward   longwave  all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwDOWN_allsky(:,:)    => null()  !< RRTMGP downward longwave  all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwUP_clrsky(:,:)      => null()  !< RRTMGP upward   longwave  clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxlwDOWN_clrsky(:,:)    => null()  !< RRTMGP downward longwave  clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswUP_allsky(:,:)      => null()  !< RRTMGP upward   shortwave all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswDOWN_allsky(:,:)    => null()  !< RRTMGP downward shortwave all-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswUP_clrsky(:,:)      => null()  !< RRTMGP upward   shortwave clr-sky flux profile
    real (kind=kind_phys), pointer      :: fluxswDOWN_clrsky(:,:)    => null()  !< RRTMGP downward shortwave clr-sky flux profile
    real (kind=kind_phys), pointer      :: sfc_emiss_byband(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: sec_diff_byband(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_nir_dir(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_nir_dif(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_uvvis_dir(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: sfc_alb_uvvis_dif(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: toa_src_lw(:,:)           => null()  !<
    real (kind=kind_phys), pointer      :: toa_src_sw(:,:)           => null()  !<
    type(proflw_type), pointer          :: flxprf_lw(:,:)            => null()  !< DDT containing RRTMGP longwave fluxes
    type(profsw_type), pointer          :: flxprf_sw(:,:)            => null()  !< DDT containing RRTMGP shortwave fluxes
    real (kind=kind_phys), pointer      :: vmr_o2(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: vmr_h2o(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_o3(:,:)               => null()  !<
    real (kind=kind_phys), pointer      :: vmr_ch4(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_n2o(:,:)              => null()  !<
    real (kind=kind_phys), pointer      :: vmr_co2(:,:)              => null()  !<

    !-- GSL drag suite
    real (kind=kind_phys), pointer      :: varss(:)           => null()  !<
    real (kind=kind_phys), pointer      :: ocss(:)            => null()  !<
    real (kind=kind_phys), pointer      :: oa4ss(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: clxss(:,:)         => null()  !<

    !-- 3D diagnostics
    integer :: rtg_ozone_index, rtg_tke_index

    !-- CCPP suite simulator
    real (kind=kind_phys), pointer      :: active_phys_tend(:,:,:) => null() ! tendencies for active physics process

    contains

      procedure :: create      => gfs_interstitial_create     !<   allocate array data
      procedure :: destroy     => gfs_interstitial_destroy    !<   deallocate array data
      procedure :: reset       => gfs_interstitial_reset      !<   reset array data

  end type GFS_interstitial_type

contains

!----------------------
! GFS_interstitial_type
!----------------------

  subroutine gfs_interstitial_create (Interstitial, ixs, ixe, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    integer,                intent(in) :: ixs, ixe
    type(GFS_control_type), intent(in) :: Model
    integer                            :: iGas
    !
    ! Allocate arrays
    allocate (Interstitial%adjsfculw_land  (ixs:ixe))
    allocate (Interstitial%adjsfculw_ice   (ixs:ixe))
    allocate (Interstitial%adjsfculw_water (ixs:ixe))
    allocate (Interstitial%adjnirbmd       (ixs:ixe))
    allocate (Interstitial%adjnirbmu       (ixs:ixe))
    allocate (Interstitial%adjnirdfd       (ixs:ixe))
    allocate (Interstitial%adjnirdfu       (ixs:ixe))
    allocate (Interstitial%adjvisbmd       (ixs:ixe))
    allocate (Interstitial%adjvisbmu       (ixs:ixe))
    allocate (Interstitial%adjvisdfu       (ixs:ixe))
    allocate (Interstitial%adjvisdfd       (ixs:ixe))
    allocate (Interstitial%aerodp          (ixs:ixe,Model%NSPC1))
    allocate (Interstitial%alb1d           (ixs:ixe))
    if (.not. Model%do_RRTMGP) then
      ! RRTMGP uses its own cloud_overlap_param
      allocate (Interstitial%alpha         (ixs:ixe,Model%levr+Model%LTP))
    end if
    allocate (Interstitial%bexp1d          (ixs:ixe))
    allocate (Interstitial%cd              (ixs:ixe))
    allocate (Interstitial%cd_ice          (ixs:ixe))
    allocate (Interstitial%cd_land         (ixs:ixe))
    allocate (Interstitial%cd_water        (ixs:ixe))
    allocate (Interstitial%cdq             (ixs:ixe))
    allocate (Interstitial%cdq_ice         (ixs:ixe))
    allocate (Interstitial%cdq_land        (ixs:ixe))
    allocate (Interstitial%cdq_water       (ixs:ixe))
    allocate (Interstitial%chh_ice         (ixs:ixe))
    allocate (Interstitial%chh_land        (ixs:ixe))
    allocate (Interstitial%chh_water       (ixs:ixe))
    allocate (Interstitial%cldf            (ixs:ixe))
    allocate (Interstitial%cldsa           (ixs:ixe,5))
    allocate (Interstitial%cldtaulw        (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%cldtausw        (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%cld1d           (ixs:ixe))
    allocate (Interstitial%clouds          (ixs:ixe,Model%levr+Model%LTP,Model%NF_CLDS))
    allocate (Interstitial%clw             (ixs:ixe,Model%levs,Model%nn))
    allocate (Interstitial%clx             (ixs:ixe,4))
    allocate (Interstitial%cmm_ice         (ixs:ixe))
    allocate (Interstitial%cmm_land        (ixs:ixe))
    allocate (Interstitial%cmm_water       (ixs:ixe))
    allocate (Interstitial%cnvc            (ixs:ixe,Model%levs))
    allocate (Interstitial%ctei_r          (ixs:ixe))
    allocate (Interstitial%ctei_rml        (ixs:ixe))
    allocate (Interstitial%cumabs          (ixs:ixe))
    allocate (Interstitial%dd_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%de_lgth         (ixs:ixe))
    allocate (Interstitial%del             (ixs:ixe,Model%levs))
    allocate (Interstitial%del_gz          (ixs:ixe,Model%levs+1))
    allocate (Interstitial%delr            (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%dlength         (ixs:ixe))
    allocate (Interstitial%dqdt            (ixs:ixe,Model%levs,Model%ntrac))
    allocate (Interstitial%dqsfc1          (ixs:ixe))
    allocate (Interstitial%drain           (ixs:ixe))
    allocate (Interstitial%dtdt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dtsfc1          (ixs:ixe))
    allocate (Interstitial%dt_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%dtzm            (ixs:ixe))
    allocate (Interstitial%dudt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dusfcg          (ixs:ixe))
    allocate (Interstitial%dusfc1          (ixs:ixe))
    allocate (Interstitial%dvdt            (ixs:ixe,Model%levs))
    allocate (Interstitial%dvsfcg          (ixs:ixe))
    allocate (Interstitial%dvsfc1          (ixs:ixe))
    allocate (Interstitial%dvdftra         (ixs:ixe,Model%levs,Model%nvdiff))
    allocate (Interstitial%dzlyr           (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%elvmax          (ixs:ixe))
    allocate (Interstitial%ep1d            (ixs:ixe))
    allocate (Interstitial%ep1d_ice        (ixs:ixe))
    allocate (Interstitial%ep1d_land       (ixs:ixe))
    allocate (Interstitial%ep1d_water      (ixs:ixe))
    allocate (Interstitial%evap_ice        (ixs:ixe))
    allocate (Interstitial%evap_land       (ixs:ixe))
    allocate (Interstitial%evap_water      (ixs:ixe))
    allocate (Interstitial%evbs            (ixs:ixe))
    allocate (Interstitial%evcw            (ixs:ixe))
    allocate (Interstitial%pah             (ixs:ixe))
    allocate (Interstitial%ecan            (ixs:ixe))
    allocate (Interstitial%etran           (ixs:ixe))
    allocate (Interstitial%edir            (ixs:ixe))
    allocate (Interstitial%faerlw          (ixs:ixe,Model%levr+Model%LTP,Model%NBDLW,Model%NF_AELW))
    allocate (Interstitial%faersw          (ixs:ixe,Model%levr+Model%LTP,Model%NBDSW,Model%NF_AESW))
    allocate (Interstitial%ffhh_ice        (ixs:ixe))
    allocate (Interstitial%ffhh_land       (ixs:ixe))
    allocate (Interstitial%ffhh_water      (ixs:ixe))
    allocate (Interstitial%fh2             (ixs:ixe))
    allocate (Interstitial%fh2_ice         (ixs:ixe))
    allocate (Interstitial%fh2_land        (ixs:ixe))
    allocate (Interstitial%fh2_water       (ixs:ixe))
    allocate (Interstitial%flag_cice       (ixs:ixe))
    allocate (Interstitial%flag_guess      (ixs:ixe))
    allocate (Interstitial%flag_iter       (ixs:ixe))
    allocate (Interstitial%flag_lakefreeze (ixs:ixe))
    allocate (Interstitial%ffmm_ice        (ixs:ixe))
    allocate (Interstitial%ffmm_land       (ixs:ixe))
    allocate (Interstitial%ffmm_water      (ixs:ixe))
    allocate (Interstitial%fm10            (ixs:ixe))
    allocate (Interstitial%fm10_ice        (ixs:ixe))
    allocate (Interstitial%fm10_land       (ixs:ixe))
    allocate (Interstitial%fm10_water      (ixs:ixe))
    allocate (Interstitial%frland          (ixs:ixe))
    allocate (Interstitial%fscav           (Model%nscav))
    allocate (Interstitial%fswtr           (Model%nscav))
    allocate (Interstitial%gabsbdlw        (ixs:ixe))
    allocate (Interstitial%gabsbdlw_ice    (ixs:ixe))
    allocate (Interstitial%gabsbdlw_land   (ixs:ixe))
    allocate (Interstitial%gabsbdlw_water  (ixs:ixe))
    allocate (Interstitial%gamma           (ixs:ixe))
    allocate (Interstitial%gamq            (ixs:ixe))
    allocate (Interstitial%gamt            (ixs:ixe))
    allocate (Interstitial%gasvmr          (ixs:ixe,Model%levr+Model%LTP,Model%NF_VGAS))
    allocate (Interstitial%gflx            (ixs:ixe))
    allocate (Interstitial%gflx_ice        (ixs:ixe))
    allocate (Interstitial%gflx_land       (ixs:ixe))
    allocate (Interstitial%gflx_water      (ixs:ixe))
    allocate (Interstitial%gwdcu           (ixs:ixe,Model%levs))
    allocate (Interstitial%gwdcv           (ixs:ixe,Model%levs))
    allocate (Interstitial%zvfun           (ixs:ixe))
    allocate (Interstitial%hffac           (ixs:ixe))
    allocate (Interstitial%hflxq           (ixs:ixe))
    allocate (Interstitial%hflx_ice        (ixs:ixe))
    allocate (Interstitial%hflx_land       (ixs:ixe))
    allocate (Interstitial%hflx_water      (ixs:ixe))
    allocate (Interstitial%htlwc           (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%htlw0           (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%htswc           (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%htsw0           (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%dry             (ixs:ixe))
    allocate (Interstitial%idxday          (ixs:ixe))
    allocate (Interstitial%icy             (ixs:ixe))
    allocate (Interstitial%lake            (ixs:ixe))
    allocate (Interstitial%ocean           (ixs:ixe))
    allocate (Interstitial%islmsk          (ixs:ixe))
    allocate (Interstitial%islmsk_cice     (ixs:ixe))
    allocate (Interstitial%wet             (ixs:ixe))
    allocate (Interstitial%kbot            (ixs:ixe))
    allocate (Interstitial%kcnv            (ixs:ixe))
    allocate (Interstitial%kinver          (ixs:ixe))
    allocate (Interstitial%kpbl            (ixs:ixe))
    allocate (Interstitial%ktop            (ixs:ixe))
    allocate (Interstitial%mbota           (ixs:ixe,3))
    allocate (Interstitial%mtopa           (ixs:ixe,3))
    allocate (Interstitial%oa4             (ixs:ixe,4))
    allocate (Interstitial%oc              (ixs:ixe))
    allocate (Interstitial%olyr            (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%plvl            (ixs:ixe,Model%levr+1+Model%LTP))
    allocate (Interstitial%plyr            (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%prnum           (ixs:ixe,Model%levs))
    allocate (Interstitial%qlyr            (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%prcpmp          (ixs:ixe))
    allocate (Interstitial%qss_ice         (ixs:ixe))
    allocate (Interstitial%qss_land        (ixs:ixe))
    allocate (Interstitial%qss_water       (ixs:ixe))
    allocate (Interstitial%raincd          (ixs:ixe))
    allocate (Interstitial%raincs          (ixs:ixe))
    allocate (Interstitial%rainmcadj       (ixs:ixe))
    allocate (Interstitial%rainp           (ixs:ixe,Model%levs))
    allocate (Interstitial%rb              (ixs:ixe))
    allocate (Interstitial%rb_ice          (ixs:ixe))
    allocate (Interstitial%rb_land         (ixs:ixe))
    allocate (Interstitial%rb_water        (ixs:ixe))
    allocate (Interstitial%rhc             (ixs:ixe,Model%levs))
    allocate (Interstitial%runoff          (ixs:ixe))
    allocate (Interstitial%save_q          (ixs:ixe,Model%levs,Model%ntrac))
    allocate (Interstitial%save_t          (ixs:ixe,Model%levs))
    allocate (Interstitial%save_tcp        (ixs:ixe,Model%levs))
    allocate (Interstitial%save_u          (ixs:ixe,Model%levs))
    allocate (Interstitial%save_v          (ixs:ixe,Model%levs))
    allocate (Interstitial%sbsno           (ixs:ixe))
    allocate (Interstitial%scmpsw          (ixs:ixe))
    allocate (Interstitial%sfcalb          (ixs:ixe,Model%NF_ALBD))
    allocate (Interstitial%sigma           (ixs:ixe))
    allocate (Interstitial%sigmaf          (ixs:ixe))
    allocate (Interstitial%sigmafrac       (ixs:ixe,Model%levs))
    allocate (Interstitial%sigmatot        (ixs:ixe,Model%levs+1))
    allocate (Interstitial%snowc           (ixs:ixe))
    allocate (Interstitial%snohf           (ixs:ixe))
    allocate (Interstitial%snowmt          (ixs:ixe))
    allocate (Interstitial%stress          (ixs:ixe))
    allocate (Interstitial%stress_ice      (ixs:ixe))
    allocate (Interstitial%stress_land     (ixs:ixe))
    allocate (Interstitial%stress_water    (ixs:ixe))
    allocate (Interstitial%theta           (ixs:ixe))
    allocate (Interstitial%tkeh            (ixs:ixe,Model%levs+1)) !Vertical turbulent kinetic energy at model layer interfaces
    allocate (Interstitial%tlvl            (ixs:ixe,Model%levr+1+Model%LTP))
    allocate (Interstitial%tlyr            (ixs:ixe,Model%levr+Model%LTP))
    allocate (Interstitial%tprcp_ice       (ixs:ixe))
    allocate (Interstitial%tprcp_land      (ixs:ixe))
    allocate (Interstitial%tprcp_water     (ixs:ixe))
    allocate (Interstitial%trans           (ixs:ixe))
    allocate (Interstitial%tseal           (ixs:ixe))
    allocate (Interstitial%tsfa            (ixs:ixe))
    allocate (Interstitial%tsfc_water      (ixs:ixe))
    allocate (Interstitial%tsfg            (ixs:ixe))
    allocate (Interstitial%tsurf_ice       (ixs:ixe))
    allocate (Interstitial%tsurf_land      (ixs:ixe))
    allocate (Interstitial%tsurf_water     (ixs:ixe))
    allocate (Interstitial%ud_mf           (ixs:ixe,Model%levs))
    allocate (Interstitial%uustar_ice      (ixs:ixe))
    allocate (Interstitial%uustar_land     (ixs:ixe))
    allocate (Interstitial%uustar_water    (ixs:ixe))
    allocate (Interstitial%vdftra          (ixs:ixe,Model%levs,Model%nvdiff))  !GJF first dimension was set as 'IX' in GFS_physics_driver
    allocate (Interstitial%vegf1d          (ixs:ixe))
    allocate (Interstitial%wcbmax          (ixs:ixe))
    allocate (Interstitial%wind            (ixs:ixe))
    allocate (Interstitial%work1           (ixs:ixe))
    allocate (Interstitial%work2           (ixs:ixe))
    allocate (Interstitial%work3           (ixs:ixe))
    allocate (Interstitial%xcosz           (ixs:ixe))
    allocate (Interstitial%xlai1d          (ixs:ixe))
    allocate (Interstitial%xmu             (ixs:ixe))
    allocate (Interstitial%z01d            (ixs:ixe))
    allocate (Interstitial%zt1d            (ixs:ixe))
    allocate (Interstitial%ztmax_ice       (ixs:ixe))
    allocate (Interstitial%ztmax_land      (ixs:ixe))
    allocate (Interstitial%ztmax_water     (ixs:ixe))

    ! RRTMGP
    if (Model%do_RRTMGP) then
       allocate (Interstitial%tracer               (ixs:ixe, Model%levs,Model%ntrac))
       allocate (Interstitial%tv_lay               (ixs:ixe, Model%levs))
       allocate (Interstitial%relhum               (ixs:ixe, Model%levs))
       allocate (Interstitial%qs_lay               (ixs:ixe, Model%levs))
       allocate (Interstitial%q_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaZ               (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaZc              (ixs:ixe, Model%levs))
       allocate (Interstitial%deltaP               (ixs:ixe, Model%levs))
       allocate (Interstitial%p_lev                (ixs:ixe, Model%levs+1))
       allocate (Interstitial%p_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%t_lev                (ixs:ixe, Model%levs+1))
       allocate (Interstitial%t_lay                (ixs:ixe, Model%levs))
       allocate (Interstitial%cloud_overlap_param  (ixs:ixe, Model%levs))
       allocate (Interstitial%precip_overlap_param (ixs:ixe, Model%levs))
       allocate (Interstitial%fluxlwUP_allsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwDOWN_allsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwUP_clrsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxlwDOWN_clrsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswUP_allsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswDOWN_allsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswUP_clrsky      (ixs:ixe, Model%levs+1))
       allocate (Interstitial%fluxswDOWN_clrsky    (ixs:ixe, Model%levs+1))
       allocate (Interstitial%aerosolslw           (ixs:ixe, Model%levs, Model%rrtmgp_nBandsLW, Model%NF_AELW))
       allocate (Interstitial%aerosolssw           (ixs:ixe, Model%levs, Model%rrtmgp_nBandsSW, Model%NF_AESW))
       allocate (Interstitial%precip_frac          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_frac         (ixs:ixe, Model%levs))
       allocate (Interstitial%cnv_cloud_overlap_param(ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_lwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_reliq        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_iwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_cnv_reice        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_lwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_reliq        (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_iwp          (ixs:ixe, Model%levs))
       allocate (Interstitial%cld_pbl_reice        (ixs:ixe, Model%levs))
       allocate (Interstitial%flxprf_lw            (ixs:ixe, Model%levs+1))
       allocate (Interstitial%flxprf_sw            (ixs:ixe, Model%levs+1))
       allocate (Interstitial%sfc_emiss_byband     (Model%rrtmgp_nBandsLW,ixs:ixe))
       allocate (Interstitial%sec_diff_byband      (Model%rrtmgp_nBandsLW,ixs:ixe))
       allocate (Interstitial%sfc_alb_nir_dir      (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_nir_dif      (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_uvvis_dir    (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%sfc_alb_uvvis_dif    (Model%rrtmgp_nBandsSW,ixs:ixe))
       allocate (Interstitial%toa_src_sw           (ixs:ixe,Model%rrtmgp_nGptsSW))
       allocate (Interstitial%toa_src_lw           (ixs:ixe,Model%rrtmgp_nGptsLW))
       allocate (Interstitial%vmr_o2               (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_h2o              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_o3               (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_ch4              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_n2o              (ixs:ixe, Model%levs))
       allocate (Interstitial%vmr_co2              (ixs:ixe, Model%levs))

    end if

! UGWP common
    allocate (Interstitial%tau_mtb         (ixs:ixe))
    allocate (Interstitial%tau_ogw         (ixs:ixe))
    allocate (Interstitial%tau_tofd        (ixs:ixe))
    allocate (Interstitial%tau_ngw         (ixs:ixe))
    allocate (Interstitial%tau_oss         (ixs:ixe))
    allocate (Interstitial%dudt_mtb        (ixs:ixe,Model%levs))
    allocate (Interstitial%dudt_tms        (ixs:ixe,Model%levs))
    allocate (Interstitial%zmtb            (ixs:ixe)           )
    allocate (Interstitial%zlwb            (ixs:ixe)           )
    allocate (Interstitial%zogw            (ixs:ixe)           )
    allocate (Interstitial%zngw            (ixs:ixe)           )

! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      allocate (Interstitial%dudt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%dvdt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%dtdt_ngw        (ixs:ixe,Model%levs))
      allocate (Interstitial%kdis_ngw        (ixs:ixe,Model%levs))
    end if

!-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       allocate (Interstitial%varss           (ixs:ixe))
       allocate (Interstitial%ocss            (ixs:ixe))
       allocate (Interstitial%oa4ss           (ixs:ixe,4))
       allocate (Interstitial%clxss           (ixs:ixe,4))
    end if
!
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       allocate (Interstitial%graupelmp  (ixs:ixe))
       allocate (Interstitial%icemp      (ixs:ixe))
       allocate (Interstitial%rainmp     (ixs:ixe))
       allocate (Interstitial%snowmp     (ixs:ixe))
    else if (Model%imp_physics == Model%imp_physics_mg) then
       allocate (Interstitial%ncgl       (ixs:ixe,Model%levs))
       allocate (Interstitial%ncpr       (ixs:ixe,Model%levs))
       allocate (Interstitial%ncps       (ixs:ixe,Model%levs))
       allocate (Interstitial%qgl        (ixs:ixe,Model%levs))
       allocate (Interstitial%qrn        (ixs:ixe,Model%levs))
       allocate (Interstitial%qsnw       (ixs:ixe,Model%levs))
       allocate (Interstitial%qlcn       (ixs:ixe,Model%levs))
       allocate (Interstitial%qicn       (ixs:ixe,Model%levs))
       allocate (Interstitial%w_upi      (ixs:ixe,Model%levs))
       allocate (Interstitial%cf_upi     (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_mfd    (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_dqldt  (ixs:ixe,Model%levs))
       allocate (Interstitial%clcn       (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_fice   (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_ndrop  (ixs:ixe,Model%levs))
       allocate (Interstitial%cnv_nice   (ixs:ixe,Model%levs))
    end if
    if (Model%lsm == Model%ilsm_noahmp) then
       allocate (Interstitial%t2mmp (ixs:ixe))
       allocate (Interstitial%q2mp  (ixs:ixe))
    end if
    !
    ! CCPP suite simulator
    if (Model%do_ccpp_suite_sim) then
       allocate (Interstitial%active_phys_tend(ixs:ixe,Model%levs,Model%physics_process(1)%nprg_active))
    endif

    ! Reset all variables
    call Interstitial%reset (Model)
    !
  end subroutine gfs_interstitial_create

  subroutine gfs_interstitial_destroy (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    type(GFS_control_type), intent(in) :: Model

    ! Allocate arrays
    deallocate (Interstitial%adjsfculw_land)
    deallocate (Interstitial%adjsfculw_ice)
    deallocate (Interstitial%adjsfculw_water)
    deallocate (Interstitial%adjnirbmd)
    deallocate (Interstitial%adjnirbmu)
    deallocate (Interstitial%adjnirdfd)
    deallocate (Interstitial%adjnirdfu)
    deallocate (Interstitial%adjvisbmd)
    deallocate (Interstitial%adjvisbmu)
    deallocate (Interstitial%adjvisdfu)
    deallocate (Interstitial%adjvisdfd)
    deallocate (Interstitial%aerodp)
    deallocate (Interstitial%alb1d)
    if (.not. Model%do_RRTMGP) then
      deallocate (Interstitial%alpha)
    end if
    deallocate (Interstitial%bexp1d)
    deallocate (Interstitial%cd)
    deallocate (Interstitial%cd_ice)
    deallocate (Interstitial%cd_land)
    deallocate (Interstitial%cd_water)
    deallocate (Interstitial%cdq)
    deallocate (Interstitial%cdq_ice)
    deallocate (Interstitial%cdq_land)
    deallocate (Interstitial%cdq_water)
    deallocate (Interstitial%chh_ice)
    deallocate (Interstitial%chh_land)
    deallocate (Interstitial%chh_water)
    deallocate (Interstitial%cldf)
    deallocate (Interstitial%cldsa)
    deallocate (Interstitial%cldtaulw)
    deallocate (Interstitial%cldtausw)
    deallocate (Interstitial%cld1d)
    deallocate (Interstitial%clouds)
    deallocate (Interstitial%clw)
    deallocate (Interstitial%clx)
    deallocate (Interstitial%cmm_ice)
    deallocate (Interstitial%cmm_land)
    deallocate (Interstitial%cmm_water)
    deallocate (Interstitial%cnvc)
    deallocate (Interstitial%ctei_r)
    deallocate (Interstitial%ctei_rml)
    deallocate (Interstitial%cumabs)
    deallocate (Interstitial%dd_mf)
    deallocate (Interstitial%de_lgth)
    deallocate (Interstitial%del)
    deallocate (Interstitial%del_gz)
    deallocate (Interstitial%delr)
    deallocate (Interstitial%dlength)
    deallocate (Interstitial%dqdt)
    deallocate (Interstitial%dqsfc1)
    deallocate (Interstitial%drain)
    deallocate (Interstitial%dtdt)
    deallocate (Interstitial%dtsfc1)
    deallocate (Interstitial%dt_mf)
    deallocate (Interstitial%dtzm)
    deallocate (Interstitial%dudt)
    deallocate (Interstitial%dusfcg)
    deallocate (Interstitial%dusfc1)
    deallocate (Interstitial%dvdt)
    deallocate (Interstitial%dvsfcg)
    deallocate (Interstitial%dvsfc1)
    deallocate (Interstitial%dvdftra)
    deallocate (Interstitial%dzlyr)
    deallocate (Interstitial%elvmax)
    deallocate (Interstitial%ep1d)
    deallocate (Interstitial%ep1d_ice)
    deallocate (Interstitial%ep1d_land)
    deallocate (Interstitial%ep1d_water)
    deallocate (Interstitial%evap_ice)
    deallocate (Interstitial%evap_land)
    deallocate (Interstitial%evap_water)
    deallocate (Interstitial%evbs)
    deallocate (Interstitial%evcw)
    deallocate (Interstitial%pah)
    deallocate (Interstitial%ecan)
    deallocate (Interstitial%etran)
    deallocate (Interstitial%edir)
    deallocate (Interstitial%faerlw)
    deallocate (Interstitial%faersw)
    deallocate (Interstitial%ffhh_ice)
    deallocate (Interstitial%ffhh_land)
    deallocate (Interstitial%ffhh_water)
    deallocate (Interstitial%fh2)
    deallocate (Interstitial%fh2_ice)
    deallocate (Interstitial%fh2_land)
    deallocate (Interstitial%fh2_water)
    deallocate (Interstitial%flag_cice)
    deallocate (Interstitial%flag_guess)
    deallocate (Interstitial%flag_iter)
    deallocate (Interstitial%flag_lakefreeze)
    deallocate (Interstitial%ffmm_ice)
    deallocate (Interstitial%ffmm_land)
    deallocate (Interstitial%ffmm_water)
    deallocate (Interstitial%fm10)
    deallocate (Interstitial%fm10_ice)
    deallocate (Interstitial%fm10_land)
    deallocate (Interstitial%fm10_water)
    deallocate (Interstitial%frland)
    deallocate (Interstitial%fscav)
    deallocate (Interstitial%fswtr)
    deallocate (Interstitial%gabsbdlw)
    deallocate (Interstitial%gabsbdlw_ice)
    deallocate (Interstitial%gabsbdlw_land)
    deallocate (Interstitial%gabsbdlw_water)
    deallocate (Interstitial%gamma)
    deallocate (Interstitial%gamq)
    deallocate (Interstitial%gamt)
    deallocate (Interstitial%gasvmr)
    deallocate (Interstitial%gflx)
    deallocate (Interstitial%gflx_ice)
    deallocate (Interstitial%gflx_land)
    deallocate (Interstitial%gflx_water)
    deallocate (Interstitial%gwdcu)
    deallocate (Interstitial%gwdcv)
    deallocate (Interstitial%zvfun)
    deallocate (Interstitial%hffac)
    deallocate (Interstitial%hflxq)
    deallocate (Interstitial%hflx_ice)
    deallocate (Interstitial%hflx_land)
    deallocate (Interstitial%hflx_water)
    deallocate (Interstitial%htlwc)
    deallocate (Interstitial%htlw0)
    deallocate (Interstitial%htswc)
    deallocate (Interstitial%htsw0)
    deallocate (Interstitial%dry)
    deallocate (Interstitial%idxday)
    deallocate (Interstitial%icy)
    deallocate (Interstitial%lake)
    deallocate (Interstitial%ocean)
    deallocate (Interstitial%islmsk)
    deallocate (Interstitial%islmsk_cice)
    deallocate (Interstitial%wet)
    deallocate (Interstitial%kbot)
    deallocate (Interstitial%kcnv)
    deallocate (Interstitial%kinver)
    deallocate (Interstitial%kpbl)
    deallocate (Interstitial%ktop)
    deallocate (Interstitial%mbota)
    deallocate (Interstitial%mtopa)
    deallocate (Interstitial%oa4)
    deallocate (Interstitial%oc)
    deallocate (Interstitial%olyr)
    deallocate (Interstitial%plvl)
    deallocate (Interstitial%plyr)
    deallocate (Interstitial%prnum)
    deallocate (Interstitial%qlyr)
    deallocate (Interstitial%prcpmp)
    deallocate (Interstitial%qss_ice)
    deallocate (Interstitial%qss_land)
    deallocate (Interstitial%qss_water)
    deallocate (Interstitial%raincd)
    deallocate (Interstitial%raincs)
    deallocate (Interstitial%rainmcadj)
    deallocate (Interstitial%rainp)
    deallocate (Interstitial%rb)
    deallocate (Interstitial%rb_ice)
    deallocate (Interstitial%rb_land)
    deallocate (Interstitial%rb_water)
    deallocate (Interstitial%rhc)
    deallocate (Interstitial%runoff)
    deallocate (Interstitial%save_q)
    deallocate (Interstitial%save_t)
    deallocate (Interstitial%save_tcp)
    deallocate (Interstitial%save_u)
    deallocate (Interstitial%save_v)
    deallocate (Interstitial%sbsno)
    deallocate (Interstitial%scmpsw)
    deallocate (Interstitial%sfcalb)
    deallocate (Interstitial%sigma)
    deallocate (Interstitial%sigmaf)
    deallocate (Interstitial%sigmafrac)
    deallocate (Interstitial%sigmatot)
    deallocate (Interstitial%snowc)
    deallocate (Interstitial%snohf)
    deallocate (Interstitial%snowmt)
    deallocate (Interstitial%stress)
    deallocate (Interstitial%stress_ice)
    deallocate (Interstitial%stress_land)
    deallocate (Interstitial%stress_water)
    deallocate (Interstitial%theta)
    deallocate (Interstitial%tkeh)
    deallocate (Interstitial%tlvl)
    deallocate (Interstitial%tlyr)
    deallocate (Interstitial%tprcp_ice)
    deallocate (Interstitial%tprcp_land)
    deallocate (Interstitial%tprcp_water)
    deallocate (Interstitial%trans)
    deallocate (Interstitial%tseal)
    deallocate (Interstitial%tsfa)
    deallocate (Interstitial%tsfc_water)
    deallocate (Interstitial%tsfg)
    deallocate (Interstitial%tsurf_ice)
    deallocate (Interstitial%tsurf_land)
    deallocate (Interstitial%tsurf_water)
    deallocate (Interstitial%ud_mf)
    deallocate (Interstitial%uustar_ice)
    deallocate (Interstitial%uustar_land)
    deallocate (Interstitial%uustar_water)
    deallocate (Interstitial%vdftra)
    deallocate (Interstitial%vegf1d)
    deallocate (Interstitial%wcbmax)
    deallocate (Interstitial%wind)
    deallocate (Interstitial%work1)
    deallocate (Interstitial%work2)
    deallocate (Interstitial%work3)
    deallocate (Interstitial%xcosz)
    deallocate (Interstitial%xlai1d)
    deallocate (Interstitial%xmu)
    deallocate (Interstitial%z01d)
    deallocate (Interstitial%zt1d)
    deallocate (Interstitial%ztmax_ice)
    deallocate (Interstitial%ztmax_land)
    deallocate (Interstitial%ztmax_water)

    ! RRTMGP
    if (Model%do_RRTMGP) then
       deallocate (Interstitial%tracer)
       deallocate (Interstitial%tv_lay)
       deallocate (Interstitial%relhum)
       deallocate (Interstitial%qs_lay)
       deallocate (Interstitial%q_lay)
       deallocate (Interstitial%deltaZ)
       deallocate (Interstitial%deltaZc)
       deallocate (Interstitial%deltaP)
       deallocate (Interstitial%p_lev)
       deallocate (Interstitial%p_lay)
       deallocate (Interstitial%t_lev)
       deallocate (Interstitial%t_lay)
       deallocate (Interstitial%cloud_overlap_param)
       deallocate (Interstitial%precip_overlap_param)
       deallocate (Interstitial%fluxlwUP_allsky)
       deallocate (Interstitial%fluxlwDOWN_allsky)
       deallocate (Interstitial%fluxlwUP_clrsky)
       deallocate (Interstitial%fluxlwDOWN_clrsky)
       deallocate (Interstitial%fluxswUP_allsky)
       deallocate (Interstitial%fluxswDOWN_allsky)
       deallocate (Interstitial%fluxswUP_clrsky)
       deallocate (Interstitial%fluxswDOWN_clrsky)
       deallocate (Interstitial%aerosolslw)
       deallocate (Interstitial%aerosolssw)
       deallocate (Interstitial%precip_frac)
       deallocate (Interstitial%cld_cnv_frac)
       deallocate (Interstitial%cnv_cloud_overlap_param)
       deallocate (Interstitial%cld_cnv_lwp)
       deallocate (Interstitial%cld_cnv_reliq)
       deallocate (Interstitial%cld_cnv_iwp)
       deallocate (Interstitial%cld_cnv_reice)
       deallocate (Interstitial%cld_pbl_lwp)
       deallocate (Interstitial%cld_pbl_reliq)
       deallocate (Interstitial%cld_pbl_iwp)
       deallocate (Interstitial%cld_pbl_reice)
       deallocate (Interstitial%flxprf_lw)
       deallocate (Interstitial%flxprf_sw)
       deallocate (Interstitial%sfc_emiss_byband)
       deallocate (Interstitial%sec_diff_byband)
       deallocate (Interstitial%sfc_alb_nir_dir)
       deallocate (Interstitial%sfc_alb_nir_dif)
       deallocate (Interstitial%sfc_alb_uvvis_dir)
       deallocate (Interstitial%sfc_alb_uvvis_dif)
       deallocate (Interstitial%toa_src_sw)
       deallocate (Interstitial%toa_src_lw)
       deallocate (Interstitial%vmr_o2)
       deallocate (Interstitial%vmr_h2o)
       deallocate (Interstitial%vmr_o3)
       deallocate (Interstitial%vmr_ch4)
       deallocate (Interstitial%vmr_n2o)
       deallocate (Interstitial%vmr_co2)
    end if

    ! UGWP common
    deallocate (Interstitial%tau_mtb)
    deallocate (Interstitial%tau_ogw)
    deallocate (Interstitial%tau_tofd)
    deallocate (Interstitial%tau_ngw)
    deallocate (Interstitial%tau_oss)
    deallocate (Interstitial%dudt_mtb)
    deallocate (Interstitial%dudt_tms)
    deallocate (Interstitial%zmtb)
    deallocate (Interstitial%zlwb)
    deallocate (Interstitial%zogw)
    deallocate (Interstitial%zngw)

    ! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      deallocate (Interstitial%dudt_ngw)
      deallocate (Interstitial%dvdt_ngw)
      deallocate (Interstitial%dtdt_ngw)
      deallocate (Interstitial%kdis_ngw)
    end if

    !-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       deallocate (Interstitial%varss)
       deallocate (Interstitial%ocss)
       deallocate (Interstitial%oa4ss)
       deallocate (Interstitial%clxss)
    end if

    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       deallocate (Interstitial%graupelmp)
       deallocate (Interstitial%icemp)
       deallocate (Interstitial%rainmp)
       deallocate (Interstitial%snowmp)
    else if (Model%imp_physics == Model%imp_physics_mg) then
       deallocate (Interstitial%ncgl)
       deallocate (Interstitial%ncpr)
       deallocate (Interstitial%ncps)
       deallocate (Interstitial%qgl)
       deallocate (Interstitial%qrn)
       deallocate (Interstitial%qsnw)
       deallocate (Interstitial%qlcn)
       deallocate (Interstitial%qicn)
       deallocate (Interstitial%w_upi)
       deallocate (Interstitial%cf_upi)
       deallocate (Interstitial%cnv_mfd)
       deallocate (Interstitial%cnv_dqldt)
       deallocate (Interstitial%clcn)
       deallocate (Interstitial%cnv_fice)
       deallocate (Interstitial%cnv_ndrop)
       deallocate (Interstitial%cnv_nice)
    end if
    if (Model%lsm == Model%ilsm_noahmp) then
       deallocate (Interstitial%t2mmp)
       deallocate (Interstitial%q2mp)
    end if
    
  end subroutine gfs_interstitial_destroy

  subroutine gfs_interstitial_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    !
    Interstitial%adjsfculw_land  = Model%clear_val
    Interstitial%adjsfculw_ice   = Model%clear_val
    Interstitial%adjsfculw_water = Model%clear_val
    Interstitial%adjnirbmd       = Model%clear_val
    Interstitial%adjnirbmu       = Model%clear_val
    Interstitial%adjnirdfd       = Model%clear_val
    Interstitial%adjnirdfu       = Model%clear_val
    Interstitial%adjvisbmd       = Model%clear_val
    Interstitial%adjvisbmu       = Model%clear_val
    Interstitial%adjvisdfu       = Model%clear_val
    Interstitial%adjvisdfd       = Model%clear_val
    Interstitial%aerodp          = Model%clear_val
    Interstitial%alb1d           = Model%clear_val
    if (.not. Model%do_RRTMGP) then
      Interstitial%alpha         = Model%clear_val
    end if
    Interstitial%bexp1d          = Model%clear_val
    Interstitial%cd              = Model%clear_val
    Interstitial%cd_ice          = Model%huge
    Interstitial%cd_land         = Model%huge
    Interstitial%cd_water        = Model%huge
    Interstitial%cdq             = Model%clear_val
    Interstitial%cdq_ice         = Model%huge
    Interstitial%cdq_land        = Model%huge
    Interstitial%cdq_water       = Model%huge
    Interstitial%chh_ice         = Model%huge
    Interstitial%chh_land        = Model%huge
    Interstitial%chh_water       = Model%huge
    Interstitial%cldf            = Model%clear_val
    Interstitial%cldsa           = Model%clear_val
    Interstitial%cldtaulw        = Model%clear_val
    Interstitial%cldtausw        = Model%clear_val
    Interstitial%cld1d           = Model%clear_val
    Interstitial%clouds          = Model%clear_val
    Interstitial%clw             = Model%clear_val
    Interstitial%clw(:,:,2)      = -999.9
    Interstitial%clx             = Model%clear_val
    Interstitial%cmm_ice         = Model%huge
    Interstitial%cmm_land        = Model%huge
    Interstitial%cmm_water       = Model%huge
    Interstitial%cnvc            = Model%clear_val
    Interstitial%ctei_r          = Model%clear_val
    Interstitial%ctei_rml        = Model%clear_val
    Interstitial%cumabs          = Model%clear_val
    Interstitial%dd_mf           = Model%clear_val
    Interstitial%de_lgth         = Model%clear_val
    Interstitial%del             = Model%clear_val
    Interstitial%del_gz          = Model%clear_val
    Interstitial%delr            = Model%clear_val
    Interstitial%dlength         = Model%clear_val
    Interstitial%dqdt            = Model%clear_val
    Interstitial%dqsfc1          = Model%clear_val
    Interstitial%drain           = Model%clear_val
    Interstitial%dtdt            = Model%clear_val
    Interstitial%dtsfc1          = Model%clear_val
    Interstitial%dt_mf           = Model%clear_val
    Interstitial%dtzm            = Model%clear_val
    Interstitial%dudt            = Model%clear_val
    Interstitial%dusfcg          = Model%clear_val
    Interstitial%dusfc1          = Model%clear_val
    Interstitial%dvdt            = Model%clear_val
    Interstitial%dvsfcg          = Model%clear_val
    Interstitial%dvsfc1          = Model%clear_val
    Interstitial%dvdftra         = Model%clear_val
    Interstitial%dzlyr           = Model%clear_val
    Interstitial%elvmax          = Model%clear_val
    Interstitial%ep1d            = Model%clear_val
    Interstitial%ep1d_ice        = Model%huge
    Interstitial%ep1d_land       = Model%huge
    Interstitial%ep1d_water      = Model%huge
    Interstitial%evap_ice        = Model%huge
    Interstitial%evap_land       = Model%huge
    Interstitial%evap_water      = Model%huge
    Interstitial%evbs            = Model%clear_val
    Interstitial%evcw            = Model%clear_val
    Interstitial%pah             = Model%clear_val
    Interstitial%ecan            = Model%clear_val
    Interstitial%etran           = Model%clear_val
    Interstitial%edir            = Model%clear_val
    Interstitial%faerlw          = Model%clear_val
    Interstitial%faersw          = Model%clear_val
    Interstitial%ffhh_ice        = Model%huge
    Interstitial%ffhh_land       = Model%huge
    Interstitial%ffhh_water      = Model%huge
    Interstitial%fh2             = Model%clear_val
    Interstitial%fh2_ice         = Model%huge
    Interstitial%fh2_land        = Model%huge
    Interstitial%fh2_water       = Model%huge
    Interstitial%flag_cice       = .false.
    Interstitial%flag_guess      = .false.
    Interstitial%flag_iter       = .true.
    Interstitial%flag_lakefreeze = .false.
    Interstitial%ffmm_ice        = Model%huge
    Interstitial%ffmm_land       = Model%huge
    Interstitial%ffmm_water      = Model%huge
    Interstitial%fm10            = Model%clear_val
    Interstitial%fm10_ice        = Model%huge
    Interstitial%fm10_land       = Model%huge
    Interstitial%fm10_water      = Model%huge
    Interstitial%frland          = Model%clear_val
    Interstitial%fscav           = Model%clear_val
    Interstitial%fswtr           = Model%clear_val
    Interstitial%gabsbdlw        = Model%clear_val
    Interstitial%gabsbdlw_ice    = Model%clear_val
    Interstitial%gabsbdlw_land   = Model%clear_val
    Interstitial%gabsbdlw_water  = Model%clear_val
    Interstitial%gamma           = Model%clear_val
    Interstitial%gamq            = Model%clear_val
    Interstitial%gamt            = Model%clear_val
    Interstitial%gasvmr          = Model%clear_val
    Interstitial%gflx            = Model%clear_val
    Interstitial%gflx_ice        = Model%clear_val
    Interstitial%gflx_land       = Model%clear_val
    Interstitial%gflx_water      = Model%clear_val
    Interstitial%gwdcu           = Model%clear_val
    Interstitial%gwdcv           = Model%clear_val
    Interstitial%zvfun           = Model%clear_val
    Interstitial%hffac           = Model%clear_val
    Interstitial%hflxq           = Model%clear_val
    Interstitial%hflx_ice        = Model%huge
    Interstitial%hflx_land       = Model%huge
    Interstitial%hflx_water      = Model%huge
    Interstitial%htlwc           = Model%clear_val
    Interstitial%htlw0           = Model%clear_val
    Interstitial%htswc           = Model%clear_val
    Interstitial%htsw0           = Model%clear_val
    Interstitial%dry             = .false.
    Interstitial%idxday          = 0
    Interstitial%icy             = .false.
    Interstitial%lake            = .false.
    Interstitial%lndp_vgf        = Model%clear_val
    Interstitial%ocean           = .false.
    Interstitial%islmsk          = 0
    Interstitial%islmsk_cice     = 0
    Interstitial%wet             = .false.
    Interstitial%kb              = 0
    Interstitial%kbot            = Model%levs
    Interstitial%kcnv            = 0
    Interstitial%kd              = 0
    Interstitial%kinver          = Model%levs
    Interstitial%kpbl            = 0
    Interstitial%kt              = 0
    Interstitial%ktop            = 1
    Interstitial%mbota           = 0
    Interstitial%mtopa           = 0
    Interstitial%nday            = 0
    Interstitial%oa4             = Model%clear_val
    Interstitial%oc              = Model%clear_val
    Interstitial%olyr            = Model%clear_val
    Interstitial%plvl            = Model%clear_val
    Interstitial%plyr            = Model%clear_val
    Interstitial%prnum           = Model%clear_val
    Interstitial%qlyr            = Model%clear_val
    Interstitial%prcpmp          = Model%clear_val
    Interstitial%qss_ice         = Model%huge
    Interstitial%qss_land        = Model%huge
    Interstitial%qss_water       = Model%huge
    Interstitial%raddt           = Model%clear_val
    Interstitial%raincd          = Model%clear_val
    Interstitial%raincs          = Model%clear_val
    Interstitial%rainmcadj       = Model%clear_val
    Interstitial%rainp           = Model%clear_val
    Interstitial%rb              = Model%clear_val
    Interstitial%rb_ice          = Model%huge
    Interstitial%rb_land         = Model%huge
    Interstitial%rb_water        = Model%huge
    Interstitial%rhc             = Model%clear_val
    Interstitial%runoff          = Model%clear_val
    Interstitial%save_q          = Model%clear_val
    Interstitial%save_t          = Model%clear_val
    Interstitial%save_tcp        = Model%clear_val
    Interstitial%save_u          = Model%clear_val
    Interstitial%save_v          = Model%clear_val
    Interstitial%sbsno           = Model%clear_val
    Interstitial%scmpsw%uvbfc    = Model%clear_val
    Interstitial%scmpsw%uvbf0    = Model%clear_val
    Interstitial%scmpsw%nirbm    = Model%clear_val
    Interstitial%scmpsw%nirdf    = Model%clear_val
    Interstitial%scmpsw%visbm    = Model%clear_val
    Interstitial%scmpsw%visdf    = Model%clear_val
    Interstitial%sfcalb          = Model%clear_val
    Interstitial%sigma           = Model%clear_val
    Interstitial%sigmaf          = Model%clear_val
    Interstitial%sigmafrac       = Model%clear_val
    Interstitial%sigmatot        = Model%clear_val
    Interstitial%snowc           = Model%clear_val
    Interstitial%snohf           = Model%clear_val
    Interstitial%snowmt          = Model%clear_val
    Interstitial%stress          = Model%clear_val
    Interstitial%stress_ice      = Model%huge
    Interstitial%stress_land     = Model%huge
    Interstitial%stress_water    = Model%huge
    Interstitial%theta           = Model%clear_val
    Interstitial%tkeh            = 0
    Interstitial%tlvl            = Model%clear_val
    Interstitial%tlyr            = Model%clear_val
    Interstitial%tprcp_ice       = Model%huge
    Interstitial%tprcp_land      = Model%huge
    Interstitial%tprcp_water     = Model%huge
    Interstitial%trans           = Model%clear_val
    Interstitial%tseal           = Model%clear_val
    Interstitial%tsfa            = Model%clear_val
    Interstitial%tsfc_water      = Model%huge
    Interstitial%tsfg            = Model%clear_val
    Interstitial%tsurf_ice       = Model%huge
    Interstitial%tsurf_land      = Model%huge
    Interstitial%tsurf_water     = Model%huge
    Interstitial%ud_mf           = Model%clear_val
    Interstitial%uustar_ice      = Model%huge
    Interstitial%uustar_land     = Model%huge
    Interstitial%uustar_water    = Model%huge
    Interstitial%vdftra          = Model%clear_val
    Interstitial%vegf1d          = Model%clear_val
    Interstitial%wcbmax          = Model%clear_val
    Interstitial%wind            = Model%huge
    Interstitial%work1           = Model%clear_val
    Interstitial%work2           = Model%clear_val
    Interstitial%work3           = Model%clear_val
    Interstitial%xcosz           = Model%clear_val
    Interstitial%xlai1d          = Model%clear_val
    Interstitial%xmu             = Model%clear_val
    Interstitial%z01d            = Model%clear_val
    Interstitial%zt1d            = Model%clear_val
    Interstitial%ztmax_ice       = Model%clear_val
    Interstitial%ztmax_land      = Model%clear_val
    Interstitial%ztmax_water     = Model%clear_val

    ! RRTMGP
    if (Model%do_RRTMGP) then
       Interstitial%tracer                      = Model%clear_val
       Interstitial%tv_lay                      = Model%clear_val
       Interstitial%relhum                      = Model%clear_val
       Interstitial%qs_lay                      = Model%clear_val
       Interstitial%q_lay                       = Model%clear_val
       Interstitial%deltaZ                      = Model%clear_val
       Interstitial%deltaZc                     = Model%clear_val
       Interstitial%deltaP                      = Model%clear_val
       Interstitial%p_lev                       = Model%clear_val
       Interstitial%p_lay                       = Model%clear_val
       Interstitial%t_lev                       = Model%clear_val
       Interstitial%t_lay                       = Model%clear_val
       Interstitial%cloud_overlap_param         = Model%clear_val
       Interstitial%precip_overlap_param        = Model%clear_val
       Interstitial%fluxlwUP_allsky             = Model%clear_val
       Interstitial%fluxlwDOWN_allsky           = Model%clear_val
       Interstitial%fluxlwUP_clrsky             = Model%clear_val
       Interstitial%fluxlwDOWN_clrsky           = Model%clear_val
       Interstitial%fluxswUP_allsky             = Model%clear_val
       Interstitial%fluxswDOWN_allsky           = Model%clear_val
       Interstitial%fluxswUP_clrsky             = Model%clear_val
       Interstitial%fluxswDOWN_clrsky           = Model%clear_val
       Interstitial%aerosolslw                  = Model%clear_val
       Interstitial%aerosolssw                  = Model%clear_val
       Interstitial%precip_frac                 = Model%clear_val
       Interstitial%cld_cnv_frac                = Model%clear_val
       Interstitial%cnv_cloud_overlap_param     = Model%clear_val
       Interstitial%cld_cnv_lwp                 = Model%clear_val
       Interstitial%cld_cnv_reliq               = Model%clear_val
       Interstitial%cld_cnv_iwp                 = Model%clear_val
       Interstitial%cld_cnv_reice               = Model%clear_val
       Interstitial%cld_pbl_lwp                 = Model%clear_val
       Interstitial%cld_pbl_reliq               = Model%clear_val
       Interstitial%cld_pbl_iwp                 = Model%clear_val
       Interstitial%cld_pbl_reice               = Model%clear_val
       Interstitial%flxprf_lw%upfxc             = Model%clear_val
       Interstitial%flxprf_lw%dnfxc             = Model%clear_val
       Interstitial%flxprf_lw%upfx0             = Model%clear_val
       Interstitial%flxprf_lw%dnfx0             = Model%clear_val
       Interstitial%flxprf_sw%upfxc             = Model%clear_val
       Interstitial%flxprf_sw%dnfxc             = Model%clear_val
       Interstitial%flxprf_sw%upfx0             = Model%clear_val
       Interstitial%flxprf_sw%dnfx0             = Model%clear_val
       Interstitial%sfc_emiss_byband            = Model%clear_val
       Interstitial%sec_diff_byband             = Model%clear_val
       Interstitial%sfc_alb_nir_dir             = Model%clear_val
       Interstitial%sfc_alb_nir_dif             = Model%clear_val
       Interstitial%sfc_alb_uvvis_dir           = Model%clear_val
       Interstitial%sfc_alb_uvvis_dif           = Model%clear_val
       Interstitial%toa_src_sw                  = Model%clear_val
       Interstitial%toa_src_lw                  = Model%clear_val
       Interstitial%vmr_o2                      = Model%clear_val
       Interstitial%vmr_h2o                     = Model%clear_val
       Interstitial%vmr_o3                      = Model%clear_val
       Interstitial%vmr_ch4                     = Model%clear_val
       Interstitial%vmr_n2o                     = Model%clear_val
       Interstitial%vmr_co2                     = Model%clear_val
    end if

    ! UGWP common
    Interstitial%tau_mtb         = Model%clear_val
    Interstitial%tau_ogw         = Model%clear_val
    Interstitial%tau_tofd        = Model%clear_val
    Interstitial%tau_ngw         = Model%clear_val
    Interstitial%tau_oss         = Model%clear_val
    Interstitial%dudt_mtb        = Model%clear_val
    Interstitial%dudt_tms        = Model%clear_val
    Interstitial%zmtb            = Model%clear_val
    Interstitial%zlwb            = Model%clear_val
    Interstitial%zogw            = Model%clear_val
    Interstitial%zngw            = Model%clear_val

! CIRES UGWP v1
    if (Model%ldiag_ugwp .or. Model%do_ugwp_v0 .or. Model%do_ugwp_v0_nst_only &
        .or. Model%do_ugwp_v1) then
      Interstitial%dudt_ngw      = Model%clear_val
      Interstitial%dvdt_ngw      = Model%clear_val
      Interstitial%dtdt_ngw      = Model%clear_val
      Interstitial%kdis_ngw      = Model%clear_val
    end if

!-- GSL drag suite
    if (Model%gwd_opt==3 .or. Model%gwd_opt==33 .or. &
        Model%gwd_opt==2 .or. Model%gwd_opt==22 ) then
       Interstitial%varss        = Model%clear_val
       Interstitial%ocss         = Model%clear_val
       Interstitial%oa4ss        = Model%clear_val
       Interstitial%clxss        = Model%clear_val
    end if
!
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson &
         .or. Model%imp_physics == Model%imp_physics_tempo .or. Model%imp_physics == Model%imp_physics_nssl &
        ) then
       Interstitial%graupelmp    = Model%clear_val
       Interstitial%icemp        = Model%clear_val
       Interstitial%rainmp       = Model%clear_val
       Interstitial%snowmp       = Model%clear_val
    else if (Model%imp_physics == Model%imp_physics_mg) then
       Interstitial%ncgl         = Model%clear_val
       Interstitial%ncpr         = Model%clear_val
       Interstitial%ncps         = Model%clear_val
       Interstitial%qgl          = Model%clear_val
       Interstitial%qrn          = Model%clear_val
       Interstitial%qsnw         = Model%clear_val
       Interstitial%qlcn         = Model%clear_val
       Interstitial%qicn         = Model%clear_val
       Interstitial%w_upi        = Model%clear_val
       Interstitial%cf_upi       = Model%clear_val
       Interstitial%cnv_mfd      = Model%clear_val
       Interstitial%cnv_dqldt    = Model%clear_val
       Interstitial%clcn         = Model%clear_val
       Interstitial%cnv_fice     = Model%clear_val
       Interstitial%cnv_ndrop    = Model%clear_val
       Interstitial%cnv_nice     = Model%clear_val
    end if
    if (Model%lsm == Model%ilsm_noahmp) then
       Interstitial%t2mmp        = Model%clear_val
       Interstitial%q2mp         = Model%clear_val
    end if

    ! DJS - Move max_hourly_reset, ext_diag_thompson_reset, and fullradar_diag to GFS_control_type?
    ! Set flag for resetting maximum hourly output fields
    Interstitial%max_hourly_reset = mod(Model%kdt-1, nint(Model%avg_max_length/Model%dtp)) == 0
    ! Use same logic in UFS to reset Thompson extended diagnostics
    Interstitial%ext_diag_thompson_reset = Interstitial%max_hourly_reset

    ! Frequency flag for computing the full radar reflectivity (water coated ice) 
    if (Model%nsfullradar_diag<0) then
      Interstitial%fullradar_diag = .true.
    else
      Interstitial%fullradar_diag = (Model%kdt == 1 .or. mod(Model%kdt, nint(Model%nsfullradar_diag/Model%dtp)) == 0)
    end if
    !

    !
    ! CCPP suite simulator
    if (Model%do_ccpp_suite_sim) then
       Interstitial%active_phys_tend = Model%clear_val
    endif

  end subroutine gfs_interstitial_reset

end module CCPP_typedefs
