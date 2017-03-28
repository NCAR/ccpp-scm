!> \file gmtb_scm_output.f90
!!  Contains output-related subroutines

module gmtb_scm_output

use netcdf
use gmtb_scm_input, only: check
use gmtb_scm_kinds, only: sp, dp, qp

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup output gmtb_scm_output
!! @{
!! Contains output-related subroutines

!> This subroutine initializes the output netCDF file, "output.nc", placed in the directory specified by the case_config file used.
subroutine output_init(output_dir, output_file, n_columns, n_levels, year, month, day, hour)
  character(len=*), intent(in)      :: output_dir !< name of the output directory
  character(len=*), intent(in)      :: output_file !< name of the output file
  integer, intent(in)               :: n_columns !< number of model columns
  integer, intent(in)               :: n_levels !< number of model levels
  integer, intent(in)               :: year, month, day, hour !<initial date


  INTEGER :: ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_i_id, dummy_id, year_id, month_id, day_id, hour_id
  !> \section output_init_alg Algorithm
  !! @{
  !> - Create the output directory if necessary using a system call.
  CALL SYSTEM('mkdir -p '//TRIM(output_dir))

  !> - Create the output file in the output directory.
  CALL CHECK(NF90_CREATE(PATH=TRIM(output_dir)//"/"//TRIM(output_file)//".nc",CMODE=NF90_CLOBBER,NCID=ncid))

  !> - Define netCDF dimensions.
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_dim",LEN=NF90_UNLIMITED,DIMID=time_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="hor_dim_layer",LEN=n_columns,DIMID=hor_dim_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_layer",LEN=n_levels,DIMID=vert_dim_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_interface",LEN=n_levels+1,DIMID=vert_dim_i_id))

  !> - Define the dimension variables.
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME="time",XTYPE=NF90_FLOAT,DIMIDS=(/ time_id /), VARID=dummy_id))

  !> - Define the pressure-related variables.
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='pres',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='pres_i',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_i_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sigma',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sigma_i',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_i_id, time_id /), &
    VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='phi',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='phi_i',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_i_id, time_id /), &
  !   VARID=dummy_id))

  !> - Define the state variables.
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='qv',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='T',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='u',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='v',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='qc',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))

  !> - Define the forcing-related variables.
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='qv_force_tend',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='T_force_tend',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='u_force_tend',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='v_force_tend',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='w_ls',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='u_g',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='v_g',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_rad_forc',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='h_advec_thil',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='h_advec_qt',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='v_advec_thil',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='v_advec_qt',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
    VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='T_s',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='pres_s',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lhf',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='shf',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='tau_u',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='tau_v',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))

  !> - Define the diagnostics variables.
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='cldcov',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_rad_heating_rate',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_rad_heating_rate',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='precip',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='rain',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='rainc',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='pwat',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_lwrad',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_swrad',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /),&
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_PBL',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_deepconv',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_shalconv',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dT_dt_micro',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dq_dt_PBL',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dq_dt_deepconv',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dq_dt_shalconv',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id, time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dq_dt_micro',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='upd_mf',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='dwn_mf',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='det_mf',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  !   VARID=dummy_id))
  ! ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='forcet',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  ! !   VARID=dummy_id))
  ! ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='forceq',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, vert_dim_id,time_id /), &
  ! !   VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='PBL_height',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_up_TOA_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_dn_TOA_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_up_TOA_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_up_sfc_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_dn_sfc_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_up_sfc_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='sw_dn_sfc_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_up_TOA_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_up_TOA_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_up_sfc_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_up_sfc_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_dn_sfc_tot',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))
  ! CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='lw_dn_sfc_clr',XTYPE=NF90_FLOAT,DIMIDS= (/ hor_dim_id, time_id /), VARID=dummy_id))

  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='init_year',XTYPE=NF90_FLOAT,VARID=year_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='init_month',XTYPE=NF90_FLOAT,VARID=month_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='init_day',XTYPE=NF90_FLOAT,VARID=day_id))
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME='init_hour',XTYPE=NF90_FLOAT,VARID=hour_id))


  !> - Close variable definition and the file.
  CALL CHECK(NF90_ENDDEF(NCID=ncid))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=year_id,VALUES=year))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=month_id,VALUES=month))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=day_id,VALUES=day))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=hour_id,VALUES=hour))
  CALL CHECK(NF90_CLOSE(NCID=ncid))

  !> @}
end subroutine output_init

!> This subroutine appends data to the "output.nc" file.
subroutine output_append(output_dir, output_file, record_number, time, pres_l, pres_i, sigma_l, sigma_i, qv, T, u, v,&
  qc, &
  u_force_tend, v_force_tend, T_force_tend, qv_force_tend, w_ls, u_g, v_g, dT_dt_rad_forc, h_advec_thil, h_advec_qt, v_advec_thil, &
  v_advec_qt, T_s, pres_s)
  character(len=*), intent(in)      :: output_dir !< name of the output directory
  character(len=*), intent(in)      :: output_file !< name of the output file
  integer, intent(in)               :: record_number !< the record counter
  real(kind=dp), intent(in)          :: time !< elapsed model time
  real(kind=dp), intent(in)          :: pres_l(:,:) !< pressure on model layer centers (Pa) (horizontal, vertical)
  real(kind=dp), intent(in)          :: pres_i(:,:) !< pressure on model layer interfaces (Pa) (horizontal, vertical)
  real(kind=dp), intent(in)          :: sigma_l(:,:) !< sigma (p/p_s) on model layer centers (Pa) (horizontal, vertical)
  real(kind=dp), intent(in)          :: sigma_i(:,:) !< sigma (p/p_s) on model layer interfaces (Pa) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: phi_l(:,:) !< geopotential on model layer centers (Pa) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: phi_i(:,:) !< geopotential on model layer interfaces (Pa) (horizontal, vertical)
  real(kind=dp), intent(in)          :: qv(:,:) !< water vapor specific humidity on model layer centers (kg/kg) (horizontal, vertical)
  real(kind=dp), intent(in)          :: T(:,:) !< absolute temperature on model layer centers (K) (horizontal, vertical)
  real(kind=dp), intent(in)          :: u(:,:) !< zonal wind on model layer centers (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: v(:,:) !< meridional wind on model layer centers (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: qc(:,:) !< suspended cloud water on model layer centers (kg/kg) (horizontal, vertical)
  real(kind=dp), intent(in)          :: u_force_tend(:,:) !< total forcing tendency for zonal momentum (m/s /s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: v_force_tend(:,:) !< total forcing tendency for meridional momentum (m/s /s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: T_force_tend(:,:) !< total forcing tendency for absolute temperature (K /s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: qv_force_tend(:,:) !< total forcing tendency for water vapor specific humidity (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: w_ls(:,:) !< large-scale vertical velocity (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: u_g(:,:) !< geostrophic zonal wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: v_g(:,:) !< geostrophic meridional wind (m/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: dT_dt_rad_forc(:,:) !< prescribed radiative heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: h_advec_thil(:,:) !< large-scale horizontal advective heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: h_advec_qt(:,:) !< large-scale horizontal advective moistening rate (kg/kg /s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: v_advec_thil(:,:) !< large-scale vertical advective heating rate (K/s) (horizontal, vertical)
  real(kind=dp), intent(in)          :: v_advec_qt(:,:) !< large-scale vertical advective moistening rate (kg/kg /s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: cldcov(:,:) !< cloud fraction (horizontal, vertical)
  real(kind=dp), intent(in)          :: T_s(:) !< surface temperature (K) (hoizontal)
  real(kind=dp), intent(in)          :: pres_s(:) !< surface pressure (Pa) (horizontal)
  ! real(kind=dp), intent(in)          :: lhf(:) !< surface latent heat flux (W/m^2) (horizontal)
  ! real(kind=dp), intent(in)          :: shf(:) !< surface sensible heat flux (W/m^2) (horizontal)
  ! real(kind=dp), intent(in)          :: tau_u(:) !< surface zonal wind stress (horizontal)
  ! real(kind=dp), intent(in)          :: tau_v(:) !< surface meridional wind stress (horizontal)
  ! real(kind=dp), intent(in)          :: swh(:,:) !< calculated SW radiative heating rate (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: hlw(:,:) !< calculated LW radiative heating rate (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: tprcp(:) !< surface precipitation (horizontal)
  ! real(kind=dp), intent(in)          :: rain(:) !< surface total rain (m/s) (horizontal)
  ! real(kind=dp), intent(in)          :: rainc(:) !< surface convective rain (m/s) (horizontal)
  ! real(kind=dp), intent(in)          :: pwat(:) !< surface rain (horizontal)
  ! real(kind=dp), intent(in)          :: dT_dt_lwrad(:,:) !< change in temperature due to lw rad scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dT_dt_swrad(:,:) !< change in temperature due to sw rad scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dT_dt_PBL(:,:) !< change in temperature due to PBL scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dT_dt_deepconv(:,:) !< change in temperature due to deep conv scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dT_dt_shalconv(:,:) !< change in temperature due to shal conv scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dT_dt_micro(:,:) !< change in temperature due to micro scheme (K/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dq_dt_PBL(:,:) !< change in moisture due to PBL scheme (kg/kg/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dq_dt_deepconv(:,:) !< change in moisture due to deep conv scheme (kg/kg/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dq_dt_shalconv(:,:) !< change in moisture due to shal conv scheme (kg/kg/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dq_dt_micro(:,:) !< change in moisture due to micro scheme (kg/kg/s) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: upd_mf(:,:) !< updraft mass flux (kg m^-2 s^-1) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: dwn_mf(:,:) !< downdraft mass flux (kg m^-2 s^-1) (horizontal, vertical)
  ! real(kind=dp), intent(in)          :: det_mf(:,:) !< detrainment mass flux (kg m^-2 s^-1) (horizontal, vertical)
  ! !real(kind=dp), intent(in)          :: forcet(:,:)
  ! !real(kind=dp), intent(in)          :: forceq(:,:)
  ! real(kind=dp), intent(in)          :: hpbl(:) !< PBL height (m) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_up_TOA_tot(:) !< total sky upward sw flux at toa (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_dn_TOA_tot(:) !< total sky downward sw flux at toa (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_up_TOA_clr(:) !< clear sky upward sw flux at toa (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_up_sfc_tot(:) !< total sky upward sw flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_dn_sfc_tot(:) !< total sky downward sw flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_up_sfc_clr(:) !< clear sky upward sw flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: sw_dn_sfc_clr(:) !< clear sky downward sw flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_up_TOA_tot(:) !< total sky upward LW flux at toa (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_up_TOA_clr(:) !< clear sky upward LW flux at toa (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_up_sfc_tot(:) !< total sky upward LW flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_up_sfc_clr(:) !< clear sky upward LW flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_dn_sfc_tot(:) !< total sky downward LW flux at sfc (\f$W/m^2\f$) (horizontal)
  ! real(kind=dp), intent(in)          :: lw_dn_sfc_clr(:) !< clear sky downward LW flux at sfc (\f$W/m^2\f$) (horizontal)


  integer :: ncid, var_id

  !> \section output_append_alg Algorithm
  !! @{

  !> - Open the file.
  CALL CHECK(NF90_OPEN(PATH=TRIM(output_dir)//"/"//TRIM(output_file)//".nc",MODE=NF90_WRITE,NCID=ncid))

  !> - Append all of the variables to the file.
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=time,START=(/ record_number /)))

  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="pres",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=pres_l(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="pres_i",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=pres_i(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sigma",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sigma_l(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sigma_i",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sigma_i(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="phi",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=phi_l(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="phi_i",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=phi_i(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="qv",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=qv(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="T",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=T(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="u",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=u(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="v",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=v(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="qc",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=qc(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="qv_force_tend",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=qv_force_tend(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="T_force_tend",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=T_force_tend(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="u_force_tend",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=u_force_tend(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="v_force_tend",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=v_force_tend(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="w_ls",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=w_ls(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="u_g",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=u_g(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="v_g",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=v_g(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_rad_forc",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_rad_forc(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="h_advec_thil",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=h_advec_thil(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="h_advec_qt",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=h_advec_qt(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="v_advec_thil",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=v_advec_thil(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="v_advec_qt",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=v_advec_qt(:,:),START=(/1,1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="T_s",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=T_s(:),START=(/1,record_number /)))
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="pres_s",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=pres_s(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lhf",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lhf(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="shf",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=shf(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="tau_u",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=tau_u(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="tau_v",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=tau_v(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="cldcov",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=cldcov(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_rad_heating_rate",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=swh(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_rad_heating_rate",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=hlw(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="precip",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=tprcp(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="rain",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=rain(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="rainc",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=rainc(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="pwat",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=pwat(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_lwrad",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_lwrad(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_swrad",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_swrad(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_PBL",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_PBL(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_deepconv",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_deepconv(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_shalconv",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_shalconv(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dT_dt_micro",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dT_dt_micro(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dq_dt_PBL",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dq_dt_PBL(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dq_dt_deepconv",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dq_dt_deepconv(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dq_dt_shalconv",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dq_dt_shalconv(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dq_dt_micro",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dq_dt_micro(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="upd_mf",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=upd_mf(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="dwn_mf",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=dwn_mf(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="det_mf",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=det_mf(:,:),START=(/1,1,record_number /)))
  ! ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="forcet",VARID=var_id))
  ! ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=forcet(:,:),START=(/1,1,record_number /)))
  ! ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="forceq",VARID=var_id))
  ! ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=forceq(:,:),START=(/1,1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="PBL_height",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=hpbl(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_up_TOA_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_up_TOA_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_dn_TOA_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_dn_TOA_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_up_TOA_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_up_TOA_clr(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_up_sfc_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_up_sfc_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_dn_sfc_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_dn_sfc_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_up_sfc_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_up_sfc_clr(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="sw_dn_sfc_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=sw_dn_sfc_clr(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_up_TOA_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_up_TOA_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_up_TOA_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_up_TOA_clr(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_up_sfc_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_up_sfc_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_up_sfc_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_up_sfc_clr(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_dn_sfc_tot",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_dn_sfc_tot(:),START=(/1,record_number /)))
  ! CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="lw_dn_sfc_clr",VARID=var_id))
  ! CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=lw_dn_sfc_clr(:),START=(/1,record_number /)))

  !> - Close the file.
  CALL CHECK(NF90_CLOSE(ncid))

  !> @}
end subroutine output_append


!> @}
!> @}
end module gmtb_scm_output
