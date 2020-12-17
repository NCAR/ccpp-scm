!> \file gmtb_scm_output.f90
!!  Contains output-related subroutines

module gmtb_scm_output

use netcdf
use NetCDF_read, only: check
use gmtb_scm_kinds, only: sp, dp, qp

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup output gmtb_scm_output
!! @{
!! Contains output-related subroutines

!> This subroutine initializes the output netCDF file, "output.nc", placed in the directory specified by the case_config file used.
subroutine output_init(scm_state, physics)
  use gmtb_scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_def, only: NetCDF_def_var
  use NetCDF_put, only: NetCDF_put_var

  type(scm_state_type), intent(in) :: scm_state
  type(physics_type),   intent(in) :: physics

  INTEGER :: i, ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_i_id, vert_dim_rad_id, dummy_id, year_id, month_id, day_id, hour_id
  character(2) :: idx
  !> \section output_init_alg Algorithm
  !! @{

  !> - Create the output file in the output directory.
  CALL CHECK(NF90_CREATE(PATH=TRIM(scm_state%output_dir)//"/"//TRIM(scm_state%output_file)//".nc",CMODE=NF90_CLOBBER,NCID=ncid))

  !> - Define netCDF dimensions.
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_dim",LEN=NF90_UNLIMITED,DIMID=time_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="hor_dim_layer",LEN=scm_state%n_cols,DIMID=hor_dim_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_layer",LEN=scm_state%n_levels,DIMID=vert_dim_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_interface",LEN=scm_state%n_levels+1,DIMID=vert_dim_i_id))
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_rad",LEN=physics%Interstitial%lmk,DIMID=vert_dim_rad_id))

  !> - Define the dimension variables.
  CALL CHECK(NF90_DEF_VAR(NCID=ncid,NAME="time",XTYPE=NF90_FLOAT,DIMIDS=(/ time_id /), VARID=dummy_id))
  CALL CHECK(NF90_PUT_ATT(NCID=ncid,VARID=dummy_id,NAME="description",VALUES="model elapsed time"))
  CALL CHECK(NF90_PUT_ATT(NCID=ncid,VARID=dummy_id,NAME="units",VALUES="s"))

  !> - Define the state variables
  CALL output_init_state(ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_i_id)
  !> - Define the forcing variables
  CALL output_init_forcing(ncid, time_id, hor_dim_id, vert_dim_id)
  
  !> - Define all diagnostic/physics variables
  CALL output_init_sfcprop(ncid, time_id, hor_dim_id, scm_state, physics)
  CALL output_init_interstitial(ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_rad_id)
  CALL output_init_radtend(ncid, time_id, hor_dim_id, vert_dim_id)
  CALL output_init_diag(ncid, time_id, hor_dim_id, vert_dim_id, physics)
  
  call NetCDF_def_var(ncid, 'init_year',  NF90_FLOAT, "model initialization year",  "year",  year_id)
  call NetCDF_def_var(ncid, 'init_month', NF90_FLOAT, "model initialization month", "month", month_id)
  call NetCDF_def_var(ncid, 'init_day',   NF90_FLOAT, "model initialization day",   "day",   day_id)
  call NetCDF_def_var(ncid, 'init_hour',  NF90_FLOAT, "model initialization hour",  "hour",  hour_id)
  
  !> - Close variable definition and the file.
  CALL CHECK(NF90_ENDDEF(NCID=ncid))
  
  call NetCDF_put_var(ncid, "init_year",  scm_state%init_year,  year_id)
  call NetCDF_put_var(ncid, "init_month", scm_state%init_month, month_id)
  call NetCDF_put_var(ncid, "init_day",   scm_state%init_day,   day_id)
  call NetCDF_put_var(ncid, "init_hour",  scm_state%init_hour,  hour_id)

  CALL CHECK(NF90_CLOSE(NCID=ncid))

  !> @}
end subroutine output_init

subroutine output_init_state(ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_i_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_i_id
  
  integer :: dummy_id
  
  !> - Define the pressure-related variables.
  call NetCDF_def_var(ncid, 'pres',    NF90_FLOAT, "pressure on model layer centers",         "Pa",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'pres_i',  NF90_FLOAT, "pressure on model layer interfaces",      "Pa",   dummy_id, (/ hor_dim_id, vert_dim_i_id, time_id /))
  call NetCDF_def_var(ncid, 'sigma',   NF90_FLOAT, "sigma (p/p_s) on model layer centers",    "none", dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'sigma_i', NF90_FLOAT, "sigma (p/p_s) on model layer interfaces", "none", dummy_id, (/ hor_dim_id, vert_dim_i_id, time_id /))
  call NetCDF_def_var(ncid, 'pres_s',  NF90_FLOAT, "surface pressure",                        "Pa",   dummy_id, (/ hor_dim_id,                time_id /))
  
  !> - Define the state variables.
  call NetCDF_def_var(ncid, 'qv', NF90_FLOAT, "water vapor specific humidity on model layer centers",                "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'T',  NF90_FLOAT, "absolute temperature on model layer centers",                         "K",       dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'u',  NF90_FLOAT, "x-wind on model layer centers",                                       "m s-1",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'v',  NF90_FLOAT, "y-wind on model layer centers",                                       "m s-1",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'ql', NF90_FLOAT, "suspended resolved liquid cloud water on model layer centers",        "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'qi', NF90_FLOAT, "suspended resolved ice cloud water on model layer centers",           "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  call NetCDF_def_var(ncid, 'qc', NF90_FLOAT, "suspended (resolved + SGS) total cloud water on model layer centers", "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_id /))
  
end subroutine output_init_state

subroutine output_init_forcing(ncid, time_id, hor_dim_id, vert_dim_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id, vert_dim_id
  
  integer :: dummy_id
  
  !> - Define the forcing-related variables.
  call NetCDF_def_var(ncid, 'qv_force_tend',  NF90_FLOAT, "total forcing tendency for water vapor specific humidity",                                 "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'T_force_tend',   NF90_FLOAT, "total forcing tendency for absolute temperature",                                          "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'u_force_tend',   NF90_FLOAT, "total forcing tendency for x-wind",                                                        "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'v_force_tend',   NF90_FLOAT, "total forcing tendency for y-wind",                                                        "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'w_ls',           NF90_FLOAT, "large-scale vertical velocity forcing",                                                    "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'u_g',            NF90_FLOAT, "geostrophic x-wind forcing",                                                               "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'v_g',            NF90_FLOAT, "geostrophic y-wind forcing",                                                               "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_rad_forc', NF90_FLOAT, "prescribed radiative heating rate forcing",                                                "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'h_advec_thil',   NF90_FLOAT, "large-scale horizontal advective heating rate for ice-liquid water potential temperature", "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'h_advec_qt',     NF90_FLOAT, "large-scale horizontal advective moistening rate for total water specific humidity",       "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'v_advec_thil',   NF90_FLOAT, "large-scale vertical advective heating rate for ice-liquid water potential temperature",   "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'v_advec_qt',     NF90_FLOAT, "large-scale vertical advective moistening rate for total water specific humidity",         "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'T_s',            NF90_FLOAT, "surface temperature forcing",                                                              "K",           dummy_id, (/ hor_dim_id,              time_id /))

end subroutine output_init_forcing

subroutine output_init_sfcprop(ncid, time_id, hor_dim_id, scm_state, physics)
  use gmtb_scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id
  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics
  
  character(2) :: idx
  integer :: i, dummy_id
  
  if (scm_state%model_ics .or. scm_state%lsm_ics) then
    if (physics%Model%lsoil > 0) then
      do i=1, physics%Model%lsoil
        write(idx,'(I2)') i
        call NetCDF_def_var(ncid, 'T_soil_'//adjustl(trim(idx)),  NF90_FLOAT, "soil temperature for layer "//idx, "K", dummy_id, (/ hor_dim_id, time_id /))
        call NetCDF_def_var(ncid, 'soil_moisture_'//adjustl(trim(idx)),  NF90_FLOAT, "soil moisture for layer "//idx, "m3 m-3", dummy_id, (/ hor_dim_id, time_id /))
      end do
    end if
  end if
  
  call NetCDF_def_var(ncid, 'lhf', NF90_FLOAT, "surface latent heat flux", "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'shf', NF90_FLOAT, "surface sensible heat flux", "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'tprcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation", "m", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'tprcp_rate_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation rate (tprcp_inst/dt)", "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  
end subroutine output_init_sfcprop

subroutine output_init_interstitial(ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_rad_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id, vert_dim_id, vert_dim_rad_id
  
  integer :: dummy_id
  
  call NetCDF_def_var(ncid, 'tau_u',  NF90_FLOAT, "surface x-wind stress", "Pa",         dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'tau_v',  NF90_FLOAT, "surface y-wind stress", "Pa",         dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'upd_mf', NF90_FLOAT, "updraft mass flux",     "kg m-2 s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dwn_mf', NF90_FLOAT, "downdraft mass flux",   "kg m-2 s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'det_mf', NF90_FLOAT, "detrainment mass flux", "kg m-2 s-1", dummy_id, (/ hor_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'sfc_up_lw_land',     NF90_FLOAT, "surface upwelling longwave flux over land fraction (valid all timesteps)",                 "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_lw_ice',      NF90_FLOAT, "surface upwelling longwave flux over ice fraction (valid all timesteps)",                  "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_lw_ocean',    NF90_FLOAT, "surface upwelling longwave flux over ocean fraction (valid all timesteps)",                "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dir_nir',  NF90_FLOAT, "surface upwelling shortwave direct near-infrared flux (valid all timesteps)",              "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dif_nir',  NF90_FLOAT, "surface upwelling shortwave diffuse near-infrared flux (valid all timesteps)",             "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dir_vis',  NF90_FLOAT, "surface upwelling shortwave direct visible and ultraviolet flux (valid all timesteps)",    "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dif_vis',  NF90_FLOAT, "surface upwelling shortwave diffuse visible and ultraviolet flux (valid all timesteps)",   "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dir_nir', NF90_FLOAT, "surface downwelling shortwave direct near-infrared flux (valid all timesteps)",            "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dif_nir', NF90_FLOAT, "surface downwelling shortwave diffuse near-infrared flux (valid all timesteps)",           "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dir_vis', NF90_FLOAT, "surface downwelling shortwave direct visible and ultraviolet flux (valid all timesteps)",  "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dif_vis', NF90_FLOAT, "surface downwelling shortwave diffuse visible and ultraviolet flux (valid all timesteps)", "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'mp_prcp_inst',   NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from microphysics scheme",       "m", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dcnv_prcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from deep convection scheme",    "m", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'scnv_prcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from shallow convection scheme", "m", dummy_id, (/ hor_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'rad_cloud_fraction', NF90_FLOAT, "instantaneous cloud fraction used in radiation",                   "fraction", dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_lwp',      NF90_FLOAT, "instantaneous cloud liquid water path used in radiation",          "g m-2",    dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_ql',     NF90_FLOAT, "instantaneous effective radius for liquid cloud used in radiation", "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_iwp',      NF90_FLOAT, "instantaneous cloud ice water path used in radiation",              "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qi',     NF90_FLOAT, "instantaneous effective radius for ice cloud used in radiation",    "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_rwp',      NF90_FLOAT, "instantaneous rain water path used in radiation",                   "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qr',     NF90_FLOAT, "instantaneous effective radius for raindrop used in radiation",     "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_swp',      NF90_FLOAT, "instantaneous snow water path used in radiation",                   "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qs',     NF90_FLOAT, "instantaneous effective radius for snowflake in radiation",         "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_id /))
  
end subroutine output_init_interstitial

subroutine output_init_radtend(ncid, time_id, hor_dim_id, vert_dim_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id, vert_dim_id
  
  integer :: dummy_id
  
  call NetCDF_def_var(ncid, 'sw_rad_heating_rate', NF90_FLOAT, "total sky shortwave radiative heating rate (radiation timesteps only)", "K s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'lw_rad_heating_rate', NF90_FLOAT, "total sky longwave radiative heating rate (radiation timesteps only)",  "K s-1", dummy_id, (/ hor_dim_id, time_id /))
  
end subroutine output_init_radtend

subroutine output_init_diag(ncid, time_id, hor_dim_id, vert_dim_id, physics)
  use gmtb_scm_type_defs, only: physics_type
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_id, hor_dim_id, vert_dim_id
  type(physics_type), intent(in) :: physics
  
  integer :: i, dummy_id
  character(2) :: idx
  
  call NetCDF_def_var(ncid, 'pwat',            NF90_FLOAT, "column precipitable water", "kg m-2", dummy_id, (/ hor_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'dT_dt_lwrad',     NF90_FLOAT, "temperature tendency due to longwave radiation scheme",           "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_swrad',     NF90_FLOAT, "temperature tendency due to shortwave radiation scheme",          "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_pbl',       NF90_FLOAT, "temperature tendency due to PBL scheme",                          "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_deepconv',  NF90_FLOAT, "temperature tendency due to deep convection scheme",              "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_shalconv',  NF90_FLOAT, "temperature tendency due to shallow convection scheme",           "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_micro',     NF90_FLOAT, "temperature tendency due to deep microphysics scheme",            "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_ogwd',      NF90_FLOAT, "temperature tendency due to orographic gravity wave drag scheme", "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_rayleigh',  NF90_FLOAT, "temperature tendency due to rayleigh damping scheme",             "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_cgwd',      NF90_FLOAT, "temperature tendency due to convective gravity wave drag scheme", "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_phys',      NF90_FLOAT, "temperature tendency due to all physics schemes",                 "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dT_dt_nonphys',   NF90_FLOAT, "temperature tendency due to all processes other than physics",    "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_pbl',       NF90_FLOAT, "moisture tendency due to PBL scheme",                             "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_deepconv',  NF90_FLOAT, "moisture tendency due to deep convection scheme",                 "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_shalconv',  NF90_FLOAT, "moisture tendency due to shallow convection scheme",              "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_micro',     NF90_FLOAT, "moisture tendency due to microphysics scheme",                    "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_phys',      NF90_FLOAT, "moisture tendency due to all physics schemes",                    "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dq_dt_nonphys',   NF90_FLOAT, "moisture tendency due to all processes other than physics",       "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_pbl',      NF90_FLOAT, "ozone tendency due to PBL scheme",                                "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_prodloss', NF90_FLOAT, "ozone tendency due to ozone production/loss",                     "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_oz',       NF90_FLOAT, "ozone tendency due to ozone",                                     "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_T',        NF90_FLOAT, "ozone tendency due to temperature",                               "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_ovhd',     NF90_FLOAT, "ozone tendency due to overhead ozone column",                     "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_phys',     NF90_FLOAT, "ozone tendency due to all physics schemes",                       "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'doz_dt_nonphys',  NF90_FLOAT, "ozone tendency due to all processes other than physics",          "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_pbl',       NF90_FLOAT, "x-wind tendency due to PBL scheme",                               "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_ogwd',      NF90_FLOAT, "x-wind tendency due to orographic GWD scheme",                    "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_deepconv',  NF90_FLOAT, "x-wind tendency due to deep convection scheme",                   "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_cgwd',      NF90_FLOAT, "x-wind tendency due to convective GWD scheme",                    "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_rayleigh',  NF90_FLOAT, "x-wind tendency due to rayleigh damping scheme",                  "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_shalconv',  NF90_FLOAT, "x-wind tendency due to shallow convection scheme",                "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_phys',      NF90_FLOAT, "x-wind tendency due to all physics schemes",                      "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'du_dt_nonphys',   NF90_FLOAT, "x-wind tendency due to all processes other than physics",         "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_pbl',       NF90_FLOAT, "y-wind tendency due to PBL scheme",                               "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_ogwd',      NF90_FLOAT, "y-wind tendency due to orographic GWD scheme",                    "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_deepconv',  NF90_FLOAT, "y-wind tendency due to deep convection scheme",                   "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_cgwd',      NF90_FLOAT, "y-wind tendency due to convective GWD scheme",                    "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_rayleigh',  NF90_FLOAT, "y-wind tendency due to rayleigh damping scheme",                  "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_shalconv',  NF90_FLOAT, "y-wind tendency due to shallow convection scheme",                "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_phys',      NF90_FLOAT, "y-wind tendency due to all physics schemes",                      "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'dv_dt_nonphys',   NF90_FLOAT, "y-wind tendency due to all processes other than physics",         "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'sfc_dwn_sw',      NF90_FLOAT, "surface downwelling shortwave flux (valid all timesteps)",                   "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw',       NF90_FLOAT, "surface upwelling shortwave flux (valid all timesteps)",                     "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_net_sw',      NF90_FLOAT, "surface net shortwave flux (downwelling - upwelling) (valid all timesteps)", "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_lw',      NF90_FLOAT, "surface downwelling longwave flux (valid all timesteps)",                    "W m-2", dummy_id, (/ hor_dim_id, time_id /))
  
  call NetCDF_def_var(ncid, 'tprcp_accum',          NF90_FLOAT, "cumulative surface liquid water equivalent thickness of total precipitation (valid over diagnostic interval)",           "m",     dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'ice_accum',            NF90_FLOAT, "cumulative surface liquid water equivalent thickness of ice precipitation (valid over diagnostic interval)",             "m",     dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'snow_accum',           NF90_FLOAT, "cumulative surface liquid water equivalent thickness of snow precipitation (valid over diagnostic interval)",            "m",     dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'graupel_accum',        NF90_FLOAT, "cumulative surface liquid water equivalent thickness of graupel precipitation (valid over diagnostic interval)",         "m",     dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'conv_prcp_accum',      NF90_FLOAT, "cumulative surface liquid water equivalent thickness of convective precipitation (valid over diagnostic interval)",      "m",     dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'tprcp_rate_accum',     NF90_FLOAT, "cumulative surface liquid water equivalent thickness of total precipitation rate (valid over diagnostic interval)",      "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'ice_rate_accum',       NF90_FLOAT, "cumulative surface liquid water equivalent thickness of ice precipitation rate (valid over diagnostic interval)",        "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'snow_rate_accum',      NF90_FLOAT, "cumulative surface liquid water equivalent thickness of snow precipitation rate (valid over diagnostic interval)",       "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'graupel_rate_accum',   NF90_FLOAT, "cumulative surface liquid water equivalent thickness of graupel precipitation rate (valid over diagnostic interval)",    "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  call NetCDF_def_var(ncid, 'conv_prcp_rate_accum', NF90_FLOAT, "cumulative surface liquid water equivalent thickness of convective precipitation rate (valid over diagnostic interval)", "m s-1", dummy_id, (/ hor_dim_id, time_id /))
  
  if (physics%Model%naux2d > 0) then
    do i=1, physics%Model%naux2d
      write(idx,'(I2)') i
      call NetCDF_def_var(ncid, 'aux2d'//adjustl(trim(idx)), NF90_FLOAT, "generic 2-D diagnostics array with index "//adjustl(trim(idx)), "unknown", dummy_id, (/ hor_dim_id, time_id /))
    end do
  end if
  
  if (physics%Model%naux3d > 0) then
    do i=1, physics%Model%naux3d
      write(idx,'(I2)') i
      call NetCDF_def_var(ncid, 'aux3d'//adjustl(trim(idx)), NF90_FLOAT, "generic 3-D diagnostics array with index "//adjustl(trim(idx)), "unknown", dummy_id, (/ hor_dim_id, vert_dim_id, time_id /))
    end do
  end if
  
end subroutine output_init_diag

!> This subroutine appends data to the "output.nc" file.
subroutine output_append(scm_state, physics)

  use gmtb_scm_type_defs, only: scm_state_type, physics_type

  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics

  real(kind=dp), allocatable :: dummy_1D(:), dummy_2d(:,:)
  
  integer :: ncid, var_id, i, j
  character(2) :: idx

  allocate(dummy_1D(scm_state%n_cols), dummy_2d(scm_state%n_cols, scm_state%n_levels))

  !> \section output_append_alg Algorithm
  !! @{

  !> - Open the file.
  CALL CHECK(NF90_OPEN(PATH=TRIM(scm_state%output_dir)//"/"//TRIM(scm_state%output_file)//".nc",MODE=NF90_WRITE,NCID=ncid))

  !> - Append all of the variables to the file.
  CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time",VARID=var_id))
  CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_out /)))

  call output_append_state(ncid, scm_state, physics)
  call output_append_forcing(ncid, scm_state)
  
  call output_append_sfcprop(ncid, scm_state, physics)
  call output_append_interstitial(ncid, scm_state, physics)
  call output_append_radtend(ncid, scm_state, physics)
  call output_append_diag(ncid, scm_state, physics)
    
  !> - Close the file.
  CALL CHECK(NF90_CLOSE(ncid))

  !> @}
end subroutine output_append

subroutine output_append_state(ncid, scm_state, physics)
  use gmtb_scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_put, only: NetCDF_put_var
  
  integer, intent(in) :: ncid
  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics
  
  call NetCDF_put_var(ncid, "pres",    scm_state%pres_l(:,:),  scm_state%itt_out)
  call NetCDF_put_var(ncid, "pres_i",  scm_state%pres_i(:,:),  scm_state%itt_out)
  call NetCDF_put_var(ncid, "sigma",   scm_state%sl(:,:),      scm_state%itt_out)
  call NetCDF_put_var(ncid, "sigma_i", scm_state%si(:,:),      scm_state%itt_out)
  call NetCDF_put_var(ncid, "pres_s",  scm_state%pres_surf(:), scm_state%itt_out)
  
  call NetCDF_put_var(ncid, "qv",      scm_state%state_tracer(:,:,scm_state%water_vapor_index,1), scm_state%itt_out)
  call NetCDF_put_var(ncid, "T",       scm_state%state_T(:,:,1), scm_state%itt_out)
  call NetCDF_put_var(ncid, "u",       scm_state%state_u(:,:,1), scm_state%itt_out)
  call NetCDF_put_var(ncid, "v",       scm_state%state_v(:,:,1), scm_state%itt_out)
  call NetCDF_put_var(ncid, "ql",      scm_state%state_tracer(:,:,scm_state%cloud_water_index,1), scm_state%itt_out)
  call NetCDF_put_var(ncid, "qi",      scm_state%state_tracer(:,:,scm_state%cloud_ice_index,1), scm_state%itt_out)
  if (physics%model%do_mynnedmf) then
    call NetCDF_put_var(ncid, "qc",    scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) + &
                                       scm_state%state_tracer(:,:,scm_state%cloud_ice_index,1)   + &
                                       physics%Tbd%QC_BL(:,:), scm_state%itt_out)
  else
    call NetCDF_put_var(ncid, "qc",    scm_state%state_tracer(:,:,scm_state%cloud_water_index,1) + &
                                       scm_state%state_tracer(:,:,scm_state%cloud_ice_index,1),    &
                                       scm_state%itt_out)
  endif
  
end subroutine output_append_state

subroutine output_append_forcing(ncid, scm_state)
    use gmtb_scm_type_defs, only: scm_state_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    
    call NetCDF_put_var(ncid, "qv_force_tend",  scm_state%qv_force_tend(:,:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "T_force_tend",   scm_state%T_force_tend(:,:),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "u_force_tend",   scm_state%u_force_tend(:,:),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "v_force_tend",   scm_state%v_force_tend(:,:),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "w_ls",           scm_state%w_ls(:,:),          scm_state%itt_out)
    call NetCDF_put_var(ncid, "u_g",            scm_state%u_g(:,:),           scm_state%itt_out)
    call NetCDF_put_var(ncid, "v_g",            scm_state%v_g(:,:),           scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_rad_forc", scm_state%dT_dt_rad(:,:),     scm_state%itt_out)
    call NetCDF_put_var(ncid, "h_advec_thil",   scm_state%h_advec_thil(:,:),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "h_advec_qt",     scm_state%h_advec_qt(:,:),    scm_state%itt_out)
    call NetCDF_put_var(ncid, "v_advec_thil",   scm_state%v_advec_thil(:,:),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "v_advec_qt",     scm_state%v_advec_qt(:,:),    scm_state%itt_out)
    call NetCDF_put_var(ncid, "T_s",            scm_state%T_surf(:),          scm_state%itt_out)
    
end subroutine output_append_forcing

subroutine output_append_sfcprop(ncid, scm_state, physics)
  use gmtb_scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_put, only: NetCDF_put_var
  use gmtb_scm_physical_constants, only: con_rd, con_fvirt, con_cp, con_hvap
  
  integer, intent(in) :: ncid
  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics
  
  character(2) :: idx
  integer :: i
  real(kind=dp), dimension(scm_state%n_cols) :: rho, shf, lhf
  
  if (scm_state%model_ics .or. scm_state%lsm_ics) then
    if (physics%Model%lsoil > 0) then
      do i=1, physics%Model%lsoil
        write(idx,'(I2)') i
        call NetCDF_put_var(ncid, 'T_soil_'//adjustl(trim(idx)),  physics%Sfcprop%stc(:,i), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_'//adjustl(trim(idx)),  physics%Sfcprop%smc(:,i), scm_state%itt_out)
      end do
    end if
  end if
  
  !convert kinematic surface fluxes to W m-2
  do i=1, scm_state%n_cols
    rho(i) = physics%Statein%pgr(i) / (con_rd*physics%Statein%tgrs(i,1)*(1.0 + con_fvirt*physics%Statein%qgrs(i,1,physics%Model%ntqv)))
    shf(i) = con_cp*rho(i)*physics%Sfcprop%hflx(i)
    lhf(i) = con_hvap*rho(i)*physics%Sfcprop%evap(i)
  end do
  call NetCDF_put_var(ncid, 'lhf', lhf, scm_state%itt_out)
  call NetCDF_put_var(ncid, 'shf', shf, scm_state%itt_out)
  call NetCDF_put_var(ncid, 'tprcp_inst', physics%Sfcprop%tprcp(:), scm_state%itt_out)
  call NetCDF_put_var(ncid, 'tprcp_rate_inst', physics%Sfcprop%tprcp(:)/scm_state%dt, scm_state%itt_out)
  
end subroutine output_append_sfcprop

subroutine output_append_interstitial(ncid, scm_state, physics)
    use gmtb_scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    call NetCDF_put_var(ncid, "tau_u",   physics%Interstitial%dusfc1(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "tau_v",   physics%Interstitial%dvsfc1(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "upd_mf",  physics%Interstitial%ud_mf(:,:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "dwn_mf",  physics%Interstitial%dd_mf(:,:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "det_mf",  physics%Interstitial%dt_mf(:,:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "sfc_up_lw_land",     physics%Interstitial%adjsfculw_land(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_lw_ice",      physics%Interstitial%adjsfculw_ice(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_lw_ocean",    physics%Interstitial%adjsfculw_ocean(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_sw_dir_nir",  physics%Interstitial%adjnirbmu(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_sw_dif_nir",  physics%Interstitial%adjnirdfu(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_sw_dir_vis",  physics%Interstitial%adjvisbmu(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_sw_dif_vis",  physics%Interstitial%adjvisdfu(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_dwn_sw_dir_nir", physics%Interstitial%adjnirbmd(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_dwn_sw_dif_nir", physics%Interstitial%adjnirdfd(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_dwn_sw_dir_vis", physics%Interstitial%adjvisbmd(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_dwn_sw_dif_vis", physics%Interstitial%adjvisdfd(:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "mp_prcp_inst",    physics%Interstitial%prcpmp(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "dcnv_prcp_inst",  physics%Interstitial%raincd(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "scnv_prcp_inst",  physics%Interstitial%raincs(:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "rad_cloud_fraction", physics%Interstitial%clouds(:,:,1), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_cloud_lwp",      physics%Interstitial%clouds(:,:,2), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_eff_rad_ql",     physics%Interstitial%clouds(:,:,3), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_cloud_iwp",      physics%Interstitial%clouds(:,:,4), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_eff_rad_qi",     physics%Interstitial%clouds(:,:,5), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_cloud_rwp",      physics%Interstitial%clouds(:,:,6), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_eff_rad_qr",     physics%Interstitial%clouds(:,:,7), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_cloud_swp",      physics%Interstitial%clouds(:,:,8), scm_state%itt_out)
    call NetCDF_put_var(ncid, "rad_eff_rad_qs",     physics%Interstitial%clouds(:,:,9), scm_state%itt_out)

end subroutine output_append_interstitial

subroutine output_append_radtend(ncid, scm_state, physics)
    use gmtb_scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    call NetCDF_put_var(ncid, "sw_rad_heating_rate",  physics%Radtend%htrsw(:,:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "lw_rad_heating_rate",  physics%Radtend%htrlw(:,:), scm_state%itt_out)
    ! TOA/SFC fluxes (on radiation timesteps)
    !physics%Diag%topfsw(i)%upfxc !sw_up_TOA_tot
    !physics%Diag%topfsw(i)%dnfxc !sw_dn_TOA_tot
    !physics%Diag%topfsw(i)%upfx0 !sw_up_TOA_clr
    !physics%Radtend%sfcfsw(i)%upfxc !sw_up_sfc_tot
    !physics%Radtend%sfcfsw(i)%dnfxc !sw_dn_sfc_tot
    !physics%Radtend%sfcfsw(i)%upfx0 !sw_up_sfc_clr
    !physics%Radtend%sfcfsw(i)%dnfx0 !sw_dn_sfc_clr
    !physics%Diag%topflw(i)%upfxc !lw_up_TOA_tot
    !physics%Diag%topflw(i)%upfx0 !lw_up_TOA_clr
    !physics%Radtend%sfcflw(i)%upfxc !lw_up_sfc_tot
    !physics%Radtend%sfcflw(i)%dnfxc !lw_dn_sfc_tot
    !physics%Radtend%sfcflw(i)%upfx0 !lw_up_sfc_clr
    !physics%Radtend%sfcflw(i)%dnfx0 !lw_dn_sfc_clr 
    
end subroutine output_append_radtend

subroutine output_append_diag(ncid, scm_state, physics)
    use gmtb_scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels, 11) :: T_tend
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels, 13) :: q_tend
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels, 8)  :: u_tend
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels, 8)  :: v_tend
    
    integer :: i,j
    character(2) :: idx
    
    real(kind=dp), dimension(scm_state%n_cols) :: temp_1d
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels) :: temp_2d
    
    !calculate tendencies from cumulative changes
    T_tend = physics%Diag%dt3dt/scm_state%dt
    q_tend = physics%Diag%dq3dt/scm_state%dt
    u_tend = physics%Diag%du3dt/scm_state%dt
    v_tend = physics%Diag%dv3dt/scm_state%dt
    
    call NetCDF_put_var(ncid, "pwat",  physics%Diag%pwat(:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "dT_dt_lwrad",     T_tend(:,:,1),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_swrad",     T_tend(:,:,2),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_pbl",       T_tend(:,:,3),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_deepconv",  T_tend(:,:,4),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_shalconv",  T_tend(:,:,5),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_micro",     T_tend(:,:,6),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_ogwd",      T_tend(:,:,7),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_rayleigh",  T_tend(:,:,8),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_cgwd",      T_tend(:,:,9),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_phys",      T_tend(:,:,10), scm_state%itt_out)
    call NetCDF_put_var(ncid, "dT_dt_nonphys",   T_tend(:,:,11), scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_pbl",       q_tend(:,:,1),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_deepconv",  q_tend(:,:,2),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_shalconv",  q_tend(:,:,3),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_micro",     q_tend(:,:,4),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_pbl",      q_tend(:,:,5),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_prodloss", q_tend(:,:,6),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_oz",       q_tend(:,:,7),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_T",        q_tend(:,:,8),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_ovhd",     q_tend(:,:,9),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_phys",      q_tend(:,:,10), scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_phys",     q_tend(:,:,11), scm_state%itt_out)
    call NetCDF_put_var(ncid, "dq_dt_nonphys",   q_tend(:,:,12), scm_state%itt_out)
    call NetCDF_put_var(ncid, "doz_dt_nonphys",  q_tend(:,:,13), scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_pbl",       u_tend(:,:,1),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_ogwd",      u_tend(:,:,2),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_deepconv",  u_tend(:,:,3),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_cgwd",      u_tend(:,:,4),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_rayleigh",  u_tend(:,:,5),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_shalconv",  u_tend(:,:,6),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_phys",      u_tend(:,:,7),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "du_dt_nonphys",   u_tend(:,:,8),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_pbl",       v_tend(:,:,1),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_ogwd",      v_tend(:,:,2),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_deepconv",  v_tend(:,:,3),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_cgwd",      v_tend(:,:,4),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_rayleigh",  v_tend(:,:,5),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_shalconv",  v_tend(:,:,6),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_phys",      v_tend(:,:,7),  scm_state%itt_out)
    call NetCDF_put_var(ncid, "dv_dt_nonphys",   v_tend(:,:,8),  scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "sfc_dwn_sw",  physics%Diag%dswsfci(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_sw",   physics%Diag%dswsfci(:) - physics%Diag%nswsfci(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_net_sw",  physics%Diag%nswsfci(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_dwn_lw",  physics%Diag%dlwsfci(:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "tprcp_accum",          physics%Diag%totprcpb(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "ice_accum",            physics%Diag%toticeb(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "snow_accum",           physics%Diag%totsnwb(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "graupel_accum",        physics%Diag%totgrpb(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "conv_prcp_accum",      physics%Diag%cnvprcpb(:), scm_state%itt_out)
    !dividing by dt only works when diagnostics are also instantaneous
    call NetCDF_put_var(ncid, "tprcp_rate_accum",     physics%Diag%totprcpb(:)/scm_state%dt, scm_state%itt_out)
    call NetCDF_put_var(ncid, "ice_rate_accum",       physics%Diag%toticeb(:)/scm_state%dt, scm_state%itt_out)
    call NetCDF_put_var(ncid, "snow_rate_accum",      physics%Diag%totsnwb(:)/scm_state%dt, scm_state%itt_out)
    call NetCDF_put_var(ncid, "graupel_rate_accum",   physics%Diag%totgrpb(:)/scm_state%dt, scm_state%itt_out)
    call NetCDF_put_var(ncid, "conv_prcp_rate_accum", physics%Diag%cnvprcpb(:)/scm_state%dt, scm_state%itt_out)
    
    if (physics%Model%naux2d > 0) then
      do j=1, physics%Model%naux2d
        write(idx,'(I2)') j
        do i=1, scm_state%n_cols
           temp_1d(i) = physics%Diag%aux2d(i,j)
        end do
        call NetCDF_put_var(ncid, "aux2d"//adjustl(trim(idx)),  temp_1d(:), scm_state%itt_out)
      end do
    end if
    
    if (physics%Model%naux3d > 0) then
      do j=1, physics%Model%naux3d
        write(idx,'(I2)') j
        do i=1, scm_state%n_cols
           temp_2d(i,:) = physics%Diag%aux3d(i,:,j)
        end do
        call NetCDF_put_var(ncid, "aux3d"//adjustl(trim(idx)),  temp_2d(:,:), scm_state%itt_out)
      end do
    end if
    
end subroutine output_append_diag

!> @}
!> @}
end module gmtb_scm_output
