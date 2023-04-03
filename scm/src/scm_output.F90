!> \file scm_output.f90
!!  Contains output-related subroutines

module scm_output

use netcdf
use NetCDF_read, only: check
use scm_kinds, only: sp, dp, qp

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup output scm_output
!! @{
!! Contains output-related subroutines

!> This subroutine initializes the output netCDF file, "output.nc", placed in the directory specified by the case_config file used.
subroutine output_init(scm_state, physics)
  use scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_def, only: NetCDF_def_var
  use NetCDF_put, only: NetCDF_put_var
  use NetCDF_read, only: missing_value

  type(scm_state_type), intent(in) :: scm_state
  type(physics_type),   intent(in) :: physics
  
  integer :: n_timesteps, n_inst, n_diag, n_swrad, n_lwrad, n_rad
  integer :: ncid, hor_dim_id, vert_dim_id, vert_dim_i_id, vert_dim_rad_id, vert_dim_soil_id, dummy_id
  integer :: time_inst_id, time_diag_id, time_swrad_id, time_lwrad_id, time_rad_id
  integer :: year_id, month_id, day_id, hour_id, min_id, time_swrad_var_id, time_lwrad_var_id, time_rad_var_id
  character(2) :: idx
  
  real(kind=dp), dimension(scm_state%n_cols) :: missing_value_1D
  real(kind=dp), dimension(scm_state%n_cols,scm_state%n_levels) :: missing_value_2D
  
  !> \section output_init_alg Algorithm
  !! @{

  !> - Create the output file in the output directory.
  CALL CHECK(NF90_CREATE(PATH=TRIM(scm_state%output_dir)//"/"//TRIM(scm_state%output_file)//".nc",CMODE=NF90_CLOBBER,NCID=ncid),'nf90_create()')

  !> - Define netCDF dimensions.
  n_timesteps = ceiling(scm_state%runtime/scm_state%dt)
  
  !output_append is hard-coded to be called once after initialization and once after the first timestep
  if(scm_state%n_itt_out > 1) then
    n_inst = n_timesteps/scm_state%n_itt_out + 1 !one added for hard-coded output after first time step
  else
    n_inst = n_timesteps + 1 !one added for hard-coded initialization output
  end if
  if(physics%Model%nszero == 1) then
    n_diag = n_timesteps/physics%Model%nszero + 1
  else
    n_diag = n_timesteps/physics%Model%nszero
  end if
  if (physics%Model%nsswr > 0) then
    n_swrad = n_timesteps/physics%Model%nsswr + 1
  else
    !if SW radiation is not called, write out missing values for one time level
    n_swrad = 1
  end if
  if (physics%Model%nslwr > 0) then
    !if LW radiation is not called, write out missing values for one time level
    n_lwrad = n_timesteps/physics%Model%nslwr + 1
  else
    n_lwrad = 1
  end if
  if (physics%Model%nslwr > 0) then
    if (physics%Model%nsswr > 0) then
      n_rad = min(n_swrad, n_lwrad)
    else
      n_rad = n_lwrad
    end if
  else
    if (physics%Model%nsswr > 0) then
      n_rad = n_swrad
    else
      n_rad = 1
    end if
  end if
  
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_inst_dim",LEN=n_inst,DIMID=time_inst_id),"nf90_def_dim(time_inst_dim)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_diag_dim",LEN=n_diag,DIMID=time_diag_id),"nf90_def_dim(time_diag_dim)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_swrad_dim",LEN=n_swrad,DIMID=time_swrad_id),"nf90_def_dim(time_swrad_dim)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_lwrad_dim",LEN=n_lwrad,DIMID=time_lwrad_id),"nf90_def_dim(time_lwrad_dim)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="time_rad_dim",LEN=n_rad,DIMID=time_rad_id),"nf90_def_dim(time_rad_dim)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="hor_dim_layer",LEN=scm_state%n_cols,DIMID=hor_dim_id),"nf90_def_dim(hor_dim_layer)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_layer",LEN=scm_state%n_levels,DIMID=vert_dim_id),"nf90_def_dim(vert_dim_layer)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_interface",LEN=scm_state%n_levels+1,DIMID=vert_dim_i_id),"nf90_def_dim(vert_dim_interface)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_rad",LEN=physics%Interstitial%lmk,DIMID=vert_dim_rad_id),"nf90_def_dim(vert_dim_rad)")
  CALL CHECK(NF90_DEF_DIM(NCID=ncid,NAME="vert_dim_soil",LEN=physics%Model%lsoil_lsm,DIMID=vert_dim_soil_id),"nf90_def_dim(vert_dim_soil)")

  !> - Define the dimension variables.
  call NetCDF_def_var(ncid, 'time_inst', NF90_FLOAT, "model elapsed time for instantaneous variables", "s", dummy_id, (/ time_inst_id /))
  call NetCDF_def_var(ncid, 'time_diag', NF90_FLOAT, "model elapsed time for diagnostic variables", "s", dummy_id, (/ time_diag_id /))
  call NetCDF_def_var(ncid, 'time_swrad', NF90_FLOAT, "model elapsed time for shortwave radiation variables", "s", time_swrad_var_id, (/ time_swrad_id /))
  call NetCDF_def_var(ncid, 'time_lwrad', NF90_FLOAT, "model elapsed time for longwave radiation variables", "s", time_lwrad_var_id, (/ time_lwrad_id /))
  call NetCDF_def_var(ncid, 'time_rad', NF90_FLOAT, "model elapsed time for either LW or SW radiation variables", "s", time_rad_var_id, (/ time_rad_id /))

  !> - Define the state variables
  CALL output_init_state(ncid, time_inst_id, hor_dim_id, vert_dim_id, vert_dim_i_id)
  !> - Define the forcing variables
  CALL output_init_forcing(ncid, time_inst_id, hor_dim_id, vert_dim_id)
  
  !> - Define all diagnostic/physics variables
  CALL output_init_sfcprop(ncid, time_inst_id, hor_dim_id, vert_dim_soil_id, scm_state, physics)
  CALL output_init_interstitial(ncid, time_inst_id, time_rad_id, hor_dim_id, vert_dim_id, vert_dim_rad_id, physics)
  CALL output_init_radtend(ncid, time_swrad_id, time_lwrad_id, hor_dim_id, vert_dim_id)
  CALL output_init_diag(ncid, time_inst_id, time_diag_id, time_rad_id, hor_dim_id, vert_dim_id, physics)
  
  call NetCDF_def_var(ncid, 'init_year',   NF90_INT, "model initialization year",   "year",   year_id)
  call NetCDF_def_var(ncid, 'init_month',  NF90_INT, "model initialization month",  "month",  month_id)
  call NetCDF_def_var(ncid, 'init_day',    NF90_INT, "model initialization day",    "day",    day_id)
  call NetCDF_def_var(ncid, 'init_hour',   NF90_INT, "model initialization hour",   "hour",   hour_id)
  call NetCDF_def_var(ncid, 'init_minute', NF90_INT, "model initialization minute", "minute", min_id)
  
  !> - Close variable definition and the file.
  CALL CHECK(NF90_ENDDEF(NCID=ncid),'nf90_enddef()')
  
  call NetCDF_put_var(ncid, "init_year",   scm_state%init_year,  year_id)
  call NetCDF_put_var(ncid, "init_month",  scm_state%init_month, month_id)
  call NetCDF_put_var(ncid, "init_day",    scm_state%init_day,   day_id)
  call NetCDF_put_var(ncid, "init_hour",   scm_state%init_hour,  hour_id)
  call NetCDF_put_var(ncid, "init_minute", scm_state%init_min,   min_id)
  
  if (physics%Model%nsswr <= 0) then
    !write out missing values at the initial time
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=time_swrad_var_id,VALUES=0.0,START=(/ 1 /)),"nf90_put_var(time_swrad_var)")
    missing_value_2D = missing_value
    call NetCDF_put_var(ncid, "sw_rad_heating_rate",  missing_value_2D, 1)
  end if
  if (physics%Model%nslwr <= 0) then
    !write out missing values at the initial time
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=time_lwrad_var_id,VALUES=0.0,START=(/ 1 /)),"nf90_put_var(time_lwrad_var)")
    missing_value_2D = missing_value
    call NetCDF_put_var(ncid, "lw_rad_heating_rate",  missing_value_2D, 1)
  end if
  if (physics%Model%nsswr <= 0 .and. physics%Model%nslwr <= 0) then
    !write out missing values at the initial time
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=time_rad_var_id,VALUES=0.0,START=(/ 1 /)),"nf90_put_var(time_rad_var)")
    missing_value_2D = missing_value
    missing_value_1D = missing_value
    call NetCDF_put_var(ncid, "rad_cloud_fraction", missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_cloud_lwp",      missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_eff_rad_ql",     missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_cloud_iwp",      missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_eff_rad_qi",     missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_cloud_rwp",      missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_eff_rad_qr",     missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_cloud_swp",      missing_value_2D, 1)
    call NetCDF_put_var(ncid, "rad_eff_rad_qs",     missing_value_2D, 1)
    call NetCDF_put_var(ncid, 'max_cloud_fraction', missing_value_1D, 1)
    call NetCDF_put_var(ncid, 'toa_total_albedo',   missing_value_1D, 1)
    call NetCDF_put_var(ncid, 'vert_int_lwp_mp',    missing_value_1D, 1)
    call NetCDF_put_var(ncid, 'vert_int_iwp_mp',    missing_value_1D, 1)
    call NetCDF_put_var(ncid, 'vert_int_lwp_cf',    missing_value_1D, 1)
    call NetCDF_put_var(ncid, 'vert_int_iwp_cf',    missing_value_1D, 1)
  end if
  
  CALL CHECK(NF90_CLOSE(NCID=ncid),"nf90_close()")

  !> @}
end subroutine output_init

subroutine output_init_state(ncid, time_inst_id, hor_dim_id, vert_dim_id, vert_dim_i_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_inst_id, hor_dim_id, vert_dim_id, vert_dim_i_id
  
  integer :: dummy_id
  
  !> - Define the pressure-related variables.
  call NetCDF_def_var(ncid, 'pres',    NF90_FLOAT, "pressure on model layer centers",         "Pa",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'pres_i',  NF90_FLOAT, "pressure on model layer interfaces",      "Pa",   dummy_id, (/ hor_dim_id, vert_dim_i_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sigma',   NF90_FLOAT, "sigma (p/p_s) on model layer centers",    "none", dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'sigma_i', NF90_FLOAT, "sigma (p/p_s) on model layer interfaces", "none", dummy_id, (/ hor_dim_id, vert_dim_i_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'pres_s',  NF90_FLOAT, "surface pressure",                        "Pa",   dummy_id, (/ hor_dim_id,                time_inst_id /))
  
  !> - Define the state variables.
  call NetCDF_def_var(ncid, 'qv', NF90_FLOAT, "water vapor specific humidity on model layer centers",                "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'T',  NF90_FLOAT, "absolute temperature on model layer centers",                         "K",       dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'u',  NF90_FLOAT, "x-wind on model layer centers",                                       "m s-1",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'v',  NF90_FLOAT, "y-wind on model layer centers",                                       "m s-1",   dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'ql', NF90_FLOAT, "suspended resolved liquid cloud water on model layer centers",        "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'qi', NF90_FLOAT, "suspended resolved ice cloud water on model layer centers",           "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  call NetCDF_def_var(ncid, 'qc', NF90_FLOAT, "suspended (resolved + SGS) total cloud water on model layer centers", "kg kg-1", dummy_id, (/ hor_dim_id, vert_dim_id,   time_inst_id /))
  
end subroutine output_init_state

subroutine output_init_forcing(ncid, time_inst_id, hor_dim_id, vert_dim_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_inst_id, hor_dim_id, vert_dim_id
  
  integer :: dummy_id
  
  !> - Define the forcing-related variables.
  call NetCDF_def_var(ncid, 'qv_force_tend',  NF90_FLOAT, "total forcing tendency for water vapor specific humidity",                                 "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'T_force_tend',   NF90_FLOAT, "total forcing tendency for absolute temperature",                                          "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'u_force_tend',   NF90_FLOAT, "total forcing tendency for x-wind",                                                        "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'v_force_tend',   NF90_FLOAT, "total forcing tendency for y-wind",                                                        "m s-2",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'w_ls',           NF90_FLOAT, "large-scale vertical velocity forcing",                                                    "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'u_g',            NF90_FLOAT, "geostrophic x-wind forcing",                                                               "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'v_g',            NF90_FLOAT, "geostrophic y-wind forcing",                                                               "m s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'dT_dt_rad_forc', NF90_FLOAT, "prescribed radiative heating rate forcing",                                                "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'h_advec_thil',   NF90_FLOAT, "large-scale horizontal advective heating rate for ice-liquid water potential temperature", "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'h_advec_qt',     NF90_FLOAT, "large-scale horizontal advective moistening rate for total water specific humidity",       "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'v_advec_thil',   NF90_FLOAT, "large-scale vertical advective heating rate for ice-liquid water potential temperature",   "K s-1",       dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'v_advec_qt',     NF90_FLOAT, "large-scale vertical advective moistening rate for total water specific humidity",         "kg kg-1 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'T_s',            NF90_FLOAT, "surface temperature forcing",                                                              "K",           dummy_id, (/ hor_dim_id,              time_inst_id /))

end subroutine output_init_forcing

subroutine output_init_sfcprop(ncid, time_inst_id, hor_dim_id, vert_dim_soil_id, scm_state, physics)
  use scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_inst_id, hor_dim_id, vert_dim_soil_id
  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics

  integer :: i, dummy_id
  
  if (scm_state%model_ics .or. scm_state%lsm_ics) then
    if (physics%Model%lsoil_lsm > 0) then
        call NetCDF_def_var(ncid, 'soil_T',  NF90_FLOAT, "soil temperature profile ", "K", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
        call NetCDF_def_var(ncid, 'soil_moisture',  NF90_FLOAT, "soil moisture profile ", "m3 m-3", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
        call NetCDF_def_var(ncid, 'soil_moisture_unfrozen',  NF90_FLOAT, "unfrozen soil moisture profile ", "m3 m-3", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
        
        call NetCDF_def_var(ncid, 'soil_moisture_total_vol_frac',  NF90_FLOAT, "total soil moisture volume fraction profile ", "", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
        call NetCDF_def_var(ncid, 'soil_moisture_unfrozen_vol_frac',  NF90_FLOAT, "unfrozen soil moisture volume fraction profile ", "", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
        call NetCDF_def_var(ncid, 'soil_moisture_frozen_vol_frac',  NF90_FLOAT, "frozen soil moisture volume fraction profile ", "", dummy_id, (/ hor_dim_id, vert_dim_soil_id, time_inst_id /))
    end if
  end if
  
  call NetCDF_def_var(ncid, 'lhf', NF90_FLOAT, "surface latent heat flux", "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'shf', NF90_FLOAT, "surface sensible heat flux", "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'tprcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation", "m", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'tprcp_rate_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation rate (tprcp_inst/dt)", "m s-1", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 't2m', NF90_FLOAT, "2-m temperature", "K", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'q2m', NF90_FLOAT, "2-m specific humidity", "kg kg-1", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'ustar', NF90_FLOAT, "surface friction velocity", "m s-1", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'tsfc', NF90_FLOAT, "surface skin temperature", "m s-1K", dummy_id, (/ hor_dim_id, time_inst_id /))
  
end subroutine output_init_sfcprop

subroutine output_init_interstitial(ncid, time_inst_id, time_rad_id, hor_dim_id, vert_dim_id, vert_dim_rad_id, physics)
  use scm_type_defs, only: physics_type
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_inst_id, time_rad_id, hor_dim_id, vert_dim_id, vert_dim_rad_id
  type(physics_type), intent(in) :: physics
  
  integer :: dummy_id
  
  call NetCDF_def_var(ncid, 'tau_u',  NF90_FLOAT, "surface x-wind stress", "Pa",         dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'tau_v',  NF90_FLOAT, "surface y-wind stress", "Pa",         dummy_id, (/ hor_dim_id, time_inst_id /))
  if (physics%model%imfdeepcnv >= 0 .or. physics%model%imfshalcnv >= 0) then
    call NetCDF_def_var(ncid, 'upd_mf', NF90_FLOAT, "updraft mass flux",     "kg m-2 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  end if
  call NetCDF_def_var(ncid, 'dwn_mf', NF90_FLOAT, "downdraft mass flux",   "kg m-2 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'det_mf', NF90_FLOAT, "detrainment mass flux", "kg m-2 s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
  
  call NetCDF_def_var(ncid, 'sfc_up_lw_land',     NF90_FLOAT, "surface upwelling longwave flux over land fraction (valid all timesteps)",                 "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_lw_ice',      NF90_FLOAT, "surface upwelling longwave flux over ice fraction (valid all timesteps)",                  "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_lw_water',    NF90_FLOAT, "surface upwelling longwave flux over water fraction (valid all timesteps)",                "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dir_nir',  NF90_FLOAT, "surface upwelling shortwave direct near-infrared flux (valid all timesteps)",              "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dif_nir',  NF90_FLOAT, "surface upwelling shortwave diffuse near-infrared flux (valid all timesteps)",             "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dir_vis',  NF90_FLOAT, "surface upwelling shortwave direct visible and ultraviolet flux (valid all timesteps)",    "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw_dif_vis',  NF90_FLOAT, "surface upwelling shortwave diffuse visible and ultraviolet flux (valid all timesteps)",   "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dir_nir', NF90_FLOAT, "surface downwelling shortwave direct near-infrared flux (valid all timesteps)",            "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dif_nir', NF90_FLOAT, "surface downwelling shortwave diffuse near-infrared flux (valid all timesteps)",           "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dir_vis', NF90_FLOAT, "surface downwelling shortwave direct visible and ultraviolet flux (valid all timesteps)",  "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_sw_dif_vis', NF90_FLOAT, "surface downwelling shortwave diffuse visible and ultraviolet flux (valid all timesteps)", "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  
  call NetCDF_def_var(ncid, 'mp_prcp_inst',   NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from microphysics scheme",       "m", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'dcnv_prcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from deep convection scheme",    "m", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'scnv_prcp_inst', NF90_FLOAT, "instantaneous surface liquid water equivalent thickness of total precipitation from shallow convection scheme", "m", dummy_id, (/ hor_dim_id, time_inst_id /))
  
  call NetCDF_def_var(ncid, 'rad_cloud_fraction', NF90_FLOAT, "instantaneous cloud fraction used in radiation (radiation timesteps only)",                   "fraction", dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_lwp',      NF90_FLOAT, "instantaneous cloud liquid water path used in radiation (radiation timesteps only)",          "g m-2",    dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_ql',     NF90_FLOAT, "instantaneous effective radius for liquid cloud used in radiation (radiation timesteps only)", "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_iwp',      NF90_FLOAT, "instantaneous cloud ice water path used in radiation (radiation timesteps only)",              "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qi',     NF90_FLOAT, "instantaneous effective radius for ice cloud used in radiation (radiation timesteps only)",    "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_rwp',      NF90_FLOAT, "instantaneous rain water path used in radiation (radiation timesteps only)",                   "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qr',     NF90_FLOAT, "instantaneous effective radius for raindrop used in radiation (radiation timesteps only)",     "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_cloud_swp',      NF90_FLOAT, "instantaneous snow water path used in radiation (radiation timesteps only)",                   "g m-2",   dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'rad_eff_rad_qs',     NF90_FLOAT, "instantaneous effective radius for snowflake in radiation (radiation timesteps only)",         "um",      dummy_id, (/ hor_dim_id, vert_dim_rad_id, time_rad_id /))
  
end subroutine output_init_interstitial

subroutine output_init_radtend(ncid, time_swrad_id, time_lwrad_id, hor_dim_id, vert_dim_id)
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_swrad_id, time_lwrad_id, hor_dim_id, vert_dim_id
  
  integer :: dummy_id
  
  call NetCDF_def_var(ncid, 'sw_rad_heating_rate', NF90_FLOAT, "total sky shortwave radiative heating rate (radiation timesteps only)", "K s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_swrad_id /))
  call NetCDF_def_var(ncid, 'lw_rad_heating_rate', NF90_FLOAT, "total sky longwave radiative heating rate (radiation timesteps only)",  "K s-1", dummy_id, (/ hor_dim_id, vert_dim_id, time_lwrad_id /))
  
end subroutine output_init_radtend

subroutine output_init_diag(ncid, time_inst_id, time_diag_id, time_rad_id, hor_dim_id, vert_dim_id, physics)
  use scm_type_defs, only: physics_type
  use NetCDF_def, only : NetCDF_def_var
  
  integer, intent(in) :: ncid, time_inst_id, time_diag_id, time_rad_id, hor_dim_id, vert_dim_id
  type(physics_type), intent(in) :: physics
  
  integer :: i, dummy_id
  character(2) :: idx
  
  call NetCDF_def_var(ncid, 'pwat',            NF90_FLOAT, "column precipitable water", "kg m-2", dummy_id, (/ hor_dim_id, time_inst_id /)) !the variable is reset every timestep in GFS_MP_generic
  
  call output_init_tendency(ncid, 'dT_dt_lwrad',     "temperature tendency due to longwave radiation scheme",        "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_longwave))
  call output_init_tendency(ncid, 'dT_dt_swrad',     "temperature tendency due to shortwave radiation scheme",       "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_shortwave))
  call output_init_tendency(ncid, 'dT_dt_pbl',       "temperature tendency due to PBL scheme",                       "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_pbl))
  call output_init_tendency(ncid, 'dT_dt_deepconv',  "temperature tendency due to deep convection scheme",           "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_dcnv))
  call output_init_tendency(ncid, 'dT_dt_shalconv',  "temperature tendency due to shallow convection scheme",        "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_scnv))
  call output_init_tendency(ncid, 'dT_dt_micro',     "temperature tendency due to microphysics scheme",              "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_mp))
  call output_init_tendency(ncid, 'dT_dt_ogwd',      "temperature tendency due to orographic GWD scheme",            "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_orographic_gwd))
  call output_init_tendency(ncid, 'dT_dt_rayleigh',  "temperature tendency due to rayleigh damping scheme",          "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_rayleigh_damping))
  call output_init_tendency(ncid, 'dT_dt_cgwd',      "temperature tendency due to convective GWD scheme",            "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_nonorographic_gwd))
  call output_init_tendency(ncid, 'dT_dt_phys',      "temperature tendency due to all physics schemes",              "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_physics))
  call output_init_tendency(ncid, 'dT_dt_nonphys',   "temperature tendency due to all processes other than physics", "K s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_non_physics))
  call output_init_tendency(ncid, 'dq_dt_pbl',       "moisture tendency due to PBL scheme",                       "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_pbl))
  call output_init_tendency(ncid, 'dq_dt_deepconv',  "moisture tendency due to deep convection scheme",           "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_dcnv))
  call output_init_tendency(ncid, 'dq_dt_shalconv',  "moisture tendency due to shallow convection scheme",        "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_scnv))
  call output_init_tendency(ncid, 'dq_dt_micro',     "moisture tendency due to microphysics scheme",              "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_mp))
  call output_init_tendency(ncid, 'dq_dt_phys',      "moisture tendency due to all physics schemes",              "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_physics))
  call output_init_tendency(ncid, 'dq_dt_nonphys',   "moisture tendency due to all processes other than physics", "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_non_physics))
  call output_init_tendency(ncid, 'doz_dt_pbl',      "ozone tendency due to PBL scheme",                          "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_pbl))
  call output_init_tendency(ncid, 'doz_dt_prodloss', "ozone tendency due to ozone production/loss",               "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_prod_loss))
  call output_init_tendency(ncid, 'doz_dt_oz',       "ozone tendency due to ozone",                               "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_ozmix))
  call output_init_tendency(ncid, 'doz_dt_T',        "ozone tendency due to temperature",                         "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_temp))
  call output_init_tendency(ncid, 'doz_dt_ovhd',     "ozone tendency due to overhead ozone column",               "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_overhead_ozone))
  call output_init_tendency(ncid, 'doz_dt_phys',     "ozone tendency due to all physics schemes",                 "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_physics))
  call output_init_tendency(ncid, 'doz_dt_nonphys',  "ozone tendency due to all processes other than physics",    "kg kg-1 s-1", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_non_physics))
  call output_init_tendency(ncid, 'du_dt_pbl',       "x-wind tendency due to PBL scheme",                       "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_pbl))
  call output_init_tendency(ncid, 'du_dt_ogwd',      "x-wind tendency due to orographic GWD scheme",            "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_orographic_gwd))
  call output_init_tendency(ncid, 'du_dt_deepconv',  "x-wind tendency due to deep convection scheme",           "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_dcnv))
  call output_init_tendency(ncid, 'du_dt_cgwd',      "x-wind tendency due to convective GWD scheme",            "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_nonorographic_gwd))
  call output_init_tendency(ncid, 'du_dt_rayleigh',  "x-wind tendency due to rayleigh damping scheme",          "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_rayleigh_damping))
  call output_init_tendency(ncid, 'du_dt_shalconv',  "x-wind tendency due to shallow convection scheme",        "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_scnv))
  call output_init_tendency(ncid, 'du_dt_phys',      "x-wind tendency due to all physics schemes",              "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_physics))
  call output_init_tendency(ncid, 'du_dt_nonphys',   "x-wind tendency due to all processes other than physics", "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_non_physics))
  call output_init_tendency(ncid, 'dv_dt_pbl',       "y-wind tendency due to PBL scheme",                       "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_pbl))
  call output_init_tendency(ncid, 'dv_dt_ogwd',      "y-wind tendency due to orographic GWD scheme",            "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_orographic_gwd))
  call output_init_tendency(ncid, 'dv_dt_deepconv',  "y-wind tendency due to deep convection scheme",           "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_dcnv))
  call output_init_tendency(ncid, 'dv_dt_cgwd',      "y-wind tendency due to convective GWD scheme",            "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_nonorographic_gwd))
  call output_init_tendency(ncid, 'dv_dt_rayleigh',  "y-wind tendency due to rayleigh damping scheme",          "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_rayleigh_damping))
  call output_init_tendency(ncid, 'dv_dt_shalconv',  "y-wind tendency due to shallow convection scheme",        "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_scnv))
  call output_init_tendency(ncid, 'dv_dt_phys',      "y-wind tendency due to all physics schemes",              "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_physics))
  call output_init_tendency(ncid, 'dv_dt_nonphys',   "y-wind tendency due to all processes other than physics", "m s-2", hor_dim_id, vert_dim_id, time_diag_id, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_non_physics))

  call NetCDF_def_var(ncid, 'sfc_dwn_sw',      NF90_FLOAT, "surface downwelling shortwave flux (valid all timesteps)",                   "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_up_sw',       NF90_FLOAT, "surface upwelling shortwave flux (valid all timesteps)",                     "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_net_sw',      NF90_FLOAT, "surface net shortwave flux (downwelling - upwelling) (valid all timesteps)", "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'sfc_dwn_lw',      NF90_FLOAT, "surface downwelling longwave flux (valid all timesteps)",                    "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'gflux',           NF90_FLOAT, "instantaneous surface ground heat flux (valid all timesteps)",               "W m-2", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'u10m',            NF90_FLOAT, "instantaneous 10-m zonal wind",                                              "m s-1", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'v10m',            NF90_FLOAT, "instantaneous 10-m meridional wind",                                         "m s-1", dummy_id, (/ hor_dim_id, time_inst_id /))
  call NetCDF_def_var(ncid, 'hpbl',            NF90_FLOAT, "instantaneous planetary boundary layer height",                              "m",     dummy_id, (/ hor_dim_id, time_inst_id /))
  
  call NetCDF_def_var(ncid, 'tprcp_accum',          NF90_FLOAT, "cumulative surface liquid water equivalent thickness of total precipitation (valid over diagnostic interval)",           "m",     dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'ice_accum',            NF90_FLOAT, "cumulative surface liquid water equivalent thickness of ice precipitation (valid over diagnostic interval)",             "m",     dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'snow_accum',           NF90_FLOAT, "cumulative surface liquid water equivalent thickness of snow precipitation (valid over diagnostic interval)",            "m",     dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'graupel_accum',        NF90_FLOAT, "cumulative surface liquid water equivalent thickness of graupel precipitation (valid over diagnostic interval)",         "m",     dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'conv_prcp_accum',      NF90_FLOAT, "cumulative surface liquid water equivalent thickness of convective precipitation (valid over diagnostic interval)",      "m",     dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'tprcp_rate_accum',     NF90_FLOAT, "cumulative surface liquid water equivalent thickness of total precipitation rate (valid over diagnostic interval)",      "m s-1", dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'ice_rate_accum',       NF90_FLOAT, "cumulative surface liquid water equivalent thickness of ice precipitation rate (valid over diagnostic interval)",        "m s-1", dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'snow_rate_accum',      NF90_FLOAT, "cumulative surface liquid water equivalent thickness of snow precipitation rate (valid over diagnostic interval)",       "m s-1", dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'graupel_rate_accum',   NF90_FLOAT, "cumulative surface liquid water equivalent thickness of graupel precipitation rate (valid over diagnostic interval)",    "m s-1", dummy_id, (/ hor_dim_id, time_diag_id /))
  call NetCDF_def_var(ncid, 'conv_prcp_rate_accum', NF90_FLOAT, "cumulative surface liquid water equivalent thickness of convective precipitation rate (valid over diagnostic interval)", "m s-1", dummy_id, (/ hor_dim_id, time_diag_id /))
  
  call NetCDF_def_var(ncid, 'max_cloud_fraction', NF90_FLOAT, "maximum cloud fraction in the column (radiation timesteps only)",                           "fraction",     dummy_id, (/ hor_dim_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'toa_total_albedo',   NF90_FLOAT, "top-of-atmosphere total sky albedo (radiation timesteps only)",                             "fraction",     dummy_id, (/ hor_dim_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'vert_int_lwp_mp',    NF90_FLOAT, "vertically-integrated liquid water path from microphysics (radiation timesteps only)",      "kg m-2",     dummy_id, (/ hor_dim_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'vert_int_iwp_mp',    NF90_FLOAT, "vertically-integrated ice water path from microphysics (radiation timesteps only)",         "kg m-2",     dummy_id, (/ hor_dim_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'vert_int_lwp_cf',    NF90_FLOAT, "vertically-integrated liquid water path from fractional clouds (radiation timesteps only)", "kg m-2",     dummy_id, (/ hor_dim_id, time_rad_id /))
  call NetCDF_def_var(ncid, 'vert_int_iwp_cf',    NF90_FLOAT, "vertically-integrated ice water path from fractional clouds (radiation timesteps only)",    "kg m-2",     dummy_id, (/ hor_dim_id, time_rad_id /))
  
  !all auxilliary diagnostics will be output every timestep (logic will need to be added to implement time averaging -- including resetting in GFS_typedefs/phys_diag_zero)
  if (physics%Model%naux2d > 0) then
    do i=1, physics%Model%naux2d
      write(idx,'(I2)') i
      call NetCDF_def_var(ncid, 'aux2d'//adjustl(trim(idx)), NF90_FLOAT, "generic 2-D diagnostics array with index "//adjustl(trim(idx)), "unknown", dummy_id, (/ hor_dim_id, time_inst_id /))
    end do
  end if
  
  if (physics%Model%naux3d > 0) then
    do i=1, physics%Model%naux3d
      write(idx,'(I2)') i
      call NetCDF_def_var(ncid, 'aux3d'//adjustl(trim(idx)), NF90_FLOAT, "generic 3-D diagnostics array with index "//adjustl(trim(idx)), "unknown", dummy_id, (/ hor_dim_id, vert_dim_id, time_inst_id /))
    end do
  end if
  
end subroutine output_init_diag

!> This subroutine appends data to the "output.nc" file.
subroutine output_append(scm_state, physics, force)

  use scm_type_defs, only: scm_state_type, physics_type

  type(scm_state_type), intent(inout) :: scm_state
  type(physics_type), intent(in) :: physics
  logical, optional, intent(in) :: force
  
  integer :: ncid, var_id
  logical :: force_inst
  
  if (PRESENT(force)) then
    force_inst = force
  else
    force_inst = .false.
  end if
  
  !> \section output_append_alg Algorithm
  !! @{
  if (.not. (mod(scm_state%itt, scm_state%n_itt_out)==0 .or. &
            physics%Model%lsswr .or. physics%Model%lslwr .or. &
            mod(scm_state%itt,physics%Model%nszero) == 0 .or. &
            force_inst)) then
    return
  end if
  !> - Open the file.
  CALL CHECK(NF90_OPEN(PATH=TRIM(scm_state%output_dir)//"/"//TRIM(scm_state%output_file)//".nc",MODE=NF90_WRITE,NCID=ncid),"nf90_open()")
  
  !> - Append all of the variables to the file.
  
  if (mod(scm_state%itt, scm_state%n_itt_out)==0 .or. force_inst) then
    scm_state%itt_out = scm_state%itt_out+1
    
    CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time_inst",VARID=var_id),"nf90_inq_varid(time_inst)")
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_out /)),"nf90_put_var(time_inst)")
    
    call output_append_state(ncid, scm_state, physics)
    call output_append_forcing(ncid, scm_state)
    call output_append_sfcprop(ncid, scm_state, physics)
    call output_append_interstitial_inst(ncid, scm_state, physics)
    call output_append_diag_inst(ncid, scm_state, physics)
  end if
  
  if(physics%Model%lslwr .or. physics%Model%lsswr) then
    if (physics%Model%lsswr) then
      scm_state%itt_swrad = scm_state%itt_swrad + 1
      CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time_swrad",VARID=var_id),"nf90_inq_varid(time_swrad)")
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_swrad /)),"nf90_put_var(time_swrad)")
    end if
    if (physics%Model%lslwr) then
      scm_state%itt_lwrad = scm_state%itt_lwrad + 1
      CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time_lwrad",VARID=var_id),"nf90_inq_varid(time_lwrad)")
      CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_lwrad /)),"nf90_put_var(time_lwrad)")
    end if
    call output_append_radtend(ncid, scm_state, physics)
    scm_state%itt_rad = scm_state%itt_rad + 1
    CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time_rad",VARID=var_id),"nf90_inq_varid(time_rad)")
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_rad /)),"nf90_put_var(time_rad)")
    call output_append_interstitial_rad(ncid, scm_state, physics)
    call output_append_diag_rad(ncid, scm_state, physics)
  end if
  
  if(mod(scm_state%itt,physics%Model%nszero) == 0) then
    scm_state%itt_diag = scm_state%itt_diag + 1
    CALL CHECK(NF90_INQ_VARID(NCID=ncid,NAME="time_diag",VARID=var_id),"nf90_inq_varid(time_diag)")
    CALL CHECK(NF90_PUT_VAR(NCID=ncid,VARID=var_id,VALUES=scm_state%model_time,START=(/ scm_state%itt_diag /)),"nf90_put_var(time_diag)")
    call output_append_diag_avg(ncid, scm_state, physics)
  end if
  
  !> - Close the file.
  CALL CHECK(NF90_CLOSE(ncid),"nf90_close()")

  !> @}
end subroutine output_append

subroutine output_append_state(ncid, scm_state, physics)
  use scm_type_defs, only: scm_state_type, physics_type
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
    use scm_type_defs, only: scm_state_type
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
  use scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_put, only: NetCDF_put_var
  use NetCDF_read, only: missing_value
  use scm_physical_constants, only: con_rd, con_fvirt, con_cp, con_hvap
  
  integer, intent(in) :: ncid
  type(scm_state_type), intent(in) :: scm_state
  type(physics_type), intent(in) :: physics
  
  integer :: i
  real(kind=dp), dimension(scm_state%n_cols) :: rho, shf, lhf
  real(kind=dp), dimension(scm_state%n_cols,physics%Model%lsoil_lsm) :: missing_value_2D
  
  if (scm_state%model_ics .or. scm_state%lsm_ics) then
    if (physics%Model%lsoil_lsm > 0) then
      missing_value_2D = missing_value
      if (physics%Model%lsm == physics%Model%lsm_noah .or. physics%Model%lsm == physics%Model%lsm_noahmp) then
        call NetCDF_put_var(ncid, 'soil_T',                          physics%Sfcprop%stc(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture',                   physics%Sfcprop%smc(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_unfrozen',          physics%Sfcprop%slc(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_total_vol_frac',    missing_value_2D, scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_unfrozen_vol_frac', missing_value_2D, scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_frozen_vol_frac',   missing_value_2D, scm_state%itt_out)
      else if (physics%Model%lsm == physics%Model%lsm_ruc) then
        call NetCDF_put_var(ncid, 'soil_T',                          physics%Sfcprop%tslb(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture',                   missing_value_2D, scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_unfrozen',          missing_value_2D, scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_total_vol_frac',    physics%Sfcprop%smois(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_unfrozen_vol_frac', physics%Sfcprop%sh2o(:,:), scm_state%itt_out)
        call NetCDF_put_var(ncid, 'soil_moisture_frozen_vol_frac',   physics%Sfcprop%keepsmfr(:,:), scm_state%itt_out)
      end if
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
  call NetCDF_put_var(ncid, 't2m', physics%Sfcprop%t2m(:), scm_state%itt_out)
  call NetCDF_put_var(ncid, 'q2m', physics%Sfcprop%q2m(:), scm_state%itt_out)
  call NetCDF_put_var(ncid, 'ustar', physics%Sfcprop%uustar(:), scm_state%itt_out)
  call NetCDF_put_var(ncid, 'tsfc', physics%Sfcprop%tsfc(:), scm_state%itt_out)
  
end subroutine output_append_sfcprop

subroutine output_append_interstitial_inst(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    call NetCDF_put_var(ncid, "tau_u",   physics%Interstitial%dusfc1(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "tau_v",   physics%Interstitial%dvsfc1(:), scm_state%itt_out)
    if (physics%model%imfdeepcnv >= 0 .or. physics%model%imfshalcnv >= 0) then
      call NetCDF_put_var(ncid, "upd_mf",  physics%Tbd%ud_mf(:,:), scm_state%itt_out)
    end if
    call NetCDF_put_var(ncid, "dwn_mf",  physics%Interstitial%dd_mf(:,:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "det_mf",  physics%Interstitial%dt_mf(:,:), scm_state%itt_out)
    
    call NetCDF_put_var(ncid, "sfc_up_lw_land",     physics%Interstitial%adjsfculw_land(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_lw_ice",      physics%Interstitial%adjsfculw_ice(:), scm_state%itt_out)
    call NetCDF_put_var(ncid, "sfc_up_lw_water",    physics%Interstitial%adjsfculw_water(:), scm_state%itt_out)
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

end subroutine output_append_interstitial_inst

subroutine output_append_interstitial_rad(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    call NetCDF_put_var(ncid, "rad_cloud_fraction", physics%Interstitial%clouds(:,:,1), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_cloud_lwp",      physics%Interstitial%clouds(:,:,2), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_eff_rad_ql",     physics%Interstitial%clouds(:,:,3), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_cloud_iwp",      physics%Interstitial%clouds(:,:,4), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_eff_rad_qi",     physics%Interstitial%clouds(:,:,5), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_cloud_rwp",      physics%Interstitial%clouds(:,:,6), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_eff_rad_qr",     physics%Interstitial%clouds(:,:,7), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_cloud_swp",      physics%Interstitial%clouds(:,:,8), scm_state%itt_rad)
    call NetCDF_put_var(ncid, "rad_eff_rad_qs",     physics%Interstitial%clouds(:,:,9), scm_state%itt_rad)

end subroutine output_append_interstitial_rad

subroutine output_append_radtend(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    if (physics%Model%lsswr) then
      call NetCDF_put_var(ncid, "sw_rad_heating_rate",  physics%Radtend%htrsw(:,:), scm_state%itt_swrad)
    end if
    if (physics%Model%lslwr) then
      call NetCDF_put_var(ncid, "lw_rad_heating_rate",  physics%Radtend%htrlw(:,:), scm_state%itt_lwrad)
    end if
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

subroutine output_append_diag_inst(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    integer :: i,j
    character(2) :: idx
    
    real(kind=dp), dimension(scm_state%n_cols) :: temp_1d
    real(kind=dp), dimension(scm_state%n_cols, scm_state%n_levels) :: temp_2d
    
    call NetCDF_put_var(ncid, "pwat",  physics%Diag%pwat(:), scm_state%itt_out)  !do not average (this variable is reset every physics timestep in GFS_MP_generic)
    
    call NetCDF_put_var(ncid, "sfc_dwn_sw",  physics%Diag%dswsfci(:), scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in dcyc2)
    call NetCDF_put_var(ncid, "sfc_up_sw",   physics%Diag%dswsfci(:) - physics%Diag%nswsfci(:), scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in dcyc2)
    call NetCDF_put_var(ncid, "sfc_net_sw",  physics%Diag%nswsfci(:), scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in dcyc2)
    call NetCDF_put_var(ncid, "sfc_dwn_lw",  physics%Diag%dlwsfci(:), scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in dcyc2)
    
    call NetCDF_put_var(ncid, "gflux",       physics%Diag%gfluxi(:), scm_state%itt_out)  !do not average (this variable is set from an interstitial variable in GFS_surface_generic_post)
    call NetCDF_put_var(ncid, "u10m",        physics%Diag%u10m(:),   scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in sfc_diag)
    call NetCDF_put_var(ncid, "v10m",        physics%Diag%v10m(:),   scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in sfc_diag)
    call NetCDF_put_var(ncid, "hpbl",        physics%Tbd%hpbl(:),    scm_state%itt_out)  !do not average (this variable is intent(out) every physics timestep in individual PBL schemes)
    
    !all auxilliary diagnostics will be output every timestep (logic will need to be added to implement time averaging [move this to output_append_diag_avg when that happens] -- including resetting in GFS_typedefs/phys_diag_zero)
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
end subroutine output_append_diag_inst

subroutine output_append_diag_avg(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
    
    real(kind=dp) :: inverse_n_diag, inverse_dt
    
    inverse_n_diag = 1.0/physics%Model%nszero
    inverse_dt = 1.0/scm_state%dt
    
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_longwave),          "dT_dt_lwrad",    scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_shortwave),         "dT_dt_swrad",    scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_pbl),               "dT_dt_pbl",      scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_dcnv),              "dT_dt_deepconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_scnv),              "dT_dt_shalconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_mp),                "dT_dt_micro",    scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_orographic_gwd),    "dT_dt_ogwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_rayleigh_damping),  "dT_dt_rayleigh", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_nonorographic_gwd), "dT_dt_cgwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_physics),           "dT_dt_phys",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_temperature,physics%Model%index_of_process_non_physics),       "dT_dt_nonphys",  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_pbl),            "dq_dt_pbl",       scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_dcnv),           "dq_dt_deepconv",  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_scnv),           "dq_dt_shalconv",  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_mp),             "dq_dt_micro",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_physics),        "dq_dt_phys",      scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntqv,physics%Model%index_of_process_non_physics),    "dq_dt_nonphys",   scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_pbl),            "doz_dt_pbl",      scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_prod_loss),      "doz_dt_prodloss", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_ozmix),          "doz_dt_oz",       scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_temp),           "doz_dt_T",        scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_overhead_ozone), "doz_dt_ovhd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_physics),        "doz_dt_phys",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(100+physics%Model%ntoz,physics%Model%index_of_process_non_physics),    "doz_dt_nonphys",  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_pbl),               "du_dt_pbl",      scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_orographic_gwd),    "du_dt_ogwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_dcnv),              "du_dt_deepconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_nonorographic_gwd), "du_dt_cgwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_rayleigh_damping),  "du_dt_rayleigh", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_scnv),              "du_dt_shalconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_physics),           "du_dt_phys",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_x_wind,physics%Model%index_of_process_non_physics),       "du_dt_nonphys",  scm_state%itt_diag, inverse_n_diag*inverse_dt)

    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_pbl),               "dv_dt_pbl",      scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_orographic_gwd),    "dv_dt_ogwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_dcnv),              "dv_dt_deepconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_nonorographic_gwd), "dv_dt_cgwd",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_rayleigh_damping),  "dv_dt_rayleigh", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_scnv),              "dv_dt_shalconv", scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_physics),           "dv_dt_phys",     scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call output_append_tendency(ncid, physics, physics%Model%dtidx(physics%Model%index_of_y_wind,physics%Model%index_of_process_non_physics),       "dv_dt_nonphys",  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    
    call NetCDF_put_var(ncid, "tprcp_accum",          physics%Diag%totprcpb(:), scm_state%itt_diag, inverse_n_diag)
    call NetCDF_put_var(ncid, "ice_accum",            physics%Diag%toticeb(:),  scm_state%itt_diag, inverse_n_diag)
    call NetCDF_put_var(ncid, "snow_accum",           physics%Diag%totsnwb(:),  scm_state%itt_diag, inverse_n_diag)
    call NetCDF_put_var(ncid, "graupel_accum",        physics%Diag%totgrpb(:),  scm_state%itt_diag, inverse_n_diag)
    call NetCDF_put_var(ncid, "conv_prcp_accum",      physics%Diag%cnvprcpb(:), scm_state%itt_diag, inverse_n_diag)
    call NetCDF_put_var(ncid, "tprcp_rate_accum",     physics%Diag%totprcpb(:), scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call NetCDF_put_var(ncid, "ice_rate_accum",       physics%Diag%toticeb(:),  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call NetCDF_put_var(ncid, "snow_rate_accum",      physics%Diag%totsnwb(:),  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call NetCDF_put_var(ncid, "graupel_rate_accum",   physics%Diag%totgrpb(:),  scm_state%itt_diag, inverse_n_diag*inverse_dt)
    call NetCDF_put_var(ncid, "conv_prcp_rate_accum", physics%Diag%cnvprcpb(:), scm_state%itt_diag, inverse_n_diag*inverse_dt)
    
end subroutine output_append_diag_avg

subroutine output_append_diag_rad(ncid, scm_state, physics)
    use scm_type_defs, only: scm_state_type, physics_type
    use NetCDF_put, only: NetCDF_put_var
    
    integer, intent(in) :: ncid
    type(scm_state_type), intent(in) :: scm_state
    type(physics_type), intent(in) :: physics
      
    call NetCDF_put_var(ncid, 'max_cloud_fraction', physics%Diag%cldfra2d(:),     scm_state%itt_rad)
    call NetCDF_put_var(ncid, 'toa_total_albedo',   physics%Diag%total_albedo(:), scm_state%itt_rad)
    call NetCDF_put_var(ncid, 'vert_int_lwp_mp',    physics%Diag%lwp_ex(:),       scm_state%itt_rad)
    call NetCDF_put_var(ncid, 'vert_int_iwp_mp',    physics%Diag%iwp_ex(:),       scm_state%itt_rad)
    call NetCDF_put_var(ncid, 'vert_int_lwp_cf',    physics%Diag%lwp_fc(:),       scm_state%itt_rad)
    call NetCDF_put_var(ncid, 'vert_int_iwp_cf',    physics%Diag%iwp_fc(:),       scm_state%itt_rad)
    
end subroutine output_append_diag_rad

subroutine output_init_tendency(ncid, label, desc, units, hor_id, vert_id, time_id, idtend)
  use NetCDF_def, only: NetCDF_def_var
  
  integer, intent(in) :: ncid, hor_id, vert_id, time_id, idtend
  character(len=*), intent(in) :: label, desc, units
  
  integer :: dummy_id
  
  if(idtend >=1) then
    call NetCDF_def_var(ncid, label, NF90_FLOAT, desc, units, dummy_id, (/ hor_id, vert_id, time_id /))
  end if
  
end subroutine output_init_tendency

subroutine output_append_tendency(ncid, physics, idtend, label, itt, mult_const)
  use scm_type_defs, only: scm_state_type, physics_type
  use NetCDF_put, only: NetCDF_put_var
  
  integer, intent(in) :: ncid, idtend, itt
  type(physics_type), intent(in) :: physics
  character(len=*), intent(in) :: label
  real(kind=dp), intent(in), optional :: mult_const
  
  if(idtend>=1) then
    call NetCDF_put_var(ncid, label, physics%Diag%dtend(:,:,idtend), itt, mult_const)
  endif
  
end subroutine output_append_tendency

!> @}
!> @}
end module scm_output
