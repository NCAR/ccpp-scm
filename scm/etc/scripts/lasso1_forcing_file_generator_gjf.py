#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import forcing_file_common as ffc
import scipy.interpolate
import gmtb_scm_plotting_routines as gspr
import datetime
import sys

reload(ffc)

#read in raw input file

nc_fid_ls  = Dataset("../../data/raw_case_input/LASSO_20160518/input_ls_forcing.nc", 'r')
nc_fid_sfc = Dataset("../../data/raw_case_input/LASSO_20160518/input_sfc_forcing.nc", 'r')
nc_fid_wrf = Dataset("../../data/raw_case_input/LASSO_20160518/wrfinput_d01.nc", 'r')

#data listed in the case instructions
lat = 36.6000 #degrees north
lon = -97.4833

time_chars = nc_fid_ls.variables['Times']
time_strings = [''.join(x) for x in time_chars]
data_date = [datetime.datetime.strptime(x, '%Y-%m-%d_%H:%M:%S') for x in time_strings]

year=str(data_date[0].year)
mon=format(data_date[0].month,'02d')
day=format(data_date[0].day,'02d')
hour=format(data_date[0].hour,'02d')

file_out="../../data/processed_case_input/LASSO_"+year+mon+day+hour+".nc"
description="GMTB SCM forcing file for the LASSO case for "+mon+"/"+day+"/"+year+" "+hour+"Z"

#height in LS forcing file
z_ls = nc_fid_ls.variables['Z_LS']

#height in WRF initial conditions file
ph = nc_fid_wrf.variables['PH']
phb = nc_fid_wrf.variables['PHB']
ph = ph[:] + phb[:]
ph = np.mean(ph, axis=(0,2,3))
z_wrf_stag = ph/9.81
z_wrf = (z_wrf_stag[1:z_wrf_stag.size]+z_wrf_stag[0:z_wrf_stag.size-1])*0.5

#pressure levels of initial conditions
pert_pressure = nc_fid_wrf.variables['P']
pert_pressure = np.mean(pert_pressure, axis=(0,2,3))
base_state_pressure = nc_fid_wrf.variables['PB']
base_state_pressure = np.mean(base_state_pressure, axis=(0,2,3))
levels = pert_pressure + base_state_pressure

# Open ozone file and get ozone data
nc_fid_o3 = Dataset("../../data/raw_case_input/mid_lat_summer_std.nc", 'r')

oz_pres = nc_fid_o3.variables['pressure']
oz_data = nc_fid_o3.variables['o3']

oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
o3 = oz_f(levels[:])
nc_fid_o3.close()

#get initial theta

theta_pert = nc_fid_wrf.variables['T']
theta_pert = np.mean(theta_pert, axis=(0,2,3))

theta = theta_pert + 300.0 #instructions say to add 300 K to theta perturbation
T = ffc.theta_to_T(theta, levels)

#get initial water species mixing ratios and convert to specific humidity
qv_mr = nc_fid_wrf.variables['QVAPOR']
qv_mr = np.mean(qv_mr, axis=(0,2,3))
qv = qv_mr/(1.0 + qv_mr) #convert to specific humidity from mixing ratio

ql_mr = nc_fid_wrf.variables['QCLOUD']
ql_mr = np.mean(ql_mr, axis=(0,2,3))
ql = ql_mr/(1.0 + ql_mr) #convert to specific humidity from mixing ratio

qi_mr = nc_fid_wrf.variables['QICE']
qi_mr = np.mean(qi_mr, axis=(0,2,3))
qi = qi_mr/(1.0 + qi_mr) #convert to specific humidity from mixing ratio

#calc qt and theta_il
qt = qv + ql + qi

theta_il = theta - ((ffc.L_v/ffc.c_p)*ql_mr + (ffc.L_s/ffc.c_p)*qi_mr)

#read in initial u and v
u_wind = nc_fid_wrf.variables['U']
u_wind = np.mean(u_wind, axis=(0,2,3))

v_wind = nc_fid_wrf.variables['V']
v_wind = np.mean(v_wind, axis=(0,2,3))

#no TKE data -- init to 0
tke = np.zeros((levels.size),dtype=float)

#read in forcing from sfc and LS file

t_surf = nc_fid_sfc.variables['PRE_TSK'][:]

sh_flux_sfc = nc_fid_sfc.variables['PRE_SH_FLX'][:]
sh_flux_sfc = sh_flux_sfc*ffc.R_dry*t_surf/(ffc.c_p*levels[0]) #convert to K m/s
lh_flux_sfc = nc_fid_sfc.variables['PRE_LH_FLX'][:]
lh_flux_sfc = lh_flux_sfc*ffc.R_dry*t_surf/(ffc.L_v*levels[0]) #convert to kg/kg m/s

#read in forcing on forcing grid and interpolate to initial conditions grid
w_ls_force_grid = nc_fid_ls.variables['W_LS']
u_ls_force_grid = nc_fid_ls.variables['U_LS']
v_ls_force_grid = nc_fid_ls.variables['V_LS']
theta_nudge_force_grid = nc_fid_ls.variables['TH_RLX']
qv_nudge_force_grid = nc_fid_ls.variables['QV_RLX']
h_advec_thil_force_grid = nc_fid_ls.variables['TH_ADV']
h_advec_qt_force_grid = nc_fid_ls.variables['QV_ADV']
w_ls = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
u_ls = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
v_ls = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
omega = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
theta_nudge = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
T_nudge = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float) #not needed
qv_nudge = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
h_advec_thil = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
h_advec_qt = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float)
v_advec_thil = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float) #not needed if w_ls is specified
v_advec_qt = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float) #not needed if w_ls is specified
rad_heating = np.zeros((z_wrf.size, w_ls_force_grid.shape[0]),dtype=float) #not needed if active rad
time = np.zeros((w_ls_force_grid.shape[0]),dtype=float)
for t in range(w_ls_force_grid.shape[0]):
    time[t] = (data_date[t] - data_date[0]).total_seconds()
    w_ls_f = scipy.interpolate.interp1d(z_ls[t,:], w_ls_force_grid[t,:])
    w_ls[:,t] = w_ls_f(z_wrf)
    omega[:,t] = ffc.w_to_omega(w_ls[:,t], levels, T)
    u_ls_f = scipy.interpolate.interp1d(z_ls[t,:], u_ls_force_grid[t,:])
    u_ls[:,t] = u_ls_f(z_wrf)
    v_ls_f = scipy.interpolate.interp1d(z_ls[t,:], v_ls_force_grid[t,:])
    v_ls[:,t] = v_ls_f(z_wrf)
    theta_nudge_f = scipy.interpolate.interp1d(z_ls[t,:], theta_nudge_force_grid[t,:])
    theta_nudge[:,t] = theta_nudge_f(z_wrf)
    qv_nudge_f = scipy.interpolate.interp1d(z_ls[t,:], qv_nudge_force_grid[t,:])
    qv_nudge[:,t] = qv_nudge_f(z_wrf)
    h_advec_thil_f = scipy.interpolate.interp1d(z_ls[t,:], h_advec_thil_force_grid[t,:])
    h_advec_thil[:,t] = h_advec_thil_f(z_wrf)
    h_advec_qt_f = scipy.interpolate.interp1d(z_ls[t,:], h_advec_qt_force_grid[t,:])
    h_advec_qt[:,t] = h_advec_qt_f(z_wrf)

#writefile_fid = Dataset('../processed_case_input/LASSO_2016051812.nc', 'w', format='NETCDF4')
#writefile_fid = Dataset('LASSO_2016061012.nc', 'w', format='NETCDF4')
writefile_fid = Dataset(file_out, 'w', format='NETCDF4')
writefile_fid.description = description #"GMTB SCM forcing file for the LASSO case for 06/10/2016 12Z"

#create groups for scalars, intitialization, and forcing

writefile_scalar_grp = writefile_fid.createGroup("scalars")
writefile_initial_grp = writefile_fid.createGroup("initial")
writefile_forcing_grp = writefile_fid.createGroup("forcing")

#create dimensions and write them out

writefile_time_dim = writefile_fid.createDimension('time', None)
writefile_time_var = writefile_fid.createVariable('time', 'f4', ('time',))
writefile_time_var[:] = time
writefile_time_var.units = 's'
writefile_time_var.description = 'elapsed time since the beginning of the simulation'

writefile_levels_dim = writefile_fid.createDimension('levels', None)
writefile_levels_var = writefile_fid.createVariable('levels', 'f4', ('levels',))
writefile_levels_var[:] = levels
writefile_levels_var.units = 'Pa'
writefile_levels_var.description = 'pressure levels'

#create variables and write them out

#scalar group

writefile_init_year_var = writefile_scalar_grp.createVariable('init_year', 'i4')
writefile_init_year_var[:] = data_date[0].year
writefile_init_year_var.units = 'years'
writefile_init_year_var.description = 'year at time of initial values'

writefile_init_month_var = writefile_scalar_grp.createVariable('init_month', 'i4')
writefile_init_month_var[:] = data_date[0].month
writefile_init_month_var.units = 'months'
writefile_init_month_var.description = 'month at time of initial values'

writefile_init_day_var = writefile_scalar_grp.createVariable('init_day', 'i4')
writefile_init_day_var[:] = data_date[0].day
writefile_init_day_var.units = 'days'
writefile_init_day_var.description = 'day at time of initial values'

writefile_init_hour_var = writefile_scalar_grp.createVariable('init_hour', 'i4')
writefile_init_hour_var[:] = data_date[0].hour
writefile_init_hour_var.units = 'hours'
writefile_init_hour_var.description = 'hour at time of initial values'

writefile_init_min_var = writefile_scalar_grp.createVariable('init_minute', 'i4')
writefile_init_min_var[:] = data_date[0].minute
writefile_init_min_var.units = 'minutes'
writefile_init_min_var.description = 'minute at time of initial values'

writefile_init_second_var = writefile_scalar_grp.createVariable('init_second', 'i4')
writefile_init_second_var[:] = data_date[0].second
writefile_init_second_var.units = 'seconds'
writefile_init_second_var.description = 'second at time of initial values'

#initial group

writefile_height_var = writefile_initial_grp.createVariable('height', 'f4', ('levels',))
writefile_height_var[:] = z_wrf
writefile_height_var.units = 'm'
writefile_height_var.description = 'physical height at pressure levels'

writefile_thetail_var = writefile_initial_grp.createVariable('thetail', 'f4', ('levels',))
writefile_thetail_var[:] = theta_il
writefile_thetail_var.units = 'K'
writefile_thetail_var.description = 'initial profile of ice-liquid water potential temperature'

writefile_qt_var = writefile_initial_grp.createVariable('qt', 'f4', ('levels',))
writefile_qt_var[:] = qt
writefile_qt_var.units = 'kg kg^-1'
writefile_qt_var.description = 'initial profile of total water specific humidity'

writefile_ql_var = writefile_initial_grp.createVariable('ql', 'f4', ('levels',))
writefile_ql_var[:] = ql
writefile_ql_var.units = 'kg kg^-1'
writefile_ql_var.description = 'initial profile of liquid water specific humidity'

writefile_qi_var = writefile_initial_grp.createVariable('qi', 'f4', ('levels',))
writefile_qi_var[:] = qi
writefile_qi_var.units = 'kg kg^-1'
writefile_qi_var.description = 'initial profile of ice water specific humidity'

writefile_u_var = writefile_initial_grp.createVariable('u', 'f4', ('levels',))
writefile_u_var[:] = u_wind
writefile_u_var.units = 'm s^-1'
writefile_u_var.description = 'initial profile of E-W horizontal wind'

writefile_v_var = writefile_initial_grp.createVariable('v', 'f4', ('levels',))
writefile_v_var[:] = v_wind
writefile_v_var.units = 'm s^-1'
writefile_v_var.description = 'initial profile of N-S horizontal wind'

writefile_tke_var = writefile_initial_grp.createVariable('tke', 'f4', ('levels',))
writefile_tke_var[:] = tke
writefile_tke_var.units = 'm^2 s^-2'
writefile_tke_var.description = 'initial profile of turbulence kinetic energy'

writefile_ozone_var = writefile_initial_grp.createVariable('ozone', 'f4', ('levels',))
writefile_ozone_var[:] = o3
writefile_ozone_var.units = 'kg kg^-1'
writefile_ozone_var.description = 'initial profile of ozone mass mixing ratio'

#forcing group

writefile_lat_var = writefile_forcing_grp.createVariable('lat', 'f4', ('time',))
writefile_lat_var[:] = lat*np.ones((time.size),dtype=float)
writefile_lat_var.units = 'degrees N'
writefile_lat_var.description = 'latitude of column'

writefile_lon_var = writefile_forcing_grp.createVariable('lon', 'f4', ('time',))
writefile_lon_var[:] = lon*np.ones((time.size),dtype=float)
writefile_lon_var.units = 'degrees E'
writefile_lon_var.description = 'longitude of column'

writefile_p_surf_var = writefile_forcing_grp.createVariable('p_surf', 'f4', ('time',))
writefile_p_surf_var[:] = levels[0]*np.ones((time.size),dtype=float)
writefile_p_surf_var.units = 'Pa'
writefile_p_surf_var.description = 'surface pressure'

writefile_T_surf_var = writefile_forcing_grp.createVariable('T_surf', 'f4', ('time',))
writefile_T_surf_var[:] = t_surf
writefile_T_surf_var.units = 'K'
writefile_T_surf_var.description = 'surface absolute temperature'

writefile_sh_flux_sfc_var = writefile_forcing_grp.createVariable('sh_flux_sfc', 'f4', ('time',))
writefile_sh_flux_sfc_var[:] = sh_flux_sfc
writefile_sh_flux_sfc_var.units = 'K m s^-1'
writefile_sh_flux_sfc_var.description = 'surface sensible heat flux'

writefile_lh_flux_sfc_var = writefile_forcing_grp.createVariable('lh_flux_sfc', 'f4', ('time',))
writefile_lh_flux_sfc_var[:] = lh_flux_sfc
writefile_lh_flux_sfc_var.units = 'kg kg^-1 m s^-1'
writefile_lh_flux_sfc_var.description = 'surface latent heat flux'

writefile_w_ls_var = writefile_forcing_grp.createVariable('w_ls', 'f4', ('levels','time',))
writefile_w_ls_var[:] = w_ls
writefile_w_ls_var.units = 'm s^-1'
writefile_w_ls_var.description = 'large scale vertical velocity'

writefile_omega_var = writefile_forcing_grp.createVariable('omega', 'f4', ('levels','time',))
writefile_omega_var[:] = omega
writefile_omega_var.units = 'Pa s^-1'
writefile_omega_var.description = 'large scale pressure vertical velocity'

writefile_u_g_var = writefile_forcing_grp.createVariable('u_g', 'f4', ('levels','time',))
writefile_u_g_var[:] = u_ls
writefile_u_g_var.units = 'm s^-1'
writefile_u_g_var.description = 'large scale geostrophic E-W wind'

writefile_v_g_var = writefile_forcing_grp.createVariable('v_g', 'f4', ('levels','time',))
writefile_v_g_var[:] = v_ls
writefile_v_g_var.units = 'm s^-1'
writefile_v_g_var.description = 'large scale geostrophic N-S wind'

writefile_u_nudge_var = writefile_forcing_grp.createVariable('u_nudge', 'f4', ('levels','time',))
writefile_u_nudge_var[:] = u_ls
writefile_u_nudge_var.units = 'm s^-1'
writefile_u_nudge_var.description = 'E-W wind to nudge toward'

writefile_v_nudge_var = writefile_forcing_grp.createVariable('v_nudge', 'f4', ('levels','time',))
writefile_v_nudge_var[:] = v_ls
writefile_v_nudge_var.units = 'm s^-1'
writefile_v_nudge_var.description = 'N-S wind to nudge toward'

writefile_T_nudge_var = writefile_forcing_grp.createVariable('T_nudge', 'f4', ('levels','time',))
writefile_T_nudge_var[:] = T_nudge
writefile_T_nudge_var.units = 'K'
writefile_T_nudge_var.description = 'absolute temperature to nudge toward'

writefile_thil_nudge_var = writefile_forcing_grp.createVariable('thil_nudge', 'f4', ('levels','time',))
writefile_thil_nudge_var[:] = theta_nudge
writefile_thil_nudge_var.units = 'K'
writefile_thil_nudge_var.description = 'potential temperature to nudge toward'

writefile_qt_nudge_var = writefile_forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
writefile_qt_nudge_var[:] = qv_nudge
writefile_qt_nudge_var.units = 'kg kg^-1'
writefile_qt_nudge_var.description = 'q_t to nudge toward'

writefile_rad_heating_var = writefile_forcing_grp.createVariable('dT_dt_rad', 'f4', ('levels','time',))
writefile_rad_heating_var[:] = rad_heating
writefile_rad_heating_var.units = 'K s^-1'
writefile_rad_heating_var.description = 'prescribed radiative heating rate'

writefile_h_advec_thil_var = writefile_forcing_grp.createVariable('h_advec_thetail', 'f4', ('levels','time',))
writefile_h_advec_thil_var[:] = h_advec_thil
writefile_h_advec_thil_var.units = 'K s^-1'
writefile_h_advec_thil_var.description = 'prescribed theta_il tendency due to horizontal advection'

writefile_v_advec_thil_var = writefile_forcing_grp.createVariable('v_advec_thetail', 'f4', ('levels','time',))
writefile_v_advec_thil_var[:] = v_advec_thil
writefile_v_advec_thil_var.units = 'K s^-1'
writefile_v_advec_thil_var.description = 'prescribed theta_il tendency due to vertical advection'

writefile_h_advec_qt_var = writefile_forcing_grp.createVariable('h_advec_qt', 'f4', ('levels','time',))
writefile_h_advec_qt_var[:] = h_advec_qt
writefile_h_advec_qt_var.units = 'kg kg^-1 s^-1'
writefile_h_advec_qt_var.description = 'prescribed q_t tendency due to horizontal advection'

writefile_v_advec_qt_var = writefile_forcing_grp.createVariable('v_advec_qt', 'f4', ('levels','time',))
writefile_v_advec_qt_var[:] = v_advec_qt
writefile_v_advec_qt_var.units = 'kg kg^-1 s^-1'
writefile_v_advec_qt_var.description = 'prescribed q_t tendency due to vertical advection'

#close processed input file
writefile_fid.close()
nc_fid_ls.close()
nc_fid_sfc.close()
nc_fid_wrf.close()
sys.exit()
