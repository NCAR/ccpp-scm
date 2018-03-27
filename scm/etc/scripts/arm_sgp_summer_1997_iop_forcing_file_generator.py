#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import forcing_file_common as ffc
import scipy.interpolate
import gmtb_scm_plotting_routines as gspr

reload(ffc)

#read in raw ASTEX input file

nc_fid = Dataset("../raw_case_input/sgp3hIOPsndgBasedV2.0_ConstrVarAnaX1.c1.19970618.000000.cdf", 'r')

#data listed in the case instructions
z_sfc = 315.0 #m above sea level
lat = 36.6000 #degrees north
lon = -97.4833

#ncdump to look at raw input file
#nc_attrs, nc_dims, nc_vars = ncdump(nc_fid, False)

#netCDF how-to
#set attributes
#file_id.setncattr(file_id.variables['variable'].ncattr(), nc_fid.variables['time'].getncattr(ncattr))
#get attributes
#file_id.variables['variable'].getncattr(index)

#get raw input variables

# day = nc_fid.variables['day'][:]
# hour = nc_fid.variables['hour'][:]
# for t in range(day.size):
#     #find time index corresponding to Jan 19, 2006 at 03Z
#     if day[t] == 19 and hour[t] == 3:
#         start_t_index = t
#         break

start_t_index = 0

time = nc_fid.variables['time_offset'][:] #number of seconds since 00Z on 6/18/1997 (file time starts at 23Z and 03 seconds)
#subtract the initial time_offset from all values to get elapsed time since the start of the simulation
time = time - time[start_t_index]
levels = nc_fid.variables['lev'][:] #pressure levels in mb
#convert levels to Pa
levels = 100.0*levels
#need to flip the column top-to-bottom
levels = np.flipud(levels)

#height = nc_fid.variables['alt'][:]
#lat = nc_fid.variables['lat'][:] #degrees south
#convert latitutde to degrees north
#lat = -1*lat
#lon = nc_fid.variables['lon'][:] #degrees east

T_abs = nc_fid.variables['Temp'][:,:,0,0] #absolute temperature (time, lev)
#swap the time and levels axes, plus flip the column top-to-bottom
T_abs = np.flipud(np.swapaxes(T_abs, 0, 1))

#calculate theta_il from absolute temperature (assuming no condensate)
thetal = np.zeros((levels.size,time.size),dtype=float)
for t in range(time.size):
    thetal[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*T_abs[:,t]

qv = nc_fid.variables['H2O_Mixing_Ratio'][:,:,0,0] #water vapor mixing ratio in g/kg (time, lev)
qt = np.zeros((levels.size,time.size),dtype=float)
qt_mr = 1.0E-3*np.flipud(np.swapaxes(qv, 0, 1)) #swap the time and levels axis, flip the column top-to-bottom, convert to kg/kg
qt = qt_mr/(1.0 + qt_mr) #convert to specific humidity from mixing ratio


#ql and tke are not specified; set to zero
ql = np.zeros((levels.size,time.size),dtype=float)
qi = np.zeros((levels.size,time.size),dtype=float)
tke = np.zeros((levels.size,time.size),dtype=float)

# ozone_mmr = nc_fid.variables['o3mmr'][:]

u_wind = nc_fid.variables['u_wind'][:,:,0,0]
u_wind = np.flipud(np.swapaxes(u_wind, 0, 1)) #swap the time and levels axis and flip column top-to-bottom


v_wind = nc_fid.variables['v_wind'][:,:,0,0]
v_wind = np.flipud(np.swapaxes(v_wind, 0, 1)) #swap the time and levels axis and flip the column top-to-bottom

w_sub = np.zeros((levels.size,time.size),dtype=float)
omega = nc_fid.variables['omega'][:,:,0,0] #vertical pressure velocity in mb/hr
omega = omega*100.0/3600.0 #convert to Pa/s
omega = np.flipud(np.swapaxes(omega, 0, 1)) #swap the time and levels axis and flip the column top-to-bottom
#convert to w
for t in range(time.size):
    w_sub[:,t] = ffc.omega_to_w(omega[:,t],levels,T_abs[:,t])

T_surf = nc_fid.variables['Ts_Air'][:,0,0]
T_surf = T_surf + 273.15 #convert to K

p_surf = nc_fid.variables['Area_Mean_Ps'][:,0,0] #mb
p_surf = p_surf*100.0 #convert to Pa

h_advec_thil = np.zeros((levels.size,time.size),dtype=float)
h_advec_T = nc_fid.variables['Horizontal_Temp_Advec'][:,:,0,0] #K/hr
h_advec_T = np.flipud(np.swapaxes(h_advec_T, 0, 1))/3600.0 #swap the time and levels axis, flip the column top-to-bottom, convert to K/s
for t in range(time.size):
    h_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*h_advec_T[:,t] #convert to potential temperature

v_advec_thil = np.zeros((levels.size,time.size),dtype=float)
v_advec_s = nc_fid.variables['Vertical_s_Advec'][:,:,0,0] #vertical advection of s includes vertical advcetion of T and adiabatic term
v_advec_T = np.flipud(np.swapaxes(v_advec_s, 0, 1))/(3600.0*ffc.c_p) #swap the time and levels axis, flip the column top-to-bottom, divide dry static energy by c_p to get T term, convert to K/s
#v_advec_T = nc_fid.variables['Vertical_T_Advec'][:,:,0,0] #K/hr
#v_advec_T = np.flipud(np.swapaxes(v_advec_T, 0, 1))/3600.0 #swap the time and levels axis, flip the column top-to-bottom, convert to K/s
for t in range(time.size):
    v_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*v_advec_T[:,t] #convert to potential temperature

dT_dt = nc_fid.variables['DT_dt'][:,:,0,0] #K/hr
dT_dt = np.flipud(np.swapaxes(dT_dt, 0, 1))/3600.0 #swap the time and levels axis, flip the column top-to-bottom, convert to K/s

# h_advec_T = h_advec_T*86400.0
# v_advec_T = v_advec_T*86400.0
# dT_dt = dT_dt*86400.0
# gspr.contour_plot_firl(time, levels, h_advec_T, np.min(h_advec_T), np.max(h_advec_T), 'h advec T', 'time', 'pressure', 'h_advec_T.eps', y_inverted = True)
# gspr.contour_plot_firl(time, levels, v_advec_T, np.min(v_advec_T), np.max(v_advec_T), 'v advec T', 'time', 'pressure', 'v_advec_T.eps', y_inverted = True)
# gspr.contour_plot_firl(time, levels, dT_dt, np.min(dT_dt), np.max(dT_dt), 'total T tend', 'time', 'pressure', 'dT_dt.eps', y_inverted = True)

h_advec_qt = np.zeros((levels.size,time.size),dtype=float)
h_advec_qt = nc_fid.variables['Horizontal_q_Advec'][:,:,0,0] #g/kg/hr
h_advec_qt = np.flipud(np.swapaxes(h_advec_qt, 0, 1))*1.0E-3/3600.0 #swap the time and levels axis, flip the column from top-to-bottom, convert to kg/kg/s
h_advec_qt = h_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

v_advec_qt = np.zeros((levels.size,time.size),dtype=float)
v_advec_qt = nc_fid.variables['Vertical_q_Advec'][:,:,0,0] #g/kg/hr
v_advec_qt = np.flipud(np.swapaxes(v_advec_qt, 0, 1))*1.0E-3/3600.0 #swap the time and levels axis, flip the column from top-to-bottom, convert to kg/kg/s
v_advec_qt = v_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

dq_dt = nc_fid.variables['Dq_dt'][:,:,0,0]
dq_dt = np.flipud(np.swapaxes(dq_dt, 0, 1))*1.0E-3/3600.0 #swap the time and levels axis, flip the column from top-to-bottom, convert to K/s

sh_flux_sfc = nc_fid.variables['SH'][:,0,0] #W/m^2
sh_flux_sfc = sh_flux_sfc*ffc.R_dry*T_surf/(ffc.c_p*p_surf) #convert to Km/s
lh_flux_sfc = nc_fid.variables['LH'][:,0,0] #W/m^2
lh_flux_sfc = lh_flux_sfc*ffc.R_dry*T_surf/(ffc.L_v*p_surf) #convert to Km/s

# h_advec_qt = h_advec_qt*86400.0
# v_advec_qt = v_advec_qt*86400.0
# dq_dt = dq_dt*86400.0
# diff = dq_dt - (h_advec_qt + v_advec_qt)
# gspr.contour_plot_firl(time, levels, h_advec_qt, np.min(h_advec_qt), np.max(h_advec_qt), 'h advec q', 'time', 'pressure', 'h_advec_q.eps', y_inverted = True)
# gspr.contour_plot_firl(time, levels, v_advec_qt, np.min(v_advec_qt), np.max(v_advec_qt), 'v advec q', 'time', 'pressure', 'v_advec_q.eps', y_inverted = True)
# gspr.contour_plot_firl(time, levels, dq_dt, np.min(dq_dt), np.max(dq_dt), 'total q tend', 'time', 'pressure', 'dq_dt.eps', y_inverted = True)
# gspr.contour_plot_firl(time, levels, diff, np.min(diff), np.max(diff), 'total q tend diff', 'time', 'pressure', 'dq_dt_diff.eps', y_inverted = True)

phi_sfc = nc_fid.variables['phis'][0,0]
#z_sfc = nc_fid.variables['alt'][:]

height = ffc.get_height_from_pres(T_abs[:,0],levels,z_sfc)

print levels, height

#the following variables are not in the astex forcing file, but are included in other cases
rad_heating = np.zeros((levels.size,time.size),dtype=float)
u_g = np.zeros((levels.size,time.size),dtype=float)
v_g = np.zeros((levels.size,time.size),dtype=float)

# Open ozone file
nc_fid_o3 = Dataset("../raw_case_input/mid_lat_summer_std.nc", 'r')

oz_pres = nc_fid_o3.variables['pressure']
oz_data = nc_fid_o3.variables['o3']

oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
o3 = oz_f(levels[1:])


#open processed input file for writing

writefile_fid = Dataset('../processed_case_input/arm_sgp_summer_1997.nc', 'w', format='NETCDF4')
writefile_fid.description = "GMTB SCM forcing file for the ARM SGP Summer of 1997 case"

#create groups for scalars, intitialization, and forcing

writefile_scalar_grp = writefile_fid.createGroup("scalars")
writefile_initial_grp = writefile_fid.createGroup("initial")
writefile_forcing_grp = writefile_fid.createGroup("forcing")

#create dimensions and write them out

writefile_time_dim = writefile_fid.createDimension('time', None)
writefile_time_var = writefile_fid.createVariable('time', 'f4', ('time',))
writefile_time_var[:] = time[start_t_index:]
writefile_time_var.units = 's'
writefile_time_var.description = 'elapsed time since the beginning of the simulation'

writefile_levels_dim = writefile_fid.createDimension('levels', None)
writefile_levels_var = writefile_fid.createVariable('levels', 'f4', ('levels',))
writefile_levels_var[:] = levels
writefile_levels_var.units = 'Pa'
writefile_levels_var.description = 'pressure levels'

#create variables and write them out

#scalar group

#initial group

writefile_height_var = writefile_initial_grp.createVariable('height', 'f4', ('levels',))
writefile_height_var[:] = height
writefile_height_var.units = 'm'
writefile_height_var.description = 'physical height at pressure levels'

writefile_thetail_var = writefile_initial_grp.createVariable('thetail', 'f4', ('levels',))
writefile_thetail_var[:] = thetal[:,start_t_index]
writefile_thetail_var.units = 'K'
writefile_thetail_var.description = 'initial profile of ice-liquid water potential temperature'

writefile_qt_var = writefile_initial_grp.createVariable('qt', 'f4', ('levels',))
writefile_qt_var[:] = qt[:,start_t_index]
writefile_qt_var.units = 'kg kg^-1'
writefile_qt_var.description = 'initial profile of total water specific humidity'

writefile_ql_var = writefile_initial_grp.createVariable('ql', 'f4', ('levels',))
writefile_ql_var[:] = ql[:,start_t_index]
writefile_ql_var.units = 'kg kg^-1'
writefile_ql_var.description = 'initial profile of liquid water specific humidity'

writefile_qi_var = writefile_initial_grp.createVariable('qi', 'f4', ('levels',))
writefile_qi_var[:] = qi[:,start_t_index]
writefile_qi_var.units = 'kg kg^-1'
writefile_qi_var.description = 'initial profile of ice water specific humidity'

writefile_u_var = writefile_initial_grp.createVariable('u', 'f4', ('levels',))
writefile_u_var[:] = u_wind[:,start_t_index]
writefile_u_var.units = 'm s^-1'
writefile_u_var.description = 'initial profile of E-W horizontal wind'

writefile_v_var = writefile_initial_grp.createVariable('v', 'f4', ('levels',))
writefile_v_var[:] = v_wind[:,start_t_index]
writefile_v_var.units = 'm s^-1'
writefile_v_var.description = 'initial profile of N-S horizontal wind'

writefile_tke_var = writefile_initial_grp.createVariable('tke', 'f4', ('levels',))
writefile_tke_var[:] = tke[:,start_t_index]
writefile_tke_var.units = 'm^2 s^-2'
writefile_tke_var.description = 'initial profile of turbulence kinetic energy'

writefile_ozone_var = writefile_initial_grp.createVariable('ozone', 'f4', ('levels',))
writefile_ozone_var[:] = o3
writefile_ozone_var.units = 'kg kg^-1'
writefile_ozone_var.description = 'initial profile of ozone mass mixing ratio'

#forcing group

writefile_lat_var = writefile_forcing_grp.createVariable('lat', 'f4', ('time',))
writefile_lat_var[:] = lat*np.ones((time.size),dtype=float)[start_t_index:]
writefile_lat_var.units = 'degrees N'
writefile_lat_var.description = 'latitude of column'

writefile_lon_var = writefile_forcing_grp.createVariable('lon', 'f4', ('time',))
writefile_lon_var[:] = lon*np.ones((time.size),dtype=float)[start_t_index:]
writefile_lon_var.units = 'degrees E'
writefile_lon_var.description = 'longitude of column'

writefile_p_surf_var = writefile_forcing_grp.createVariable('p_surf', 'f4', ('time',))
writefile_p_surf_var[:] = p_surf[start_t_index:]
writefile_p_surf_var.units = 'Pa'
writefile_p_surf_var.description = 'surface pressure'

writefile_T_surf_var = writefile_forcing_grp.createVariable('T_surf', 'f4', ('time',))
writefile_T_surf_var[:] = T_surf[start_t_index:]
writefile_T_surf_var.units = 'K'
writefile_T_surf_var.description = 'surface absolute temperature'

writefile_sh_flux_sfc_var = writefile_forcing_grp.createVariable('sh_flux_sfc', 'f4', ('time',))
writefile_sh_flux_sfc_var[:] = sh_flux_sfc[start_t_index:]
writefile_sh_flux_sfc_var.units = 'K m s^-1'
writefile_sh_flux_sfc_var.description = 'surface sensible heat flux'

writefile_lh_flux_sfc_var = writefile_forcing_grp.createVariable('lh_flux_sfc', 'f4', ('time',))
writefile_lh_flux_sfc_var[:] = lh_flux_sfc[start_t_index:]
writefile_lh_flux_sfc_var.units = 'kg kg^-1 m s^-1'
writefile_lh_flux_sfc_var.description = 'surface latent heat flux'

writefile_w_ls_var = writefile_forcing_grp.createVariable('w_ls', 'f4', ('levels','time',))
writefile_w_ls_var[:] = w_sub[:,start_t_index:]
writefile_w_ls_var.units = 'm s^-1'
writefile_w_ls_var.description = 'large scale vertical velocity'

writefile_omega_var = writefile_forcing_grp.createVariable('omega', 'f4', ('levels','time',))
writefile_omega_var[:] = omega[:,start_t_index:]
writefile_omega_var.units = 'Pa s^-1'
writefile_omega_var.description = 'large scale pressure vertical velocity'

writefile_u_g_var = writefile_forcing_grp.createVariable('u_g', 'f4', ('levels','time',))
writefile_u_g_var[:] = u_g[:,start_t_index:]
writefile_u_g_var.units = 'm s^-1'
writefile_u_g_var.description = 'large scale geostrophic E-W wind'

writefile_v_g_var = writefile_forcing_grp.createVariable('v_g', 'f4', ('levels','time',))
writefile_v_g_var[:] = v_g[:,start_t_index:]
writefile_v_g_var.units = 'm s^-1'
writefile_v_g_var.description = 'large scale geostrophic N-S wind'

writefile_u_nudge_var = writefile_forcing_grp.createVariable('u_nudge', 'f4', ('levels','time',))
writefile_u_nudge_var[:] = u_wind[:,start_t_index:]
writefile_u_nudge_var.units = 'm s^-1'
writefile_u_nudge_var.description = 'E-W wind to nudge toward'

writefile_v_nudge_var = writefile_forcing_grp.createVariable('v_nudge', 'f4', ('levels','time',))
writefile_v_nudge_var[:] = v_wind[:,start_t_index:]
writefile_v_nudge_var.units = 'm s^-1'
writefile_v_nudge_var.description = 'N-S wind to nudge toward'

writefile_T_nudge_var = writefile_forcing_grp.createVariable('T_nudge', 'f4', ('levels','time',))
writefile_T_nudge_var[:] = T_abs[:,start_t_index:]
writefile_T_nudge_var.units = 'K'
writefile_T_nudge_var.description = 'absolute temperature to nudge toward'

writefile_thil_nudge_var = writefile_forcing_grp.createVariable('thil_nudge', 'f4', ('levels','time',))
writefile_thil_nudge_var[:] = thetal[:,start_t_index:]
writefile_thil_nudge_var.units = 'K'
writefile_thil_nudge_var.description = 'potential temperature to nudge toward'

writefile_qt_nudge_var = writefile_forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
writefile_qt_nudge_var[:] = qt[:,start_t_index:]
writefile_qt_nudge_var.units = 'kg kg^-1'
writefile_qt_nudge_var.description = 'q_t to nudge toward'

writefile_rad_heating_var = writefile_forcing_grp.createVariable('dT_dt_rad', 'f4', ('levels','time',))
writefile_rad_heating_var[:] = rad_heating[:,start_t_index:]
writefile_rad_heating_var.units = 'K s^-1'
writefile_rad_heating_var.description = 'prescribed radiative heating rate'

writefile_h_advec_thil_var = writefile_forcing_grp.createVariable('h_advec_thetail', 'f4', ('levels','time',))
writefile_h_advec_thil_var[:] = h_advec_thil[:,start_t_index:]
writefile_h_advec_thil_var.units = 'K s^-1'
writefile_h_advec_thil_var.description = 'prescribed theta_il tendency due to horizontal advection'

writefile_v_advec_thil_var = writefile_forcing_grp.createVariable('v_advec_thetail', 'f4', ('levels','time',))
writefile_v_advec_thil_var[:] = v_advec_thil[:,start_t_index:]
writefile_v_advec_thil_var.units = 'K s^-1'
writefile_v_advec_thil_var.description = 'prescribed theta_il tendency due to vertical advection'

writefile_h_advec_qt_var = writefile_forcing_grp.createVariable('h_advec_qt', 'f4', ('levels','time',))
writefile_h_advec_qt_var[:] = h_advec_qt[:,start_t_index:]
writefile_h_advec_qt_var.units = 'kg kg^-1 s^-1'
writefile_h_advec_qt_var.description = 'prescribed q_t tendency due to horizontal advection'

writefile_v_advec_qt_var = writefile_forcing_grp.createVariable('v_advec_qt', 'f4', ('levels','time',))
writefile_v_advec_qt_var[:] = v_advec_qt[:,start_t_index:]
writefile_v_advec_qt_var.units = 'kg kg^-1 s^-1'
writefile_v_advec_qt_var.description = 'prescribed q_t tendency due to vertical advection'


#close processed input file
writefile_fid.close()

#close raw input file
nc_fid.close()
