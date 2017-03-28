#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import forcing_file_common as ffc

#read in raw ASTEX input file

nc_fid = Dataset("../raw_case_input/astex_input_v5.nc", 'r')

#ncdump to look at raw input file
#nc_attrs, nc_dims, nc_vars = ncdump(nc_fid, False)

#netCDF how-to
#set attributes
#file_id.setncattr(file_id.variables['variable'].ncattr(), nc_fid.variables['time'].getncattr(ncattr))
#get attributes
#file_id.variables['variable'].getncattr(index)

#get raw input variables

time = nc_fid.variables['time'][:]
levels = nc_fid.variables['lev'][:]
height = nc_fid.variables['height'][:]
lat = nc_fid.variables['lat'][:]
lon = nc_fid.variables['lon'][:]

thetal = nc_fid.variables['thetal'][:]
qt = nc_fid.variables['qt'][:]
ql = nc_fid.variables['ql'][:]
tke = nc_fid.variables['tke'][:]
ozone_mmr = nc_fid.variables['o3mmr'][:]
u_wind = nc_fid.variables['u'][:]
v_wind = nc_fid.variables['v'][:]

w_sub = nc_fid.variables['w'][:]
u_g = nc_fid.variables['ug'][:]
v_g = nc_fid.variables['vg'][:]

T_surf = nc_fid.variables['Tg'][:]
p_surf = nc_fid.variables['Ps'][:]
u_nudge = nc_fid.variables['ufa'][:]
v_nudge = nc_fid.variables['vfa'][:]

#the following variables are not in the astex forcing file, but are included in other cases
rad_heating = np.zeros((levels.size,time.size),dtype=float)
h_advec_thil = np.zeros((levels.size,time.size),dtype=float)
v_advec_thil = np.zeros((levels.size,time.size),dtype=float)
h_advec_qt = np.zeros((levels.size,time.size),dtype=float)
v_advec_qt = np.zeros((levels.size,time.size),dtype=float)
u_nudge = np.zeros((levels.size,time.size),dtype=float)
v_nudge = np.zeros((levels.size,time.size),dtype=float)
T_nudge = np.zeros((levels.size,time.size),dtype=float)
thil_nudge = np.zeros((levels.size,time.size),dtype=float)
qt_nudge = np.zeros((levels.size,time.size),dtype=float)


#open processed input file for writing

writefile_fid = Dataset('../processed_case_input/astex.nc', 'w', format='NETCDF4')
writefile_fid.description = "GMTB SCM forcing file for ASTEX case"

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

#initial group

writefile_height_var = writefile_initial_grp.createVariable('height', 'f4', ('levels',))
writefile_height_var[:] = height
writefile_height_var.units = 'm'
writefile_height_var.description = 'physical height at pressure levels'

writefile_thetail_var = writefile_initial_grp.createVariable('thetail', 'f4', ('levels',))
writefile_thetail_var[:] = thetal
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
writefile_qi_var[:] = np.zeros(len(levels))
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
writefile_ozone_var[:] = ozone_mmr
writefile_ozone_var.units = 'kg kg^-1'
writefile_ozone_var.description = 'initial profile of ozone mass mixing ratio'

#forcing group

writefile_lat_var = writefile_forcing_grp.createVariable('lat', 'f4', ('time',))
writefile_lat_var[:] = lat
writefile_lat_var.units = 'degrees N'
writefile_lat_var.description = 'latitude of column'

writefile_lon_var = writefile_forcing_grp.createVariable('lon', 'f4', ('time',))
writefile_lon_var[:] = lon
writefile_lon_var.units = 'degrees E'
writefile_lon_var.description = 'longitude of column'

writefile_p_surf_var = writefile_forcing_grp.createVariable('p_surf', 'f4', ('time',))
writefile_p_surf_var[:] = p_surf
writefile_p_surf_var.units = 'Pa'
writefile_p_surf_var.description = 'surface pressure'

writefile_T_surf_var = writefile_forcing_grp.createVariable('T_surf', 'f4', ('time',))
writefile_T_surf_var[:] = T_surf
writefile_T_surf_var.units = 'K'
writefile_T_surf_var.description = 'surface absolute temperature'

writefile_w_ls_var = writefile_forcing_grp.createVariable('w_ls', 'f4', ('levels','time',))
writefile_w_ls_var[:] = w_sub
writefile_w_ls_var.units = 'm s^-1'
writefile_w_ls_var.description = 'large scale vertical velocity'

#convert to w
omega = np.zeros((levels.size,time.size),dtype=float)
T_abs = ffc.theta_to_T(thetal,levels)
for t in range(time.size):
    omega[:,t] = ffc.w_to_omega(w_sub[:,t],levels,T_abs)

writefile_omega_var = writefile_forcing_grp.createVariable('omega', 'f4', ('levels','time',))
writefile_omega_var[:] = omega
writefile_omega_var.units = 'Pa s^-1'
writefile_omega_var.description = 'large scale pressure velocity'

#ASTEX input file has one value per time for the following variables; need to expand the dimenions of the array to have levels and times
ug_tile = np.tile(u_g.transpose(),(len(levels),1))
writefile_u_g_var = writefile_forcing_grp.createVariable('u_g', 'f4', ('levels','time',))
writefile_u_g_var[:] = ug_tile
writefile_u_g_var.units = 'm s^-1'
writefile_u_g_var.description = 'large scale geostrophic E-W wind'

vg_tile = np.tile(v_g.transpose(),(len(levels),1))
writefile_v_g_var = writefile_forcing_grp.createVariable('v_g', 'f4', ('levels','time',))
writefile_v_g_var[:] = vg_tile
writefile_v_g_var.units = 'm s^-1'
writefile_v_g_var.description = 'large scale geostrophic N-S wind'

writefile_u_nudge_var = writefile_forcing_grp.createVariable('u_nudge', 'f4', ('levels','time',))
writefile_u_nudge_var[:] = u_nudge
writefile_u_nudge_var.units = 'm s^-1'
writefile_u_nudge_var.description = 'E-W wind to nudge toward'

writefile_v_nudge_var = writefile_forcing_grp.createVariable('v_nudge', 'f4', ('levels','time',))
writefile_v_nudge_var[:] = v_nudge
writefile_v_nudge_var.units = 'm s^-1'
writefile_v_nudge_var.description = 'N-S wind to nudge toward'

writefile_T_nudge_var = writefile_forcing_grp.createVariable('T_nudge', 'f4', ('levels','time',))
writefile_T_nudge_var[:] = T_nudge
writefile_T_nudge_var.units = 'K'
writefile_T_nudge_var.description = 'absolute temperature to nudge toward'

writefile_thil_nudge_var = writefile_forcing_grp.createVariable('thil_nudge', 'f4', ('levels','time',))
writefile_thil_nudge_var[:] = thil_nudge
writefile_thil_nudge_var.units = 'K'
writefile_thil_nudge_var.description = 'potential temperature to nudge toward'

writefile_qt_nudge_var = writefile_forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
writefile_qt_nudge_var[:] = qt_nudge
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

#close raw input file
nc_fid.close()
