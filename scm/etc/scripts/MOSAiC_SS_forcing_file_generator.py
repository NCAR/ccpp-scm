#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import forcing_file_common as ffc
import scipy.interpolate
import scm_plotting_routines as spr

#reload(ffc)

#read in raw input file

nc_fid = Dataset("../../data/raw_case_input/MOSAiC_2Mar20200Z_raw.nc", 'r')

#ncdump to look at raw input file
#nc_attrs, nc_dims, nc_vars = ncdump(nc_fid, False)

#netCDF how-to
#set attributes
#file_id.setncattr(file_id.variables['variable'].ncattr(), nc_fid.variables['time'].getncattr(ncattr))
#get attributes
#file_id.variables['variable'].getncattr(index)

#get raw input variables

day = nc_fid.variables['day'][:]
hour = nc_fid.variables['hour'][:]
for t in range(day.size):
    #find time index corresponding to October 29, 2019 at 0Z
    if day[t] == 2 and hour[t] == 0:
        start_t_index = t
        break

time = nc_fid.variables['time_offset'][:] #number of seconds since 00Z on 1/17/2006 (starts at 03Z)
#subtract the initial time_offset from all values to get elapsed time since the start of the simulation
time = time - time[start_t_index]
levels = nc_fid.variables['levels'][:] #pressure levels in mb
#convert levels to Pa
levels = 100.0*levels
ice_thickness = np.zeros((2),dtype=float)
#height = nc_fid.variables['alt'][:]
lat = nc_fid.variables['lat'][:] #degrees north
#convert latitutde to degrees north
lon = nc_fid.variables['lon'][:] #degrees east
#use upstream T
T_abs = nc_fid.variables['Tu'][:] #absolute temperature (time, lev)
T_abs = np.swapaxes(T_abs, 0, 1)
thetailu = nc_fid.variables['thetailu'][:] #theta_il (time, lev)
thetailu = np.swapaxes(thetailu, 0, 1)
#use unpstream theta instead of thetail
thetail = nc_fid.variables['thetail'][:] #theta (time, lev)
thetail = np.swapaxes(thetail, 0, 1)
#calculate theta_il from absolute temperature (assuming no condensate)
#thetal = np.zeros((levels.size,time.size),dtype=float)
#for t in range(time.size):
#    thetal[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*T_abs[:,t]
qv = nc_fid.variables['qu'][:] #water vapor mixing ratio in kg/kg (time, lev)
qv = np.swapaxes(qv, 0, 1)
qt = np.zeros((levels.size,time.size),dtype=float)
qt_mr = nc_fid.variables['qt'][:] #water vapor mixing ratio in kg/kg (time, lev)
qt_mr = np.swapaxes(qt_mr, 0, 1) #swap the time and levels axis
qt = qt_mr/(1.0 + qt_mr) #convert to specific humidity from mixing ratio
qtu = np.zeros((levels.size,time.size),dtype=float)
#use upstream qv instead of qt
qtu_mr = nc_fid.variables['qtu'][:] #water vapor mixing ratio in kg/kg (time, lev)
qtu_mr = np.swapaxes(qtu_mr, 0, 1) #swap the time and levels axis
qtu = qtu_mr/(1.0 + qtu_mr) #convert to specific humidity from mixing ratio

#ql and tke are not specified; set to zero
#ql = np.zeros((levels.size,time.size),dtype=float)
#qi = np.zeros((levels.size,time.size),dtype=float)
ql = nc_fid.variables['ql'][:] #ql (time, lev)
ql = np.swapaxes(ql, 0, 1)
qi = nc_fid.variables['qi'][:] #ql (time, lev)
qi = np.swapaxes(qi, 0, 1)
tke = np.zeros((levels.size,time.size),dtype=float)
# ozone_mmr = nc_fid.variables['o3mmr'][:]
u_wind = nc_fid.variables['u'][:]
u_wind = np.swapaxes(u_wind, 0, 1) #swap the time and levels axis
v_wind = nc_fid.variables['v'][:]
v_wind = np.swapaxes(v_wind, 0, 1) #swap the time and levels axis

#w_sub = np.zeros((levels.size,time.size),dtype=float)
omega = nc_fid.variables['omega'][:] #vertical pressure velocity in Pa/s
omega = np.swapaxes(omega, 0, 1) #swap the time and levels axis
w_sub = nc_fid.variables['w'][:] #vertical velocity in m/s
w_sub = np.swapaxes(w_sub, 0, 1) #swap the time and levels axis
#convert to w
#for t in range(time.size):
#    w_sub[:,t] = ffc.omega_to_w(omega[:,t],levels,T_abs[:,t])

T_surf = nc_fid.variables['T_skin'][:]
#T_surf = T_surf + 273.15 #convert to K
#T_surf = (29 + 273.15)*np.ones((time.size),dtype=float) #forcing instructions specify time-invariant 29 deg C.
p_surf = nc_fid.variables['p_srf'][:] #Pa
#p_surf = p_surf*100.0 #convert to Pa

#h_advec_thil = np.zeros((levels.size,time.size),dtype=float)
h_advec_thil = nc_fid.variables['h_advec_thetail'][:] #K/s
h_advec_thil = np.swapaxes(h_advec_thil, 0, 1) #swap the time and levels axis
#for t in range(time.size):
#    h_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*h_advec_T[:,t] #convert to potential temperature

#v_advec_thil = np.zeros((levels.size,time.size),dtype=float)
v_advec_thil = nc_fid.variables['v_advec_thetail'][:] #K/s
v_advec_thil = np.swapaxes(v_advec_thil, 0, 1) #swap the time and levels axis
#for t in range(time.size):
#    v_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*v_advec_T[:,t] #convert to potential temperature

dT_dt = np.zeros((levels.size,time.size),dtype=float)

# h_advec_T = h_advec_T*86400.0
# v_advec_T = v_advec_T*86400.0
# dT_dt = dT_dt*86400.0
# spr.contour_plot_firl(time, levels, h_advec_T, np.min(h_advec_T), np.max(h_advec_T), 'h advec T', 'time', 'pressure', 'h_advec_T.eps', y_inverted = True)
# spr.contour_plot_firl(time, levels, v_advec_T, np.min(v_advec_T), np.max(v_advec_T), 'v advec T', 'time', 'pressure', 'v_advec_T.eps', y_inverted = True)
# spr.contour_plot_firl(time, levels, dT_dt, np.min(dT_dt), np.max(dT_dt), 'total T tend', 'time', 'pressure', 'dT_dt.eps', y_inverted = True)

#h_advec_qt = np.zeros((levels.size,time.size),dtype=float)
h_advec_qt = nc_fid.variables['h_advec_qt'][:] #kg/kg/s
h_advec_qt = np.swapaxes(h_advec_qt, 0, 1) #swap the time and levels axis
h_advec_qt = h_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

h_advec_qi = np.zeros((levels.size,time.size),dtype=float)
h_advec_ql = np.zeros((levels.size,time.size),dtype=float)

#v_advec_qt = np.zeros((levels.size,time.size),dtype=float)
v_advec_qt = nc_fid.variables['v_advec_qt'][:] #kg/kg/s
v_advec_qt = np.swapaxes(v_advec_qt, 0, 1) #swap the time and levels axis
v_advec_qt = v_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

dq_dt = np.zeros((levels.size,time.size),dtype=float)
#dq_dt = nc_fid.variables['dqdt'][:]
#dq_dt = np.swapaxes(dq_dt, 0, 1)*1.0E-3/3600.0 #swap the time and levels axis, convert to K/s

#h_advec_qt = h_advec_qt*86400.0
#v_advec_qt = v_advec_qt*86400.0
#dq_dt = dq_dt*86400.0
#spr.contour_plot_firl(time, levels, h_advec_qt, np.min(h_advec_qt), np.max(h_advec_qt), 'h advec q', 'time', 'pressure', 'h_advec_q.eps', y_inverted = True)
#spr.contour_plot_firl(time, levels, v_advec_qt, np.min(v_advec_qt), np.max(v_advec_qt), 'v advec q', 'time', 'pressure', 'v_advec_q.eps', y_inverted = True)
#spr.contour_plot_firl(time, levels, dq_dt, np.min(dq_dt), np.max(dq_dt), 'total q tend', 'time', 'pressure', 'dq_dt.eps', y_inverted = True)

#phi_sfc = nc_fid.variables['phis'][:]
#z_sfc = nc_fid.variables['alt'][:]
height = nc_fid.variables['height'][:]
#height = ffc.get_height_from_pres(T_abs[:,0],levels,z_sfc)

#the following variables are not in this forcing file, but are included in other cases
soil_depth = np.zeros(4,dtype=float)
soil_depth[:] = [0.1, 0.4, 1.0, 2.0]

#rad_heating = nc_fid.variables['dT_dt_rad'][:] #K/s
#rad_heating = np.swapaxes(rad_heating, 0, 1) #swap the time and levels axis
rad_heating = np.zeros((levels.size,time.size),dtype=float)
u_g = np.zeros((levels.size,time.size),dtype=float)
v_g = np.zeros((levels.size,time.size),dtype=float)

ozone = np.zeros((levels.size,time.size),dtype=float)
ozone = ozone + 2.e-8

tiice = np.zeros((ice_thickness.size),dtype=float)
##tiice[0] = T_surf[start_t_index] 
tiice[0] = 271.35 + .75*(T_surf[start_t_index] - 271.35)
tiice[1] = 271.35 + .25*(T_surf[start_t_index] - 271.35)
stc      = tiice[0]
smc      = 0.33
slc      = 0.33
hice     = 0.3
slmsk    = 2.0    #sea ice
tsfco    = T_surf[start_t_index]
weasd    = 200.0 #mm water equivalent snow depth
fice     = 1.0
tisfc    = tsfco
snwdph   = 2.e-4
tg3      = 271.35
zorl     = 15.0
alvsf    = 0.06
alnsf    = 0.06
alvwf    = 0.06
alnwf    = 0.06
facsf    = 0.0
facwf    = 0.0
vegfrac  = 0.0
canopy   = 0.0
vegtyp   = 0
soiltyp  = 0
scolor   = 1
uustar   = 0.3828793
shdmin   = 0.0
shdmax   = 0.0
slopetyp = 1
snoalb   = 0.0
sncovr   = 1.0

# Open ozone file
#f = open('../../data/raw_case_input/twpice_CRM_ozone.txt', 'r')

# Read and ignore header lines
#header1 = f.readline()

#oz_pres = []
#oz_data = []
# Loop over lines and extract variables of interest
#for line in f:
#    line = line.strip()
#    columns = line.split()
#    oz_pres.append(float(columns[1]))
#    oz_data.append(float(columns[2]))

#f.close()

#oz_pres = 100.0*np.array(oz_pres)
#oz_data = np.array(oz_data)
#oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
#ozone_ppb = oz_f(levels[1:].tolist())
#ozone_ppb = np.insert(ozone_ppb, 0, oz_data[0])
#ozone_mmr = ozone_ppb*1.0E-9

#
#open processed input file for writing

writefile_fid = Dataset('../../data/processed_case_input/MOSAiC.nc', 'w', format='NETCDF4')
writefile_fid.description = "CCPP SCM forcing file for MOSAiC case"

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

writefile_soil_depth_dim = writefile_fid.createDimension('soil_depth', None)
writefile_soil_depth_var = writefile_fid.createVariable('soil_depth', 'f4', ('soil_depth',))
writefile_soil_depth_var[:] = soil_depth[:]
writefile_soil_depth_var.units = 'm'
writefile_soil_depth_var.description = 'soil depth'

writefile_ice_thickness_dim = writefile_fid.createDimension('ice_thickness', None)
writefile_ice_thickness_var = writefile_fid.createVariable('ice_thickness', 'f4', ('ice_thickness',))
writefile_ice_thickness_var[:] = ice_thickness
writefile_ice_thickness_var.units = 'm'
writefile_ice_thickness_var.description = 'depth of ice layers'

#create variables and write them out

#scalar group

writefile_lat_var = writefile_scalar_grp.createVariable('lat', 'f4')
writefile_lat_var[:] = lat
writefile_lat_var.units = 'degrees N'
writefile_lat_var.description = 'latitude of column'

writefile_lon_var = writefile_scalar_grp.createVariable('lon', 'f4')
writefile_lon_var[:] = lon
writefile_lon_var.units = 'degrees E'
writefile_lon_var.description = 'longitude of column'

writefile_hice_var = writefile_scalar_grp.createVariable('hice', 'f4')
writefile_hice_var[:] = hice
writefile_hice_var.units = 'm'
writefile_hice_var.description = 'sea ice thickness'

writefile_slmsk_var = writefile_scalar_grp.createVariable('slmsk', 'f4')
writefile_slmsk_var[:] = slmsk
writefile_slmsk_var.units = ''
writefile_slmsk_var.description = 'land-sea-ice mask'

writefile_tsfco_var = writefile_scalar_grp.createVariable('tsfco', 'f4')
writefile_tsfco_var[:] = tsfco
writefile_tsfco_var.units = 'm'
writefile_tsfco_var.description = 'sea ice surface skin temperature'

writefile_weasd_var = writefile_scalar_grp.createVariable('weasd', 'f4')
writefile_weasd_var[:] = weasd
writefile_weasd_var.units = 'mm'
writefile_weasd_var.description = 'water equivalent accumulated snow depth'

writefile_fice_var = writefile_scalar_grp.createVariable('fice', 'f4')
writefile_fice_var[:] = fice
writefile_fice_var.units = '1'
writefile_fice_var.description = 'ice fraction'

writefile_tisfc_var = writefile_scalar_grp.createVariable('tisfc', 'f4')
writefile_tisfc_var[:] = tisfc
writefile_tisfc_var.units = 'K'
writefile_tisfc_var.description = 'ice surface temperature'

writefile_snwdph_var = writefile_scalar_grp.createVariable('snwdph', 'f4')
writefile_snwdph_var[:] = snwdph
writefile_snwdph_var.units = 'mm'
writefile_snwdph_var.description = 'water equivalent snow depth'

writefile_tg3_var = writefile_scalar_grp.createVariable('tg3', 'f4')
writefile_tg3_var[:] = tg3
writefile_tg3_var.units = 'K'
writefile_tg3_var.description = 'deep soil temperature'

writefile_zorl_var = writefile_scalar_grp.createVariable('zorl', 'f4')
writefile_zorl_var[:] = zorl
writefile_zorl_var.units = 'cm'
writefile_zorl_var.description = 'composite surface roughness length'

writefile_zorll_var = writefile_scalar_grp.createVariable('zorll', 'f4')
writefile_zorll_var[:] = zorl
writefile_zorll_var.units = 'cm'
writefile_zorll_var.description = 'surface roughness length over land'

writefile_zorlw_var = writefile_scalar_grp.createVariable('zorlw', 'f4')
writefile_zorlw_var[:] = zorl
writefile_zorlw_var.units = 'cm'
writefile_zorlw_var.description = 'surface roughness length over ocean'

writefile_zorli_var = writefile_scalar_grp.createVariable('zorli', 'f4')
writefile_zorli_var[:] = zorl
writefile_zorli_var.units = 'cm'
writefile_zorli_var.description = 'surface roughness length over ice'

writefile_alvsf_var = writefile_scalar_grp.createVariable('alvsf', 'f4')
writefile_alvsf_var[:] = alvsf
writefile_alvsf_var.units = ''
writefile_alvsf_var.description = '60 degree vis albedo with strong cosz dependency'

writefile_alnsf_var = writefile_scalar_grp.createVariable('alnsf', 'f4')
writefile_alnsf_var[:] = alnsf
writefile_alnsf_var.units = ''
writefile_alnsf_var.description = '60 degree nir albedo with strong cosz dependency'

writefile_alvwf_var = writefile_scalar_grp.createVariable('alvwf', 'f4')
writefile_alvwf_var[:] = alvwf
writefile_alvwf_var.units = ''
writefile_alvwf_var.description = '60 degree vis albedo with weak cosz dependency'

writefile_alnwf_var = writefile_scalar_grp.createVariable('alnwf', 'f4')
writefile_alnwf_var[:] = alnwf
writefile_alnwf_var.units = ''
writefile_alnwf_var.description = '60 degree nir albedo with weak cosz dependency'

writefile_facsf_var = writefile_scalar_grp.createVariable('facsf', 'f4')
writefile_facsf_var[:] = facsf
writefile_facsf_var.units = ''
writefile_facsf_var.description = 'fractional coverage with strong cosz dependency'

writefile_facwf_var = writefile_scalar_grp.createVariable('facwf', 'f4')
writefile_facwf_var[:] = facwf
writefile_facwf_var.units = ''
writefile_facwf_var.description = 'fractional coverage with weak cosz dependency'

writefile_vegfrac_var = writefile_scalar_grp.createVariable('vegfrac', 'f4')
writefile_vegfrac_var[:] = vegfrac
writefile_vegfrac_var.units = '1'
writefile_vegfrac_var.description = 'vegetation fraction'

writefile_canopy_var = writefile_scalar_grp.createVariable('canopy', 'f4')
writefile_canopy_var[:] = canopy
writefile_canopy_var.units = 'kg m-2'
writefile_canopy_var.description = 'amount of water stored in camopy'

writefile_vegtyp_var = writefile_scalar_grp.createVariable('vegtyp', 'f4')
writefile_vegtyp_var[:] = vegtyp
writefile_vegtyp_var.units = ''
writefile_vegtyp_var.description = 'vegetation type 1-12'

writefile_soiltyp_var = writefile_scalar_grp.createVariable('soiltyp', 'f4')
writefile_soiltyp_var[:] = soiltyp
writefile_soiltyp_var.units = ''
writefile_soiltyp_var.description = 'soil type 1-12'

writefile_scolor_var = writefile_scalar_grp.createVariable('scolor', 'f4')
writefile_scolor_var[:] = scolor
writefile_scolor_var.units = ''
writefile_scolor_var.description = 'soil color'

writefile_uustar_var = writefile_scalar_grp.createVariable('uustar', 'f4')
writefile_uustar_var[:] = uustar
writefile_uustar_var.units = 'm s-1'
writefile_uustar_var.description = 'friction velocity'

writefile_shdmin_var = writefile_scalar_grp.createVariable('shdmin', 'f4')
writefile_shdmin_var[:] = shdmin
writefile_shdmin_var.units = '1'
writefile_shdmin_var.description = 'minimum vegetation fraction'

writefile_shdmax_var = writefile_scalar_grp.createVariable('shdmax', 'f4')
writefile_shdmax_var[:] = shdmax
writefile_shdmax_var.units = '1'
writefile_shdmax_var.description = 'maximum vegetation fraction'

writefile_slopetyp_var = writefile_scalar_grp.createVariable('slopetyp', 'f4')
writefile_slopetyp_var[:] = slopetyp
writefile_slopetyp_var.units = ''
writefile_slopetyp_var.description = 'slope type 1-9'

writefile_snoalb_var = writefile_scalar_grp.createVariable('snoalb', 'f4')
writefile_snoalb_var[:] = snoalb
writefile_snoalb_var.units = '1'
writefile_snoalb_var.description = 'maximum snow albedo'

writefile_sncovr_var = writefile_scalar_grp.createVariable('sncovr', 'f4')
writefile_sncovr_var[:] = sncovr
writefile_sncovr_var.units = '1'
writefile_sncovr_var.description = 'surface snow area fraction'

#initial group

writefile_height_var = writefile_initial_grp.createVariable('height', 'f4', ('levels',))
writefile_height_var[:] = height
writefile_height_var.units = 'm'
writefile_height_var.description = 'physical height at pressure levels'

writefile_tiice_var = writefile_initial_grp.createVariable('tiice', 'f4', ('ice_thickness',))
writefile_tiice_var[:] = tiice
writefile_tiice_var.units = 'K'
writefile_tiice_var.description = 'initial profile of sea ice internal temperature'

writefile_stc_var = writefile_initial_grp.createVariable('stc', 'f4', ('soil_depth',))
writefile_stc_var[:] = stc
writefile_stc_var.units = 'K'
writefile_stc_var.description = 'initial profile of sea ice internal temperature'

writefile_smc_var = writefile_initial_grp.createVariable('smc', 'f4', ('soil_depth',))
writefile_smc_var[:] = smc
writefile_smc_var.units = 'm3 m-3'
writefile_smc_var.description = 'initial profile of soil moisture'

writefile_slc_var = writefile_initial_grp.createVariable('slc', 'f4', ('soil_depth',))
writefile_slc_var[:] = slc
writefile_slc_var.units = 'm3 m-3'
writefile_slc_var.description = 'initial profile of soil liquid water'

writefile_thetail_var = writefile_initial_grp.createVariable('thetail', 'f4', ('levels',))
writefile_thetail_var[:] = thetail[:,start_t_index]
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
writefile_ozone_var[:] = ozone[:,start_t_index]
writefile_ozone_var.units = 'kg kg^-1'
writefile_ozone_var.description = 'initial profile of ozone mass mixing ratio'

#forcing group

writefile_p_surf_var = writefile_forcing_grp.createVariable('p_surf', 'f4', ('time',))
writefile_p_surf_var[:] = p_surf[start_t_index:]
writefile_p_surf_var.units = 'Pa'
writefile_p_surf_var.description = 'surface pressure'

writefile_T_surf_var = writefile_forcing_grp.createVariable('T_surf', 'f4', ('time',))
writefile_T_surf_var[:] = T_surf[start_t_index:]
writefile_T_surf_var.units = 'K'
writefile_T_surf_var.description = 'surface absolute temperature'

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
writefile_thil_nudge_var[:] = thetailu[:,start_t_index:]
writefile_thil_nudge_var.units = 'K'
writefile_thil_nudge_var.description = 'potential temperature to nudge toward'

writefile_qt_nudge_var = writefile_forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
writefile_qt_nudge_var[:] = qtu[:,start_t_index:]
writefile_qt_nudge_var.units = 'kg kg^-1'
writefile_qt_nudge_var.description = 'q_t to nudge toward'

writefile_qi_nudge_var = writefile_forcing_grp.createVariable('qi_nudge', 'f4', ('levels','time',))
writefile_qi_nudge_var[:] = qi[:,start_t_index:]
writefile_qi_nudge_var.units = 'kg kg^-1'
writefile_qi_nudge_var.description = 'q_i to nudge toward'

writefile_ql_nudge_var = writefile_forcing_grp.createVariable('ql_nudge', 'f4', ('levels','time',))
writefile_ql_nudge_var[:] = ql[:,start_t_index:]
writefile_ql_nudge_var.units = 'kg kg^-1'
writefile_ql_nudge_var.description = 'q_l to nudge toward'

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

writefile_h_advec_qi_var = writefile_forcing_grp.createVariable('h_advec_qi', 'f4', ('levels','time',))
writefile_h_advec_qi_var[:] = h_advec_qi[:,start_t_index:]
writefile_h_advec_qi_var.units = 'kg kg^-1 s^-1'
writefile_h_advec_qi_var.description = 'prescribed q_i tendency due to horizontal advection'

writefile_h_advec_ql_var = writefile_forcing_grp.createVariable('h_advec_ql', 'f4', ('levels','time',))
writefile_h_advec_ql_var[:] = h_advec_ql[:,start_t_index:]
writefile_h_advec_ql_var.units = 'kg kg^-1 s^-1'
writefile_h_advec_ql_var.description = 'prescribed q_l tendency due to horizontal advection'

writefile_v_advec_qt_var = writefile_forcing_grp.createVariable('v_advec_qt', 'f4', ('levels','time',))
writefile_v_advec_qt_var[:] = v_advec_qt[:,start_t_index:]
writefile_v_advec_qt_var.units = 'kg kg^-1 s^-1'
writefile_v_advec_qt_var.description = 'prescribed q_t tendency due to vertical advection'


#close processed input file
writefile_fid.close()

#close raw input file
nc_fid.close()
