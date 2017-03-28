#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import forcing_file_common as ffc
import scipy.interpolate
import gmtb_scm_plotting_routines as gspr

def generate_forcing_file(input_file, output_file, z_sfc):
    #ensemble forcing files don't have the 'alt' variable; must pass in from best estimate forcing


    #read in raw ASTEX input file

    nc_fid = Dataset(input_file, 'r')

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
        #find time index corresponding to Jan 19, 2006 at 03Z
        if day[t] == 19 and hour[t] == 3:
            start_t_index = t
            break

    time = nc_fid.variables['time_offset'][:] #number of seconds since 00Z on 1/17/2006 (starts at 03Z)
    #subtract the initial time_offset from all values to get elapsed time since the start of the simulation
    time = time - time[start_t_index]
    levels = nc_fid.variables['lev'][:] #pressure levels in mb
    #convert levels to Pa
    levels = 100.0*levels
    #height = nc_fid.variables['alt'][:]
    lat = nc_fid.variables['lat'][:] #degrees south
    #convert latitutde to degrees north
    lat = -1*lat
    lon = nc_fid.variables['lon'][:] #degrees east
    T_abs = nc_fid.variables['T'][:] #absolute temperature (time, lev)
    T_abs = np.swapaxes(T_abs, 0, 1)
    #calculate theta_il from absolute temperature (assuming no condensate)
    thetal = np.zeros((levels.size,time.size),dtype=float)
    for t in range(time.size):
        thetal[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*T_abs[:,t]
    qv = nc_fid.variables['q'][:] #water vapor mixing ratio in g/kg (time, lev)
    qt = np.zeros((levels.size,time.size),dtype=float)
    qt_mr = 1.0E-3*np.swapaxes(qv, 0, 1) #swap the time and levels axis, convert to kg/kg
    qt = qt_mr/(1.0 + qt_mr) #convert to specific humidity from mixing ratio

    #ql and tke are not specified; set to zero
    ql = np.zeros((levels.size,time.size),dtype=float)
    qi = np.zeros((levels.size,time.size),dtype=float)
    tke = np.zeros((levels.size,time.size),dtype=float)
    # ozone_mmr = nc_fid.variables['o3mmr'][:]
    u_wind = nc_fid.variables['u'][:]
    u_wind = np.swapaxes(u_wind, 0, 1) #swap the time and levels axis
    v_wind = nc_fid.variables['v'][:]
    v_wind = np.swapaxes(v_wind, 0, 1) #swap the time and levels axis

    w_sub = np.zeros((levels.size,time.size),dtype=float)
    omega = nc_fid.variables['omega'][:] #vertical pressure velocity in mb/hr
    omega = omega*100.0/3600.0 #convert to Pa/s
    omega = np.swapaxes(omega, 0, 1) #swap the time and levels axis
    #convert to w
    for t in range(time.size):
        w_sub[:,t] = ffc.omega_to_w(omega[:,t],levels,T_abs[:,t])

    #T_surf = nc_fid.variables['T_skin'][:]
    #T_surf = T_surf + 273.15 #convert to K
    T_surf = (29 + 273.15)*np.ones((time.size),dtype=float) #forcing instructions specify time-invariant 29 deg C.
    p_surf = nc_fid.variables['p_srf_aver'][:] #mb
    p_surf = p_surf*100.0 #convert to Pa

    h_advec_thil = np.zeros((levels.size,time.size),dtype=float)
    h_advec_T = nc_fid.variables['T_adv_h'][:] #K/hr
    h_advec_T = np.swapaxes(h_advec_T, 0, 1)/3600.0 #swap the time and levels axis, convert to K/s
    for t in range(time.size):
        h_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*h_advec_T[:,t] #convert to potential temperature

    v_advec_thil = np.zeros((levels.size,time.size),dtype=float)
    v_advec_T = nc_fid.variables['T_adv_v'][:] #K/hr
    v_advec_T = np.swapaxes(v_advec_T, 0, 1)/3600.0 #swap the time and levels axis, convert to K/s
    for t in range(time.size):
        v_advec_thil[:,t] = (ffc.p0/levels)**(ffc.R_dry/ffc.c_p)*v_advec_T[:,t] #convert to potential temperature

    dT_dt = nc_fid.variables['dTdt'][:] #K/hr
    dT_dt = np.swapaxes(dT_dt, 0, 1)/3600.0 #swap the time and levels axis, convert to K/s

    # h_advec_T = h_advec_T*86400.0
    # v_advec_T = v_advec_T*86400.0
    # dT_dt = dT_dt*86400.0
    # gspr.contour_plot_firl(time, levels, h_advec_T, np.min(h_advec_T), np.max(h_advec_T), 'h advec T', 'time', 'pressure', 'h_advec_T.eps', y_inverted = True)
    # gspr.contour_plot_firl(time, levels, v_advec_T, np.min(v_advec_T), np.max(v_advec_T), 'v advec T', 'time', 'pressure', 'v_advec_T.eps', y_inverted = True)
    # gspr.contour_plot_firl(time, levels, dT_dt, np.min(dT_dt), np.max(dT_dt), 'total T tend', 'time', 'pressure', 'dT_dt.eps', y_inverted = True)

    h_advec_qt = np.zeros((levels.size,time.size),dtype=float)
    h_advec_qt = nc_fid.variables['q_adv_h'][:] #g/kg/hr
    h_advec_qt = np.swapaxes(h_advec_qt, 0, 1)*1.0E-3/3600.0 #swap the time and levels axis, convert to kg/kg/s
    h_advec_qt = h_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

    v_advec_qt = np.zeros((levels.size,time.size),dtype=float)
    v_advec_qt = nc_fid.variables['q_adv_v'][:] #g/kg/hr
    v_advec_qt = np.swapaxes(v_advec_qt, 0, 1)*1.0E-3/3600.0 #swap the time and levels axis, convert to kg/kg/s
    v_advec_qt = v_advec_qt/(1.0 + qt_mr)**2 #convert to specific humidity

    dq_dt = nc_fid.variables['dqdt'][:]
    dq_dt = np.swapaxes(dq_dt, 0, 1)*1.0E-3/3600.0 #swap the time and levels axis, convert to K/s

    # h_advec_qt = h_advec_qt*86400.0
    # v_advec_qt = v_advec_qt*86400.0
    # dq_dt = dq_dt*86400.0
    #
    # gspr.contour_plot_firl(time, levels, h_advec_qt, np.min(h_advec_qt), np.max(h_advec_qt), 'h advec q', 'time', 'pressure', 'h_advec_q_ens.eps', y_inverted = True)
    # gspr.contour_plot_firl(time, levels, v_advec_qt, np.min(v_advec_qt), np.max(v_advec_qt), 'v advec q', 'time', 'pressure', 'v_advec_q_ens.eps', y_inverted = True)
    # gspr.contour_plot_firl(time, levels, dq_dt, np.min(dq_dt), np.max(dq_dt), 'total q tend', 'time', 'pressure', 'dq_dt_ens.eps', y_inverted = True)



    phi_sfc = nc_fid.variables['phis'][:]
    #z_sfc = nc_fid.variables['alt'][:]

    height = ffc.get_height_from_pres(T_abs[:,0],levels,z_sfc)

    #the following variables are not in the astex forcing file, but are included in other cases
    rad_heating = np.zeros((levels.size,time.size),dtype=float)
    u_g = np.zeros((levels.size,time.size),dtype=float)
    v_g = np.zeros((levels.size,time.size),dtype=float)

    # Open ozone file
    f = open('../raw_case_input/twpice_CRM_ozone.txt', 'r')

    # Read and ignore header lines
    header1 = f.readline()

    oz_pres = []
    oz_data = []
    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
        oz_pres.append(float(columns[1]))
        oz_data.append(float(columns[2]))

    f.close()

    oz_pres = 100.0*np.array(oz_pres)
    oz_data = np.array(oz_data)
    oz_f = scipy.interpolate.interp1d(oz_pres, oz_data)
    ozone_ppb = oz_f(levels[1:])
    ozone_ppb = np.insert(ozone_ppb, 0, oz_data[0])
    ozone_mmr = ozone_ppb*1.0E-9

    #
    #open processed input file for writing

    writefile_fid = Dataset(output_file, 'w', format='NETCDF4')
    writefile_fid.description = "GMTB SCM forcing file for TWP-ICE case"

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
    writefile_ozone_var[:] = ozone_mmr
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

#subroutine for printing progress to the command line
def print_progress(n_complete, n_total):
    print str(n_complete) + ' of ' + str(n_total) + ' complete: (' + str(100.0*n_complete/float(n_total)) + '%)'

reload(ffc)
num_files = 100

nc_fid = Dataset('../raw_case_input/twp180iopsndgvarana_v2.1_C3.c1.20060117.000000.cdf', 'r')
z_sfc = nc_fid.variables['alt'][:]
nc_fid.close()

for i in range(num_files):
    infile =  '../raw_case_input/twpice_ensemble/twpice_p' + str(i).zfill(2) + '.cdf'
    outfile =  '../processed_case_input/twpice_ens' + str(i).zfill(2) + '.nc'

    generate_forcing_file(infile, outfile, z_sfc)
    print_progress(i+1, num_files)
