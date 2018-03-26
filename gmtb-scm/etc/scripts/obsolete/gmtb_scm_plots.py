#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gmtb_scm_plotting_routines as gspr
import bisect

plot_dir = "../plots/"

reload(gspr)

nc_fid = Dataset("../output_astex/output.nc", 'r')
LES_fid = Dataset("../comparison_data/astex_LES.nc", 'r')

#read in SCM data
time = nc_fid.variables['time'][:]
pres_l = nc_fid.variables['pres'][:]
pres_i = nc_fid.variables['pres_i'][:]
sigma_l = nc_fid.variables['sigma'][:]
sigma_i = nc_fid.variables['sigma_i'][:]
phi_l = nc_fid.variables['phi'][:]
phi_i = nc_fid.variables['phi_i'][:]
qv = nc_fid.variables['qv'][:]
T = nc_fid.variables['T'][:]
u = nc_fid.variables['u'][:]
v = nc_fid.variables['v'][:]
qc = nc_fid.variables['qc'][:]
qv_force_tend = nc_fid.variables['qv_force_tend'][:]
T_force_tend = nc_fid.variables['T_force_tend'][:]
u_force_tend = nc_fid.variables['u_force_tend'][:]
v_force_tend = nc_fid.variables['v_force_tend'][:]
w_ls = nc_fid.variables['w_ls'][:]
u_g = nc_fid.variables['u_g'][:]
v_g = nc_fid.variables['v_g'][:]
dT_dt_rad_forc = nc_fid.variables['dT_dt_rad_forc'][:]
h_advec_thil = nc_fid.variables['h_advec_thil'][:]
h_advec_qt = nc_fid.variables['h_advec_qt'][:]
T_s = nc_fid.variables['T_s'][:]
pres_s = nc_fid.variables['pres_s'][:]
lhf = nc_fid.variables['lhf'][:]
shf = nc_fid.variables['shf'][:]
tau_u = nc_fid.variables['tau_u'][:]
tau_v = nc_fid.variables['tau_v'][:]
cldcov = nc_fid.variables['cldcov'][:]
sw_rad_heating_rate = nc_fid.variables['sw_rad_heating_rate'][:]
lw_rad_heating_rate = nc_fid.variables['lw_rad_heating_rate'][:]
precip = nc_fid.variables['precip'][:]
rain = nc_fid.variables['rain'][:]
dT_dt_lwrad = nc_fid.variables['dT_dt_lwrad'][:]
dT_dt_swrad = nc_fid.variables['dT_dt_swrad'][:]
dT_dt_PBL = nc_fid.variables['dT_dt_PBL'][:]
dT_dt_deepconv = nc_fid.variables['dT_dt_deepconv'][:]
dT_dt_shalconv = nc_fid.variables['dT_dt_shalconv'][:]
dT_dt_micro = nc_fid.variables['dT_dt_micro'][:]
dq_dt_PBL = nc_fid.variables['dq_dt_PBL'][:]
dq_dt_deepconv = nc_fid.variables['dq_dt_deepconv'][:]
dq_dt_shalconv = nc_fid.variables['dq_dt_shalconv'][:]
dq_dt_micro = nc_fid.variables['dq_dt_micro'][:]

#read in LES data
LES_zt = LES_fid.variables['zt'][:]
LES_time = LES_fid.variables['time'][:]
LES_snapshot = LES_fid.variables['snapshot'][:]
LES_u = LES_fid.variables['u'][:]
LES_v = LES_fid.variables['v'][:]
LES_thil = LES_fid.variables['thil'][:]
LES_qt = LES_fid.variables['qt'][:]
LES_ql = LES_fid.variables['ql'][:]
LES_qr = LES_fid.variables['qr'][:]
LES_cld_frc = LES_fid.variables['cld_frc'][:]
LES_z_base = LES_fid.variables['z_base'][:]
LES_z_top = LES_fid.variables['z_top'][:]
LES_shf = LES_fid.variables['shf'][:]
LES_lhf = LES_fid.variables['lhf'][:]
LES_lwp = LES_fid.variables['lwp'][:]
LES_sfc_p = LES_fid.variables['sfc_p'][:]

LES_time_hours = LES_time/3600.0
LES_rho = pres_i[0,0,0]/(287.0*T_s[0])
LES_lhf = LES_lhf*LES_rho*2.5E6
LES_shf = LES_shf*LES_rho*1004.0

#adjust units
time_hours = time/3600.0

pres_l = 1.0E-2*pres_l
pres_i = 1.0E-2*pres_i
pres_s = 1.0E-2*pres_s
qv = 1.0E3*qv
qc = 1.0E3*qc
qt = qv + qc
thil = (1000.0/pres_l)**(287.0/1004.0)*T - 2.5E6/1004.0*(qc*1.0E-3)
qv_force_tend = 1.0E3*qv_force_tend*86400.0 # now in g/kg/day
T_force_tend = T_force_tend*86400.0 #now in K/day
u_force_tend = u_force_tend*86400.0 #now in m/s/day
v_force_tend = v_force_tend*86400.0 #now in m/s/day
sw_rad_heating_rate = sw_rad_heating_rate*86400.0
lw_rad_heating_rate = lw_rad_heating_rate*86400.0

dT_dt_lwrad = dT_dt_lwrad*86400 #was in units of K/s , now it is K/day
dT_dt_swrad = dT_dt_swrad*86400
dT_dt_PBL = dT_dt_PBL*86400
dT_dt_deepconv = dT_dt_deepconv*86400
dT_dt_shalconv = dT_dt_shalconv*86400
dT_dt_micro = dT_dt_micro*86400
dq_dt_PBL = dq_dt_PBL*86400*1.0E3 #now in g/kg/day
dq_dt_deepconv = dq_dt_deepconv*86400*1.0E3
dq_dt_shalconv = dq_dt_shalconv*86400*1.0E3
dq_dt_micro = dq_dt_micro*86400*1.0E3

dT_dt_lwrad_mean = np.mean(dT_dt_lwrad, (0,2))
dT_dt_swrad_mean = np.mean(dT_dt_swrad, (0,2))
dT_dt_PBL_mean = np.mean(dT_dt_PBL, (0,2))
dT_dt_deepconv_mean = np.mean(dT_dt_deepconv, (0,2))
dT_dt_shalconv_mean = np.mean(dT_dt_shalconv, (0,2))
dT_dt_micro_mean = np.mean(dT_dt_micro, (0,2))
T_force_tend_mean = np.mean(T_force_tend, (0,2))

dq_dt_PBL_mean = np.mean(dq_dt_PBL, (0,2))
dq_dt_deepconv_mean = np.mean(dq_dt_deepconv, (0,2))
dq_dt_shalconv_mean = np.mean(dq_dt_shalconv, (0,2))
dq_dt_micro_mean = np.mean(dq_dt_micro, (0,2))
qv_force_tend_mean = np.mean(qv_force_tend, (0,2))

h_advec_qt = 1.0E3*h_advec_qt

plot_pres_top = 600.0

#global plot options
time_start = 0.0
time_end = time_hours[-1]
time_ticks = 5

y_inverted_val = True
y_log_val = True
y_lim_val = [plot_pres_top, np.amax(pres_l)]
x_ticks_val = [time_start, time_end, time_ticks]
y_ticks_val = [600, 1000, 5]


### need to interpolate data to a regular pressure/time grid; pressure levels are potentially changing... ###

### for now, pressure levels are constant...

pressure_axis = pres_l[0,:,0]
phi_axis = phi_l[1,:,0]
plot_pres_top_index = pressure_axis.size - bisect.bisect_left(pressure_axis[::-1], plot_pres_top)

gspr.plot_profile_multi(phi_axis/9.8, [T_force_tend_mean, dT_dt_lwrad_mean, dT_dt_swrad_mean, dT_dt_PBL_mean, dT_dt_deepconv_mean, dT_dt_shalconv_mean, dT_dt_micro_mean], ['FORCE','LW','SW','PBL','DEEP','SHAL','MICRO'], r'$\frac{\partial T}{\partial t} \left(\frac{K}{day}\right)$', 'Z (m)', plot_dir + 'dT_.eps', y_lim = [0.0, 3000.0])
gspr.plot_profile_multi(phi_axis/9.8, [qv_force_tend_mean, dq_dt_PBL_mean, dq_dt_deepconv_mean, dq_dt_shalconv_mean, dq_dt_micro_mean], ['FORCE','PBL','DEEP','SHAL','MICRO'], r'$\frac{\partial q_v}{\partial t} \left(\frac{g/kg}{day}\right)$ ', 'Z (m)', plot_dir + 'dq_.eps', y_lim = [0.0, 3000.0])

### SCM comparison to LES ###
# get time indices that correspond to the LES snapshots
t_index = []
for t in LES_snapshot:
    t_index.append(np.where(time_hours == t))

for t in range(len(t_index)):
    gspr.plot_profile_compare(phi_axis/9.8, u[t_index[t],:,0].flatten(), LES_u[:,t,:], LES_zt, r'$u(m/s)$', 'Z (m)', plot_dir + 'u_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0])
    gspr.plot_profile_compare(phi_axis/9.8, v[t_index[t],:,0].flatten(), LES_v[:,t,:], LES_zt, r'$v(m/s)$', 'Z (m)', plot_dir + 'v_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0])
    gspr.plot_profile_compare(phi_axis/9.8, qt[t_index[t],:,0].flatten(), LES_qt[:,t,:], LES_zt, r'$q_t(g/kg)$', 'Z (m)', plot_dir + 'qt_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0])
    gspr.plot_profile_compare(phi_axis/9.8, thil[t_index[t],:,0].flatten(), LES_thil[:,t,:], LES_zt, r'$\theta_l (K)$', 'Z (m)', plot_dir + 'thil_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0], x_lim=[280.0, 320.0])
    gspr.plot_profile_compare(phi_axis/9.8, qc[t_index[t],:,0].flatten(), LES_ql[:,t,:], LES_zt, r'$q_c(g/kg)$', 'Z (m)', plot_dir + 'qc_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0])
    gspr.plot_profile_compare(phi_axis/9.8, cldcov[t_index[t],:,0].flatten(), LES_cld_frc[:,t,:], LES_zt, 'cloud fraction', 'Z (m)', plot_dir + 'cld_frc_comp_' + str(int(LES_snapshot[t])) + 'hr.eps', y_lim = [0.0, 3000.0])

gspr.plot_time_series_compare(time_hours, lhf, LES_time_hours, LES_lhf, 'time (h)', 'lhf', plot_dir + 'lhf_comp.eps')
gspr.plot_time_series_compare(time_hours, shf, LES_time_hours, LES_shf, 'time (h)', 'shf', plot_dir + 'shf_comp.eps')

#gspr.plot_profile_compare(phi_axis/9.8, qv[t_index[0],:,0].flatten(), LES_qt[:,0,:], LES_zt, r'$q_v(g/kg)$', 'Z (m)', plot_dir + 'test_qv_comp_8hr.eps', y_lim = [0.0, 3000.0])
#gspr.plot_profile_compare(phi_axis/9.8, qv[t_index[1],:,0].flatten(), LES_qt[:,1,:], LES_zt, r'$q_v(g/kg)$', 'Z (m)', plot_dir + 'test_qv_comp_19hr.eps', y_lim = [0.0, 3000.0])
#gspr.plot_profile_compare(phi_axis/9.8, qv[t_index[2],:,0].flatten(), LES_qt[:,2,:], LES_zt, r'$q_v(g/kg)$', 'Z (m)', plot_dir + 'test_qv_comp_36hr.eps', y_lim = [0.0, 3000.0])



#forcing

gspr.plot_time_series(time_hours, T_s, 'time (h)', r'$T_s (K)$', plot_dir + 'T_surf_ts.eps')
gspr.plot_time_series(time_hours, pres_s, 'time (h)', r'$p_s (hPa)$', plot_dir + 'p_surf_ts.eps')

gspr.plot_time_series(time_hours, lhf, 'time (h)', 'lhf', plot_dir + 'lhf.eps')
gspr.plot_time_series(time_hours, shf, 'time (h)', 'shf', plot_dir + 'shf.eps')
gspr.plot_time_series(time_hours, tau_u, 'time (h)', r'$\tau_u$', plot_dir + 'tau_u.eps')
gspr.plot_time_series(time_hours, tau_v, 'time (h)', r'$\tau_v$', plot_dir + 'tau_v.eps')
gspr.plot_time_series(time_hours, precip, 'time (h)', 'precip', plot_dir + 'precip.eps')
gspr.plot_time_series(time_hours, rain, 'time (h)', 'rain', plot_dir + 'rain.eps')

plot_info = [[qv_force_tend, r'$q_v$ forcing tendency (g/kg/day)',plot_dir + 'qv_force_tend_cs.eps'],[T_force_tend, r'$T$ forcing tendency (K/day)',plot_dir + 'T_force_tend_cs.eps'],[u_force_tend, r'$u$ forcing tendency (m/s/day)',plot_dir + 'u_force_tend_cs.eps'],[v_force_tend, r'$v$ forcing tendency (m/s/day)',plot_dir + 'v_force_tend_cs.eps'],[w_ls, r'$w$ large scale (m/s)',plot_dir + 'w_ls_cs.eps']]
for plot in plot_info:
    gspr.contour_plot_firl(time_hours, pressure_axis, np.transpose(plot[0][:,:,0]), np.amin(plot[0][:,0:plot_pres_top_index,0]), np.amax(plot[0][:,0:plot_pres_top_index,0]), plot[1], 'time (h)', 'p (hPa)', plot[2], xticks=x_ticks_val, yticks=y_ticks_val, y_inverted=y_inverted_val, y_log = y_log_val, y_lim = y_lim_val)

#state vars

#pressure as vertical coordinate
#plot_info = [[qv, r'$q_v(g/kg)$',plot_dir + 'qv_cs.eps'],[qc, r'$q_c(g/kg)$',plot_dir + 'qc_cs.eps'],[T, r'$T(K)$',plot_dir + 'T_cs.eps'],[u, r'$u(m/s)$',plot_dir + 'u_cs.eps'],[v, r'$v(m/s)$',plot_dir + 'v_cs.eps'],[cldcov, 'cloud fraction',plot_dir + 'cld_frc_cs.eps']]

#for plot in plot_info:
#    gspr.contour_plot_firl(time_hours, pressure_axis, np.transpose(plot[0][:,:,0]), np.amin(plot[0][:,0:plot_pres_top_index,0]), np.amax(plot[0][:,0:plot_pres_top_index,0]), plot[1], 'time (h)', 'p (hPa)', plot[2], xticks=x_ticks_val, yticks=y_ticks_val, y_inverted=y_inverted_val, y_log = y_log_val, y_lim = y_lim_val)

#geopotential height as vertical coordinate
plot_info = [[qv, r'$q_v(g/kg)$',plot_dir + 'qv_cs_phi.eps'],
    [qc, r'$q_c(g/kg)$',plot_dir + 'qc_cs_phi.eps'],
    [T, r'$T(K)$',plot_dir + 'T_cs_phi.eps'],
    [u, r'$u(m/s)$',plot_dir + 'u_cs_phi.eps'],
    [v, r'$v(m/s)$',plot_dir + 'v_cs_phi.eps'],
    [cldcov, 'cloud fraction',plot_dir + 'cld_frc_phi_cs.eps'],
    [sw_rad_heating_rate + lw_rad_heating_rate, 'rad. heating rate',plot_dir + 'rad_heating_rate_phi_cs.eps'],
    [T_force_tend, 'dT dt forcing',plot_dir + 'dT_dt_forcing_phi_cs.eps'],
    [dT_dt_lwrad, 'dT dt lwrad',plot_dir + 'dT_dt_lwrad_phi_cs.eps'],
    [dT_dt_swrad, 'dT dt swrad',plot_dir + 'dT_dt_swrad_phi_cs.eps'],
    [dT_dt_PBL, 'dT dt PBL',plot_dir + 'dT_dt_PBL_phi_cs.eps'],
    [dT_dt_deepconv, 'dT dt deepconv',plot_dir + 'dT_dt_deepconv_phi_cs.eps'],
    [dT_dt_shalconv, 'dT dt shalconv',plot_dir + 'dT_dt_shalconv_phi_cs.eps'],
    [dT_dt_micro, 'dT dt micro',plot_dir + 'dT_dt_micro_phi_cs.eps'],
    [qv_force_tend, 'dq dt forcing',plot_dir + 'dq_dt_forcing_phi_cs.eps'],
    [dq_dt_PBL, 'dq dt PBL',plot_dir + 'dq_dt_PBL_phi_cs.eps'],
    [dq_dt_deepconv, 'dq dt deepconv',plot_dir + 'dq_dt_deepconv_phi_cs.eps'],
    [dq_dt_shalconv, 'dq dt shalconv',plot_dir + 'dq_dt_shalconv_phi_cs.eps'],
    [dq_dt_micro, 'dq dt micro',plot_dir + 'dq_dt_micro_phi_cs.eps']]

for plot in plot_info:
    gspr.contour_plot_firl(time_hours, phi_axis/9.8, np.transpose(plot[0][:,:,0]), np.amin(plot[0][:,0:plot_pres_top_index,0]), np.amax(plot[0][:,0:plot_pres_top_index,0]), plot[1], 'time (h)', 'Z (m)', plot[2], xticks=x_ticks_val, y_lim = [0.0, 4000.0])


#gspr.plot_profile(pres_l[time_step], qv_force_tend[time_step], r'$q_v(g/kg)$', 'p (hPa)', plot_dir + 'test_qv_force_tend.eps', y_inverted=True, y_log = True, y_lim = [100, max(pres_l[time_step])], yticks=[100, 1000, 10])
#
#gspr.plot_profile(pres_l[time_step], qv[time_step], r'$q_v(g/kg)$', 'p (hPa)', 'test_qv.eps', y_inverted=True, y_log = True, y_lim = [100, max(pres_l[time_step])], yticks=[100, 1000, 10])
#gspr.plot_profile(pres_l[time_step], T[time_step], r'$T(K)$', 'p (hPa)', 'test_T.eps', y_inverted=True, y_log = True, y_lim = [100, max(pres_l[time_step])], yticks=[100, 1000, 10])
#
#print time.shape, pres_l[0,:,0].shape, qv[:,:,0].shape
#gspr.contour_plot_firl(time, pres_l[0,:,0], np.transpose(qv[:,:,0]), 0.0, np.amax(qv), 0, time[-1], 10, 10000, 110000, 10, 'test', 'time', 'p', 'test_qv_cs.eps', y_inverted=True, y_log = True, y_lim = [600, np.amax(pres_l)], yticks=[600, 1000, 5])
#gspr.contour_plot_firl(time, pres_l[0,:,0], np.transpose(qc[:,:,0]), 0.0, np.amax(qc), 0, time[-1], 10, 10000, 110000, 10, 'test', 'time', 'p', 'test_qc_cs.eps', y_inverted=True, y_log = True, y_lim = [600, max(pres_l[time_step])], yticks=[600, 1000, 5])



nc_fid.close()
LES_fid.close()
