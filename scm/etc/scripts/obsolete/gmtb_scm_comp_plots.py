#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import gmtb_scm_plotting_routines as gspr
import bisect
import argparse
import datetime
import os
import scipy.interpolate
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('output_paths', help='paths to output files', nargs='+')
parser.add_argument('--labels', help='names with which to label the output files', nargs='+')
parser.add_argument('--LES_file', help='optional file holding LES info', nargs=1)
parser.add_argument('--obs_file', help='optional file holding observations', nargs=1)
parser.add_argument('--out_dir', help='optional directory to place plots', nargs=1, default='../plots')

args = parser.parse_args()

#check for proper number of labels associated with output files
if args.labels and (len(args.output_paths) != len(args.labels)):
    print 'This script expected exactly ',len(args.output_paths), 'labels associated with the output paths',args.output_paths,' Exiting...'
    quit()
elif not args.labels:
    labels = []
    for i in range(len(args.output_paths)):
        labels.append('exp_'+str(i+1))
else:
    labels = args.labels

plot_dir = args.out_dir[0]

#open the output files for reading
year = []
month = []
day = []
hour = []
time = []
date = []
pres_l = []
pres_i = []
sigma_l = []
sigma_i = []
phi_l = []
phi_i = []
qv = []
T = []
u = []
v = []
qc = []
qv_force_tend = []
T_force_tend = []
u_force_tend = []
v_force_tend = []
w_ls = []
u_g = []
v_g = []
dT_dt_rad_forc = []
h_advec_thil = []
h_advec_qt = []
T_s = []
pres_s = []
lhf = []
shf = []
tau_u = []
tau_v = []
cldcov = []
sw_rad_heating_rate = []
lw_rad_heating_rate = []
precip = []
rain = []
pwat = []
dT_dt_lwrad = []
dT_dt_swrad = []
dT_dt_PBL = []
dT_dt_deepconv = []
dT_dt_shalconv = []
dT_dt_micro = []
dq_dt_PBL = []
dq_dt_deepconv = []
dq_dt_shalconv = []
dq_dt_micro = []
upd_mf = []
dwn_mf = []
det_mf = []
active_start_index = []
active_end_index = []
suppressed_start_index = []
suppressed_end_index = []
for i in range(len(args.output_paths)):
    nc_fid = Dataset(args.output_paths[i], 'r')

    year.append(nc_fid.variables['init_year'][:])
    month.append(nc_fid.variables['init_month'][:])
    day.append(nc_fid.variables['init_day'][:])
    hour.append(nc_fid.variables['init_hour'][:])

    time.append(nc_fid.variables['time'][:])
    pres_l.append(nc_fid.variables['pres'][:])
    pres_i.append(nc_fid.variables['pres_i'][:])
    sigma_l.append(nc_fid.variables['sigma'][:])
    sigma_i.append(nc_fid.variables['sigma_i'][:])
    phi_l.append(nc_fid.variables['phi'][:])
    phi_i.append(nc_fid.variables['phi_i'][:])
    qv.append(nc_fid.variables['qv'][:])
    T.append(nc_fid.variables['T'][:])
    u.append(nc_fid.variables['u'][:])
    v.append(nc_fid.variables['v'][:])
    qc.append(nc_fid.variables['qc'][:])
    qv_force_tend.append(nc_fid.variables['qv_force_tend'][:]*86400.0*1.0E3)
    T_force_tend.append(nc_fid.variables['T_force_tend'][:]*86400.0)
    u_force_tend.append(nc_fid.variables['u_force_tend'][:]*86400.0)
    v_force_tend.append(nc_fid.variables['v_force_tend'][:]*86400.0)
    w_ls.append(nc_fid.variables['w_ls'][:])
    u_g.append(nc_fid.variables['u_g'][:])
    v_g.append(nc_fid.variables['v_g'][:])
    dT_dt_rad_forc.append(nc_fid.variables['dT_dt_rad_forc'][:])
    h_advec_thil.append(nc_fid.variables['h_advec_thil'][:])
    h_advec_qt.append(nc_fid.variables['h_advec_qt'][:])
    T_s.append(nc_fid.variables['T_s'][:])
    pres_s.append(nc_fid.variables['pres_s'][:])
    lhf.append(nc_fid.variables['lhf'][:])
    shf.append(nc_fid.variables['shf'][:])
    tau_u.append(nc_fid.variables['tau_u'][:])
    tau_v.append(nc_fid.variables['tau_v'][:])
    cldcov.append(nc_fid.variables['cldcov'][:])
    sw_rad_heating_rate.append(nc_fid.variables['sw_rad_heating_rate'][:])
    lw_rad_heating_rate.append(nc_fid.variables['lw_rad_heating_rate'][:])
    precip.append(nc_fid.variables['precip'][:]*3600.0) #convert to mm/hr
    pwat.append(nc_fid.variables['pwat'][:]/(1.0E3)*100.0) #convert to cm
    dT_dt_lwrad.append(nc_fid.variables['dT_dt_lwrad'][:]*86400.0)
    dT_dt_swrad.append(nc_fid.variables['dT_dt_swrad'][:]*86400.0)
    dT_dt_PBL.append(nc_fid.variables['dT_dt_PBL'][:]*86400.0)
    dT_dt_deepconv.append(nc_fid.variables['dT_dt_deepconv'][:]*86400.0)
    dT_dt_shalconv.append(nc_fid.variables['dT_dt_shalconv'][:]*86400.0)
    dT_dt_micro.append(nc_fid.variables['dT_dt_micro'][:]*86400.0)
    dq_dt_PBL.append(nc_fid.variables['dq_dt_PBL'][:]*86400.0*1.0E3)
    dq_dt_deepconv.append(nc_fid.variables['dq_dt_deepconv'][:]*86400.0*1.0E3)
    dq_dt_shalconv.append(nc_fid.variables['dq_dt_shalconv'][:]*86400.0*1.0E3)
    dq_dt_micro.append(nc_fid.variables['dq_dt_micro'][:]*86400.0*1.0E3)
    upd_mf.append(nc_fid.variables['upd_mf'][:])
    dwn_mf.append(nc_fid.variables['dwn_mf'][:])
    det_mf.append(nc_fid.variables['det_mf'][:])
    nc_fid.variables['rain'][:]
    initial_date = datetime.datetime(year[i], month[i], day[i], hour[i], 0, 0, 0)

    #convert times to datetime objects starting from initial date
    date.append(np.array([initial_date + datetime.timedelta(seconds=s) for s in np.int_(time[-1])]))

    #find the indices corresponding to the start and end times of the "active" and suppressed periods
    active_start_index.append(np.where(date[i] == datetime.datetime(2006, 1, 20, 0))[0])
    active_end_index.append(np.where(date[i] == datetime.datetime(2006, 1, 25, 12))[0])
    suppressed_start_index.append(np.where(date[i] == datetime.datetime(2006, 1, 28, 0))[0])
    suppressed_end_index.append(np.where(date[i] == datetime.datetime(2006, 2, 2, 12))[0])

    nc_fid.close()

#open the obs file for reading
if(args.obs_file):
    obs_fid = Dataset(args.obs_file[0], 'r')

    obs_year = obs_fid.variables['year'][:]
    obs_month = obs_fid.variables['month'][:]
    obs_day = obs_fid.variables['day'][:]
    obs_hour = obs_fid.variables['hour'][:]
    obs_time_offset = obs_fid.variables['time_offset'][:]

    obs_date = []
    for i in range(obs_hour.size):
        obs_date.append(datetime.datetime(obs_year[i], obs_month[i], obs_day[i], obs_hour[i], 0, 0, 0))
    obs_date = np.array(obs_date)

    #find the indices corresponding to the start and end times of the "active" and suppressed periods
    obs_active_start_index = np.where(obs_date == datetime.datetime(2006, 1, 20, 0))[0]
    obs_active_end_index = np.where(obs_date == datetime.datetime(2006, 1, 25, 12))[0]
    obs_suppressed_start_index = np.where(obs_date == datetime.datetime(2006, 1, 28, 0))[0]
    obs_suppressed_end_index = np.where(obs_date == datetime.datetime(2006, 2, 2, 12))[0]

    #find the index corresponding to the start of the simulations
    obs_start_index = np.where(obs_date == date[0][0])[0]
    obs_time_offset = obs_time_offset - obs_time_offset[obs_start_index]

    obs_levels = obs_fid.variables['lev'][:] #pressure levels in mb

    obs_cld = obs_fid.variables['cld'][:]/100.0
    obs_T = obs_fid.variables['T'][:]
    obs_q = obs_fid.variables['q'][:]
    obs_u = obs_fid.variables['u'][:]
    obs_v = obs_fid.variables['v'][:]
    obs_precip = obs_fid.variables['prec_srf'][:]
    obs_shf = obs_fid.variables['SH'][:]
    obs_lhf = obs_fid.variables['LH'][:]
    obs_pwat = obs_fid.variables['PW'][:]


    obs_cld_active_profile = np.mean(obs_cld[obs_active_start_index:obs_active_end_index], axis=0)
    obs_cld_suppressed_profile = np.mean(obs_cld[obs_suppressed_start_index:obs_suppressed_end_index], axis=0)
    obs_T_active_profile = np.mean(obs_T[obs_active_start_index:obs_active_end_index], axis=0)
    obs_T_suppressed_profile = np.mean(obs_T[obs_suppressed_start_index:obs_suppressed_end_index], axis=0)
    obs_q_active_profile = np.mean(obs_q[obs_active_start_index:obs_active_end_index], axis=0)
    obs_q_suppressed_profile = np.mean(obs_q[obs_suppressed_start_index:obs_suppressed_end_index], axis=0)

    #obs_cld_active_profile_f = scipy.interpolate.interp1d(obs_levels*100.0, obs_cld_active_profile)

    #ozone_ppb = oz_f(levels[1:])
    #ozone_ppb = np.insert(ozone_ppb, 0, oz_data[0])

#plots
try:
    os.makedirs(plot_dir)
except OSError:
    if not os.path.isdir(plot_dir):
        raise

active_pressure_axis = np.mean(pres_l[0][active_start_index[0]:active_end_index[0],:,0], axis=0)*1.0E-2
suppressed_pressure_axis = np.mean(pres_l[0][suppressed_start_index[0]:suppressed_end_index[0],:,0], axis=0)*1.0E-2

#obs_cld_active_profile_interp = obs_cld_active_profile_f(active_pressure_axis[:-2])

#print obs_cld_active_profile_interp, obs_cld_active_profile_interp.shape

mean_active_cld_profiles = []
mean_suppressed_cld_profiles = []
mean_active_T_profiles = []
mean_suppressed_T_profiles = []
mean_active_q_profiles = []
mean_suppressed_q_profiles = []
mean_active_upd_mf_profiles = []
mean_suppressed_upd_mf_profiles = []
mean_active_dwn_mf_profiles = []
mean_suppressed_dwn_mf_profiles = []
mean_active_det_mf_profiles = []
mean_suppressed_det_mf_profiles = []
mean_active_dT_dt_lwrad_profiles = []
mean_suppressed_dT_dt_lwrad_profiles = []
mean_active_dT_dt_swrad_profiles = []
mean_suppressed_dT_dt_swrad_profiles = []
mean_active_dT_dt_PBL_profiles = []
mean_suppressed_dT_dt_PBL_profiles = []
mean_active_dT_dt_deepconv_profiles = []
mean_suppressed_dT_dt_deepconv_profiles = []
mean_active_dT_dt_shalconv_profiles = []
mean_suppressed_dT_dt_shalconv_profiles = []
mean_active_dT_dt_micro_profiles = []
mean_suppressed_dT_dt_micro_profiles = []
mean_active_dq_dt_PBL_profiles = []
mean_suppressed_dq_dt_PBL_profiles = []
mean_active_dq_dt_deepconv_profiles = []
mean_suppressed_dq_dt_deepconv_profiles = []
mean_active_dq_dt_shalconv_profiles = []
mean_suppressed_dq_dt_shalconv_profiles = []
mean_active_dq_dt_micro_profiles = []
mean_suppressed_dq_dt_micro_profiles = []
mean_active_T_force_tend_profiles = []
mean_suppressed_T_force_tend_profiles = []
mean_active_q_force_tend_profiles = []
mean_suppressed_q_force_tend_profiles = []
for i in range(len(args.output_paths)):
    mean_active_cld_profiles.append(np.mean(cldcov[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_cld_profiles.append(np.mean(cldcov[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_T_profiles.append(np.mean(T[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_T_profiles.append(np.mean(T[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_q_profiles.append(1.0E3*np.mean(qv[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_q_profiles.append(1.0E3*np.mean(qv[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_upd_mf_profiles.append(np.mean(upd_mf[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_upd_mf_profiles.append(np.mean(upd_mf[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dwn_mf_profiles.append(np.mean(dwn_mf[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dwn_mf_profiles.append(np.mean(dwn_mf[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_det_mf_profiles.append(np.mean(det_mf[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_det_mf_profiles.append(np.mean(det_mf[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_lwrad_profiles.append(np.mean(dT_dt_lwrad[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_lwrad_profiles.append(np.mean(dT_dt_lwrad[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_swrad_profiles.append(np.mean(dT_dt_swrad[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_swrad_profiles.append(np.mean(dT_dt_swrad[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_PBL_profiles.append(np.mean(dT_dt_PBL[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_PBL_profiles.append(np.mean(dT_dt_PBL[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_deepconv_profiles.append(np.mean(dT_dt_deepconv[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_deepconv_profiles.append(np.mean(dT_dt_deepconv[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_shalconv_profiles.append(np.mean(dT_dt_shalconv[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_shalconv_profiles.append(np.mean(dT_dt_shalconv[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dT_dt_micro_profiles.append(np.mean(dT_dt_micro[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dT_dt_micro_profiles.append(np.mean(dT_dt_micro[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dq_dt_PBL_profiles.append(np.mean(dq_dt_PBL[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dq_dt_PBL_profiles.append(np.mean(dq_dt_PBL[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dq_dt_deepconv_profiles.append(np.mean(dq_dt_deepconv[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dq_dt_deepconv_profiles.append(np.mean(dq_dt_deepconv[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dq_dt_shalconv_profiles.append(np.mean(dq_dt_shalconv[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dq_dt_shalconv_profiles.append(np.mean(dq_dt_shalconv[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_dq_dt_micro_profiles.append(np.mean(dq_dt_micro[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_dq_dt_micro_profiles.append(np.mean(dq_dt_micro[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_T_force_tend_profiles.append(np.mean(T_force_tend[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_T_force_tend_profiles.append(np.mean(T_force_tend[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))
    mean_active_q_force_tend_profiles.append(np.mean(qv_force_tend[i][active_start_index[i]:active_end_index[i],:,0], axis=0))
    mean_suppressed_q_force_tend_profiles.append(np.mean(qv_force_tend[i][suppressed_start_index[i]:suppressed_end_index[i],:,0], axis=0))

#profile controls
plot_pres_top = 100.0
y_inverted_val = True
y_log_val = False
y_lim_val = [plot_pres_top, np.amax(active_pressure_axis)]
y_ticks_val = [plot_pres_top, 1000, 10]

time_start = 0.0
time_end = time[0][-1]/3600.0
time_ticks = 10
x_ticks_val = [time_start, time_end, time_ticks]

gspr.plot_profile_multi(active_pressure_axis, mean_active_cld_profiles, labels, 'cloud fraction', 'p (hPa)', plot_dir + '/comp_active_cld_frc.eps', obs_z = obs_levels, obs_values=obs_cld_active_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_cld_profiles, labels, 'cloud fraction', 'p (hPa)', plot_dir + '/comp_suppressed_cld_frc.eps', obs_z = obs_levels, obs_values=obs_cld_suppressed_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, mean_active_T_profiles, labels, r'T (K)', 'p (hPa)', plot_dir + '/comp_active_T.eps', obs_z = obs_levels, obs_values=obs_T_active_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_T_profiles, labels, r'$T$ (K)', 'p (hPa)', plot_dir + '/comp_suppressed_T.eps', obs_z = obs_levels, obs_values=obs_T_suppressed_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, mean_active_q_profiles, labels, r'$q_v$ (g/kg)', 'p (hPa)', plot_dir + '/comp_active_q_v.eps', obs_z = obs_levels, obs_values=obs_q_active_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_q_profiles, labels, r'$q_v$ (g/kg)', 'p (hPa)', plot_dir + '/comp_suppressed_q_v.eps', obs_z = obs_levels, obs_values=obs_q_suppressed_profile, y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, mean_active_upd_mf_profiles, labels, 'updraft mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_active_upd_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_upd_mf_profiles, labels, 'updraft mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_suppressed_upd_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, mean_active_dwn_mf_profiles, labels, 'downdraft mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_active_dwn_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_dwn_mf_profiles, labels, 'downdraft mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_suppressed_dwn_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, mean_active_det_mf_profiles, labels, 'detrainment mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_active_det_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, mean_suppressed_det_mf_profiles, labels, 'detrainment mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_suppressed_det_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_upd_mf_profiles,mean_active_dwn_mf_profiles,mean_active_det_mf_profiles], [['upd','dwn','det'],labels], 'conv. mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_active_conv_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_upd_mf_profiles,mean_suppressed_dwn_mf_profiles,mean_suppressed_det_mf_profiles], [['upd','dwn','det'],labels], 'conv. mass flux' + r'$(kg/m^2/s)$', 'p (hPa)', plot_dir + '/comp_suppressed_conv_mf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_T_force_tend_profiles, mean_active_dT_dt_PBL_profiles,mean_active_dT_dt_deepconv_profiles,mean_active_dT_dt_shalconv_profiles, mean_active_dT_dt_micro_profiles, mean_active_dT_dt_lwrad_profiles,mean_active_dT_dt_swrad_profiles,], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_active_dT.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_T_force_tend_profiles,mean_suppressed_dT_dt_PBL_profiles,mean_suppressed_dT_dt_deepconv_profiles,mean_suppressed_dT_dt_shalconv_profiles, mean_suppressed_dT_dt_micro_profiles,mean_suppressed_dT_dt_lwrad_profiles,mean_suppressed_dT_dt_swrad_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dT.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_q_force_tend_profiles,mean_active_dq_dt_PBL_profiles,mean_active_dq_dt_deepconv_profiles,mean_active_dq_dt_shalconv_profiles, mean_active_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_active_dq.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_q_force_tend_profiles,mean_suppressed_dq_dt_PBL_profiles,mean_suppressed_dq_dt_deepconv_profiles,mean_suppressed_dq_dt_shalconv_profiles, mean_suppressed_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dq.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_T_force_tend_profiles, mean_active_dT_dt_PBL_profiles,mean_active_dT_dt_deepconv_profiles,mean_active_dT_dt_shalconv_profiles, mean_active_dT_dt_micro_profiles, mean_active_dT_dt_lwrad_profiles,mean_active_dT_dt_swrad_profiles,], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_active_dT_ctl.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=0)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_T_force_tend_profiles, mean_active_dT_dt_PBL_profiles,mean_active_dT_dt_deepconv_profiles,mean_active_dT_dt_shalconv_profiles, mean_active_dT_dt_micro_profiles, mean_active_dT_dt_lwrad_profiles,mean_active_dT_dt_swrad_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_active_dT_gf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=1)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_T_force_tend_profiles,mean_suppressed_dT_dt_PBL_profiles,mean_suppressed_dT_dt_deepconv_profiles,mean_suppressed_dT_dt_shalconv_profiles, mean_suppressed_dT_dt_micro_profiles,mean_suppressed_dT_dt_lwrad_profiles,mean_suppressed_dT_dt_swrad_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dT_ctl.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=0)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_T_force_tend_profiles,mean_suppressed_dT_dt_PBL_profiles,mean_suppressed_dT_dt_deepconv_profiles,mean_suppressed_dT_dt_shalconv_profiles, mean_suppressed_dT_dt_micro_profiles,mean_suppressed_dT_dt_lwrad_profiles,mean_suppressed_dT_dt_swrad_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP','LW','SW'],labels], r'$\frac{\partial T}{\partial t}$' + r'$(K/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dT_gf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=1)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_q_force_tend_profiles,mean_active_dq_dt_PBL_profiles,mean_active_dq_dt_deepconv_profiles,mean_active_dq_dt_shalconv_profiles, mean_active_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_active_dq_ctl.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=0)
gspr.plot_profile_multi(active_pressure_axis, [mean_active_q_force_tend_profiles,mean_active_dq_dt_PBL_profiles,mean_active_dq_dt_deepconv_profiles,mean_active_dq_dt_shalconv_profiles, mean_active_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_active_dq_gf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=1)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_q_force_tend_profiles,mean_suppressed_dq_dt_PBL_profiles,mean_suppressed_dq_dt_deepconv_profiles,mean_suppressed_dq_dt_shalconv_profiles, mean_suppressed_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dq_ctl.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=0)
gspr.plot_profile_multi(suppressed_pressure_axis, [mean_suppressed_q_force_tend_profiles,mean_suppressed_dq_dt_PBL_profiles,mean_suppressed_dq_dt_deepconv_profiles,mean_suppressed_dq_dt_shalconv_profiles, mean_suppressed_dq_dt_micro_profiles], [['FORCE','PBL','Deep Con.','Shal Con.','MP'],labels], r'$\frac{\partial q}{\partial t}$' + r'$(g/kg/day)$', 'p (hPa)', plot_dir + '/comp_suppressed_dq_gf.eps', y_lim=y_lim_val, y_inverted=y_inverted_val, y_log=y_log_val, yticks=y_ticks_val, freeze_axis=1)

#time-series controls
precip_RM = []
for i in range(len(args.output_paths)):
    precip_RM.append(pd.rolling_mean(precip[i][:,0], 36, center=False))

precip_active_time_series = []
precip_suppressed_time_series = []
lhf_active_time_series = []
lhf_suppressed_time_series = []
shf_active_time_series = []
shf_suppressed_time_series = []
pwat_active_time_series = []
pwat_suppressed_time_series = []
for i in range(len(args.output_paths)):
    precip_active_time_series.append(precip_RM[i][active_start_index[i]:active_end_index[i]])
    precip_suppressed_time_series.append(precip_RM[i][suppressed_start_index[i]:suppressed_end_index[i]])
    lhf_active_time_series.append(lhf[i][active_start_index[i]:active_end_index[i]])
    lhf_suppressed_time_series.append(lhf[i][suppressed_start_index[i]:suppressed_end_index[i]])
    shf_active_time_series.append(shf[i][active_start_index[i]:active_end_index[i]])
    shf_suppressed_time_series.append(shf[i][suppressed_start_index[i]:suppressed_end_index[i]])
    pwat_active_time_series.append(pwat[i][active_start_index[i]:active_end_index[i]])
    pwat_suppressed_time_series.append(pwat[i][suppressed_start_index[i]:suppressed_end_index[i]])


gspr.plot_time_series_multi(time[0][active_start_index[0]:active_end_index[0]]/3600.0, precip_active_time_series, labels, 'time (h)', 'precip (mm/hr)', plot_dir + '/comp_active_precip.eps', obs_time = obs_time_offset[obs_active_start_index:obs_active_end_index]/3600.0, obs_values = obs_precip[obs_active_start_index:obs_active_end_index])
gspr.plot_time_series_multi(time[0][suppressed_start_index[0]:suppressed_end_index[0]]/3600.0, precip_suppressed_time_series, labels, 'time (h)', 'precip (mm/hr)', plot_dir + '/comp_suppressed_precip.eps', obs_time = obs_time_offset[obs_suppressed_start_index:obs_suppressed_end_index]/3600.0, obs_values = obs_precip[obs_suppressed_start_index:obs_suppressed_end_index])
gspr.plot_time_series_multi(time[0][active_start_index[0]:active_end_index[0]]/3600.0, lhf_active_time_series , labels, 'time (h)', 'latent heat flux (W/m2)', plot_dir + '/comp_active_lhf.eps', obs_time = obs_time_offset[obs_active_start_index:obs_active_end_index]/3600.0, obs_values = obs_lhf[obs_active_start_index:obs_active_end_index])
gspr.plot_time_series_multi(time[0][suppressed_start_index[0]:suppressed_end_index[0]]/3600.0, lhf_suppressed_time_series , labels, 'time (h)', 'latent heat flux (W/m2)', plot_dir + '/comp_suppressed_lhf.eps', obs_time = obs_time_offset[obs_suppressed_start_index:obs_suppressed_end_index]/3600.0, obs_values = obs_lhf[obs_suppressed_start_index:obs_suppressed_end_index])
gspr.plot_time_series_multi(time[0][active_start_index[0]:active_end_index[0]]/3600.0, shf_active_time_series , labels, 'time (h)', 'sensible heat flux (W/m2)', plot_dir + '/comp_active_shf.eps', obs_time = obs_time_offset[obs_active_start_index:obs_active_end_index]/3600.0, obs_values = obs_shf[obs_active_start_index:obs_active_end_index])
gspr.plot_time_series_multi(time[0][suppressed_start_index[0]:suppressed_end_index[0]]/3600.0, shf_suppressed_time_series , labels, 'time (h)', 'sensible heat flux (W/m2)', plot_dir + '/comp_suppressed_shf.eps', obs_time = obs_time_offset[obs_suppressed_start_index:obs_suppressed_end_index]/3600.0, obs_values = obs_shf[obs_suppressed_start_index:obs_suppressed_end_index])
gspr.plot_time_series_multi(time[0][active_start_index[0]:active_end_index[0]]/3600.0, pwat_active_time_series , labels, 'time (h)', 'Precipitable Water (cm)', plot_dir + '/comp_active_pwat.eps', obs_time = obs_time_offset[obs_active_start_index:obs_active_end_index]/3600.0, obs_values = obs_pwat[obs_active_start_index:obs_active_end_index])
gspr.plot_time_series_multi(time[0][suppressed_start_index[0]:suppressed_end_index[0]]/3600.0, pwat_suppressed_time_series , labels, 'time (h)', 'Precipitable Water (cm)', plot_dir + '/comp_suppressed_pwat.eps', obs_time = obs_time_offset[obs_suppressed_start_index:obs_suppressed_end_index]/3600.0, obs_values = obs_pwat[obs_suppressed_start_index:obs_suppressed_end_index])

#contour plots
#model output
plot_info = [[cldcov[0], 'cloud fraction',plot_dir + '/ctl_cld_frc_cs.eps'],[cldcov[1], 'cloud fraction',plot_dir + '/gf_cld_frc_cs.eps'],[qc[0], 'qc',plot_dir + '/ctl_qc_cs.eps'],[qc[1], 'qc',plot_dir + '/gf_qc_cs.eps']] #[qv, r'$q_v(g/kg)$',plot_dir + 'qv_cs.eps'],[qc, r'$q_c(g/kg)$',plot_dir + 'qc_cs.eps'],[T, r'$T(K)$',plot_dir + 'T_cs.eps'],[u, r'$u(m/s)$',plot_dir + 'u_cs.eps'],[v, r'$v(m/s)$',plot_dir + 'v_cs.eps'],

pressure_axis = np.mean(pres_l[0][:,:,0], axis=0)*1.0E-2
plot_pres_top_index = pressure_axis.size - bisect.bisect_left(pressure_axis[::-1], plot_pres_top)

for plot in plot_info:
    gspr.contour_plot_firl(time[0]/3600.0, pressure_axis, np.transpose(plot[0][:,:,0]), np.amin(plot[0][:,0:plot_pres_top_index,0]), np.amax(plot[0][:,0:plot_pres_top_index,0]), plot[1], 'time (h)', 'p (hPa)', plot[2], xticks=x_ticks_val, yticks=y_ticks_val, y_inverted=y_inverted_val, y_log = y_log_val, y_lim = y_lim_val)

#observations

plot_info = [[obs_cld, 'cloud fraction',plot_dir + '/obs_cld_frc_cs.eps']] #[qv, r'$q_v(g/kg)$',plot_dir + 'qv_cs.eps'],[qc, r'$q_c(g/kg)$',plot_dir + 'qc_cs.eps'],[T, r'$T(K)$',plot_dir + 'T_cs.eps'],[u, r'$u(m/s)$',plot_dir + 'u_cs.eps'],[v, r'$v(m/s)$',plot_dir + 'v_cs.eps'],

pressure_axis = obs_levels
plot_pres_top_index = pressure_axis.size - bisect.bisect_left(pressure_axis[::-1], plot_pres_top)

for plot in plot_info:
    gspr.contour_plot_firl(obs_time_offset/3600.0, pressure_axis, np.transpose(plot[0][:,:]), np.amin(plot[0][:,0:plot_pres_top_index]), np.amax(plot[0][:,0:plot_pres_top_index]), plot[1], 'time (h)', 'p (hPa)', plot[2], xticks=x_ticks_val, yticks=y_ticks_val, y_inverted=y_inverted_val, y_log = y_log_val, y_lim = y_lim_val)
