#!/usr/bin/env python
from __future__ import print_function

import sys
PYTHON2 = sys.version_info[0] < 3
if PYTHON2:
    from imp import reload
else:
    from importlib import reload

from configobj import ConfigObj, flatten_errors
from validate import Validator
import argparse
import glob
import os
from netCDF4 import Dataset
import datetime
import numpy as np
import f90nml
import scm_plotting_routines as spr
import pandas as pd
import bisect
import scm_read_obs as sro

Rd = 287.0
Rv = 461.0
g = 9.81
missing_value = -999
missing_soil_levels = 4

plot_ext = '.pdf' #.pdf, .eps, .ps, .png (.png is fastest, but raster)

reload(spr)
reload(sro)

try:
  pd.plotting.register_matplotlib_converters()
except (AttributeError):
  print("Warning: The version of the pandas package you are using may lead to Future Warnings being generated. These can be ignored for now.")

#subroutine for printing progress to the command line
def print_progress(n_complete, n_total):
    print(str(n_complete) + ' of ' + str(n_total) + ' complete: (' + str(100.0*n_complete/float(n_total)) + '%)')

def replace_fill_with_nan(nc_ds, var_name, var, group, time_diag, pres_l, dataset):
    try:
        raw_data = nc_ds.variables[var_name][:]
        raw_data[raw_data == nc_fid.variables[var_name]._FillValue] = np.nan
        var.append(nc_fid.variables[var_name][:])
    except KeyError:
        print('{0} is not in the output file {1}'.format(var_name, dataset))
        print('Missing variables are replaced with {0}'.format(0))
        var.append(np.zeros((len(time_diag[-1]),pres_l[-1].shape[1],pres_l[-1].shape[2])))
    group.append(var_name)
    return[var,group]

#set up command line argument parser to read in name of config file to use
parser = argparse.ArgumentParser()
parser.add_argument('config', help='configuration file for CCPP SCM analysis', nargs=1)
parser.add_argument('-d', '--docker', help='include if scm is being run in a docker container to mount volumes', action='store_true', default=False)

args = parser.parse_args()

#read in the configuration file specified on the command line (check against the configspec.ini file for validation)
config = ConfigObj(args.config[0],configspec="configspec.ini")
validator = Validator()
results = config.validate(validator)

#standardized output for config file validation errors (quit if error)
if results != True:
    for (section_list, key, _) in flatten_errors(config, results):
        if key is not None:
            print('The "%s" key in the section "%s" failed validation' % (key, ', '.join(section_list)))
        else:
            print('The following section was missing:%s ' % ', '.join(section_list))
    print('Since CCPP SCM analysis configuration file did not pass validation, the program will quit.')
    quit()

#get the config file variables
scm_datasets = config['scm_datasets']
scm_datasets_labels = config['scm_datasets_labels']
obs_file = config['obs_file']
obs_compare = config['obs_compare']
plot_dir = config['plot_dir']
time_slices = config['time_slices']
plot_ind_datasets = config['plot_ind_datasets']
time_series_resample = config['time_series_resample']
skill_scores_val = False
bias_val = True

plots = config['plots']
profiles_mean = plots['profiles_mean']
time_series = plots['time_series']
contours = plots['contours']
profiles_mean_multi = plots['profiles_mean_multi']
time_series_multi = plots['time_series_multi']

#keep track of progress
num_base_plots = len(time_slices)*(len(profiles_mean['vars']) + len(profiles_mean_multi) + len(time_series['vars']) + len(time_series_multi) + len(contours['vars']))
num_total_plots = 0
if(plot_ind_datasets):
    num_total_plots = len(scm_datasets)*num_base_plots
if(len(scm_datasets) > 1):
    num_total_plots += num_base_plots - len(time_slices)*len(contours['vars'])

num_plots_completed = 0

#make sure that each time slice is a list of length 5 -- year, month, day, hour, min (append zero for minute if necessary)
for time_slice in time_slices:
    if (len(time_slices[time_slice]['start']) < 5):
        for i in range(5 - len(time_slices[time_slice]['start'])):
            time_slices[time_slice]['start'].append(0)
    if (len(time_slices[time_slice]['end']) < 5):
        for i in range(5 - len(time_slices[time_slice]['end'])):
            time_slices[time_slice]['end'].append(0)

#perform any special checks on the config file data
if len(scm_datasets) != len(scm_datasets_labels):
    print('The number of scm_datasets must match the number of scm_datasets_labels. Quitting...')
    print('scm_datasets = ',scm_datasets)
    print('scm_datasets_labels = ',scm_datasets_labels)
    quit()

#if running in a Docker container, the output is being copied to a different directory within the container 
#and the plots should go into that same (home) directory in order for the volume to be correctly mounted and the host to see them.
if args.docker:
    plot_dir = '/home/'+plot_dir
    for f in scm_datasets:
        f = '/home/'+f

#read in the case name from the case_config namelist (just use first dataset dir namelist)
i = scm_datasets[0].rfind('/')
dir = scm_datasets[0][:i]
for filename in glob.glob(os.path.join(dir, '*.nml')):
    nml = f90nml.read(filename)
    if nml['case_config']:
        case_name = nml['case_config']['case_name']

#initialize lists for output variables
year = []
month = []
day = []
hour = []
minute = []
time_inst = []
time_diag = []
time_swrad = []
time_lwrad = []
time_rad = []
date_inst = []
date_diag = []
date_swrad = []
date_lwrad = []
date_rad = []

pres_l = []
pres_i = []
sigma_l = []
sigma_i = []
pres_s = []
phi_l = []
phi_i = []
qv = []
T = []
u = []
v = []
qc = []
ql = []
qi = []

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
v_advec_thil = []
v_advec_qt = []
T_s = []

soil_T = []
soil_moisture = []
soil_moisture_unfrozen = []
lhf = []
shf = []
tprcp_inst = []
tprcp_rate_inst = []
t2m = []
q2m = []
ustar = []
tsfc = []

tau_u = []
tau_v = []
upd_mf = []
dwn_mf = []
det_mf = []
sfc_up_lw_land     = []
sfc_up_lw_ice      = []
sfc_up_lw_water    = []
sfc_up_sw_dir_nir  = []
sfc_up_sw_dif_nir  = []
sfc_up_sw_dir_vis  = []
sfc_up_sw_dif_vis  = []
sfc_dwn_sw_dir_nir = []
sfc_dwn_sw_dif_nir = []
sfc_dwn_sw_dir_vis = []
sfc_dwn_sw_dif_vis = []
mp_prcp_inst       = []
dcnv_prcp_inst     = []
scnv_prcp_inst     = []
rad_cloud_fraction = []
rad_cloud_lwp      = []
rad_eff_rad_ql     = []
rad_cloud_iwp      = []
rad_eff_rad_qi     = []
rad_cloud_rwp      = []
rad_eff_rad_qr     = []
rad_cloud_swp      = []
rad_eff_rad_qs     = []

sw_rad_heating_rate = []
lw_rad_heating_rate = []

pwat       = []
sfc_dwn_sw = []
sfc_up_sw  = []
sfc_net_sw = []
sfc_dwn_lw = []
gflux = []
u10m = []
v10m = []
hpbl = []
sfc_rad_net_land   = []
sfc_rad_net_ice    = []
sfc_rad_net_water  = []


dT_dt_lwrad = []
dT_dt_swrad = []
dT_dt_pbl = []
dT_dt_deepconv = []
dT_dt_shalconv = []
dT_dt_micro = []
dT_dt_conv = []
dT_dt_ogwd = []
dT_dt_rayleigh = []
dT_dt_cgwd = []
dT_dt_phys = []
dT_dt_nonphys = []
dq_dt_pbl = []
dq_dt_deepconv = []
dq_dt_shalconv = []
dq_dt_micro = []
dq_dt_conv = []
dq_dt_phys = []
dq_dt_nonphys = []
doz_dt_pbl = []
doz_dt_prodloss = []
doz_dt_oz = []
doz_dt_T = []
doz_dt_ovhd = []
doz_dt_phys = []
doz_dt_nonphys = []
du_dt_pbl = []
du_dt_ogwd = []
du_dt_deepconv = []
du_dt_cgwd = []
du_dt_rayleigh = []
du_dt_shalconv = []
du_dt_conv = []
du_dt_phys = []
du_dt_nonphys = []
dv_dt_pbl = []
dv_dt_ogwd = []
dv_dt_deepconv = []
dv_dt_cgwd = []
dv_dt_rayleigh = []
dv_dt_shalconv = []
dv_dt_conv = []
dv_dt_phys = []
dv_dt_nonphys = []

tprcp_accum = []
ice_accum = []
snow_accum = []
graupel_accum = []
conv_prcp_accum = []
tprcp_rate_accum = []
ice_rate_accum = []
snow_rate_accum = []
graupel_rate_accum = []
conv_prcp_rate_accum = []

max_cloud_fraction = []
toa_total_albedo = []
vert_int_lwp_mp = []
vert_int_iwp_mp = []
vert_int_lwp_cf = []
vert_int_iwp_cf = []

sw_up_TOA_tot = []
sw_dn_TOA_tot = []
sw_up_TOA_clr = []
sw_up_sfc_tot = []
sw_dn_sfc_tot = []
sw_up_sfc_clr = []
sw_dn_sfc_clr = []
lw_up_TOA_tot = []
lw_up_TOA_clr = []
lw_up_sfc_tot = []
lw_up_sfc_clr = []
lw_dn_sfc_tot = []
lw_dn_sfc_clr = []
rh = []
rh_500 = []

inst_time_group = []
diag_time_group = []
swrad_time_group = []
lwrad_time_group = []
rad_time_group = []

time_slice_indices_inst = []
time_slice_indices_diag = []
time_slice_indices_swrad = []
time_slice_indices_lwrad = []
time_slice_indices_rad = []
time_slice_labels = []

for i in range(len(scm_datasets)):
    nc_fid = Dataset(scm_datasets[i], 'r')
    nc_fid.set_auto_mask(False)

    #3D vars have dimensions (time, levels, horizontal)

    year.append(nc_fid.variables['init_year'][:])
    month.append(nc_fid.variables['init_month'][:])
    day.append(nc_fid.variables['init_day'][:])
    hour.append(nc_fid.variables['init_hour'][:])
    minute.append(nc_fid.variables['init_minute'][:])

    time_inst.append(nc_fid.variables['time_inst'][:])
    time_diag.append(nc_fid.variables['time_diag'][:])
    time_swrad.append(nc_fid.variables['time_swrad'][:])
    time_lwrad.append(nc_fid.variables['time_lwrad'][:])
    time_rad.append(nc_fid.variables['time_rad'][:])
    
    pres_l.append(nc_fid.variables['pres'][:])
    inst_time_group.append('pres_l')
    
    pres_i.append(nc_fid.variables['pres_i'][:])
    inst_time_group.append('pres_l')
    
    sigma_l.append(nc_fid.variables['sigma'][:])
    inst_time_group.append('sigma')
    
    sigma_i.append(nc_fid.variables['sigma_i'][:])
    inst_time_group.append('sigma_i')
    
    pres_s.append(nc_fid.variables['pres_s'][:])
    inst_time_group.append('pres_s')
    
    #phi_l.append(nc_fid.variables['phi'][:])
    #inst_time_group.append('phi_l')
    
    #phi_i.append(nc_fid.variables['phi_i'][:])
    #inst_time_group.append('phi_i')
    
    qv.append(nc_fid.variables['qv'][:])
    inst_time_group.append('qv')
    
    T.append(nc_fid.variables['T'][:])
    inst_time_group.append('T')
    
    u.append(nc_fid.variables['u'][:])
    inst_time_group.append('u')
    
    v.append(nc_fid.variables['v'][:])
    inst_time_group.append('v')
    
    qc.append(nc_fid.variables['qc'][:])
    inst_time_group.append('qc')
    
    ql.append(nc_fid.variables['ql'][:])
    inst_time_group.append('ql')
    
    qi.append(nc_fid.variables['qi'][:])
    inst_time_group.append('qi')
    
    qv_force_tend.append(nc_fid.variables['qv_force_tend'][:])
    inst_time_group.append('qv_force_tend')
    
    T_force_tend.append(nc_fid.variables['T_force_tend'][:])
    inst_time_group.append('T_force_tend')
    
    u_force_tend.append(nc_fid.variables['u_force_tend'][:])
    inst_time_group.append('u_force_tend')
    
    v_force_tend.append(nc_fid.variables['v_force_tend'][:])
    inst_time_group.append('v_force_tend')
    
    w_ls.append(nc_fid.variables['w_ls'][:])
    inst_time_group.append('w_ls')
    
    u_g.append(nc_fid.variables['u_g'][:])
    inst_time_group.append('u_g')
    
    v_g.append(nc_fid.variables['v_g'][:])
    inst_time_group.append('v_g')
    
    dT_dt_rad_forc.append(nc_fid.variables['dT_dt_rad_forc'][:])
    inst_time_group.append('dT_dt_rad_forc')
    
    h_advec_thil.append(nc_fid.variables['h_advec_thil'][:])
    inst_time_group.append('h_advec_thil')
    
    h_advec_qt.append(nc_fid.variables['h_advec_qt'][:])
    inst_time_group.append('h_advec_qt')
    
    v_advec_thil.append(nc_fid.variables['v_advec_thil'][:])
    inst_time_group.append('v_advec_thil')
    
    v_advec_qt.append(nc_fid.variables['v_advec_qt'][:])
    inst_time_group.append('v_advec_qt')
    
    T_s.append(nc_fid.variables['T_s'][:])
    inst_time_group.append('T_s')
    
    try:
        soil_T.append(nc_fid.variables['soil_T'][:])
        soil_moisture.append(nc_fid.variables['soil_moisture'][:])
        soil_moisture_unfrozen.append(nc_fid.variables['soil_moisture_unfrozen'][:])
    except KeyError:
        print('soil_T, soil_moisture, and/or soil_moisture_frozen are not in the output file {0}'.format(scm_datasets[i]))
        print('Missing variables are replaced with {0}'.format(missing_value))
        soil_T.append(missing_value*np.ones((len(time_inst[-1]),missing_soil_levels)))
        soil_moisture.append(missing_value*np.ones((len(time_inst[-1]),missing_soil_levels)))
        soil_moisture_unfrozen.append(missing_value*np.ones((len(time_inst[-1]),missing_soil_levels)))
    inst_time_group.extend(('soil_T','soil_moisture','soil_moisture_unfrozen'))
    
    lhf.append(nc_fid.variables['lhf'][:])
    inst_time_group.append('lhf')
    
    shf.append(nc_fid.variables['shf'][:])
    inst_time_group.append('shf')
    
    tprcp_inst.append(nc_fid.variables['tprcp_inst'][:])
    inst_time_group.append('tprcp_inst')
    
    tprcp_rate_inst.append(nc_fid.variables['tprcp_rate_inst'][:])
    inst_time_group.append('tprcp_rate_inst')
    
    t2m.append(nc_fid.variables['t2m'][:])
    inst_time_group.append('t2m')
    
    q2m.append(nc_fid.variables['q2m'][:])
    inst_time_group.append('q2m')
    
    ustar.append(nc_fid.variables['ustar'][:])
    inst_time_group.append('ustar')
    
    tsfc.append(nc_fid.variables['tsfc'][:])
    inst_time_group.append('tsfc')
    
    tau_u.append(nc_fid.variables['tau_u'][:])
    inst_time_group.append('tau_u')
    
    tau_v.append(nc_fid.variables['tau_v'][:])
    inst_time_group.append('tau_v')
    
    try:
        upd_mf.append(nc_fid.variables['upd_mf'][:])
    except KeyError:
        print('upd_mf is not in the output file {0}'.format(scm_datasets[i]))
        print('Missing variables are replaced with {0}'.format(missing_value))
        upd_mf.append(missing_value*np.ones((len(time_inst[-1]),pres_l[-1].shape[1],pres_l[-1].shape[2])))
    inst_time_group.append('upd_mf')
    
    dwn_mf.append(nc_fid.variables['dwn_mf'][:])
    inst_time_group.append('dwn_mf')
    
    det_mf.append(nc_fid.variables['det_mf'][:])
    inst_time_group.append('det_mf')
    
    sfc_up_lw_land.append(nc_fid.variables['sfc_up_lw_land'][:])
    inst_time_group.append('sfc_up_lw_land')
    
    sfc_up_lw_ice.append(nc_fid.variables['sfc_up_lw_ice'][:])
    inst_time_group.append('sfc_up_lw_ice')
    
    sfc_up_lw_water.append(nc_fid.variables['sfc_up_lw_water'][:])
    inst_time_group.append('sfc_up_lw_water')
    
    sfc_up_sw_dir_nir.append(nc_fid.variables['sfc_up_sw_dir_nir'][:])
    inst_time_group.append('sfc_up_sw_dir_nir')
    
    sfc_up_sw_dif_nir.append(nc_fid.variables['sfc_up_sw_dif_nir'][:])
    inst_time_group.append('sfc_up_sw_dif_nir')
    
    sfc_up_sw_dir_vis.append(nc_fid.variables['sfc_up_sw_dir_vis'][:])
    inst_time_group.append('sfc_up_sw_dir_vis')
    
    sfc_up_sw_dif_vis.append(nc_fid.variables['sfc_up_sw_dif_vis'][:])
    inst_time_group.append('sfc_up_sw_dif_vis')
    
    sfc_dwn_sw_dir_nir.append(nc_fid.variables['sfc_dwn_sw_dir_nir'][:])
    inst_time_group.append('sfc_dwn_sw_dir_nir')
    
    sfc_dwn_sw_dif_nir.append(nc_fid.variables['sfc_dwn_sw_dif_nir'][:])
    inst_time_group.append('sfc_dwn_sw_dif_nir')
    
    sfc_dwn_sw_dir_vis.append(nc_fid.variables['sfc_dwn_sw_dir_vis'][:])
    inst_time_group.append('sfc_dwn_sw_dir_vis')
    
    sfc_dwn_sw_dif_vis.append(nc_fid.variables['sfc_dwn_sw_dif_vis'][:])
    inst_time_group.append('sfc_dwn_sw_dif_vis')
    
    mp_prcp_inst.append(nc_fid.variables['mp_prcp_inst'][:])
    inst_time_group.append('mp_prcp_inst')
    
    dcnv_prcp_inst.append(nc_fid.variables['dcnv_prcp_inst'][:])
    inst_time_group.append('dcnv_prcp_inst')
    
    scnv_prcp_inst.append(nc_fid.variables['scnv_prcp_inst'][:])
    inst_time_group.append('scnv_prcp_inst')
    
    rad_cloud_fraction.append(nc_fid.variables['rad_cloud_fraction'][:])
    rad_time_group.append('rad_cloud_fraction')
    
    rad_cloud_lwp.append(nc_fid.variables['rad_cloud_lwp'][:])
    rad_time_group.append('rad_cloud_lwp')
    
    rad_eff_rad_ql.append(nc_fid.variables['rad_eff_rad_ql'][:])
    rad_time_group.append('rad_eff_rad_ql')
    
    rad_cloud_iwp.append(nc_fid.variables['rad_cloud_iwp'][:])
    rad_time_group.append('rad_cloud_iwp')
    
    rad_eff_rad_qi.append(nc_fid.variables['rad_eff_rad_qi'][:])
    rad_time_group.append('rad_eff_rad_qi')
    
    rad_cloud_rwp.append(nc_fid.variables['rad_cloud_rwp'][:])
    rad_time_group.append('rad_cloud_rwp')
    
    rad_eff_rad_qr.append(nc_fid.variables['rad_eff_rad_qr'][:])
    rad_time_group.append('rad_eff_rad_qr')
    
    rad_cloud_swp.append(nc_fid.variables['rad_cloud_swp'][:])
    rad_time_group.append('rad_cloud_swp')
    
    rad_eff_rad_qs.append(nc_fid.variables['rad_eff_rad_qs'][:])
    rad_time_group.append('rad_eff_rad_qs')
    
    sw_rad_heating_rate.append(nc_fid.variables['sw_rad_heating_rate'][:])
    swrad_time_group.append('sw_rad_heating_rate')
    
    lw_rad_heating_rate.append(nc_fid.variables['lw_rad_heating_rate'][:])
    lwrad_time_group.append('lw_rad_heating_rate')
    
    pwat.append(nc_fid.variables['pwat'][:]) #convert to cm
    inst_time_group.append('pwat')
    
    sfc_dwn_sw.append(nc_fid.variables['sfc_dwn_sw'][:])
    inst_time_group.append('sfc_dwn_sw')
    
    sfc_up_sw.append(nc_fid.variables['sfc_up_sw'][:])
    inst_time_group.append('sfc_up_sw')
    
    sfc_net_sw.append(nc_fid.variables['sfc_net_sw'][:])
    inst_time_group.append('sfc_net_sw')
    
    sfc_dwn_lw.append(nc_fid.variables['sfc_dwn_lw'][:])
    inst_time_group.append('sfc_dwn_lw')
    
    gflux.append(nc_fid.variables['gflux'][:])
    inst_time_group.append('gflux')
    
    u10m.append(nc_fid.variables['u10m'][:])
    inst_time_group.append('u10m')
    
    v10m.append(nc_fid.variables['v10m'][:])
    inst_time_group.append('v10m')
    
    hpbl.append(nc_fid.variables['hpbl'][:])
    inst_time_group.append('hpbl')
    
    [dT_dt_lwrad, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_lwrad', dT_dt_lwrad, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dT_dt_swrad, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_swrad', dT_dt_swrad, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_pbl, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_pbl', dT_dt_pbl, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_deepconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_deepconv', dT_dt_deepconv, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_shalconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_shalconv', dT_dt_shalconv, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_micro, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_micro', dT_dt_micro, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    dT_dt_conv.append(dT_dt_deepconv[-1] + dT_dt_shalconv[-1])
    diag_time_group.append('dT_dt_conv')
    
    [dT_dt_ogwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_ogwd', dT_dt_ogwd, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dT_dt_rayleigh, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_rayleigh', dT_dt_rayleigh, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_cgwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_cgwd', dT_dt_cgwd, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_phys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_phys', dT_dt_phys, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dT_dt_nonphys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dT_dt_nonphys', dT_dt_nonphys, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dq_dt_pbl, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_pbl', dq_dt_pbl, diag_time_group, time_diag, pres_l, scm_datasets[i])    
        
    [dq_dt_deepconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_deepconv', dq_dt_deepconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dq_dt_shalconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_shalconv', dq_dt_shalconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dq_dt_micro, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_micro', dq_dt_micro, diag_time_group, time_diag, pres_l, scm_datasets[i])    
        
    dq_dt_conv.append(dq_dt_deepconv[-1] + dq_dt_shalconv[-1])
    diag_time_group.append('dq_dt_conv')
    
    [dq_dt_phys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_phys', dq_dt_phys, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dq_dt_nonphys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dq_dt_nonphys', dq_dt_nonphys, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_pbl, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_pbl', doz_dt_pbl, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_prodloss, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_prodloss', doz_dt_prodloss, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_oz, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_oz', doz_dt_oz, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_T, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_T', doz_dt_T, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_ovhd, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_ovhd', doz_dt_ovhd, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [doz_dt_phys, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_phys', doz_dt_phys, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [doz_dt_nonphys, diag_time_group] = replace_fill_with_nan(nc_fid, 'doz_dt_nonphys', doz_dt_nonphys, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [du_dt_pbl, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_pbl', du_dt_pbl, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [du_dt_ogwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_ogwd', du_dt_ogwd, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [du_dt_deepconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_deepconv', du_dt_deepconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [du_dt_cgwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_cgwd', du_dt_cgwd, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [du_dt_rayleigh, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_rayleigh', du_dt_rayleigh, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [du_dt_shalconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_shalconv', du_dt_shalconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    du_dt_conv.append(du_dt_deepconv[-1] + du_dt_shalconv[-1])
    diag_time_group.append('du_dt_conv')
    
    [du_dt_phys, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_phys', du_dt_phys, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [du_dt_nonphys, diag_time_group] = replace_fill_with_nan(nc_fid, 'du_dt_nonphys', du_dt_nonphys, diag_time_group, time_diag, pres_l, scm_datasets[i])
        
    [dv_dt_pbl, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_pbl', dv_dt_pbl, diag_time_group, time_diag, pres_l, scm_datasets[i])    
    
    [dv_dt_ogwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_ogwd', dv_dt_ogwd, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dv_dt_deepconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_deepconv', dv_dt_deepconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dv_dt_cgwd, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_cgwd', dv_dt_cgwd, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dv_dt_rayleigh, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_rayleigh', dv_dt_rayleigh, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dv_dt_shalconv, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_shalconv', dv_dt_shalconv, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    dv_dt_conv.append(dv_dt_deepconv[-1] + dv_dt_shalconv[-1])
    diag_time_group.append('dv_dt_conv')
    
    [dv_dt_phys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_phys', dv_dt_phys, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    [dv_dt_nonphys, diag_time_group] = replace_fill_with_nan(nc_fid, 'dv_dt_nonphys', dv_dt_nonphys, diag_time_group, time_diag, pres_l, scm_datasets[i])
    
    tprcp_accum.append(nc_fid.variables['tprcp_accum'][:])
    diag_time_group.append('tprcp_accum')
    
    ice_accum.append(nc_fid.variables['ice_accum'][:])
    diag_time_group.append('ice_accum')
    
    snow_accum.append(nc_fid.variables['snow_accum'][:])
    diag_time_group.append('snow_accum')
    
    graupel_accum.append(nc_fid.variables['graupel_accum'][:])
    diag_time_group.append('graupel_accum')
    
    conv_prcp_accum.append(nc_fid.variables['conv_prcp_accum'][:])
    diag_time_group.append('conv_prcp_accum')
    
    tprcp_rate_accum.append(nc_fid.variables['tprcp_rate_accum'][:])
    diag_time_group.append('tprcp_rate_accum')
    
    ice_rate_accum.append(nc_fid.variables['ice_rate_accum'][:])
    diag_time_group.append('ice_rate_accum')
    
    snow_rate_accum.append(nc_fid.variables['snow_rate_accum'][:])
    diag_time_group.append('snow_rate_accum')
    
    graupel_rate_accum.append(nc_fid.variables['graupel_rate_accum'][:])
    diag_time_group.append('graupel_rate_accum')
    
    conv_prcp_rate_accum.append(nc_fid.variables['conv_prcp_rate_accum'][:])
    diag_time_group.append('conv_prcp_rate_accum')
    
    sfc_rad_net_land.append((sfc_dwn_sw[-1] - sfc_up_sw[-1]) + (sfc_dwn_lw[-1] - sfc_up_lw_land[-1]))
    inst_time_group.append('sfc_rad_net_land')
    
    sfc_rad_net_ice.append((sfc_dwn_sw[-1] - sfc_up_sw[-1]) + (sfc_dwn_lw[-1] - sfc_up_lw_ice[-1]))
    inst_time_group.append('sfc_rad_net_ice')
    
    sfc_rad_net_water.append((sfc_dwn_sw[-1] - sfc_up_sw[-1]) + (sfc_dwn_lw[-1] - sfc_up_lw_water[-1]))
    inst_time_group.append('sfc_rad_net_water')
    
    max_cloud_fraction.append(nc_fid.variables['max_cloud_fraction'][:])
    rad_time_group.append('max_cloud_fraction')
    
    toa_total_albedo.append(nc_fid.variables['toa_total_albedo'][:])
    rad_time_group.append('toa_total_albedo')
    
    vert_int_lwp_mp.append(nc_fid.variables['vert_int_lwp_mp'][:])
    rad_time_group.append('vert_int_lwp_mp')
    
    vert_int_iwp_mp.append(nc_fid.variables['vert_int_iwp_mp'][:])
    rad_time_group.append('vert_int_iwp_mp')
    
    vert_int_lwp_cf.append(nc_fid.variables['vert_int_lwp_cf'][:])
    rad_time_group.append('vert_int_lwp_cf')
    
    vert_int_iwp_cf.append(nc_fid.variables['vert_int_iwp_cf'][:])
    rad_time_group.append('vert_int_iwp_cf')
    
    # sw_up_TOA_tot.append(nc_fid.variables['sw_up_TOA_tot'][:])
    # sw_dn_TOA_tot.append(nc_fid.variables['sw_dn_TOA_tot'][:])
    # sw_up_TOA_clr.append(nc_fid.variables['sw_up_TOA_clr'][:])
    # sw_up_sfc_tot.append(nc_fid.variables['sw_up_sfc_tot'][:])
    # sw_dn_sfc_tot.append(nc_fid.variables['sw_dn_sfc_tot'][:])
    # sw_up_sfc_clr.append(nc_fid.variables['sw_up_sfc_clr'][:])
    # sw_dn_sfc_clr.append(nc_fid.variables['sw_dn_sfc_clr'][:])
    # lw_up_TOA_tot.append(nc_fid.variables['lw_up_TOA_tot'][:])
    # lw_up_TOA_clr.append(nc_fid.variables['lw_up_TOA_clr'][:])
    # lw_up_sfc_tot.append(nc_fid.variables['lw_up_sfc_tot'][:])
    # lw_up_sfc_clr.append(nc_fid.variables['lw_up_sfc_clr'][:])
    # lw_dn_sfc_tot.append(nc_fid.variables['lw_dn_sfc_tot'][:])
    # lw_dn_sfc_clr.append(nc_fid.variables['lw_dn_sfc_clr'][:])

    initial_date = datetime.datetime(year[i], month[i], day[i], hour[i], minute[i], 0, 0)
    
    #convert times to datetime objects starting from initial date
    date_inst.append(np.array([initial_date + datetime.timedelta(seconds=int(s)) for s in time_inst[-1]]))
    date_diag.append(np.array([initial_date + datetime.timedelta(seconds=int(s)) for s in time_diag[-1]]))
    date_swrad.append(np.array([initial_date + datetime.timedelta(seconds=int(s)) for s in time_swrad[-1]]))
    date_lwrad.append(np.array([initial_date + datetime.timedelta(seconds=int(s)) for s in time_lwrad[-1]]))
    date_rad.append(np.array([initial_date + datetime.timedelta(seconds=int(s)) for s in time_rad[-1]]))
    nc_fid.close()

    #calculate diagnostic values from model output
    # e_s = 6.1078*np.exp(17.2693882*(T[-1] - 273.16)/(T[-1] - 35.86))*100.0 #Tetens formula produces e_s in mb (convert to Pa)
    # e = qv[-1]*pres_l[-1]/(qv[-1] + (Rd/Rv)*(1.0 - qv[-1])) #compute vapor pressure from specific humidity
    # rh.append(np.clip(e/e_s, 0.0, 1.0))
    #
    # rh_500_kj = np.zeros((pres_l[-1].shape[0],pres_l[-1].shape[2]))
    # lwp_kj = np.zeros((pres_l[-1].shape[0],pres_l[-1].shape[2]))
    # for j in range(pres_l[-1].shape[0]): #loop over times
    #     for k in range(pres_l[-1].shape[2]): #loop over hor. index
    #         index_500 = np.where(pres_l[-1][j,:,k]*0.01 < 500.0)[0][0]
    #         lifrac = (pres_l[-1][j,index_500-1,k] - 50000.0)/(pres_l[-1][j,index_500-1,k] - pres_l[-1][j,index_500,k])
    #         rh_500_kj[j,k] = rh[-1][j,index_500-1,k] + lifrac*(rh[-1][j,index_500,k] - rh[-1][j,index_500-1,k])
    #         #print index_500, pres_l[-1][j,index_500,k], pres_l[-1][j,index_500-1,k], rh_500_kj, rh[-1][j,index_500,k], rh[-1][j,index_500-1,k]
    #         temp_lwp = 0.0
    #         for l in range(pres_l[-1].shape[1]):
    #             temp_lwp += qc[-1][j,l,k]*(pres_i[-1][j,l,k]-pres_i[-1][j,l+1,k])/g
    #         lwp_kj[j,k] = temp_lwp
    # rh_500.append(rh_500_kj)
    # lwp.append(lwp_kj)
    #
    # rad_net_srf.append((sw_dn_sfc_tot[-1] - sw_up_sfc_tot[-1]) + (lw_dn_sfc_tot[-1] - lw_up_sfc_tot[-1]))

#only keep unique elements in time group lists
inst_time_group = list(set(inst_time_group))
diag_time_group = list(set(diag_time_group))
swrad_time_group = list(set(swrad_time_group))
lwrad_time_group = list(set(lwrad_time_group))
rad_time_group = list(set(rad_time_group))

time_h_inst = [x/3600.0 for x in time_inst]
time_h_diag = [x/3600.0 for x in time_diag]
time_h_swrad = [x/3600.0 for x in time_swrad]
time_h_lwrad = [x/3600.0 for x in time_lwrad]
time_h_rad = [x/3600.0 for x in time_rad]

#find the indices corresponding to the start and end times of the time slices defined in the config file
for time_slice in time_slices:
    time_slice_labels.append(time_slice)
    start_date = datetime.datetime(time_slices[time_slice]['start'][0], time_slices[time_slice]['start'][1],time_slices[time_slice]['start'][2], time_slices[time_slice]['start'][3],time_slices[time_slice]['start'][4])
    end_date = datetime.datetime(time_slices[time_slice]['end'][0], time_slices[time_slice]['end'][1],time_slices[time_slice]['end'][2], time_slices[time_slice]['end'][3],time_slices[time_slice]['end'][4])
    
    valid_inst_indices = np.where((date_inst[0] >= start_date) & (date_inst[0] <= end_date))
    start_date_index_inst = valid_inst_indices[0][0]
    end_date_index_inst = valid_inst_indices[0][-1]
    
    valid_diag_indices = np.where((date_diag[0] >= start_date) & (date_diag[0] <= end_date))
    start_date_index_diag = valid_diag_indices[0][0]
    end_date_index_diag = valid_diag_indices[0][-1]
    
    valid_swrad_indices = np.where((date_swrad[0] >= start_date) & (date_swrad[0] <= end_date))
    
    if (len(valid_swrad_indices[0]) > 1):
        start_date_index_swrad = valid_swrad_indices[0][0]
        end_date_index_swrad = valid_swrad_indices[0][-1]
        time_slice_indices_swrad.append([start_date_index_swrad, end_date_index_swrad])
    else:
        start_date_index_swrad = valid_swrad_indices[0][0]
        end_date_index_swrad = valid_swrad_indices[0][-1]
        time_slice_indices_swrad.append([start_date_index_swrad, end_date_index_swrad+1])    
    
    valid_lwrad_indices = np.where((date_lwrad[0] >= start_date) & (date_lwrad[0] <= end_date))
    if (len(valid_lwrad_indices[0]) > 1):
        start_date_index_lwrad = valid_lwrad_indices[0][0]
        end_date_index_lwrad = valid_lwrad_indices[0][-1]
        time_slice_indices_lwrad.append([start_date_index_lwrad, end_date_index_lwrad])
    else:
        start_date_index_lwrad = valid_lwrad_indices[0][0]
        end_date_index_lwrad = valid_lwrad_indices[0][-1]
        time_slice_indices_lwrad.append([start_date_index_lwrad, end_date_index_lwrad+1])
    
    valid_rad_indices = np.where((date_rad[0] >= start_date) & (date_rad[0] <= end_date))
        
    if (len(valid_rad_indices[0]) > 1):
        start_date_index_rad = valid_rad_indices[0][0]
        end_date_index_rad = valid_rad_indices[0][-1]
        time_slice_indices_rad.append([start_date_index_rad, end_date_index_rad])
    else:
        start_date_index_rad = valid_rad_indices[0][0]
        end_date_index_rad = valid_rad_indices[0][-1]
        time_slice_indices_rad.append([start_date_index_rad, end_date_index_rad+1])
    
    time_slice_indices_inst.append([start_date_index_inst, end_date_index_inst])
    time_slice_indices_diag.append([start_date_index_diag, end_date_index_diag])

#fill the obs_dict by calling the appropriate observation file read routine
if(obs_compare and obs_file):
    if(case_name.strip() == 'twpice'):
        obs_dict = sro.read_twpice_obs(obs_file, time_slices, date_inst)
    elif(case_name.strip() == 'arm_sgp_summer_1997_A'):
        obs_dict = sro.read_arm_sgp_summer_1997_obs(obs_file, time_slices, date_inst)
    elif('LASSO' in case_name.strip()):
        obs_dict = sro.read_LASSO_obs(obs_file, time_slices, date_inst)
    elif('gabls3' in case_name.strip()):
        obs_dict = sro.read_gabls3_obs(obs_file, time_slices, date_inst)
    elif('MOSAiC' in case_name.strip()):
        obs_dict = sro.read_MOSAiC_obs(obs_file, time_slices, date_inst)

try:
    os.makedirs(plot_dir)
except OSError:
    if not os.path.isdir(plot_dir):
        raise

if(profiles_mean['vert_axis'] in locals()):
    if(obs_compare):
        obs_vert_axis = np.array(obs_dict[profiles_mean['vert_axis']])
    vert_axis_label_pm = profiles_mean['vert_axis_label']
    y_inverted_val_pm = profiles_mean['y_inverted']
    y_log_val_pm = profiles_mean['y_log']
    y_min_option_pm = profiles_mean['y_min_option']
    y_max_option_pm = profiles_mean['y_max_option']

if(contours['vert_axis'] in locals()):
    vert_axis_label_c = contours['vert_axis_label']
    y_inverted_val_c = contours['y_inverted']
    y_log_val_c = contours['y_log']
    y_min_option_c = contours['y_min_option']
    y_max_option_c = contours['y_max_option']


#make plots for each dataset individually (colors should stay the same for each dataset [using color_index keyword])
if(plot_ind_datasets):
    for i in range(len(scm_datasets)):
        #loop through the time slices
        for j in range(len(time_slice_labels)):
            ind_dir = plot_dir + scm_datasets_labels[i] + '/' + time_slice_labels[j]

            #make the directory for the current dataset
            try:
                os.makedirs(ind_dir)
            except OSError:
                if not os.path.isdir(ind_dir):
                    raise



            ### Mean Profiles ###

            #check if the specified vertical axis data exists and create the vertical axis for the mean_profile plots
            if(profiles_mean['vert_axis'] in locals()):
                vert_axis_data = np.array(locals()[profiles_mean['vert_axis']][i])
                #this is not technically correct, but it won't blow up
                vert_axis = np.mean(vert_axis_data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:], (0,2))
                if y_min_option_pm == 'min':
                    y_min_val = np.amin(vert_axis)
                elif (y_min_option_pm == 'max'):
                    y_min_val = np.amax(vert_axis)
                else:
                    y_min_val = profiles_mean['y_min']
                if(y_max_option_pm == 'min'):
                    y_max_val = np.amin(vert_axis)
                elif(y_max_option_pm == 'max'):
                    y_max_val = np.amax(vert_axis)
                else:
                    y_max_val = profiles_mean['y_max']
                y_lim_val = [y_min_val, y_max_val]

                #plot mean profiles
                for k in range(len(profiles_mean['vars'])):
                    #get the python variable associated with the vars listed in the config file
                    if(profiles_mean['vars'][k] in locals()):
                        data = np.array(locals()[profiles_mean['vars'][k]][i])
                        if (profiles_mean['vars'][k] in inst_time_group):
                            #print("{} in time group inst".format(profiles_mean['vars'][k]))
                            data_time_slice = data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:]
                        elif (profiles_mean['vars'][k] in diag_time_group):
                            #print("{} in time group diag".format(profiles_mean['vars'][k]))
                            data_time_slice = data[time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:,:]
                        elif (profiles_mean['vars'][k] in swrad_time_group):
                            #print("{} in time group swrad".format(profiles_mean['vars'][k]))
                            data_time_slice = data[time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:,:]
                        elif (profiles_mean['vars'][k] in lwrad_time_group):
                            #print("{} in time group lwrad".format(profiles_mean['vars'][k]))
                            data_time_slice = data[time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:,:]
                        elif (profiles_mean['vars'][k] in rad_time_group):
                            #print("{} in time group rad".format(profiles_mean['vars'][k]))
                            data_time_slice = data[time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:,:]
                        else:
                            print("{} not found in any time groups".format(profiles_mean['vars'][k]))
                            exit()
                        label = profiles_mean['vars_labels'][k]
                        if profiles_mean['conversion_factor']:
                            conversion_factor = profiles_mean['conversion_factor'][k]
                        else:
                            conversion_factor = 1.0
                        
                        #mean profile is obtained by averaging over dimensions 0 and 2
                        mean_data = np.mean(data_time_slice, (0,2))

                        #spr.plot_profile(vert_axis, mean_data, label, vert_axis_label, ind_dir + '/profiles_mean_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val, y_log=y_log_val, y_lim=y_lim_val)
                        if(obs_compare and profiles_mean['vars'][k] in obs_dict):
                            obs_data = np.array(obs_dict[profiles_mean['vars'][k]])
                            obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1],:]

                            obs_mean_data = np.mean(obs_data_time_slice, (0))

                            spr.plot_profile_multi(vert_axis, [mean_data], [scm_datasets_labels[i]], label, vert_axis_label_pm, ind_dir + '/profiles_mean_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, obs_z=obs_vert_axis, obs_values=obs_mean_data, line_type='color', color_index=i, conversion_factor=conversion_factor)

                        else:
                            spr.plot_profile_multi(vert_axis, [mean_data], [scm_datasets_labels[i]], label, vert_axis_label_pm, ind_dir + '/profiles_mean_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, line_type='color', color_index=i, conversion_factor=conversion_factor)

                    else:
                        print('The variable ' + profiles_mean['vars'][k] + ' found in ' + args.config[0] + ' in the profiles_mean section is invalid.')
                    num_plots_completed += 1
                    print_progress(num_plots_completed, num_total_plots)


                for multiplot in profiles_mean_multi:
                    #multiplot_vars = profiles_mean_multi[multiplot]['vars']

                    #check if all the vars exist
                    all_vars_exist = True
                    for l in range(len(profiles_mean_multi[multiplot]['vars'])):
                        if (profiles_mean_multi[multiplot]['vars'][l] not in locals()):
                            all_vars_exist = False
                            print('The variable ' + profiles_mean_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' in the ' + multiplot + ' section of profiles_mean_multi is invalid.')

                    if all_vars_exist:
                        if profiles_mean_multi[multiplot]['conversion_factor']:
                            conversion_factor = profiles_mean_multi[multiplot]['conversion_factor']
                        else:
                            conversion_factor = 1.0
                            
                        data_list = []
                        for l in range(len(profiles_mean_multi[multiplot]['vars'])):
                            data = np.array(locals()[profiles_mean_multi[multiplot]['vars'][l]])
                            if (profiles_mean_multi[multiplot]['vars'][l] in inst_time_group):
                                #print("{} in time group inst".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in diag_time_group):
                                #print("{} in time group diag".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in swrad_time_group):
                                #print("{} in time group swrad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in lwrad_time_group):
                                #print("{} in time group lwrad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in rad_time_group):
                                #print("{} in time group rad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:,:]
                            else:
                                print("{} not found in any time groups".format(profiles_mean_multi[multiplot]['vars'][l]))
                                exit()
                            data_list.append(np.mean(data_time_slice, (0,2)))

                        spr.plot_profile_multi(vert_axis, data_list, profiles_mean_multi[multiplot]['vars_labels'], profiles_mean_multi[multiplot]['x_label'], vert_axis_label_pm, ind_dir + '/profiles_mean_multi_' + multiplot + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, line_type='style', color_index=i, conversion_factor=conversion_factor)

                    num_plots_completed += 1
                    print_progress(num_plots_completed, num_total_plots)
            else:
                print('The variable ' + profiles_mean['vert_axis'] + ' found in ' + args.config[0] + ' in the profiles_mean section is invalid.')

            ### Time-Series ###
            for k in range(len(time_series['vars'])):
                if(time_series['vars'][k] in locals()):
                    if(time_series['levels']):                        
                        if (time_series['levels'][k] == -999):
                            data = np.array(locals()[time_series['vars'][k]][i])
                            plot_name = time_series['vars'][k]
                        else:
                            #minus one to account for 0-based indexing (also, 0 passed in is not working)
                            try:
                                data = np.array(locals()[time_series['vars'][k]][i])[:,time_series['levels'][k]-1]
                                plot_name = time_series['vars'][k] + "_L" + str(time_series['levels'][k])
                            except (IndexError):
                                print('The variable ' + time_series['vars'][k] + ' found in ' + args.config[0] + ' was given an invalid vertical level index: ' + str(time_series['levels'][k]))
                                continue
                    else:
                        data = np.array(locals()[time_series['vars'][k]][i])
                        plot_name = time_series['vars'][k]
                    if (time_series['vars'][k] in inst_time_group):
                        #print("{} in time group inst".format(time_series['vars'][k]))
                        data_time_slice = data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                        time_h_slice = time_h_inst[i][time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                        
                        if time_series_resample:
                            data_delta_seconds = time_inst[i][time_slice_indices_inst[j][1]] - time_inst[0][time_slice_indices_inst[j][1]-1]
                            #create date range for the model data
                            data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                            data_time_slice_periods = time_slice_indices_inst[j][1] - time_slice_indices_inst[j][0]
                            data_date_range = pd.date_range(start=date_inst[i][time_slice_indices_inst[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same

                    elif (time_series['vars'][k] in diag_time_group):
                        #print("{} in time group diag".format(time_series['vars'][k]))
                        data_time_slice = data[time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                        time_h_slice = time_h_diag[i][time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                        
                        if time_series_resample:
                            data_delta_seconds = time_diag[i][time_slice_indices_diag[j][1]] - time_diag[0][time_slice_indices_diag[j][1]-1]
                            #create date range for the model data
                            data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                            data_time_slice_periods = time_slice_indices_diag[j][1] - time_slice_indices_diag[j][0]
                            data_date_range = pd.date_range(start=date_diag[i][time_slice_indices_diag[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                    elif (time_series['vars'][k] in swrad_time_group):
                        #print("{} in time group swrad".format(time_series['vars'][k]))
                        data_time_slice = data[time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                        time_h_slice = time_h_swrad[i][time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                        
                        if time_series_resample:
                            data_delta_seconds = time_swrad[i][time_slice_indices_swrad[j][1]] - time_swrad[0][time_slice_indices_swrad[j][1]-1]
                            #create date range for the model data
                            data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                            data_time_slice_periods = time_slice_indices_swrad[j][1] - time_slice_indices_swrad[j][0]
                            data_date_range = pd.date_range(start=date_swrad[i][time_slice_indices_swrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                    elif (time_series['vars'][k] in lwrad_time_group):
                        #print("{} in time group lwrad".format(time_series['vars'][k]))
                        data_time_slice = data[time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                        time_h_slice = time_h_lwrad[i][time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                        
                        if time_series_resample:
                            data_delta_seconds = time_lwrad[i][time_slice_indices_lwrad[j][1]] - time_lwrad[0][time_slice_indices_lwrad[j][1]-1]
                            #create date range for the model data
                            data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                            data_time_slice_periods = time_slice_indices_lwrad[j][1] - time_slice_indices_lwrad[j][0]
                            data_date_range = pd.date_range(start=date_lwrad[i][time_slice_indices_lwrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                    elif (time_series['vars'][k] in rad_time_group):
                        #print("{} in time group rad".format(time_series['vars'][k]))
                        data_time_slice = data[time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                        time_h_slice = time_h_rad[i][time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                        
                        if time_series_resample:
                            data_delta_seconds = time_rad[i][time_slice_indices_rad[j][1]] - time_rad[0][time_slice_indices_rad[j][1]-1]
                            #create date range for the model data
                            data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                            data_time_slice_periods = time_slice_indices_rad[j][1] - time_slice_indices_rad[j][0]
                            data_date_range = pd.date_range(start=date_rad[i][time_slice_indices_rad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                    else:
                        print("{} not found in any time groups".format(time_series['vars'][k]))
                        exit()

                    label = time_series['vars_labels'][k]
                    if time_series['conversion_factor']:
                        conversion_factor = time_series['conversion_factor'][k]
                    else:
                        conversion_factor = 1.0

                    if(obs_compare and time_series['vars'][k] in obs_dict):
                        #get the corresponding obs data
                        obs_data = np.array(obs_dict[time_series['vars'][k]])

                        #isolate the subset of the obs data for the current time slice
                        obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]

                        #print obs_dict['time_slice_indices'][j][0], obs_dict['time_slice_indices'][j][1]

                        if time_series_resample:
                            #determine whether obs data frequency matches model output frequency
                            obs_delta_seconds = obs_dict['time'][obs_dict['time_slice_indices'][j][1]] - obs_dict['time'][obs_dict['time_slice_indices'][j][1]-1]
                            if(obs_delta_seconds > data_delta_seconds):
                                #need to downsample the model data to the obs data frequency

                                #create date range for the obs data
                                obs_data_dateoffset = pd.DateOffset(seconds=int(obs_delta_seconds))
                                obs_time_slice_periods = obs_dict['time_slice_indices'][j][1] - obs_dict['time_slice_indices'][j][0]
                                obs_date_range = pd.date_range(start=obs_dict['date'][obs_dict['time_slice_indices'][j][0]], periods=obs_time_slice_periods, freq=obs_data_dateoffset)

                                #print data_date_range, obs_date_range

                                resample_string = str(int(obs_delta_seconds)) + 'S'
                                data_time_slice_series = pd.Series(data_time_slice[:,0], index = data_date_range)
                                data_time_slice_series_rs = data_time_slice_series.resample(resample_string).mean()

                                #print obs_data_time_slice.shape, obs_date_range.shape, data_time_slice_series_rs.shape

                                spr.plot_time_series_multi(obs_date_range, [data_time_slice_series_rs], [scm_datasets_labels[i]], 'date', label, ind_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_date_range, obs_values = obs_data_time_slice, line_type='color', color_index=i, conversion_factor=conversion_factor)
                            elif(obs_delta_seconds < data_delta_seconds):
                                print('The case where observations are more frequent than model output has not been implmented yet... ')
                            else:
                                obs_time_time_slice = obs_time_h[obs_time_slice_indices[j][0]:obs_time_slice_indices[j][1]]
                                spr.plot_time_series_multi(time_h_slice, [data_time_slice], [scm_datasets_labels[i]], 'time (h)', label, ind_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, line_type='color', color_index=i, conversion_factor=conversion_factor)
                        else:
                            obs_time_time_slice = obs_dict['time_h'][obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                            spr.plot_time_series_multi(time_h_slice, [data_time_slice], [scm_datasets_labels[i]], 'time (h)', label, ind_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, line_type='color', color_index=i, conversion_factor=conversion_factor)
                    else:
                        spr.plot_time_series_multi(time_h_slice, [data_time_slice], [scm_datasets_labels[i]], 'time (h)', label, ind_dir + '/time_series_' + plot_name + plot_ext, line_type='color', color_index=i, conversion_factor=conversion_factor)
                else:
                    print('The variable ' + time_series['vars'][k] + ' found in ' + args.config[0] + ' in the time_series section is invalid.')
                num_plots_completed += 1
                print_progress(num_plots_completed, num_total_plots)

            for multiplot in time_series_multi:
                #check if all the vars exist

                all_vars_exist = True
                for l in range(len(time_series_multi[multiplot]['vars'])):
                    if (time_series_multi[multiplot]['vars'][l] not in locals()):
                        all_vars_exist = False
                        print('The variable ' + time_series_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' in the ' + multiplot + ' section of time_series_multi is invalid.')

                if all_vars_exist:
                    if time_series_multi[multiplot]['conversion_factor']:
                        conversion_factor = time_series_multi[multiplot]['conversion_factor']
                    else:
                        conversion_factor = 1.0
                        
                    data_list = []
                    for l in range(len(time_series_multi[multiplot]['vars'])):
                        if(time_series_multi[multiplot]['levels']):                        
                            if (time_series_multi[multiplot]['levels'][l] == -999):
                                data = np.array(locals()[time_series_multi[multiplot]['vars'][l]][i])
                                #plot_name = time_series['vars'][k]
                            else:
                                #minus one to account for 0-based indexing (also, 0 passed in is not working)
                                try:
                                    data = np.array(locals()[time_series_multi[multiplot]['vars'][l]][i])[:,time_series_multi[multiplot]['levels'][l]-1]
                                    #plot_name = time_series['vars'][k] + "_L" + str(time_series['levels'][k])
                                except (IndexError):
                                    print('The variable ' + time_series_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' was given an invalid vertical level index: ' + str(time_series_multi[multiplot]['levels'][l]))
                                    continue
                        else:
                            data = np.array(locals()[time_series_multi[multiplot]['vars'][l]][i])

                        if (time_series_multi[multiplot]['vars'][l] in inst_time_group):
                            #print("{} in time group inst".format(time_series_multi[multiplot]['vars'][l]))
                            data_time_slice = data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                            time_h_slice = time_h_inst[i][time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_inst[i][time_slice_indices_inst[j][1]] - time_inst[0][time_slice_indices_inst[j][1]-1]
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_inst[j][1] - time_slice_indices_inst[j][0]
                                data_date_range = pd.date_range(start=date_inst[i][time_slice_indices_inst[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same

                        elif (time_series_multi[multiplot]['vars'][l] in diag_time_group):
                            #print("{} in time group diag".format(time_series_multi[multiplot]['vars'][l]))
                            data_time_slice = data[time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                            time_h_slice = time_h_diag[i][time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_diag[i][time_slice_indices_diag[j][1]] - time_diag[0][time_slice_indices_diag[j][1]-1]
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_diag[j][1] - time_slice_indices_diag[j][0]
                                data_date_range = pd.date_range(start=date_diag[i][time_slice_indices_diag[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                        elif (time_series_multi[multiplot]['vars'][l] in swrad_time_group):
                            #print("{} in time group swrad".format(time_series_multi[multiplot]['vars'][l]))
                            data_time_slice = data[time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                            time_h_slice = time_h_swrad[i][time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_swrad[i][time_slice_indices_swrad[j][1]] - time_swrad[0][time_slice_indices_swrad[j][1]-1]
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_swrad[j][1] - time_slice_indices_swrad[j][0]
                                data_date_range = pd.date_range(start=date_swrad[i][time_slice_indices_swrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                        elif (time_series_multi[multiplot]['vars'][l] in lwrad_time_group):
                            #print("{} in time group lwrad".format(time_series_multi[multiplot]['vars'][l]))
                            data_time_slice = data[time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                            time_h_slice = time_h_lwrad[i][time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_lwrad[i][time_slice_indices_lwrad[j][1]] - time_lwrad[0][time_slice_indices_lwrad[j][1]-1]
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_lwrad[j][1] - time_slice_indices_lwrad[j][0]
                                data_date_range = pd.date_range(start=date_lwrad[i][time_slice_indices_lwrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                        elif (time_series_multi[multiplot]['vars'][l] in rad_time_group):
                            #print("{} in time group rad".format(time_series_multi[multiplot]['vars'][l]))
                            data_time_slice = data[time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                            time_h_slice = time_h_rad[i][time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_rad[i][time_slice_indices_rad[j][1]] - time_rad[0][time_slice_indices_rad[j][1]-1]
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_rad[j][1] - time_slice_indices_rad[j][0]
                                data_date_range = pd.date_range(start=date_rad[i][time_slice_indices_rad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                        else:
                            print("{} not found in any time groups".format(time_series_multi[multiplot]['vars'][l]))
                            exit()
                        data_list.append(data_time_slice)

                    if(obs_compare and time_series_multi[multiplot]['obs_var'] in obs_dict):
                        #get the corresponding obs data
                        obs_data = np.array(obs_dict[time_series_multi[multiplot]['obs_var']])

                        #isolate the subset of the obs data for the current time slice
                        obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                        if time_series_resample:
                            #determine whether obs data frequency matches model output frequency
                            obs_delta_seconds = obs_dict['time'][obs_dict['time_slice_indices'][j][1]] - obs_dict['time'][obs_dict['time_slice_indices'][j][1]-1]
                            if(obs_delta_seconds > data_delta_seconds):
                                #need to downsample the model data to the obs data frequency

                                #create date range for the obs data
                                obs_data_dateoffset = pd.DateOffset(seconds=int(obs_delta_seconds))
                                obs_time_slice_periods = obs_dict['time_slice_indices'][j][1] - obs_dict['time_slice_indices'][j][0]
                                obs_date_range = pd.date_range(start=obs_dict['date'][obs_dict['time_slice_indices'][j][0]], periods=obs_time_slice_periods, freq=obs_data_dateoffset)

                                resample_string = str(int(obs_delta_seconds)) + 'S'
                                data_list_rs = []
                                for l in range(len(data_list)):
                                    data_time_slice_series = pd.Series(data_list[l][:,0], index = data_date_range)
                                    data_time_slice_series_rs = data_time_slice_series.resample(resample_string).mean()
                                    data_list_rs.append(data_time_slice_series_rs)

                                spr.plot_time_series_multi(obs_date_range, data_list_rs, time_series_multi[multiplot]['vars_labels'], 'date', time_series_multi[multiplot]['y_label'], ind_dir + '/time_series_multi_' + multiplot + plot_ext, obs_time = obs_date_range, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'], line_type='style', color_index=i, conversion_factor=conversion_factor)
                            elif(obs_delta_seconds < data_delta_seconds):
                                print('The case where observations are more frequent than model output has not been implmented yet... ')
                            else:
                                obs_time_time_slice = obs_time_h[obs_time_slice_indices[j][0]:obs_time_slice_indices[j][1]]
                                spr.plot_time_series_multi(time_h_slice, data_list, time_series_multi[multiplot]['vars_labels'], 'time (h)', time_series_multi[multiplot]['y_label'], ind_dir + '/time_series_multi_' + multiplot + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'], line_type='style', color_index=i, conversion_factor=conversion_factor)
                        else:
                            obs_time_time_slice = obs_dict['time_h'][obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                            spr.plot_time_series_multi(time_h_slice, data_list, time_series_multi[multiplot]['vars_labels'], 'time (h)', time_series_multi[multiplot]['y_label'], ind_dir + '/time_series_multi_' + multiplot + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'], line_type='style', color_index=i, conversion_factor=conversion_factor)
                    else:
                        spr.plot_time_series_multi(time_h_slice, data_list, time_series_multi[multiplot]['vars_labels'], 'time (h)', time_series_multi[multiplot]['y_label'], ind_dir + '/time_series_multi_' + multiplot + plot_ext, line_type='style', color_index=i, conversion_factor=conversion_factor)

                num_plots_completed += 1
                print_progress(num_plots_completed, num_total_plots)

            ### Contour Plots ###
            x_ticks_num = contours['x_ticks_num']
            y_ticks_num = contours['y_ticks_num']

            if(contours['vert_axis'] in locals()):
                vert_axis_data = np.array(locals()[contours['vert_axis']][i])
                #this is not technically correct, but it won't blow up
                vert_axis = np.mean(vert_axis_data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:], (0,2))
                if y_min_option_c == 'min':
                    y_min_val = np.amin(vert_axis)
                elif (y_min_option_c == 'max'):
                    y_min_val = np.amax(vert_axis)
                else:
                    y_min_val = contours['y_min']
                if(y_max_option_c == 'min'):
                    y_max_val = np.amin(vert_axis)
                elif(y_max_option_c == 'max'):
                    y_max_val = np.amax(vert_axis)
                else:
                    y_max_val = contours['y_max']
                y_lim_val = [y_min_val, y_max_val]

                vert_axis_top_index = vert_axis.size - bisect.bisect_left(vert_axis[::-1], y_min_val)
                y_ticks_val = [y_min_val, y_max_val, y_ticks_num]

                for k in range(len(contours['vars'])):
                    if(contours['vars'][k] in locals()):
                        data = np.array(locals()[contours['vars'][k]][i])
                        if (contours['vars'][k] in inst_time_group):
                            #print("{} in time group inst".format(contours['vars'][k]))
                            x_ticks_val = [time_h_inst[i][time_slice_indices_inst[j][0]], time_h_inst[i][time_slice_indices_inst[j][1]], x_ticks_num]
                            data_time_slice = data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:]
                            time_h_slice = time_h_inst[i][time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                        elif (contours['vars'][k] in diag_time_group):
                            #print("{} in time group diag".format(contours['vars'][k]))
                            x_ticks_val = [time_h_diag[i][time_slice_indices_diag[j][0]], time_h_diag[i][time_slice_indices_diag[j][1]], x_ticks_num]
                            data_time_slice = data[time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:,:]
                            time_h_slice = time_h_diag[i][time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                        elif (contours['vars'][k] in swrad_time_group):
                            #print("{} in time group swrad".format(contours['vars'][k]))
                            x_ticks_val = [time_h_swrad[i][time_slice_indices_swrad[j][0]], time_h_swrad[i][time_slice_indices_swrad[j][1]], x_ticks_num]
                            data_time_slice = data[time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:,:]
                            time_h_slice = time_h_swrad[i][time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                        elif (contours['vars'][k] in lwrad_time_group):
                            #print("{} in time group lwrad".format(contours['vars'][k]))
                            x_ticks_val = [time_h_lwrad[i][time_slice_indices_lwrad[j][0]], time_h_lwrad[i][time_slice_indices_lwrad[j][1]], x_ticks_num]
                            data_time_slice = data[time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:,:]
                            time_h_slice = time_h_lwrad[i][time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                        elif (contours['vars'][k] in rad_time_group):
                            #print("{} in time group rad".format(contours['vars'][k]))
                            x_ticks_val = [time_h_rad[i][time_slice_indices_rad[j][0]], time_h_rad[i][time_slice_indices_rad[j][1]], x_ticks_num]
                            data_time_slice = data[time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:,:]
                            time_h_slice = time_h_rad[i][time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                        else:
                            print("{} not found in any time groups".format(contours['vars'][k]))
                            exit()
                        
                        label = contours['vars_labels'][k]
                        
                        if contours['conversion_factor']:
                            conversion_factor = contours['conversion_factor'][k]
                        else:
                            conversion_factor = 1.0
                        
                        spr.contour_plot_firl(time_h_slice, vert_axis, np.transpose(data_time_slice[:,:,0]), np.amin(data_time_slice[:,0:vert_axis_top_index,0]), np.amax(data_time_slice[:,0:vert_axis_top_index,0]), label, 'time (h)', vert_axis_label_c, ind_dir + '/contour_' + contours['vars'][k] + plot_ext, xticks=x_ticks_val, yticks=y_ticks_val, y_inverted=y_inverted_val_c, y_log = y_log_val_c, y_lim = y_lim_val, conversion_factor=conversion_factor)
                    else:
                        print('The variable ' + contours['vars'][k] + ' found in ' + args.config[0] + ' in the contours section is invalid.')

                    num_plots_completed += 1
                    print_progress(num_plots_completed, num_total_plots)
            else:
                print('The variable ' + contours['vert_axis'] + ' found in ' + args.config[0] + ' in the contours section is invalid.')


#make comparison plots
if(len(scm_datasets) > 1):
    #loop through the time slices
    for j in range(len(time_slice_labels)):
        comp_dir = plot_dir + 'comp/' + time_slice_labels[j]

        #make the directory for the current dataset
        try:
            os.makedirs(comp_dir)
        except OSError:
            if not os.path.isdir(comp_dir):
                raise

        #check if the specified vertical axis data exists and create the vertical axis for the mean_profile plots
        if(profiles_mean['vert_axis'] in locals()):
            vert_axis_data = np.array(locals()[profiles_mean['vert_axis']][0])
            #not technically correct, but doesn't blow up for now
            vert_axis = np.mean(vert_axis_data[time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:], (0,2))
            if y_min_option_pm == 'min':
                y_min_val = np.amin(vert_axis)
            elif (y_min_option_pm == 'max'):
                y_min_val = np.amax(vert_axis)
            else:
                y_min_val = profiles_mean['y_min']
            if(y_max_option_pm == 'min'):
                y_max_val = np.amin(vert_axis)
            elif(y_max_option_pm == 'max'):
                y_max_val = np.amax(vert_axis)
            else:
                y_max_val = profiles_mean['y_max']
            y_lim_val = [y_min_val, y_max_val]
            #plot mean profiles
            for k in range(len(profiles_mean['vars'])):
                #get the python variable associated with the vars listed in the config file
                if(profiles_mean['vars'][k] in locals()):
                    data = np.array(locals()[profiles_mean['vars'][k]])
                    if (profiles_mean['vars'][k] in inst_time_group):
                        #print("{} in time group inst".format(profiles_mean['vars'][k]))
                        data_time_slice = data[:,time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:]
                    elif (profiles_mean['vars'][k] in diag_time_group):
                        #print("{} in time group diag".format(profiles_mean['vars'][k]))
                        data_time_slice = data[:,time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:,:]
                    elif (profiles_mean['vars'][k] in swrad_time_group):
                        #print("{} in time group swrad".format(profiles_mean['vars'][k]))
                        data_time_slice = data[:,time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:,:]
                    elif (profiles_mean['vars'][k] in lwrad_time_group):
                        #print("{} in time group lwrad".format(profiles_mean['vars'][k]))
                        data_time_slice = data[:,time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:,:]
                    elif (profiles_mean['vars'][k] in rad_time_group):
                        #print("{} in time group rad".format(profiles_mean['vars'][k]))
                        data_time_slice = data[:,time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:,:]
                    else:
                        print("{} not found in any time groups".format(profiles_mean['vars'][k]))
                        exit()
                    label = profiles_mean['vars_labels'][k]
                    if profiles_mean['conversion_factor']:
                        conversion_factor = profiles_mean['conversion_factor'][k]
                    else:
                        conversion_factor = 1.0
                    
                    #mean profile is obtained by averaging over dimensions 0 and 2
                    mean_data = []
                    for i in range(len(scm_datasets)):
                        mean_data.append(np.mean(data_time_slice[i,:,:,:], (0,2)))

                    if(obs_compare and profiles_mean['vars'][k] in obs_dict):
                        obs_data = np.array(obs_dict[profiles_mean['vars'][k]])
                        obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1],:]

                        obs_mean_data = np.mean(obs_data_time_slice, (0))

                        spr.plot_profile_multi(vert_axis, mean_data, scm_datasets_labels, label, vert_axis_label_pm, comp_dir + '/profiles_mean_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, obs_z=obs_vert_axis, obs_values=obs_mean_data, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)

                        if(bias_val):
                            bias_data = []
                            for i in range(len(scm_datasets)):
                                interp_values = np.flipud(np.interp(np.flipud(obs_vert_axis), np.flipud(vert_axis), np.flipud(mean_data[i])))
                                bias_data.append(interp_values - obs_mean_data)
                            spr.plot_profile_multi(obs_vert_axis, bias_data, scm_datasets_labels, label + ' bias', vert_axis_label_pm, comp_dir + '/profiles_bias_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, line_type='color', zero_line=True, conversion_factor=conversion_factor)

                    else:
                        spr.plot_profile_multi(vert_axis, mean_data, scm_datasets_labels, label, vert_axis_label_pm, comp_dir + '/profiles_mean_' + profiles_mean['vars'][k] + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)
                else:
                    print('The variable ' + profiles_mean['vars'][k] + ' found in ' + args.config[0] + ' in the profiles_mean section is invalid.')
                num_plots_completed += 1
                print_progress(num_plots_completed, num_total_plots)

            for multiplot in profiles_mean_multi:
                #multiplot_vars = profiles_mean_multi[multiplot]['vars']

                #check if all the vars exist
                all_vars_exist = True
                for l in range(len(profiles_mean_multi[multiplot]['vars'])):
                    if (profiles_mean_multi[multiplot]['vars'][l] not in locals()):
                        all_vars_exist = False
                        print('The variable ' + profiles_mean_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' in the ' + multiplot + ' section of profiles_mean_multi is invalid.')

                if all_vars_exist:
                    if profiles_mean_multi[multiplot]['conversion_factor']:
                        conversion_factor = profiles_mean_multi[multiplot]['conversion_factor']
                    else:
                        conversion_factor = 1.0

                    data_list_of_list = []
                    for l in range(len(profiles_mean_multi[multiplot]['vars'])):
                        data = np.array(locals()[profiles_mean_multi[multiplot]['vars'][l]])
                        data_list = []
                        for i in range(len(scm_datasets)):
                            if (profiles_mean_multi[multiplot]['vars'][l] in inst_time_group):
                                #print("{} in time group inst".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in diag_time_group):
                                #print("{} in time group diag".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in swrad_time_group):
                                #print("{} in time group swrad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in lwrad_time_group):
                                #print("{} in time group lwrad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:,:]
                            elif (profiles_mean_multi[multiplot]['vars'][l] in rad_time_group):
                                #print("{} in time group rad".format(profiles_mean_multi[multiplot]['vars'][l]))
                                data_time_slice = data[i,time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:,:]
                            else:
                                print("{} not found in any time groups".format(profiles_mean_multi[multiplot]['vars'][l]))
                                exit()
                            data_list.append(np.mean(data_time_slice, (0,2)))
                        data_list_of_list.append(data_list)
                    spr.plot_profile_multi(vert_axis, data_list_of_list, [profiles_mean_multi[multiplot]['vars_labels'],scm_datasets_labels], profiles_mean_multi[multiplot]['x_label'], vert_axis_label_pm, comp_dir + '/profiles_mean_multi_' + multiplot + plot_ext, y_inverted=y_inverted_val_pm, y_log=y_log_val_pm, y_lim=y_lim_val, line_type='style', conversion_factor=conversion_factor)

                num_plots_completed += 1
                print_progress(num_plots_completed, num_total_plots)
        else:
            print('The variable ' + profiles_mean['vert_axis'] + ' found in ' + args.config[0] + ' in the profiles_mean section is invalid.')

        ### Time-Series ###
        for k in range(len(time_series['vars'])):
            if(time_series['vars'][k] in locals()):
                if(time_series['levels']):                        
                    if (time_series['levels'][k] == -999):
                        data = np.array(locals()[time_series['vars'][k]])
                        plot_name = time_series['vars'][k]
                    else:
                        #minus one to account for 0-based indexing (also, 0 passed in is not working)
                        try:
                            data = np.array(locals()[time_series['vars'][k]])[:,:,time_series['levels'][k]-1]
                            plot_name = time_series['vars'][k] + "_L" + str(time_series['levels'][k])
                        except (IndexError):
                            print('The variable ' + time_series['vars'][k] + ' found in ' + args.config[0] + ' was given an invalid vertical level index: ' + str(time_series['levels'][k]))
                            continue
                else:
                    data = np.array(locals()[time_series['vars'][k]])
                    plot_name = time_series['vars'][k]
                
                if (time_series['vars'][k] in inst_time_group):
                    #print("{} in time group inst".format(time_series['vars'][k]))
                    data_time_slice = data[:,time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1],:]
                    time_h_slice = time_h_inst[0][time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                    
                    if time_series_resample:
                        data_delta_seconds = time_inst[0][time_slice_indices_inst[j][1]] - time_inst[0][time_slice_indices_inst[j][1]-1]
                        #create date range for the model data
                        data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                        data_time_slice_periods = time_slice_indices_inst[j][1] - time_slice_indices_inst[j][0]
                        data_date_range = pd.date_range(start=date_inst[0][time_slice_indices_inst[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same

                elif (time_series['vars'][k] in diag_time_group):
                    #print("{} in time group diag".format(time_series['vars'][k]))
                    data_time_slice = data[:,time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1],:]
                    time_h_slice = time_h_diag[0][time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                    
                    if time_series_resample:
                        data_delta_seconds = time_diag[0][time_slice_indices_diag[j][1]] - time_diag[0][time_slice_indices_diag[j][1]-1]
                        #create date range for the model data
                        data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                        data_time_slice_periods = time_slice_indices_diag[j][1] - time_slice_indices_diag[j][0]
                        data_date_range = pd.date_range(start=date_diag[0][time_slice_indices_diag[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                elif (time_series['vars'][k] in swrad_time_group):
                    #print("{} in time group swrad".format(time_series['vars'][k]))
                    data_time_slice = data[:,time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1],:]
                    time_h_slice = time_h_swrad[0][time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                    
                    if time_series_resample:
                        data_delta_seconds = time_swrad[0][time_slice_indices_swrad[j][1]] - time_swrad[0][time_slice_indices_swrad[j][1]-1]
                        #create date range for the model data
                        data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                        data_time_slice_periods = time_slice_indices_swrad[j][1] - time_slice_indices_swrad[j][0]
                        data_date_range = pd.date_range(start=date_swrad[0][time_slice_indices_swrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                elif (time_series['vars'][k] in lwrad_time_group):
                    #print("{} in time group lwrad".format(time_series['vars'][k]))
                    data_time_slice = data[:,time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1],:]
                    time_h_slice = time_h_lwrad[0][time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                    
                    if time_series_resample:
                        data_delta_seconds = time_lwrad[0][time_slice_indices_lwrad[j][1]] - time_lwrad[0][time_slice_indices_lwrad[j][1]-1]
                        #create date range for the model data
                        data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                        data_time_slice_periods = time_slice_indices_lwrad[j][1] - time_slice_indices_lwrad[j][0]
                        data_date_range = pd.date_range(start=date_lwrad[0][time_slice_indices_lwrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                elif (time_series['vars'][k] in rad_time_group):
                    #print("{} in time group rad".format(time_series['vars'][k]))
                    data_time_slice = data[:,time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1],:]
                    time_h_slice = time_h_rad[0][time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                    
                    if time_series_resample:
                        data_delta_seconds = time_rad[0][time_slice_indices_rad[j][1]] - time_rad[0][time_slice_indices_rad[j][1]-1]
                        #create date range for the model data
                        data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                        data_time_slice_periods = time_slice_indices_rad[j][1] - time_slice_indices_rad[j][0]
                        data_date_range = pd.date_range(start=date_rad[0][time_slice_indices_rad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset) #assumes dates for all model datasets are the same
                else:
                    print("{} not found in any time groups".format(time_series['vars'][k]))
                    exit()
                label = time_series['vars_labels'][k]
                
                if time_series['conversion_factor']:
                    conversion_factor = time_series['conversion_factor'][k]
                else:
                    conversion_factor = 1.0
                
                if(obs_compare and time_series['vars'][k] in obs_dict):
                    #get the corresponding obs data
                    obs_data = np.array(obs_dict[time_series['vars'][k]])

                    #isolate the subset of the obs data for the current time slice
                    obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]

                    if time_series_resample:
                        #determine whether obs data frequency matches model output frequency
                        
                        obs_delta_seconds = obs_dict['time'][obs_dict['time_slice_indices'][j][1]] - obs_dict['time'][obs_dict['time_slice_indices'][j][1]-1]
                        if(obs_delta_seconds > data_delta_seconds):
                            #need to downsample the model data to the obs data frequency

                            #create date range for the obs data
                            obs_data_dateoffset = pd.DateOffset(seconds=int(obs_delta_seconds))
                            obs_time_slice_periods = obs_dict['time_slice_indices'][j][1] - obs_dict['time_slice_indices'][j][0]
                            obs_date_range = pd.date_range(start=obs_dict['date'][obs_dict['time_slice_indices'][j][0]], periods=obs_time_slice_periods, freq=obs_data_dateoffset)

                            resample_string = str(int(obs_delta_seconds)) + 'S'

                            #build list of resampled model datasets for this time slice
                            data_time_slice_series_rs = []
                            for i in range(len(scm_datasets)):
                                data_time_slice_series = pd.Series(data_time_slice[i,:,0], index = data_date_range)
                                data_time_slice_series_rs.append(data_time_slice_series.resample(resample_string).mean())

                            spr.plot_time_series_multi(obs_date_range, data_time_slice_series_rs, scm_datasets_labels, 'date', label, comp_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_date_range, obs_values = obs_data_time_slice, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)
                        elif(obs_delta_seconds < data_delta_seconds):
                            print('The case where observations are more frequent than model output has not been implmented yet... ')
                        else:
                            obs_time_time_slice = obs_dict['time_h'][obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                            spr.plot_time_series_multi(time_h_slice, data_time_slice, scm_datasets_labels, 'time (h)', label, comp_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)
                    else:
                        obs_time_time_slice = obs_dict['time_h'][obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                        spr.plot_time_series_multi(time_h_slice, data_time_slice, scm_datasets_labels, 'time (h)', label, comp_dir + '/time_series_' + plot_name + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)
                else:
                    spr.plot_time_series_multi(time_h_slice, data_time_slice, scm_datasets_labels, 'time (h)', label, comp_dir + '/time_series_' + plot_name + plot_ext, line_type='color',skill_scores=skill_scores_val, conversion_factor=conversion_factor)
            else:
                print('The variable ' + time_series['vars'][k] + ' found in ' + args.config[0] + ' in the time_series section is invalid.')
            num_plots_completed += 1
            print_progress(num_plots_completed, num_total_plots)

        for multiplot in time_series_multi:
            #check if all the vars exist

            all_vars_exist = True
            for l in range(len(time_series_multi[multiplot]['vars'])):
                if (time_series_multi[multiplot]['vars'][l] not in locals()):
                    all_vars_exist = False
                    print('The variable ' + time_series_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' in the ' + multiplot + ' section of time_series_multi is invalid.')

            if all_vars_exist:
                if time_series_multi[multiplot]['conversion_factor']:
                    conversion_factor = time_series_multi[multiplot]['conversion_factor']
                else:
                    conversion_factor = 1.0
                data_list_of_list = []
                for l in range(len(time_series_multi[multiplot]['vars'])):
                    if(time_series_multi[multiplot]['levels']):                        
                        if (time_series_multi[multiplot]['levels'][l] == -999):
                            data = np.array(locals()[time_series_multi[multiplot]['vars'][l]])
                            #plot_name = time_series['vars'][k]
                        else:
                            #minus one to account for 0-based indexing (also, 0 passed in is not working)
                            try:
                                data = np.array(locals()[time_series_multi[multiplot]['vars'][l]])[:,:,time_series_multi[multiplot]['levels'][l]-1]
                                #plot_name = time_series['vars'][k] + "_L" + str(time_series['levels'][k])
                            except (IndexError):
                                print('The variable ' + time_series_multi[multiplot]['vars'][l] + ' found in ' + args.config[0] + ' was given an invalid vertical level index: ' + str(time_series_multi[multiplot]['levels'][l]))
                                continue
                    else:
                        data = np.array(locals()[time_series_multi[multiplot]['vars'][l]])
                    data_list = []
                    for i in range(len(scm_datasets)):
                        if (time_series_multi[multiplot]['vars'][l] in inst_time_group):
                            data_time_slice = data[i,time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                            time_h_slice = time_h_inst[0][time_slice_indices_inst[j][0]:time_slice_indices_inst[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_inst[i][time_slice_indices_inst[j][1]] - time_inst[i][time_slice_indices_inst[j][1]-1]
                                
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_inst[j][1] - time_slice_indices_inst[j][0]
                                data_date_range = pd.date_range(start=date_inst[0][time_slice_indices_inst[j][0]], periods=data_time_slice_periods, freq=data_dateoffset)
                                
                        elif(time_series_multi[multiplot]['vars'][l] in diag_time_group):
                            data_time_slice = data[i,time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                            time_h_slice = time_h_diag[0][time_slice_indices_diag[j][0]:time_slice_indices_diag[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_diag[i][time_slice_indices_diag[j][1]] - time_diag[i][time_slice_indices_diag[j][1]-1]
                                
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_diag[j][1] - time_slice_indices_diag[j][0]
                                data_date_range = pd.date_range(start=date_diag[0][time_slice_indices_diag[j][0]], periods=data_time_slice_periods, freq=data_dateoffset)
                        elif(time_series_multi[multiplot]['vars'][l] in swrad_time_group):
                            data_time_slice = data[i,time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                            time_h_slice = time_h_swrad[0][time_slice_indices_swrad[j][0]:time_slice_indices_swrad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_swrad[i][time_slice_indices_swrad[j][1]] - time_swrad[i][time_slice_indices_swrad[j][1]-1]
                                
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_swrad[j][1] - time_slice_indices_swrad[j][0]
                                data_date_range = pd.date_range(start=date_swrad[0][time_slice_indices_swrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset)
                        elif(time_series_multi[multiplot]['vars'][l] in lwrad_time_group):
                            data_time_slice = data[i,time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                            time_h_slice = time_h_lwrad[0][time_slice_indices_lwrad[j][0]:time_slice_indices_lwrad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_lwrad[i][time_slice_indices_lwrad[j][1]] - time_lwrad[i][time_slice_indices_lwrad[j][1]-1]
                                
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_lwrad[j][1] - time_slice_indices_lwrad[j][0]
                                data_date_range = pd.date_range(start=date_lwrad[0][time_slice_indices_lwrad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset)
                        elif(time_series_multi[multiplot]['vars'][l] in rad_time_group):
                            data_time_slice = data[i,time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                            time_h_slice = time_h_rad[0][time_slice_indices_rad[j][0]:time_slice_indices_rad[j][1]]
                            
                            if time_series_resample:
                                data_delta_seconds = time_rad[i][time_slice_indices_rad[j][1]] - time_rad[i][time_slice_indices_rad[j][1]-1]
                                
                                #create date range for the model data
                                data_dateoffset = pd.DateOffset(seconds=int(data_delta_seconds))
                                data_time_slice_periods = time_slice_indices_rad[j][1] - time_slice_indices_rad[j][0]
                                data_date_range = pd.date_range(start=date_rad[0][time_slice_indices_rad[j][0]], periods=data_time_slice_periods, freq=data_dateoffset)
                        data_list.append(data_time_slice)
                    data_list_of_list.append(data_list)

                if(obs_compare and time_series_multi[multiplot]['obs_var'] in obs_dict):
                    #get the corresponding obs data
                    obs_data = np.array(obs_dict[time_series_multi[multiplot]['obs_var']])

                    #isolate the subset of the obs data for the current time slice
                    obs_data_time_slice = obs_data[obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                    if time_series_resample:
                        #determine whether obs data frequency matches model output frequency
                        
                        obs_delta_seconds = obs_dict['time'][obs_dict['time_slice_indices'][j][1]] - obs_dict['time'][obs_dict['time_slice_indices'][j][1]-1]
                        if(obs_delta_seconds > data_delta_seconds):
                            #need to downsample the model data to the obs data frequency

                            #create date range for the obs data
                            obs_data_dateoffset = pd.DateOffset(seconds=int(obs_delta_seconds))
                            obs_time_slice_periods = obs_dict['time_slice_indices'][j][1] - obs_dict['time_slice_indices'][j][0]
                            obs_date_range = pd.date_range(start=obs_dict['date'][obs_dict['time_slice_indices'][j][0]], periods=obs_time_slice_periods, freq=obs_data_dateoffset)

                            resample_string = str(int(obs_delta_seconds)) + 'S'
                            data_list_of_list_rs = []
                            for l in range(len(data_list_of_list)):
                                data_list_rs = []
                                for i in range(len(scm_datasets)):
                                    data_time_slice_series = pd.Series(data_list_of_list[l][i][:,0], index = data_date_range)
                                    data_time_slice_series_rs = data_time_slice_series.resample(resample_string).mean()
                                    data_list_rs.append(data_time_slice_series_rs)
                                data_list_of_list_rs.append(data_list_rs)

                            spr.plot_time_series_multi(obs_date_range, data_list_of_list_rs, [time_series_multi[multiplot]['vars_labels'],scm_datasets_labels], 'date', time_series_multi[multiplot]['y_label'], comp_dir + '/time_series_multi_' + multiplot + plot_ext, conversion_factor=conversion_factor)#, obs_time = obs_date_range, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'])
                        elif(obs_delta_seconds < data_delta_seconds):
                            print('The case where observations are more frequent than model output has not been implmented yet... ')
                        else:
                            obs_time_time_slice = obs_time_h[obs_time_slice_indices[j][0]:obs_time_slice_indices[j][1]]
                            spr.plot_time_series_multi(time_h_slice, data_list_of_list, [time_series_multi[multiplot]['vars_labels'],scm_datasets_labels], 'time (h)', time_series_multi[multiplot]['y_label'], comp_dir + '/time_series_multi_' + multiplot + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'], conversion_factor=conversion_factor)
                    else:
                        obs_time_time_slice = obs_dict['time_h'][obs_dict['time_slice_indices'][j][0]:obs_dict['time_slice_indices'][j][1]]
                        spr.plot_time_series_multi(time_h_slice, data_list_of_list, [time_series_multi[multiplot]['vars_labels'],scm_datasets_labels], 'time (h)', time_series_multi[multiplot]['y_label'], comp_dir + '/time_series_multi_' + multiplot + plot_ext, obs_time = obs_time_time_slice, obs_values = obs_data_time_slice, obs_label = time_series_multi[multiplot]['obs_var_label'], conversion_factor=conversion_factor)
                else:
                    spr.plot_time_series_multi(time_h_slice, data_list_of_list, [time_series_multi[multiplot]['vars_labels'],scm_datasets_labels], 'time (h)', time_series_multi[multiplot]['y_label'], comp_dir + '/time_series_multi_' + multiplot + plot_ext, conversion_factor=conversion_factor)

            num_plots_completed += 1
            print_progress(num_plots_completed, num_total_plots)
