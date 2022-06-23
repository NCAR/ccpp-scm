#!/usr/bin/env python

from netCDF4 import Dataset
import datetime
import numpy as np
import math
import forcing_file_common as ffc

def read_twpice_obs(obs_file, time_slices, date):
  obs_time_slice_indices = []

  obs_fid = Dataset(obs_file, 'r')

  obs_year = obs_fid.variables['year'][:]
  obs_month = obs_fid.variables['month'][:]
  obs_day = obs_fid.variables['day'][:]
  obs_hour = obs_fid.variables['hour'][:]
  obs_time = obs_fid.variables['time_offset'][:]

  obs_date = []
  for i in range(obs_hour.size):
      obs_date.append(datetime.datetime(obs_year[i], obs_month[i], obs_day[i], obs_hour[i], 0, 0, 0))
  obs_date = np.array(obs_date)

  for time_slice in time_slices:
      start_date = datetime.datetime(time_slices[time_slice]['start'][0], time_slices[time_slice]['start'][1],time_slices[time_slice]['start'][2], time_slices[time_slice]['start'][3], time_slices[time_slice]['start'][4])
      end_date = datetime.datetime(time_slices[time_slice]['end'][0], time_slices[time_slice]['end'][1],time_slices[time_slice]['end'][2], time_slices[time_slice]['end'][3], time_slices[time_slice]['end'][4])
      start_date_index = np.where(obs_date == start_date)[0][0]
      end_date_index = np.where(obs_date == end_date)[0][0]
      obs_time_slice_indices.append([start_date_index, end_date_index])

  #find the index corresponding to the start of the simulations
  obs_start_index = np.where(obs_date == date[0][0])[0]
  obs_time = obs_time - obs_time[obs_start_index]

  obs_pres_l = obs_fid.variables['lev'][:]*100.0 #pressure levels in mb

  obs_cld = obs_fid.variables['cld'][:]/100.0
  obs_T = obs_fid.variables['T'][:]
  obs_q = obs_fid.variables['q'][:]/1000.0
  obs_u = obs_fid.variables['u'][:]
  obs_v = obs_fid.variables['v'][:]
  obs_precip = obs_fid.variables['prec_srf'][:]/3.6E7 #convert from mm/hr to m/s
  obs_shf = obs_fid.variables['SH'][:]
  obs_lhf = obs_fid.variables['LH'][:]
  obs_pwat = obs_fid.variables['PW'][:]*10.0 #convert from cm to kg/m2
  obs_lw_net_toa = obs_fid.variables['lw_net_toa'][:]
  obs_rad_net_srf = obs_fid.variables['rad_net_srf'][:]
  obs_sw_dn_toa = obs_fid.variables['sw_dn_toa'][:]
  obs_sw_dn_srf = obs_fid.variables['sw_dn_srf'][:]
  obs_lw_dn_srf = obs_fid.variables['lw_dn_srf'][:]
  obs_lwp = obs_fid.variables['LWP'][:]*10.0 #convert from cm to kg/m2
  #obs_T_forcing = obs_fid.variables['dTdt'][:]*24.0 #convert from K/hour to K/day
  #obs_q_forcing = obs_fid.variables['dqdt'][:]*24.0 #convert from g/kg/hour to g/kg/day
  obs_h_advec_T = obs_fid.variables['T_adv_h'][:]*24.0
  obs_h_advec_q = obs_fid.variables['q_adv_h'][:]*24.0
  obs_v_advec_T = obs_fid.variables['T_adv_v'][:]*24.0
  obs_v_advec_q = obs_fid.variables['q_adv_v'][:]*24.0

  obs_T_forcing = obs_h_advec_T + obs_v_advec_T
  obs_q_forcing = obs_h_advec_q + obs_v_advec_q

  obs_time_h = obs_time/3600.0

  Rd = 287.0
  Rv = 461.0

  e_s = 6.1078*np.exp(17.2693882*(obs_T - 273.16)/(obs_T - 35.86))*100.0 #Tetens formula produces e_s in mb (convert to Pa)
  e = obs_q*obs_pres_l/(obs_q + (Rd/Rv)*(1.0 - obs_q)) #compute vapor pressure from specific humidity
  obs_rh = np.clip(e/e_s, 0.0, 1.0)

  obs_rh_500 = np.zeros(obs_rh.shape[0])
  index_500 = np.where(obs_pres_l[:]*0.01 < 500.0)[0][0]
  lifrac = (obs_pres_l[index_500-1] - 50000.0)/(obs_pres_l[index_500-1] - obs_pres_l[index_500])
  for j in range(obs_rh.shape[0]): #loop over times
    obs_rh_500[j] = obs_rh[j,index_500-1] + lifrac*(obs_rh[j,index_500] - obs_rh[j,index_500-1])
          #print index_500, pres_l[-1][j,index_500,k], pres_l[-1][j,index_500-1,k], rh_500_kj, rh[-1][j,index_500,k], rh[-1][j,index_500-1,k]

  return_dict = {'year': obs_year, 'month': obs_month, 'day': obs_day, 'hour': obs_hour,
    'time': obs_time, 'date': obs_date, 'time_slice_indices': obs_time_slice_indices,
    'pres_l': obs_pres_l, 'cld': obs_cld, 'T': obs_T, 'q': obs_q, 'u': obs_u, 'v': obs_v,
    'pwat': obs_pwat, 'time_h': obs_time_h,
    'tprcp_rate_accum': obs_precip, 'qv': obs_q, 'rh': obs_rh, 'rh_500': obs_rh_500,
    'lw_up_TOA_tot': obs_lw_net_toa, 'rad_net_srf': obs_rad_net_srf, 'sw_dn_TOA_tot': obs_sw_dn_toa,
    'lw_dn_sfc_tot': obs_lw_dn_srf, 'sw_dn_sfc_tot': obs_sw_dn_srf, 'lwp': obs_lwp,
    'T_force_tend': obs_T_forcing, 'qv_force_tend': obs_q_forcing}

  obs_fid.close()

  return return_dict

def read_arm_sgp_summer_1997_obs(obs_file, time_slices, date):
  obs_time_slice_indices = []

  obs_fid = Dataset(obs_file, 'r')

  obs_year = obs_fid.variables['Year'][:]
  obs_month = obs_fid.variables['Month'][:]
  obs_day = obs_fid.variables['Day'][:]
  #obs_hour = obs_fid.variables['hour'][:]
  obs_time = obs_fid.variables['time_offset'][:]

  #this file doesn't have the hour variable - calculate from the time offset (seconds from 00Z on 6/18/1997)
  obs_hour = (((obs_time - 3)/3600.0)%24).astype(int)

  obs_date = []
  for i in range(obs_hour.size):
      obs_date.append(datetime.datetime(obs_year[i], obs_month[i], obs_day[i], obs_hour[i], 0, 0, 0))
  obs_date = np.array(obs_date)

  for time_slice in time_slices:
    start_date = datetime.datetime(time_slices[time_slice]['start'][0], time_slices[time_slice]['start'][1],time_slices[time_slice]['start'][2], time_slices[time_slice]['start'][3], time_slices[time_slice]['start'][4])
    end_date = datetime.datetime(time_slices[time_slice]['end'][0], time_slices[time_slice]['end'][1],time_slices[time_slice]['end'][2], time_slices[time_slice]['end'][3], time_slices[time_slice]['end'][4])
    start_date_index = np.where(obs_date == start_date)[0][0]
    end_date_index = np.where(obs_date == end_date)[0][0]
    obs_time_slice_indices.append([start_date_index, end_date_index])
    #print start_date, end_date, start_date_index, end_date_index, obs_date[start_date_index], obs_date[end_date_index]

  #find the index corresponding to the start of the simulations
  obs_start_index = np.where(obs_date == date[0][0])[0]
  obs_time = obs_time - obs_time[obs_start_index]


  obs_pres_l = np.flipud(obs_fid.variables['lev'][:])*100.0 #pressure levels in mb

  obs_cld = np.fliplr(obs_fid.variables['ARSCL_Cld'][:,:,0,0])/100.0
  obs_T = np.fliplr(obs_fid.variables['Temp'][:,:,0,0])
  obs_q = np.fliplr(obs_fid.variables['H2O_Mixing_Ratio'][:,:,0,0]/1000.0)
  obs_u = np.fliplr(obs_fid.variables['u_wind'][:,:,0,0])
  obs_v = np.fliplr(obs_fid.variables['v_wind'][:,:,0,0])
  obs_precip = obs_fid.variables['Prec'][:,0,0]
  # obs_shf = obs_fid.variables['SH'][:]
  # obs_lhf = obs_fid.variables['LH'][:]
  # obs_pwat = obs_fid.variables['PW'][:]
  # obs_lw_net_toa = obs_fid.variables['lw_net_toa'][:]
  # obs_rad_net_srf = obs_fid.variables['rad_net_srf'][:]
  # obs_sw_dn_toa = obs_fid.variables['sw_dn_toa'][:]
  # obs_sw_dn_srf = obs_fid.variables['sw_dn_srf'][:]
  # obs_lw_dn_srf = obs_fid.variables['lw_dn_srf'][:]
  # obs_lwp = obs_fid.variables['LWP'][:]*10.0 #convert from cm to kg/m2
  # #obs_T_forcing = obs_fid.variables['dTdt'][:]*24.0 #convert from K/hour to K/day
  # #obs_q_forcing = obs_fid.variables['dqdt'][:]*24.0 #convert from g/kg/hour to g/kg/day
  # obs_h_advec_T = obs_fid.variables['T_adv_h'][:]*24.0
  # obs_h_advec_q = obs_fid.variables['q_adv_h'][:]*24.0
  # obs_v_advec_T = obs_fid.variables['T_adv_v'][:]*24.0
  # obs_v_advec_q = obs_fid.variables['q_adv_v'][:]*24.0
  #
  # obs_T_forcing = obs_h_advec_T + obs_v_advec_T
  # obs_q_forcing = obs_h_advec_q + obs_v_advec_q
  #
  # obs_time_h = obs_time/3600.0
  #
  # Rd = 287.0
  # Rv = 461.0
  #
  # e_s = 6.1078*np.exp(17.2693882*(obs_T - 273.16)/(obs_T - 35.86))*100.0 #Tetens formula produces e_s in mb (convert to Pa)
  # e = obs_q*obs_pres_l/(obs_q + (Rd/Rv)*(1.0 - obs_q)) #compute vapor pressure from specific humidity
  # obs_rh = np.clip(e/e_s, 0.0, 1.0)
  #
  # obs_rh_500 = np.zeros(obs_rh.shape[0])
  # index_500 = np.where(obs_pres_l[:]*0.01 < 500.0)[0][0]
  # lifrac = (obs_pres_l[index_500-1] - 50000.0)/(obs_pres_l[index_500-1] - obs_pres_l[index_500])
  # for j in range(obs_rh.shape[0]): #loop over times
  #   obs_rh_500[j] = obs_rh[j,index_500-1] + lifrac*(obs_rh[j,index_500] - obs_rh[j,index_500-1])
  #         #print index_500, pres_l[-1][j,index_500,k], pres_l[-1][j,index_500-1,k], rh_500_kj, rh[-1][j,index_500,k], rh[-1][j,index_500-1,k]

  return_dict = {'year': obs_year, 'month': obs_month, 'day': obs_day, 'hour': obs_hour,
    'time': obs_time, 'date': obs_date, 'time_slice_indices': obs_time_slice_indices,
    'pres_l': obs_pres_l, 'cld': obs_cld, 'T': obs_T, 'qv': obs_q, 'u': obs_u, 'v': obs_v,
    'precip': obs_precip}#, 'shf': obs_shf, 'lhf': obs_lhf, 'pwat': obs_pwat, 'time_h': obs_time_h,
    # 'rain': obs_precip, 'rainc': obs_precip, 'qv': obs_q, 'rh': obs_rh, 'rh_500': obs_rh_500,
    # 'lw_up_TOA_tot': obs_lw_net_toa, 'rad_net_srf': obs_rad_net_srf, 'sw_dn_TOA_tot': obs_sw_dn_toa,
    # 'lw_dn_sfc_tot': obs_lw_dn_srf, 'sw_dn_sfc_tot': obs_sw_dn_srf, 'lwp': obs_lwp,
    # 'T_force_tend': obs_T_forcing, 'qv_force_tend': obs_q_forcing}

  # return_dict = {'year': obs_year, 'month': obs_month, 'day': obs_day, 'hour': obs_hour,
  #   'time': obs_time, 'date': obs_date, 'time_slice_indices': obs_time_slice_indices,
  #   'pres_l': obs_pres_l, 'cld': obs_cld, 'T': obs_T, 'q': obs_q, 'u': obs_u, 'v': obs_v,
  #   'precip': obs_precip, 'shf': obs_shf, 'lhf': obs_lhf, 'pwat': obs_pwat, 'time_h': obs_time_h,
  #   'rain': obs_precip, 'rainc': obs_precip, 'qv': obs_q, 'rh': obs_rh, 'rh_500': obs_rh_500,
  #   'lw_up_TOA_tot': obs_lw_net_toa, 'rad_net_srf': obs_rad_net_srf, 'sw_dn_TOA_tot': obs_sw_dn_toa,
  #   'lw_dn_sfc_tot': obs_lw_dn_srf, 'sw_dn_sfc_tot': obs_sw_dn_srf, 'lwp': obs_lwp,
  #   'T_force_tend': obs_T_forcing, 'qv_force_tend': obs_q_forcing}

  obs_fid.close()

  return return_dict

def read_LASSO_obs(obs_file, time_slices, date):
    obs_time_slice_indices = []

    obs_fid = Dataset(obs_file, 'r')
    obs_fid.set_auto_mask(False)
        
    obs_time = obs_fid.variables['time_offset'][:]
    obs_datetime_string = obs_fid.getncattr('output_start_datetime')

    #get initial date from global file attribute
    obs_init_datetime = datetime.datetime.strptime(obs_datetime_string, '%Y%m%d.%H%M%S %Z')

    obs_date = []
    for i in range(obs_time.size):
        obs_date.append(obs_init_datetime + datetime.timedelta(seconds = obs_time[i]))
    obs_date = np.array(obs_date)

    for time_slice in time_slices:
      start_date = datetime.datetime(time_slices[time_slice]['start'][0], time_slices[time_slice]['start'][1],time_slices[time_slice]['start'][2], time_slices[time_slice]['start'][3], time_slices[time_slice]['start'][4])
      end_date = datetime.datetime(time_slices[time_slice]['end'][0], time_slices[time_slice]['end'][1],time_slices[time_slice]['end'][2], time_slices[time_slice]['end'][3], time_slices[time_slice]['end'][4])
      start_date_index = np.where(obs_date == start_date)[0][0]
      end_date_index = np.where(obs_date == end_date)[0][0]
      obs_time_slice_indices.append([start_date_index, end_date_index])
      #print start_date, end_date, start_date_index, end_date_index, obs_date[start_date_index], obs_date[end_date_index]

    #find the index corresponding to the start of the simulations
    obs_start_index = np.where(obs_date == date[0][0])[0]
    obs_time = obs_time - obs_time[obs_start_index]

    #pressure stored in kPa (pressure is 0 for initial time)
    obs_pres = obs_fid.variables['bar_pres'][1:,:]*1.0E3

    #get average pressure levels
    obs_pres_l = np.mean(obs_pres[:,:], (0))

    obs_theta = obs_fid.variables['potential_temp'][1:,:]
    #print obs_theta[0,:]
    obs_T = (obs_pres/ffc.p0)**(ffc.R_dry/ffc.c_p)*obs_theta

    obs_q = obs_fid.variables['water_vapor_mixing_ratio'][1:,:]*1.0E-3
    #print obs_q[0,:]
    obs_cld = obs_fid.variables['cloud_fraction'][1:,:]

    #print obs_cld[0,:]

    return_dict = {'time': obs_time, 'date': obs_date, 'time_slice_indices': obs_time_slice_indices, 'pres_l': obs_pres_l,
        'T': obs_T, 'qv': obs_q, 'cld': obs_cld}

    obs_fid.close()

    return return_dict

def read_gabls3_obs(obs_file, time_slices, date):
    
    con_g = 9.81
    con_rd = 287.0
    con_cp = 1004.0
    con_vir = 0.61
    p0 = 100000.0
    
    obs_time_slice_indices = []

    obs_fid = Dataset(obs_file, 'r')
    obs_fid.set_auto_mask(False)
    
    obs_time_hours_elapsed = obs_fid.variables['time'][:]
    obs_date_ints = obs_fid.variables['date'][:]
    n_times = len(obs_time_hours_elapsed)
    
    obs_time_hours = np.mod(obs_time_hours_elapsed,24.0*np.ones(n_times))
        
    obs_date = []
    for i in range(obs_time_hours.size):
        obs_year = int(str(obs_date_ints[i])[0:4])
        obs_month = int(str(obs_date_ints[i])[4:6])
        obs_day = int(str(obs_date_ints[i])[6:])
        obs_hour = int(math.floor(obs_time_hours[i]))
        obs_minutes = int((obs_time_hours[i] - obs_hour)*60.0)
        obs_seconds = int(((obs_time_hours[i] - obs_hour) - obs_minutes/60.0)*3600.0)
        #disregard milliseconds
        obs_date.append(datetime.datetime(obs_year, obs_month, obs_day, obs_hour, obs_minutes, obs_seconds, 0))
    obs_date = np.array(obs_date)
    
    for time_slice in time_slices:
      start_date = datetime.datetime(time_slices[time_slice]['start'][0], time_slices[time_slice]['start'][1],time_slices[time_slice]['start'][2], time_slices[time_slice]['start'][3], time_slices[time_slice]['start'][4])
      end_date = datetime.datetime(time_slices[time_slice]['end'][0], time_slices[time_slice]['end'][1],time_slices[time_slice]['end'][2], time_slices[time_slice]['end'][3], time_slices[time_slice]['end'][4])
      start_date_index = np.where(obs_date == start_date)[0][0]
      try:
          end_date_index = np.where(obs_date == end_date)[0][0]
      except IndexError:
          end_date_index = len(obs_date) - 1
      obs_time_slice_indices.append([start_date_index, end_date_index])
    
    obs_start_index = np.where(obs_date == date[0][0])[0]
    obs_time = 3600.0*(obs_time_hours_elapsed - obs_time_hours_elapsed[obs_start_index])
    
    obs_zt = np.fliplr(obs_fid.variables['zt'][:])
    obs_zf = np.fliplr(obs_fid.variables['zf'][:])
    obs_t = np.fliplr(obs_fid.variables['t'][:])
    obs_t_fill_value = obs_fid.variables['t']._FillValue
    obs_q = np.fliplr(obs_fid.variables['q'][:])
    obs_q_fill_value = obs_fid.variables['q']._FillValue
    obs_th = np.fliplr(obs_fid.variables['th'][:])
    
    obs_hpbl = obs_fid.variables['hpbl'][:]
    obs_tsk = obs_fid.variables['tsk'][:]
    obs_shf = obs_fid.variables['shf'][:]
    obs_lhf = obs_fid.variables['lhf'][:]
    obs_lw_up = obs_fid.variables['lup'][:]
    obs_lw_dn = obs_fid.variables['ldw'][:]
    obs_sw_up = obs_fid.variables['qup'][:]
    obs_sw_dn = obs_fid.variables['qdw'][:]
    obs_gflux = obs_fid.variables['g'][:]
    obs_t2m = obs_fid.variables['t2m'][:]
    obs_q2m = obs_fid.variables['q2m'][:]
    obs_ustar = obs_fid.variables['ustar'][:]
    obs_u10m = obs_fid.variables['u10m'][:]
    obs_v10m = obs_fid.variables['v10m'][:]
    
    obs_fid.close()
    
    p_surf = 102440.0 #case specifications
    
    #find where T has valid values
    good_t_indices = np.where(obs_t[0,:] != obs_t_fill_value)[0]
    
    #assume missing values always appear at same vertical indices
    good_t = obs_t[:,good_t_indices]
    good_zf = obs_zf[:,good_t_indices]
    
    #good_q_indices = np.where(obs_q[0,:] != obs_q_fill_value)[0]
    good_q = obs_q[:,good_t_indices]
    good_th = obs_th[:,good_t_indices]
        
    #p_good_lev = np.zeros(len(good_t_indices))
    #calculate p from hydrostatic equation, surface pressure, and T and q
    #p_good_lev[0] = p_surf*np.exp(-con_g/(con_rd*good_t[0]*(1.0 + con_vir*good_q[0]))*good_zf[0])
    #for k in range(1,len(good_t_indices)):
    #    p_good_lev[k] = p_good_lev[k-1]*np.exp((-con_g/(con_rd*0.5*(good_t[k-1]+good_t[k])*(1.0 + con_vir*0.5*(good_q[k-1]+good_q[k])))*(good_zf[k]-good_zf[k-1])))
    
    #calculate p from temperature and potential temperature
    
    p_good_lev_from_th = p0/(good_th[0,:]/good_t[0,:])**(con_cp/con_rd)
    
    obs_sfc_rad_net = (obs_sw_dn - obs_sw_up) + (obs_lw_dn - obs_lw_up)
    
    return_dict = {'time': obs_time, 'date': obs_date, 'time_slice_indices': obs_time_slice_indices, 'pres_l': p_good_lev_from_th,
        'T': good_t, 'qv': good_q, 'shf': obs_shf, 'lhf': obs_lhf, 'time_h': obs_time/3600.0, 'sfc_up_lw_land': obs_lw_up, 'sfc_dwn_lw': obs_lw_dn,
        'sfc_dwn_sw': obs_sw_dn, 'sfc_up_sw': obs_sw_up, 'sfc_rad_net_land': obs_sfc_rad_net, 'gflux': -1*obs_gflux, 't2m':obs_t2m, 'q2m':obs_q2m, 
        'ustar':obs_ustar,'u10m':obs_u10m, 'v10m':obs_v10m, 'hpbl':obs_hpbl, 'tsfc':obs_tsk}
    
    return return_dict