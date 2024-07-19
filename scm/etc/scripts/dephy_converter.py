#!/usr/bin/env python

import argparse
import logging
import f90nml
import os
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta


###############################################################################
# Global settings                                                             #
###############################################################################

# Path to the directory containing processed case input files
CASE_NML_DIR = '../case_config'

# Path to the directory containing processed case input files
PROCESSED_CASE_DIR = '../../data/processed_case_input'

# For developers: set logging level to DEBUG for additional output
LOGLEVEL = logging.DEBUG
#LOGLEVEL = logging.INFO

DEFAULT_MISSING_VALUE = -9999.0
DEFAULT_NUDGING_TIMESCALE = 7200.0 #s

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-n',   '--case_name',       help='name of case',         required=True)
parser.add_argument('-a',   '--use_area',        help='use column_area namelist attribute as forcing_scale',  action='store_true')


########################################################################################
#
########################################################################################
def parse_arguments():
    """Parse command line arguments"""
    args           = parser.parse_args()
    case_name      = args.case_name
    use_area       = args.use_area

        
    return (case_name, use_area)

########################################################################################
#
########################################################################################
def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)

class Case_Data(object):
    def __init__(self, name, missing_value, time, levels, lat, lon, height, theta_il, qt, ql, qi, u, v, tke, ozone, \
                 p_surf, T_surf, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, \
                 dT_dt_rad, h_advec_thil, v_advec_thil, h_advec_qt, v_advec_qt, sh_flux_sfc, lh_flux_sfc):
        self._name = name
        self._missing_value = missing_value
        self._time = time
        self._levels = levels
        self._lat = lat
        self._lon = lon
        self._height = height
        self._theta_il = theta_il
        self._qt = qt
        self._ql = ql
        self._qi = qi
        self._u = u
        self._v = v
        self._tke = tke
        self._ozone = ozone
        self._p_surf = p_surf
        self._T_surf = T_surf
        self._w_ls = w_ls
        self._omega = omega
        self._u_g = u_g
        self._v_g = v_g
        self._u_nudge = u_nudge
        self._v_nudge = v_nudge
        self._T_nudge = T_nudge
        self._thil_nudge = thil_nudge
        self._qt_nudge = qt_nudge
        self._dT_dt_rad = dT_dt_rad
        self._h_advec_thil = h_advec_thil
        self._v_advec_thil = v_advec_thil
        self._h_advec_qt = h_advec_qt
        self._v_advec_qt = v_advec_qt
        self._sh_flux_sfc = sh_flux_sfc
        self._lh_flux_sfc = lh_flux_sfc
    
    def __repr__(self):
        return "Case_Data \n Globals: \n" \
                "  name:% s \n" \
                "  missing_value:% s \n" \
                " Dimensions: \n" \
                "  time:% s \n" \
                "  levels: % s \n" \
                " Scalars: \n" \
                "  lat: % s \n" \
                "  lon: % s \n" \
                " Initial: \n" \
                "  height: % s \n" \
                "  theta_il: % s \n" \
                "  qt: % s \n" \
                "  ql: % s \n" \
                "  qi: % s \n" \
                "  u: % s \n" \
                "  v: % s \n" \
                "  tke: %s \n" \
                "  ozone: %s \n" \
                " Forcing: \n" \
                "  p_surf: %s \n" \
                "  T_surf: %s \n" \
                "  w_ls (time avg): %s \n" \
                "  omega (time avg): %s \n" \
                "  u_g (time avg): %s \n" \
                "  v_g (time avg): %s \n" \
                "  u_nudge (time avg): %s \n" \
                "  v_nudge (time avg): %s \n" \
                "  T_nudge (time avg): %s \n" \
                "  thil_nudge (time avg): %s \n" \
                "  qt_nudge (time avg): %s \n" \
                "  dT_dt_rad (time avg): %s \n" \
                "  h_advec_thil (time avg): %s \n" \
                "  v_advec_thil (time avg): %s \n" \
                "  h_advec_qt (time avg): %s \n" \
                "  v_advec_qt (time avg): %s \n" \
                "  sh_flux_sfc: %s \n" \
                "  lh_flux_sfc: %s \n" \
                  % (self._name, self._missing_value, self._time, self._levels,
                     self._lat, self._lon, self._height, self._theta_il,
                     self._qt, self._ql, self._qi, self._u, self._v,
                     self._tke, self._ozone, 
                     self._p_surf, self._T_surf, np.mean(self._w_ls, axis=1),
                     np.mean(self._omega, axis=1), np.mean(self._u_g, axis=1),
                     np.mean(self._v_g, axis=1), np.mean(self._u_nudge, axis=1),
                     np.mean(self._v_nudge, axis=1), np.mean(self._T_nudge, axis=1),
                     np.mean(self._thil_nudge, axis=1),np.mean(self._qt_nudge, axis=1),
                     np.mean(self._dT_dt_rad, axis=1),
                     np.mean(self._h_advec_thil, axis=1),np.mean(self._v_advec_thil, axis=1),
                     np.mean(self._h_advec_qt, axis=1),np.mean(self._v_advec_qt, axis=1),
                     self._sh_flux_sfc, self._lh_flux_sfc)

def get_case_nml(case_name):
    """Returns case configuration Fortran namelist"""
    
    filename = os.path.join(CASE_NML_DIR, case_name + '.nml')
    
    print(filename)
    
    error = False
    nml = ''
    if (os.path.exists(filename)):
        nml = f90nml.read(filename)
    else:
        error = True
    
    return (nml, error)

def get_case_data(case_name):
    """Returns proprietery CCPP SCM case data in NetCDF Dataset format"""
    
    #TODO: need to handle LSM ICs
    
    filename = os.path.join(PROCESSED_CASE_DIR, case_name + '.nc')
    
    error = False
    try:
        nc_fid = Dataset(filename , 'r')
    except:
        error = True
        
    if (not error):
        #read global variables
        try:
            missing_value = nc_fid.getncattr('missing_value')
        except:
            missing_value = DEFAULT_MISSING_VALUE
        
        time   = nc_fid.variables['time'][:]
        levels = nc_fid.variables['levels'][:]
        
        #read variables from scalar group
        scalars_grp = nc_fid.groups['scalars']
        lat = scalars_grp.variables['lat'][:]
        lon = scalars_grp.variables['lon'][:]
        
        #read variables from initial group
        initial_grp = nc_fid.groups['initial']
        height = initial_grp.variables['height'][:]
        theta_il = initial_grp.variables['thetail'][:]
        qt = initial_grp.variables['qt'][:]
        ql = initial_grp.variables['ql'][:]
        qi = initial_grp.variables['qi'][:]
        u = initial_grp.variables['u'][:]
        v = initial_grp.variables['v'][:]
        tke = initial_grp.variables['tke'][:]
        ozone = initial_grp.variables['ozone'][:]
        
        #read variables from forcing group
        forcing_grp = nc_fid.groups['forcing']
        p_surf = forcing_grp.variables['p_surf'][:]
        T_surf = forcing_grp.variables['T_surf'][:]
        w_ls = forcing_grp.variables['w_ls'][:]
        omega = forcing_grp.variables['omega'][:]
        u_g = forcing_grp.variables['u_g'][:]
        v_g = forcing_grp.variables['v_g'][:]
        u_nudge = forcing_grp.variables['u_nudge'][:]
        v_nudge = forcing_grp.variables['v_nudge'][:]
        T_nudge = forcing_grp.variables['T_nudge'][:]
        thil_nudge = forcing_grp.variables['thil_nudge'][:]
        qt_nudge = forcing_grp.variables['qt_nudge'][:]
        dT_dt_rad = forcing_grp.variables['dT_dt_rad'][:]
        h_advec_thil = forcing_grp.variables['h_advec_thetail'][:]
        v_advec_thil = forcing_grp.variables['v_advec_thetail'][:]
        h_advec_qt = forcing_grp.variables['h_advec_qt'][:]
        v_advec_qt = forcing_grp.variables['v_advec_qt'][:]
        try:
            sh_flux_sfc = forcing_grp.variables['sh_flux_sfc'][:]
        except KeyError:
            sh_flux_sfc = ''    
        try:
            lh_flux_sfc = forcing_grp.variables['lh_flux_sfc'][:]
        except KeyError:
            lh_flux_sfc = ''
        
        nc_fid.close()
    
    case_data = Case_Data(case_name, missing_value, time, levels, lat, lon, 
                          height, theta_il, qt, ql, qi, u, v, tke, ozone,
                          p_surf, T_surf, w_ls, omega, u_g, v_g,
                          u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge,
                          dT_dt_rad, h_advec_thil, v_advec_thil, h_advec_qt,
                          v_advec_qt, sh_flux_sfc, lh_flux_sfc)
    
    return(case_data, error)

def write_SCM_case_file(case_nml, case_data, use_area):
    """Write all data to a netCDF file in the DEPHY-SCM format"""
    
    #TODO: need to handle LSM ICs
    
    # Working types
    wp = np.float64
    wi = np.int32
    
    # Local switches
    forcing_on  = 1
    forcing_off = 0
    
    nml_filename = os.path.join(CASE_NML_DIR, case_nml['case_config']['case_name'] + '.nml')
    
    # Output file
    com = 'mkdir -p ' + PROCESSED_CASE_DIR
    logging.info(com)
    os.system(com)
    fileOUT = os.path.join(PROCESSED_CASE_DIR, case_nml['case_config']['case_name'] + '_dephy' + '_SCM_driver.nc')
    
    nc_file = Dataset(fileOUT, 'w', format='NETCDF3_CLASSIC')
    nc_file.description = "Case data for {} from CCPP SCM".format(case_nml['case_config']['case_name'])

    nc_file.missing_value   = case_data._missing_value
    
    #not all namelists will have minutes, set to 0 if nml doesn't have
    try:
        minute = case_nml['case_config']['minute']
    except KeyError:
        minute = 0    
    
    start_date = datetime(case_nml['case_config']['year'],case_nml['case_config']['month'],case_nml['case_config']['day'],case_nml['case_config']['hour'],minute,0)
    start_date_string = start_date.strftime("%Y-%m-%d %H:%M:%S")
    runtime = case_nml['case_config']['runtime']
    delta = timedelta(seconds=runtime)
    end_date = start_date + delta
    end_date_string   = end_date.strftime("%Y-%m-%d %H:%M:%S")
    loc_string  = str(case_data._lon) + "E" + str(case_data._lat) + "N"
    case_string = case_nml['case_config']['case_name'] + '_' + start_date_string + '_' + loc_string
    
    logging.debug('Case string: {}'.format(case_string))
    logging.debug('Case start date: {}'.format(start_date))
    logging.debug('Case duration: {}'.format(delta))
    logging.debug('Case end date: {}'.format(end_date))
    
    if (case_nml['case_config']['sfc_type'] == 0):
        surface_string = 'ocean'
    elif (case_nml['case_config']['sfc_type'] == 1):
        surface_string = 'land'
    elif (case_nml['case_config']['sfc_type'] == 2):
        surface_string = 'ice'
    
    
    #DEPHY v1 format specifies the global attributes in this order. Some attributes are rewritten below after the order is established in the file.
    nc_file.case              = case_string
    nc_file.title             = 'Forcing and Initial Conditions for ' + case_string
    nc_file.reference         = 'https://dtcenter.org/sites/default/files/paragraph/scm-ccpp-guide-v6-0-0.pdf'
    nc_file.author            = 'Grant J. Firl and Dustin Swales'
    nc_file.version           = 'Created on ' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    nc_file.format_version    = 'DEPHY SCM format version 1'
    nc_file.modifications     = ''
    nc_file.script            = os.path.basename(__file__)
    nc_file.comment           = 'converted from ' + case_nml['case_config']['case_name'] + '.nc'
    nc_file.start_date        = start_date_string
    nc_file.end_date          = end_date_string
    
    if (use_area and case_nml['case_config']['column_area']):
        nc_file.forcing_scale     = case_nml['case_config']['column_area']
    else:    
        nc_file.forcing_scale     = -1
        
    nc_file.adv_ta            = forcing_off
    nc_file.adv_qv            = forcing_off
    nc_file.adv_ua            = forcing_off  #no mom_forcing_type implemented this (providing pre-calculated advective terms for u)
    nc_file.adv_va            = forcing_off  #no mom_forcing_type implemented this (providing pre-calculated advective terms for v)
    nc_file.adv_theta         = forcing_off
    nc_file.adv_thetal        = forcing_off
    nc_file.adv_qt            = forcing_off
    nc_file.adv_rv            = forcing_off 
    nc_file.adv_rt            = forcing_off
    nc_file.radiation         = "on"         #not implemented in CCPP SCM - controlled by CCPP SDF and/or namelist
    nc_file.forc_wap          = forcing_off
    nc_file.forc_wa           = forcing_off
    nc_file.forc_geo          = forcing_off
    nc_file.nudging_ua        = forcing_off
    nc_file.nudging_va        = forcing_off
    nc_file.nudging_ta        = forcing_off
    nc_file.nudging_theta     = forcing_off
    nc_file.nudging_thetal    = forcing_off
    nc_file.nudging_qv        = forcing_off
    nc_file.nudging_qt        = forcing_off
    nc_file.nudging_rv        = forcing_off
    nc_file.nudging_rt        = forcing_off
    nc_file.zh_nudging_ta     = forcing_off
    nc_file.zh_nudging_theta  = forcing_off
    nc_file.zh_nudging_thetal = forcing_off
    nc_file.zh_nudging_qv     = forcing_off
    nc_file.zh_nudging_qt     = forcing_off
    nc_file.zh_nudging_rv     = forcing_off
    nc_file.zh_nudging_rt     = forcing_off
    nc_file.zh_nudging_ua     = forcing_off
    nc_file.zh_nudging_va     = forcing_off
    nc_file.pa_nudging_ta     = forcing_off
    nc_file.pa_nudging_theta  = forcing_off
    nc_file.pa_nudging_thetal = forcing_off
    nc_file.pa_nudging_qv     = forcing_off
    nc_file.pa_nudging_qt     = forcing_off
    nc_file.pa_nudging_rv     = forcing_off
    nc_file.pa_nudging_rt     = forcing_off
    nc_file.pa_nudging_ua     = forcing_off
    nc_file.pa_nudging_va     = forcing_off
    #
    nc_file.surface_type      = surface_string
    nc_file.surface_forcing_temp     = 'none'
    nc_file.surface_forcing_moisture = 'none'
    nc_file.surface_forcing_wind     = 'none'
    
    #rewrite forc_wa, forc_wap, forc_geo, nudging_ua, nudging_va depending on mom_forcing_type provided in case_config nml    
    if (case_nml['case_config']['mom_forcing_type'] == 2):
        #CCPP SCM proprietery forcing interprets mom_forcing_type = 2 as calculating vertical advective terms from provided vertical velocity AND applying geostrophic winds
        
        #CCPP SCM proprietery cases could have either w or omega available (or both); use omega by default?
        w_ls_avail = True if np.any(case_data._w_ls[:,:]) else False
        omega_avail = True if np.any(case_data._omega[:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using vertical velocity (through mom_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
            #logging.critical(message)
            #raise Exception(message)
        
        geostrophic_avail = True if (np.any(case_data._u_g[:,:]) or np.any(case_data._v_g[:,:])) else False
        if geostrophic_avail:
            nc_file.forc_geo = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using geostrophic winds (through mom_forcing = 2), but neither u_g or v_g have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_geo = forcing_off
            #logging.critical(message)
            #raise Exception(message)
        nc_file.nudging_ua = forcing_off
        nc_file.nudging_va = forcing_off
                        
    elif (case_nml['case_config']['mom_forcing_type'] == 3):
        #CCPP SCM proprietery forcing interprets mom_forcing_type = 3 as calculating momentum forcing as nudging toward u and v profiles (only)
        
        nc_file.forc_wa  = forcing_off
        nc_file.forc_wap = forcing_off
        nc_file.forc_geo = forcing_off
        
        u_nudge_avail = True if np.any(case_data._u_nudge[:,:]) else False
        v_nudge_avail = True if np.any(case_data._v_nudge[:,:]) else False
        relax_time_avail = True if case_nml['case_config']['relax_time'] else False
        if relax_time_avail:
            relax_time = case_nml['case_config']['relax_time']
        else:
            relax_time = DEFAULT_NUDGING_TIMESCALE
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using nudging (through mom_forcing = 3), but relax_time was not provided -- using default value of {1}s '.format(nml_filename, DEFAULT_NUDGING_TIMESCALE)
            logging.info(message)
        
        nc_file.nudging_ua = forcing_on*relax_time
        nc_file.nudging_va = forcing_on*relax_time
                
    if (case_nml['case_config']['thermo_forcing_type'] == 1):
        #total advective forcing + radiative heating
        nc_file.adv_thetal        = forcing_on
        nc_file.adv_qt            = forcing_on
        #nc_file.radiation         = 'tend'  #radiation isn't turned off in the CCPP SCM through the forcing
    elif (case_nml['case_config']['thermo_forcing_type'] == 2):
        #horizontal advective forcing + (radiative heating) + vertical velocity
        nc_file.adv_thetal        = forcing_on
        nc_file.adv_qt            = forcing_on
        #nc_file.radiation         = 'tend'  #radiation isn't turned off in the CCPP SCM through the forcing
        
        w_ls_avail = True if np.any(case_data._w_ls[:,:]) else False
        omega_avail = True if np.any(case_data._omega[:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using vertical velocity (through thermo_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
            #logging.critical(message)
            #raise Exception(message)
            
    elif (case_nml['case_config']['thermo_forcing_type'] == 3):
        #nudging + vertical velocity
        
        T_nudge_avail = True if np.any(case_data._T_nudge[:,:]) else False
        qt_nudge_avail = True if np.any(case_data._qt_nudge[:,:]) else False
        relax_time_avail = True if case_nml['case_config']['relax_time'] else False
        if relax_time_avail:
            relax_time = case_nml['case_config']['relax_time']
        else:
            relax_time = DEFAULT_NUDGING_TIMESCALE
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using nudging (through thermo_forcing = 3), but relax_time was not provided -- using default value of {1}s '.format(nml_filename, DEFAULT_NUDGING_TIMESCALE)
            logging.info(message)
        
        nc_file.nudging_ta = forcing_on*relax_time
        nc_file.nudging_qt = forcing_on*relax_time
        
        w_ls_avail = True if np.any(case_data._w_ls[:,:]) else False
        omega_avail = True if np.any(case_data._omega[:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using vertical velocity (through thermo_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
            #logging.critical(message)
            #raise Exception(message)
    
    time_dim   = nc_file.createDimension('time', case_data._time.shape[0])
    timei_dim  = nc_file.createDimension('t0',    1)
    lev_dim    = nc_file.createDimension('lev',  case_data._levels.shape[0])
    
    #
    timei_var                    = nc_file.createVariable('t0', wp, ('t0'))
    timei_var.units              = 'seconds since ' + start_date_string
    timei_var.standard_name      = 'Initial time'
    timei_var.calendar           = 'gregorian'
    timei_var[:]                 = 0.0
    #
    timef_var                    = nc_file.createVariable('time', wp, ('time'))
    timef_var.units              = 'seconds since ' + start_date_string
    timef_var.standard_name      = 'Forcing time'
    timef_var.calendar           = 'gregorian'
    timef_var[:]                 = case_data._time[:]
    #
    lev_var                      = nc_file.createVariable('lev', wp, ('lev'))
    lev_var.units                = 'Pa'
    lev_var.standard_name        = 'pressure'
    lev_var[:]                   = case_data._levels[:]
    
    #
    lon_var                      = nc_file.createVariable('lon', wp, ('time'))
    lon_var.units                = 'degrees_east'
    lon_var.standard_name        = 'longitude'
    lon_var[:]                   = case_data._lon

    #
    lat_var                      = nc_file.createVariable('lat', wp, ('time'))
    lat_var.units                = 'degrees_north'
    lat_var.standard_name        = 'latitude'
    lat_var[:]                   = case_data._lat
    
    #
    thetal_var                   = nc_file.createVariable('thetal', wp, ('t0','lev'))
    thetal_var.units             = 'K'
    thetal_var.standard_name     = 'air_liquid_potential_temperature'
    thetal_var[:]                = case_data._theta_il[:]
    
    #
    qt_var                       = nc_file.createVariable('qt', wp, ('t0','lev'))
    qt_var.units                 = 'kg kg-1'
    qt_var.standard_name         = 'mass_fraction_of_water_in_air'
    qt_var[:]                    = case_data._qt[:]
    
    #
    u_var                        = nc_file.createVariable('ua', wp, ('t0','lev'))
    u_var.units                  = 'm s-1'
    u_var.standard_name          = 'eastward_wind'
    u_var[:]                     = case_data._u[:]
    
    #
    v_var                        = nc_file.createVariable('va', wp, ('t0','lev'))
    v_var.units                  = 'm s-1'
    v_var.standard_name          = 'northward_wind'
    v_var[:]                     = case_data._v[:]
    
    #
    p_var                        = nc_file.createVariable('pa', wp, ('t0','lev'))
    p_var.units                  = 'Pa'
    p_var.standard_name          = 'air_pressure'
    p_var[:]                     = case_data._levels[:]
    
    #
    z_var                        = nc_file.createVariable('zh', wp, ('t0','lev'))
    z_var.units                  = 'm'
    z_var.standard_name          = 'height'
    z_var[:]                     = case_data._height[:]
    
    #
    ps_var                       = nc_file.createVariable('ps', wp, ('t0'))
    ps_var.units                 = 'Pa'
    ps_var.standard_name         = 'surface_air_pressure'
    ps_var[:]                    = case_data._p_surf[0]
    
    #
    ql_var                       = nc_file.createVariable('ql', wp, ('t0','lev'))
    ql_var.units                 = 'kg kg-1'
    ql_var.standard_name         = 'mass_fraction_of_cloud_liquid_water_in_air'
    ql_var[:]                    = case_data._ql[:]
    
    #
    qi_var                       = nc_file.createVariable('qi', wp, ('t0','lev'))
    qi_var.units                 = 'kg kg-1'
    qi_var.standard_name         = 'mass_fraction_of_cloud_ice_water_in_air'
    qi_var[:]                    = case_data._qi[:]
    
    #
    tke_var                      = nc_file.createVariable('tke', wp, ('t0','lev'))
    tke_var.units                = 'm2 s-2'
    tke_var.standard_name        = 'specific_turbulent_kinetic_energy'
    tke_var[:]                   = case_data._tke[:]
    
    #
    ozone_var                    = nc_file.createVariable('o3', wp, ('t0','lev'))
    ozone_var.units              = 'kg kg-1'
    ozone_var.standard_name      = 'mole_fraction_of_ozone_in_air'
    ozone_var[:]                 = case_data._ozone[:]
    
    ps_forc_var                  = nc_file.createVariable('ps_forc', wp, ('time'))
    ps_forc_var.units            = 'Pa'
    ps_forc_var.standard_name    = 'forcing_surface_air_pressure'
    ps_forc_var[:]               = case_data._p_surf[:]
    
    pa_forc_var                  = nc_file.createVariable('pa_forc', wp, ('time','lev'))
    pa_forc_var.units            = 'Pa'
    pa_forc_var.standard_name    = 'air_pressure_forcing'
    pa_forc_var[:]               = case_data._levels[:]
    
    zh_forc_var                  = nc_file.createVariable('zh_forc', wp, ('time','lev'))
    zh_forc_var.units            = 'm'
    zh_forc_var.standard_name    = 'height_forcing'
    zh_forc_var[:]               = case_data._height[:]
    
    if (nc_file.adv_ta == forcing_on):
        message = 'adv_ta is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnta_adv_var                    = nc_file.createVariable('tnta_adv', wp, ('time','lev'))
        # tnta_adv_var.units              = 'K s-1'
        # tnta_adv_var.standard_name      = 'tendency_of_air_temperature_due_to_advection'
        # tnta_adv_var[:]                 = np.swapaxes(case_data._tnta_adv[:],0,1)
        
    if (nc_file.adv_qv == forcing_on):
        message = 'adv_qv is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnqv_adv_var                    = nc_file.createVariable('tnqv_adv', wp, ('time','lev'))
        # tnqv_adv_var.units              = 'kg kg-1 s-1'
        # tnqv_adv_var.standard_name      = 'tendency_of_specific_humidity_due_to_advection'
        # tnqv_adv_var[:]                 = np.swapaxes(case_data._tnqv_adv[:],0,1)
    
    if (nc_file.adv_ua == forcing_on):
        message = 'adv_ua is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnua_adv_var                    = nc_file.createVariable('tnua_adv', wp, ('time','lev'))
        # tnua_adv_var.units              = 'm s-2'
        # tnua_adv_var.standard_name      = 'tendency_of_eastward_wind_due_to_advection'
        # tnua_adv_var[:]                 = np.swapaxes(case_data._tnua_adv[:],0,1)
    
    if (nc_file.adv_va == forcing_on):
        message = 'adv_va is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnva_adv_var                    = nc_file.createVariable('tnva_adv', wp, ('time','lev'))
        # tnva_adv_var.units              = 'm s-2'
        # tnva_adv_var.standard_name      = 'tendency_of_northward_wind_due_to_advection'
        # tnva_adv_var[:]                 = np.swapaxes(case_data._tnva_adv[:],0,1)
    
    if (nc_file.adv_theta == forcing_on):
        message = 'adv_theta is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tntheta_adv_var                    = nc_file.createVariable('tntheta_adv', wp, ('time','lev'))
        # tntheta_adv_var.units              = 'K s-1'
        # tntheta_adv_var.standard_name      = 'tendency_of_air_potential_temperature_due_to_advection'
        # tntheta_adv_var[:]                 = np.swapaxes(case_data._tntheta_adv[:],0,1)
    
    if (nc_file.adv_thetal == forcing_on):
        tnthetal_adv_var                    = nc_file.createVariable('tnthetal_adv', wp, ('time','lev'))
        tnthetal_adv_var.units              = 'K s-1'
        tnthetal_adv_var.standard_name      = 'tendency_of_air_liquid_potential_temperature_due_to_advection'
        if (nc_file.forc_wap == forcing_on or nc_file.forc_wa == forcing_on):
            tnthetal_adv_var[:]                 = np.swapaxes(case_data._h_advec_thil[:],0,1)
        else:
            tnthetal_adv_var[:]                 = np.swapaxes(case_data._h_advec_thil[:] + case_data._v_advec_thil[:],0,1)
    
    if (nc_file.adv_qt == forcing_on):
        tnqt_adv_var                    = nc_file.createVariable('tnqt_adv', wp, ('time','lev'))
        tnqt_adv_var.units              = 'kg kg-1 s-1'
        tnqt_adv_var.standard_name      = 'tendency_of_mass_fraction_of_water_in_air_due_to_advection'
        if (nc_file.forc_wap == forcing_on or nc_file.forc_wa == forcing_on):
            tnqt_adv_var[:]                 = np.swapaxes(case_data._h_advec_qt[:],0,1)
        else:
            tnqt_adv_var[:]                 = np.swapaxes(case_data._h_advec_qt[:] + case_data._v_advec_qt[:],0,1)
    
    if (nc_file.adv_rv == forcing_on):
        message = 'adv_rv is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnrv_adv_var                    = nc_file.createVariable('tnrv_adv', wp, ('time','lev'))
        # tnrv_adv_var.units              = 'kg kg-1 s-1'
        # tnrv_adv_var.standard_name      = 'tendency_of_humidity_mixing_ratio_due_to_advection'
        # tnrv_adv_var[:]                 = np.swapaxes(case_data._tnrv_adv[:],0,1)
    
    if (nc_file.adv_rt == forcing_on):
        message = 'adv_rt is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical()
        raise Exception(message)
        # tnrt_adv_var                    = nc_file.createVariable('tnrt_adv', wp, ('time','lev'))
        # tnrt_adv_var.units              = 'kg kg-1 s-1'
        # tnrt_adv_var.standard_name      = 'tendency_of_water_mixing_ratio_due_to_advection'
        # tnrt_adv_var[:]                 = np.swapaxes(case_data._tnrt_adv[:],0,1)
    
    if (nc_file.forc_wap == forcing_on):
        wap_var                    = nc_file.createVariable('wap', wp, ('time','lev'))
        wap_var.units              = 'Pa s-1'
        wap_var.standard_name      = 'lagrangian_tendency_of_air_pressure'
        wap_var[:]                 = np.swapaxes(case_data._omega[:],0,1)
    elif (nc_file.forc_wa == forcing_on):
        wa_var                    = nc_file.createVariable('wa', wp, ('time','lev'))
        wa_var.units              = 'm s-1'
        wa_var.standard_name      = 'upward_air_velocity'
        wa_var[:]                 = np.swapaxes(case_data._w_ls[:],0,1)
        
    if (nc_file.forc_geo == forcing_on):
        ug_var                    = nc_file.createVariable('ug', wp, ('time','lev'))
        ug_var.units              = 'm s-1'
        ug_var.standard_name      = 'geostrophic_eastward_wind'
        ug_var[:]                 = np.swapaxes(case_data._u_g[:],0,1)
        
        vg_var                    = nc_file.createVariable('vg', wp, ('time','lev'))
        vg_var.units              = 'm s-1'
        vg_var.standard_name      = 'geostrophic_northward_wind'
        vg_var[:]                 = np.swapaxes(case_data._v_g[:],0,1)
        
    if (nc_file.nudging_ua != forcing_off):
        ua_nud_var                    = nc_file.createVariable('ua_nud', wp, ('time','lev'))
        ua_nud_var.units              = 'm s-1'
        ua_nud_var.standard_name      = 'nudging_eastward_wind'
        ua_nud_var[:]                 = np.swapaxes(case_data._u_nudge[:],0,1)
    
    if (nc_file.nudging_va != forcing_off):
        va_nud_var                    = nc_file.createVariable('va_nud', wp, ('time','lev'))
        va_nud_var.units              = 'm s-1'
        va_nud_var.standard_name      = 'nudging_northward_wind'
        va_nud_var[:]                 = np.swapaxes(case_data._v_nudge[:],0,1)
    
    if (nc_file.nudging_ta != forcing_off):
        ta_nud_var                    = nc_file.createVariable('ta_nud', wp, ('time','lev'))
        ta_nud_var.units              = 'K'
        ta_nud_var.standard_name      = 'nudging_air_temperature'
        ta_nud_var[:]                 = np.swapaxes(case_data._T_nudge[:],0,1)
    
    if (nc_file.nudging_qt != forcing_off):
        qt_nud_var                    = nc_file.createVariable('qt_nud', wp, ('time','lev'))
        qt_nud_var.units              = 'kg kg-1'
        qt_nud_var.standard_name      = 'nudging_mass_fraction_of_water_in_air'
        qt_nud_var[:]                 = np.swapaxes(case_data._qt_nudge[:],0,1)
    
    if (case_nml['case_config']['sfc_flux_spec']):
        nc_file.surface_forcing_temp = 'kinematic'
        nc_file.surface_forcing_moisture = 'kinematic'
        nc_file.surface_forcing_wind = 'z0'
        
        wpthetap_s_var                  = nc_file.createVariable('wpthetap_s', wp, ('time'))
        wpthetap_s_var.units            = 'K m s-1'
        wpthetap_s_var.standard_name    = 'surface_upward_potential_temperature_flux'
        wpthetap_s_var[:]               = case_data._sh_flux_sfc[:]
        
        wpqtp_s_var                  = nc_file.createVariable('wpqvp_s', wp, ('time'))
        wpqtp_s_var.units            = 'kg kg-1 m s-1'
        wpqtp_s_var.standard_name    = 'surface_upward_specific_humidity_flux'
        wpqtp_s_var[:]               = case_data._lh_flux_sfc[:]
        
        z0_var                  = nc_file.createVariable('z0', wp, ('time'))
        z0_var.units            = 'm'
        z0_var.standard_name    = 'surface_roughness_length_for_momentum_in_air'
        z0_var[:]               = case_nml['case_config']['sfc_roughness_length_cm']*1.0E-2
        
        
        ts_var                  = nc_file.createVariable('ts_forc', wp, ('time'))
        ts_var.units            = 'K'
        ts_var.standard_name    = 'forcing_surface_temperature'
        if np.any(case_data._T_surf[:]):
            ts_var[:]               = case_data._T_surf[:]
        else:
            ts_var[:]               = case_data._missing_value
        
    else:
        nc_file.surface_forcing_temp = 'ts'
        nc_file.surface_forcing_wind = 'z0'
        
        z0_var                  = nc_file.createVariable('z0', wp, ('time'))
        z0_var.units            = 'm'
        z0_var.standard_name    = 'surface_roughness_length_for_momentum_in_air'
        z0_var[:]               = case_nml['case_config']['sfc_roughness_length_cm']*1.0E-2
        
        ts_var                  = nc_file.createVariable('ts_forc', wp, ('time'))
        ts_var.units            = 'K'
        ts_var.standard_name    = 'forcing_surface_temperature'
        ts_var[:]               = case_data._T_surf[:]
    
    nc_file.close()
    
    return(fileOUT)

def write_SCM_nml_file(case_nml):
    filename = os.path.join(CASE_NML_DIR, case_nml['case_config']['case_name'] + '_dephy.nml')
    
    #Go through existing case namelist and only add necessary items to new DEPHY-based namelist
    
    #add _dephy to case (temporary - to differentiate from old format case)
    int_dict = {'case_name':case_nml['case_config']['case_name']+'_dephy',
                'input_type':1}
    
    nml_keys = case_nml['case_config'].todict().keys()
    if ('npz_type' in nml_keys):
        int_dict['npz_type'] = case_nml['case_config']['npz_type']
        if int_dict['npz_type'] == 'input' and 'vert_coord_file' in nml_keys:
            int_dict['vert_coord_file'] = case_nml['case_config']['vert_coord_file']
    
    if ('dt' in nml_keys):
        int_dict['dt'] = case_nml['case_config']['dt']
    
    #runtime is in netCDF file
    
    if ('output_dir' in nml_keys):
        int_dict['output_dir'] = case_nml['case_config']['output_dir']
        
    if ('model_ics' in nml_keys):
        int_dict['model_ics'] = case_nml['case_config']['model_ics']
    
    if ('lsm_ics' in nml_keys):
        int_dict['lsm_ics'] = case_nml['case_config']['lsm_ics']
    
    if ('do_spinup' in nml_keys):
        int_dict['do_spinup'] = case_nml['case_config']['do_spinup']
    
    if ('spinup_timesteps' in nml_keys):
        int_dict['spinup_timesteps'] = case_nml['case_config']['spinup_timesteps']
    
    if ('C_RES' in nml_keys):
        int_dict['C_RES'] = case_nml['case_config']['C_RES']
    
    #relax_time is in netCDF file
    
    #sfc_type is in netCDF file
    
    #sfc_flux_spec is in netCDF file
    
    if ('sfc_roughness_length_cm' in nml_keys):
        int_dict['sfc_roughness_length_cm'] = case_nml['case_config']['sfc_roughness_length_cm']
    
    if ('reference_profile_choice' in nml_keys):
        int_dict['reference_profile_choice'] = case_nml['case_config']['reference_profile_choice']
    
    if ('column_area' in nml_keys):
        int_dict['column_area'] = case_nml['case_config']['column_area']
    
    nml_dict = {'case_config':int_dict}
    
    nml = f90nml.namelist.Namelist(nml_dict)
    
    #print(nml)
    nml.write(filename, force=True)
    
    return(filename)

########################################################################################
#
########################################################################################
def write_SCM_case_file_UFS(state, surface, oro, forcing, case, date, vertical_forcing):
    """Write all data to a netCDF file in the DEPHY-SCM format"""

    # Working types
    wp = np.float64
    wi = np.int32

    # Local switches
    forcing_on  = 1
    forcing_off = 0

    # Output file
    com = 'mkdir -p ' + PROCESSED_CASE_DIR
    print(com)
    os.system(com)
    fileOUT = os.path.join(PROCESSED_CASE_DIR, case + '_SCM_driver.nc')

    nc_file = Dataset(fileOUT, 'w', format='NETCDF3_CLASSIC')
    nc_file.description = "FV3GFS model profile input (UFS forcings)"

    nc_file.missing_value   = missing_value

    start_date = datetime(date["year"],date["month"],date["day"],date["hour"],date["minute"],date["second"])

    #
    # Create surface type string (Saved as GLOBAL attribute)
    #
    if surface["slmsk"] > 1.5:
        surface_string = 'ice'
    elif surface["slmsk"] > 0.5:
        surface_string = 'land'
    else:
        surface_string = 'ocean'

    #
    # Global file attributes.
    #
    runtime           = timedelta(seconds=forcing['time'][-1])
    end_date          = start_date + runtime
    end_date_string   = end_date.strftime("%Y-%m-%d %H:%M:%S")
    start_date_string = start_date.strftime("%Y-%m-%d %H:%M:%S")
    #
    loc_string  = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string
    #
    nc_file.case              = case_string
    nc_file.title             = 'Forcing and Initial Conditions for ' + case_string
    nc_file.reference         = 'https://dtcenter.org/sites/default/files/paragraph/scm-ccpp-guide-v6-0-0.pdf'
    nc_file.author            = 'Grant J. Firl and Dustin Swales'
    nc_file.version           = 'Created on ' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    nc_file.format_version    = 'DEPHY SCM format version 1'
    nc_file.modifications     = ''
    nc_file.script            = os.path.basename(__file__)
    nc_file.comment           = ''
    nc_file.start_date        = start_date_string
    nc_file.end_date          = end_date_string
    nc_file.forcing_scale     = -1
    nc_file.radiation         = "off"
    nc_file.adv_ta            = forcing_off
    nc_file.adv_qv            = forcing_off
    nc_file.adv_ua            = forcing_off
    nc_file.adv_va            = forcing_off
    nc_file.adv_theta         = forcing_off
    nc_file.adv_thetal        = forcing_off
    nc_file.adv_qt            = forcing_off
    nc_file.adv_rv            = forcing_off
    nc_file.adv_rt            = forcing_off
    if (vertical_forcing == 2):
        nc_file.forc_wa       = forcing_on
    else:
        nc_file.forc_wa       = forcing_off
    nc_file.forc_wap          = forcing_off
    nc_file.forc_geo          = forcing_off
    nc_file.nudging_ua        = forcing_off
    nc_file.nudging_va        = forcing_off
    nc_file.nudging_ta        = forcing_off
    nc_file.nudging_theta     = forcing_off
    nc_file.nudging_thetal    = forcing_off
    nc_file.nudging_qv        = forcing_off
    nc_file.nudging_qt        = forcing_off
    nc_file.nudging_rv        = forcing_off
    nc_file.nudging_rt        = forcing_off
    nc_file.zh_nudging_ta     = forcing_off
    nc_file.zh_nudging_theta  = forcing_off
    nc_file.zh_nudging_thetal = forcing_off
    nc_file.zh_nudging_qv     = forcing_off
    nc_file.zh_nudging_qt     = forcing_off
    nc_file.zh_nudging_rv     = forcing_off
    nc_file.zh_nudging_rt     = forcing_off
    nc_file.zh_nudging_ua     = forcing_off
    nc_file.zh_nudging_va     = forcing_off
    nc_file.pa_nudging_ta     = forcing_off
    nc_file.pa_nudging_theta  = forcing_off
    nc_file.pa_nudging_thetal = forcing_off
    nc_file.pa_nudging_qv     = forcing_off
    nc_file.pa_nudging_qt     = forcing_off
    nc_file.pa_nudging_rv     = forcing_off
    nc_file.pa_nudging_rt     = forcing_off
    nc_file.pa_nudging_ua     = forcing_off
    nc_file.pa_nudging_va     = forcing_off
    #
    nc_file.surface_type      = surface_string
    #
    nc_file.adv_ta            = forcing_on
    nc_file.adv_qv            = forcing_on
    nc_file.adv_ua            = forcing_on
    nc_file.adv_va            = forcing_on
    #
    nc_file.surface_forcing_temp     = 'none'
    nc_file.surface_forcing_moisture = 'none'
    nc_file.surface_forcing_wind     = 'none'
    nc_file.surface_forcing_lsm      = 'none' #'noah' #'noahmp' #'ruc'
    nc_file.surface_forcing_lsm  = 'lsm'
    # Set file dimension
    time_dim   = nc_file.createDimension('time', len(forcing['time']))
    timei_dim  = nc_file.createDimension('t0',    1)
    lev_dim    = nc_file.createDimension('lev',   state["nlevs"])
    soil_dim   = nc_file.createDimension('nsoil', len(surface["stc"]))
    snow_dim   = nc_file.createDimension('nsnow', len(surface["snicexy"]))
    nslsnw_dim = nc_file.createDimension('nsoil_plus_nsnow',len(surface["snicexy"]) + len(surface["stc"]))
    ice_dim    = nc_file.createDimension('nice',  len(surface["tiice"]))
    
    #
    timei_var                    = nc_file.createVariable('t0', wp, ('t0'))
    timei_var.units              = 'seconds since ' + start_date_string
    timei_var.standard_name      = 'Initial time'
    timei_var.calendar           = 'gregorian'
    timei_var[:]                 = 0.0
    #
    timef_var                    = nc_file.createVariable('time', wp, ('time'))
    timef_var.units              = 'seconds since ' + start_date_string
    timef_var.standard_name      = 'Forcing time'
    timef_var.calendar           = 'gregorian'
    timef_var[:]                 = forcing['time']
    #
    lev_var                      = nc_file.createVariable('lev', wp, ('lev'))
    lev_var.units                = 'm'
    lev_var.standard_name        = 'height'
    lev_var[:]                   = 0.0

    #
    lon_var                      = nc_file.createVariable('lon', wp, ('time'))
    lon_var.units                = 'degrees_east'
    lon_var.standard_name        = 'longitude'
    lon_var[:]                   = surface["lon"]

    #
    lat_var                      = nc_file.createVariable('lat', wp, ('time'))
    lat_var.units                = 'degrees_north'
    lat_var.standard_name        = 'latitude'
    lat_var[:]                   = surface["lat"]

    #
    soil_depth_var               = nc_file.createVariable('soil_depth', wp, ('nsoil'))
    soil_depth_var.units         = 'm'
    soil_depth_var.standard_name = 'depth of bottom of soil layers'
    soil_depth_var[:]            = [0.1,0.4,1.0,2.0]
    #
    theta_oro                    = nc_file.createVariable('theta_oro',wp, ('t0'))
    theta_oro.units              = "deg"
    theta_oro.standard_name      = "angle with respect to east of maximum subgrid orographic variations"
    theta_oro[:]                 = oro["theta"]
    #
    z0_var                       = nc_file.createVariable('zorl', wp, ('time'))
    z0_var.units                 =  "cm"
    z0_var.standard_name         = 'surface_roughness_length_for_momentum_in_air'
    z0_var[:]                    = surface["z0"]
    #
    zorlw_var                  = nc_file.createVariable('zorlw', wp, ('t0'))
    zorlw_var.units            = "cm"
    zorlw_var.standard_name    = "surface roughness length over ocean"
    zorlw_var[:]               = surface["z0"]
    #
    zorll_var                  = nc_file.createVariable('zorll', wp, ('t0'))
    zorll_var.units            = "cm"
    zorll_var.standard_name    = "surface roughness length over land"
    zorll_var[:]               = surface["zorll"]
    #
    zorli_var                  = nc_file.createVariable('zorli', wp, ('t0'))
    zorli_var.units            = "cm"
    zorli_var.standard_name    = "surface roughness length over ice"
    zorli_var[:]               = surface["zorli"]
    #
    zorlwav_var                = nc_file.createVariable('zorlwav', wp, ('time'))
    zorlwav_var.units          =  "cm"
    zorlwav_var.standard_name  = 'surface_roughness_length_from_wave_model'
    zorlwav_var[:]             = surface["zorlw"]

    #
    # Variables to be output to SCM input file. Only fields that come directly from forcing, 
    # surface, state, and oro. Fields that get renamed are done above.
    #
    dict = {}
    dict.update(date)
    dict.update(surface)
    dict.update(state)
    dict.update(oro)
    dict.update(forcing)

    ########################################################################################
    #
    # Dictonary format:
    # {"name": "", "type", "dimd": (), "units": "", "desc": ""}
    #
    ######################################################################################## 
    var_dict = [{"name": "orog",         "type":wp, "dimd": ('t0'         ),    "units": "m",             "desc": "surface_altitude"},\
                {"name": "zh",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m",             "desc": "height"},\
                {"name": "pa",           "type":wp, "dimd": ('t0',   'lev'),    "units": "Pa",            "desc": "air_pressure"}, \
                {"name": "ta",           "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_temperature"}, \
                {"name": "theta",        "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_potential_temperature"}, \
                {"name": "thetal",       "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_liquid_potential_temperature"}, \
                {"name": "rv",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "humidity_mixing_ratio"}, \
                {"name": "rl",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "cloud_liquid_water_mixing_ratio"}, \
                {"name": "ri",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "cloud_ice_water_mixing_ratio"}, \
                {"name": "rt",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "water_mixing_ratio"}, \
                {"name": "qv",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "specific_humidity"}, \
                {"name": "ql",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_cloud_liquid_water_in_air"}, \
                {"name": "qi",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_cloud_ice_water_in_air", "default_value": 0.0}, \
                {"name": "qt",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_water_in_air"}, \
                {"name": "hur",          "type":wp, "dimd": ('t0',   'lev'),    "units": "%",             "desc": "relative_humidity"}, \
                {"name": "tke",          "type":wp, "dimd": ('t0',   'lev'),    "units": "m2 s-2",        "desc": "specific_turbulen_kinetic_energy", "default_value": 0.0}, \
                {"name": "ua",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "eastward_wind"}, \
                {"name": "va",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "northward_wind"}, \
                {"name": "ts",           "type":wp, "dimd": ('t0'         ),    "units": "K",             "desc": "surface_temperature"},\
                {"name": "tskin",        "type":wp, "dimd": ('t0'         ),    "units": "K",             "desc": "surface_skin_temperature"}, \
                {"name": "ps",           "type":wp, "dimd": ('t0'         ),    "units": "Pa",            "desc": "surface_air_pressure"}, \
                {"name": "beta",         "type":wp, "dimd": ('t0'         ),    "units": "m",             "desc": "soil_water_stress_factor"}, \
                {"name": "mrsos",        "type":wp, "dimd": ('t0'         ),    "units": "kg m-2",        "desc": "mass_content_of_water_in_soil_layer"}, \
                {"name": "o3",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mole_fraction_of_ozone_in_air"}, \
                {"name": "sza",          "type":wp, "dimd": ('t0'         ),    "units": "degree",        "desc": "solar_zenith_angle"}, \
                {"name": "io",           "type":wp, "dimd": ('t0'         ),    "units": "W m-2",         "desc": "solar_irradiance"}, \
                {"name": "alb",          "type":wp, "dimd": ('t0'         ),    "units": "1",             "desc": "surface_albedo"}, \
                {"name": "emis",         "type":wp, "dimd": ('t0'         ),    "units": "1",             "desc": "surface_longwave_emissivity"}, \
                {"name": "slmsk",        "type":wp, "dimd": ('t0'         ),    "units": "none",          "desc": "land_sea_ice_mask"}]
    #
    var_frc  = [{"name": "zh_forc",      "type":wp, "dimd": ('time', 'lev'),    "units": "m",             "desc": "height_forcing","default_value": 1.},\
                {"name": "pa_forc",      "type":wp, "dimd": ('time', 'lev'),    "units": "Pa",            "desc": "air_pressure_forcing"}, \
                {"name": "wa",           "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "upward_air_velocity"}, \
                {"name": "wap",          "type":wp, "dimd": ('time', 'lev'),    "units": "Pa s-1",        "desc": "lagrangian_tendency_of_air_pressure"}, \
                {"name": "ug",           "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "geostrophic_eastward_wind"}, \
                {"name": "vg",           "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "geostrophic_northward_wind"}, \
                {"name": "tnua_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "m s-2",         "desc": "tendency_of_eastward_wind_due_to_advection"},\
                {"name": "tnva_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "m s-2",         "desc": "tendency_of_northward_wind_due_to_advection"}, \
                {"name": "tnta_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_air_temperature_due_to_advection"}, \
                {"name": "tntheta_adv",  "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_air_potential_temperature_due_to_advection"}, \
                {"name": "tnthetal_adv", "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_air_liquid_potential_temperature_due_to_advection"}, \
                {"name": "tnqv_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "kg kg-1 s-1",   "desc": "tendency_of_specific_humidity_due_to_advection"},\
                {"name": "tnqt_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "kg kg-1 s-1",   "desc": "tendency_of_mass_fraction_of_water_in_air_due_to_advection"},\
                {"name": "tnrv_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "kg kg-1 s-1",   "desc": "tendency_of_humidity_mixing_ratio_due_to_advection"},\
                {"name": "tnrt_adv",     "type":wp, "dimd": ('time', 'lev'),    "units": "kg kg-1 s-1",   "desc": "tendency_of_water_mixing_ratio_due_to_advection"},\
                {"name": "tnta_rad",     "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_air_temperature_due_to_radiative_heating"}, \
                {"name": "tntheta_rad",  "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_potential_air_temperature_due_to_radiative_heating"}, \
                {"name": "tnthetal_rad", "type":wp, "dimd": ('time', 'lev'),    "units": "K s-1",         "desc": "tendency_of_air_liquid_potential_temperature_due_to_radiative_heating"}, \
                {"name": "ta_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "K",             "desc": "nudging_air_temperature"}, \
                {"name": "theta_nud",    "type":wp, "dimd": ('time', 'lev'),    "units": "K",             "desc": "nudging_air_potential_temperature"}, \
                {"name": "thetal_nud",   "type":wp, "dimd": ('time', 'lev'),    "units": "K",             "desc": "nudging_air_liquid_potential_temperature"}, \
                {"name": "qt_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "kg kg-1",       "desc": "nudging_mass_fraction_of_water_in_air"}, \
                {"name": "rv_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "nudging_humidity_mixing_ratio"}, \
                {"name": "rt_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "nudging_water_mixing_ratio"}, \
                {"name": "ua_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "nudging_eastward_wind"}, \
                {"name": "va_nud",       "type":wp, "dimd": ('time', 'lev'),    "units": "m s-1",         "desc": "nudging_northward_wind"}, \
                {"name": "hfss",         "type":wp, "dimd": ('time'       ),    "units": "W m-2",         "desc": "surface_upward_sensible_heat_flux"}, \
                {"name": "hfls",         "type":wp, "dimd": ('time'       ),    "units": "W m-2",         "desc": "surface_upward_latent_heat_flux"}, \
                {"name": "wpthetap_s",   "type":wp, "dimd": ('time'       ),    "units": "K m s-1",       "desc": "surface_upward_potential_temperature_flux"}, \
                {"name": "wpqvp_s",      "type":wp, "dimd": ('time'       ),    "units": "m s-1",         "desc": "surface_upward_specific_humidity_flux"}, \
                {"name": "wpqtp_s",      "type":wp, "dimd": ('time'       ),    "units": "m s-1",         "desc": "surface_upward_water_mass_fraction_flux"}, \
                {"name": "wprvp_s",      "type":wp, "dimd": ('time'       ),    "units": "m s-1",         "desc": "surface_upward_humidity_mixing_ratio_flux"}, \
                {"name": "wprtp_s",      "type":wp, "dimd": ('time'       ),    "units": "m s-1",         "desc": "surface_upward_water_mixing_ratio_flux"}, \
                {"name": "ts_forc",      "type":wp, "dimd": ('time'       ),    "units": "K",             "desc": "forcing_surface_temperature"},\
                {"name": "ps_forc",      "type":wp, "dimd": ('time'       ),    "units": "Pa",            "desc": "forcing_surface_air_pressure"},\
                {"name": "uustar",       "type":wp, "dimd": ('time'       ),    "units": "m s-1",         "desc": "surface_friction_velocity"}, \
                {"name": "z0h",          "type":wp, "dimd": ('time'       ),    "units": "m",             "desc": "surface_roughness_length_for_heat_in_air"}, \
                {"name": "z0q",          "type":wp, "dimd": ('time'       ),    "units": "m",             "desc": "surface_roughness_length_for_humidity_in_air"}, \
                {"name": "mrsos_forc",   "type":wp, "dimd": ('time'       ),    "units": "kg m-2",        "desc": "forcing_mass_content_of_water_in_soil_layer"}]

    #
    var_oro  = [{"name": "area",         "type":wp, "dimd": ('t0'),             "units": "m2",      "desc": "grid_cell_area"},\
                {"name": "stddev",       "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "standard deviation of subgrid orography"}, \
                {"name": "convexity",    "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "convexity of subgrid orography"}, \
                {"name": "oa1",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "assymetry of subgrid orography 1"}, \
                {"name": "oa2",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "assymetry of subgrid orography 2"}, \
                {"name": "oa3",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "assymetry of subgrid orography 3"}, \
                {"name": "oa4",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "assymetry of subgrid orography 4"}, \
                {"name": "ol1",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of grid box with subgrid orography higher than critical height 1"}, \
                {"name": "ol2",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of grid box with subgrid orography higher than critical height 2"}, \
                {"name": "ol3",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of grid box with subgrid orography higher than critical height 3"}, \
                {"name": "ol4",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of grid box with subgrid orography higher than critical height 4"}, \
                {"name": "sigma",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "slope of subgrid orography"}, \
                {"name": "gamma",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "anisotropy of subgrid orography"}, \
                {"name": "elvmax",       "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "maximum of subgrid orography"}, \
                {"name": "oro",          "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "orography"}, \
                {"name": "oro_uf",       "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "unfiltered orography"}, \
                {"name": "landfrac",     "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of horizontal grid area occupied by land"}, \
                {"name": "lakefrac",     "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fraction of horizontal grid area occupied by lake", "default_value":0}, \
                {"name": "lakedepth",    "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "lake depth", "default_value":0}]
    #
    var_nsst = [{"name": "tref",         "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "sea surface reference temperature for NSST"}, \
                {"name": "z_c",          "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "sub-layer cooling thickness for NSST"}, \
                {"name": "c_0",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "coefficient 1 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "c_d",          "type":wp, "dimd": ('t0'),             "units": "nonw",    "desc": "coefficient 2 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "w_0",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "coefficient 3 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "w_d",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "coefficient 4 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "xt",           "type":wp, "dimd": ('t0'),             "units": "K m",     "desc": "heat content in diurnal thermocline layer for NSST"}, \
                {"name": "xs",           "type":wp, "dimd": ('t0'),             "units": "ppt m",   "desc": "salinity content in diurnal thermocline layer for NSST"}, \
                {"name": "xu",           "type":wp, "dimd": ('t0'),             "units": "m2 s-1",  "desc": "u-current in diurnal thermocline layer for NSST"}, \
                {"name": "xv",           "type":wp, "dimd": ('t0'),             "units": "m2 s-1",  "desc": "v-current in diurnal thermocline layer for NSST"}, \
                {"name": "xz",           "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "thickness of diurnal thermocline layer for NSST"}, \
                {"name": "zm"   ,        "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "thickness of ocean mixed layer for NSST"}, \
                {"name": "xtts",         "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST"},\
                {"name": "xzts",         "type":wp, "dimd": ('t0'),             "units": "m K-1",   "desc": "sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST"}, \
                {"name": "d_conv",       "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "thickness of free convection layer for NSST"}, \
                {"name": "ifd",          "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "index to start DTM run for NSST"}, \
                {"name": "dt_cool",      "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "sub-layer cooling amount for NSST"}, \
                {"name": "qrain",        "type":wp, "dimd": ('t0'),             "units": "W m-2",   "desc": "sensible heat due to rainfall for NSST"}]
    #
    var_frgd = [{"name": "tiice",        "type":wp, "dimd": ('t0','nice'),      "units": "K",       "desc": "sea ice internal temperature"}]
    #
    var_noah = [{"name": "vegsrc",       "type":wi, "dimd": ('t0'),             "units": "none",    "desc": "vegetation source (1-2)", "default_value": 1}, \
                {"name": "tsfco",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "sea/skin/ice surface temperature"}, \
                {"name": "weasd",        "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "water equivalent accumulated snow depth"}, \
                {"name": "tg3",          "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "deep soil temperature"}, \
                {"name": "alvsf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "60 degree vis albedo with strong cosz dependency"}, \
                {"name": "alnsf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "60 degree nir albedo with strong cosz dependency"}, \
                {"name": "alvwf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "60 degree vis albedo with weak cosz dependency"}, \
                {"name": "alnwf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "60 degree nir albedo with weak cosz dependency"}, \
                {"name": "facsf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fractional coverage with strong cosz dependency"}, \
                {"name": "facwf",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "fractional coverage with weak cosz dependency"}, \
                {"name": "vegfrac",      "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "vegetation fraction"}, \
                {"name": "canopy",       "type":wp, "dimd": ('t0'),             "units": "kg m-2",  "desc": "amount of water stored in canopy"}, \
                {"name": "f10m",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "ratio of sigma level 1 wind and 10m wind"}, \
                {"name": "t2m",          "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "2-meter absolute temperature"}, \
                {"name": "q2m",          "type":wp, "dimd": ('t0'),             "units": "kg kg-1", "desc": "2-meter specific humidity"}, \
                {"name": "vegtyp",       "type":wi, "dimd": ('t0'),             "units": "none",    "desc": "vegetation type (1-12)"}, \
                {"name": "soiltyp",      "type":wi, "dimd": ('t0'),             "units": "none",    "desc": "soil type (1-12)"}, \
                {"name": "scolor",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "soil color"}, \
                {"name": "ffmm",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "Monin-Obukhov similarity function for momentum"}, \
                {"name": "ffhh",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "Monin-Obukhov similarity function for heat"}, \
                {"name": "hice",         "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "sea ice thickness"}, \
                {"name": "fice",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "ice fraction"}, \
                {"name": "tisfc",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "ice surface temperature"}, \
                {"name": "tprcp",        "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "instantaneous total precipitation amount"}, \
                {"name": "srflag",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "snow/rain flag for precipitation"}, \
                {"name": "snowd",        "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "water equivalent snow depth"}, \
                {"name": "shdmin",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "minimum vegetation fraction"}, \
                {"name": "shdmax",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "maximum vegetation fraction"}, \
                {"name": "slopetyp",     "type":wi, "dimd": ('t0'),             "units": "none",    "desc": "slope type (1-9)"}, \
                {"name": "snoalb",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "maximum snow albedo"}, \
                {"name": "sncovr",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "surface snow area fraction"}, \
                {"name": "tsfcl",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "surface skin temperature over land"}, \
                {"name": "stc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "K",       "desc": "initial profile of soil liquid moisture"}, \
                {"name": "smc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "kg",      "desc": "initial profile of soil moisture"}, \
                {"name": "slc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "kg",      "desc": "initial profile of soil temperature"}]
    #
    var_noahmp=[{"name": "tvxy",         "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "vegetation temperature for NoahMP"}, \
                {"name": "tgxy",         "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "ground temperature for NoahMP"}, \
                {"name": "tahxy",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "canopy air temperature for NoahMP"}, \
                {"name": "canicexy",     "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "canopy intercepted ice mass for NoahMP"}, \
                {"name": "canliqxy",     "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "canopy intercepted liquid water for NoahMP"}, \
                {"name": "eahxy",        "type":wp, "dimd": ('t0'),             "units": "Pa",      "desc": "canopy air vapor pressure for NoahMP"}, \
                {"name": "cmxy",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "surface drag coefficient for momentum for NoahMP"}, \
                {"name": "chxy",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "surface exchange coeff heat & moisture for NoahMP"}, \
                {"name": "fwetxy",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "area fraction of canopy that is wetted/snowed for NoahMP"}, \
                {"name": "sneqvoxy",     "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "snow mass at previous time step for NoahMP"}, \
                {"name": "alboldxy",     "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "snow albedo at previous time step for NoahMP"}, \
                {"name": "qsnowxy",      "type":wp, "dimd": ('t0'),             "units": "mm s-1",  "desc": "snow precipitation rate at surface for NoahMP"}, \
                {"name": "wslakexy",     "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "lake water storage for NoahMP"}, \
                {"name": "taussxy",      "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "non-dimensional snow age for NoahMP"}, \
                {"name": "waxy",         "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "water storage in aquifer for NoahMP"}, \
                {"name": "wtxy",         "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "water storage in aquifer and saturated soil for NoahMP"}, \
                {"name": "zwtxy",        "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "water table depth for NoahMP"}, \
                {"name": "xlaixy",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "leaf area index for NoahMP"}, \
                {"name": "xsaixy",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "stem area index for NoahMP"}, \
                {"name": "lfmassxy",     "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "leaf mass for NoahMP"}, \
                {"name": "stmassxy",     "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "stem mass for NoahMP"}, \
                {"name": "rtmassxy",     "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "fine root mass for NoahMP"}, \
                {"name": "woodxy",       "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "wood mass including woody roots for NoahMP"}, \
                {"name": "stblcpxy",     "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "stable carbon in deep soil for NoahMP"}, \
                {"name": "fastcpxy",     "type":wp, "dimd": ('t0'),             "units": "g m-2",   "desc": "short-lived carbon in shallow soil for NoahMP"}, \
                {"name": "smcwtdxy",     "type":wp, "dimd": ('t0'),             "units": "m3 m-3",  "desc": "soil water content between the bottom of the soil and the water table for NoahMP"}, \
                {"name": "deeprechxy",   "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "recharge to or from the water table when deep for NoahMP"}, \
                {"name": "rechxy",       "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "recharge to or from the water table when shallow for NoahMP"}, \
                {"name": "snowxy",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "number of snow layers for NoahMP"}, \
                {"name": "snicexy",      "type":wp, "dimd": ('t0','nsnow'),     "units": "mm",      "desc": "initial profile of snow layer ice"}, \
                {"name": "snliqxy",      "type":wp, "dimd": ('t0','nsnow'),     "units": "mm",      "desc": "initial profile of snow layer liquid"}, \
                {"name": "tsnoxy",       "type":wp, "dimd": ('t0','nsnow'),     "units": "K",       "desc": "initial profile of snow layer temperature"}, \
                {"name": "smoiseq",      "type":wp, "dimd": ('t0','nsoil'),     "units": "m3 m-3",  "desc": "initial profile of equilibrium soil water content"}, \
                {"name": "zsnsoxy",      "type":wp, "dimd": ('t0','nsoil_plus_nsnow'), "units": "m","desc": "layer bottom depth from snow surface"}]
    #
    var_ruc  = [{"name": "wetness",          "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "normalized soil wetness for RUC LSM"}, \
                {"name": "lai",              "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "leaf area index for RUC LSM"}, \
                {"name": "clw_surf_land",    "type":wp, "dimd": ('t0'),         "units": "kg kg-1", "desc": "cloud condensed water mixing ratio at surface over land for RUC LSM"},\
                {"name": "clw_surf_ice",     "type":wp, "dimd": ('t0'),         "units": "kg kg-1", "desc": "cloud condensed water mixing ratio at surface over ice for RUC LSM"},\
                {"name": "qwv_surf_land",    "type":wp, "dimd": ('t0'),         "units": "kg kg-1", "desc": "water vapor mixing ratio at surface over land for RUC LSM"},\
                {"name": "qwv_surf_ice",     "type":wp, "dimd": ('t0'),         "units": "kg kg-1", "desc": "water vapor mixing ratio at surface over ice for RUC LSM"},\
                {"name": "tsnow_land",       "type":wp, "dimd": ('t0'),         "units": "K",       "desc": "snow temperature at the bottom of the first snow layer over land for RUC LSM"},\
                {"name": "tsnow_ice",        "type":wp, "dimd": ('t0'),         "units": "K",       "desc": "snow temperature at the bottom of the first snow layer over ice for RUC LSM"},\
                {"name": "snowfall_acc_land","type":wp, "dimd": ('t0'),         "units": "kg m-2",  "desc": "run-total snow accumulation on the ground over land for RUC LSM"},\
                {"name": "snowfall_acc_ice", "type":wp, "dimd": ('t0'),         "units": "kg m-2",  "desc": "run-total snow accumulation on the ground over ice for RUC LSM"},\
                {"name": "sfalb_lnd",        "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "surface albedo over land for RUC LSM"},\
                {"name": "sfalb_lnd_bck",    "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "surface snow-free albedo over land for RUC LSM"},\
                {"name": "sfalb_ice",        "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "surface albedo over ice for RUC LSM"},\
                {"name": "emis_lnd",         "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "surface emissivity over land for RUC LSM"},\
                {"name": "emis_ice",         "type":wp, "dimd": ('t0'),         "units": "none",    "desc": "surface emissivity over ice for RUC LSM"}, \
                {"name": "tslb",             "type":wp, "dimd": ('t0','nsoil'), "units": "K",       "desc": "soil temperature for RUC LSM"}, \
                {"name": "smois",            "type":wp, "dimd": ('t0','nsoil'), "units": "none",    "desc": "volume fraction of soil moisture for RUC LSM"}, \
                {"name": "sh2o",             "type":wp, "dimd": ('t0','nsoil'), "units": "none",    "desc": "volume fraction of unfrozen soil moisture for RUC LSM"}, \
                {"name": "smfr",             "type":wp, "dimd": ('t0','nsoil'), "units": "none",    "desc": "volume fraction of frozen soil moisture for RUC LSM"},
                {"name": "flfr",             "type":wp, "dimd": ('t0','nsoil'), "units": "none",    "desc": "flag for frozen soil physics for RUC LSM"}]

    #
    var_dict.extend(var_frc)
    var_dict.extend(var_oro)
    var_dict.extend(var_nsst)
    var_dict.extend(var_frgd)
    var_dict.extend(var_ruc)
    var_dict.extend(var_noah)
    var_dict.extend(var_noahmp)

    #
    for var in var_dict:
        if (var["name"] in dict):
            var_temp               = nc_file.createVariable(var["name"], var["type"], var["dimd"])
            var_temp.units         = var["units"]
            var_temp.standard_name = var["desc"]
            var_temp[:]            = dict[var["name"]]
        elif "default_value" in var:
            var_temp               = nc_file.createVariable(var["name"], var["type"], var["dimd"])
            var_temp.units         = var["units"]
            var_temp.standard_name = var["desc"]
            var_temp[:]            = var["default_value"]
        if "override" in var:
            var_temp[:]            = var["default_value"]
    #
    # Close file
    #
    nc_file.close()

    return(fileOUT)

########################################################################################
#
########################################################################################    
def main():
    setup_logging()
    
    #read in arguments
    (case_name, use_area) = parse_arguments()
    
    (case_nml, error) = get_case_nml(case_name)
    if (error):
        message = 'The directory {0} does not contain a config file for case {1}'.format(CASE_NML_DIR, case_name)
        logging.critical(message)
        raise Exception(message)
    else:
        logging.info(case_nml)
        
    (case_data, error) = get_case_data(case_name)
    if (error):
        message = 'The directory {0} does not contain a data file for case {1}'.format(PROCESSED_CASE_DIR, case_name)
        logging.critical(message)
        raise Exception(message)
    else:
        logging.debug(case_data)
    
    fileOUT = write_SCM_case_file(case_nml, case_data, use_area)
    logging.debug("Created {}".format(fileOUT))
    
    write_SCM_nml_file(case_nml)
    #logging.debug("Created {}".format(nmlOUT))
    
if __name__ == '__main__':
    main()
