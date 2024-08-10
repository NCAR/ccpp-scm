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
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

DEFAULT_MISSING_VALUE = -9999.0
DEFAULT_NUDGING_TIMESCALE = 7200.0 #s

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-n',   '--case_name',       help='name of case',         required=True)
parser.add_argument('-a',   '--use_area',        help='use column_area namelist attribute as forcing_scale',  action='store_true')
parser.add_argument('-d',   '--debug',           help='enable debugging output', action='store_true')


########################################################################################
#
########################################################################################
def parse_arguments():
    """Parse command line arguments"""
    args           = parser.parse_args()
        
    return (args.case_name, args.use_area, args.debug)

########################################################################################
#
########################################################################################
def setup_logging(debug):
    """Sets up the logging module."""
    if debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

class Case_Data(object):
    def __init__(self, name, missing_value, time, levels, soil_depth, lat, lon, slmsk, vegsrc, vegtyp, soiltyp, \
                 scolor, slopetyp, tsfco, vegfrac, shdmin, shdmax, canopy, hice, fice, tisfc, snowd, snoalb, tg3, \
                 uustar, alvsf, alnsf, alvwf, alnwf, facsf, facwf, weasd, f10m, t2m, q2m, ffmm, ffhh, \
                 tprcp, srflag, sncovr, tsfcl, zorl, zorll, zorli, zorlw, zorlwav, tvxy, tgxy, tahxy, canicexy, canliqxy, eahxy, \
                 cmxy, chxy, fwetxy, sneqvoxy, alboldxy, qsnowxy, wslakexy, taussxy, waxy, wtxy, zwtxy, xlaixy, xsaixy, \
                 lfmassxy, stmassxy, rtmassxy, woodxy, stblcpxy, fastcpxy, smcwtdxy, deeprechxy, rechxy, snowxy, \
                 wetness, clw_surf_land, clw_surf_ice, qwv_surf_land, qwv_surf_ice, tsnow_land, tsnow_ice, \
                 snowfallac_land, snowfallac_ice, sncovr_ice, sfalb_lnd, sfalb_lnd_bck, emis_ice, lai, area, \
                 stddev, convexity, oa1, oa2, oa3, oa4, ol1, ol2, ol3, ol4, theta_oro, gamma, sigma, elvmax, \
                 oro, oro_uf, landfrac, lakefrac, lakedepth, tref, z_c, c_0, c_d, w_0, w_d, xt, xs, xu, xv, xz, zm, \
                 xtts, xzts, d_conv, ifd, dt_cool, qrain, height, theta_il, t,\
                 qt, ql, qi, u, v, tke, ozone, stc, smc, slc, snicexy, snliqxy, tsnoxy, smoiseq, zsnsoxy, tiice, \
                 tslb, smois, sh2o, smfr, flfr, \
                 p_surf, T_surf, w_ls, omega, u_g, v_g, u_nudge, v_nudge, T_nudge, thil_nudge, qt_nudge, \
                 dT_dt_rad, h_advec_thil, v_advec_thil, h_advec_qt, v_advec_qt, sh_flux_sfc, lh_flux_sfc):
        self._name = name
        self._missing_value = missing_value
        #dimensions
        self._time = time
        self._levels = levels
        self._soil_depth = soil_depth
        #scalars
        self._lat = lat
        self._lon = lon
        #NOAH LSM scalars
        self._slmsk = slmsk
        self._vegsrc = vegsrc
        self._vegtyp = vegtyp
        self._soiltyp = soiltyp
        self._scolor = scolor
        self._slopetyp = slopetyp
        self._tsfco = tsfco
        self._vegfrac = vegfrac
        self._shdmin = shdmin
        self._shdmax = shdmax
        self._canopy = canopy
        self._hice = hice
        self._fice = fice
        self._tisfc = tisfc
        self._snowd = snowd
        self._snoalb = snoalb
        self._tg3 = tg3
        self._uustar = uustar
        self._alvsf = alvsf
        self._alnsf = alnsf
        self._alvwf = alvwf
        self._alnwf = alnwf
        self._facsf = facsf
        self._facwf = facwf
        self._weasd = weasd
        self._f10m = f10m
        self._t2m = t2m
        self._q2m = q2m
        self._ffmm = ffmm
        self._ffhh = ffhh
        self._tprcp = tprcp
        self._srflag = srflag
        self._sncovr = sncovr
        self._tsfcl = tsfcl
        self._zorl = zorl
        self._zorll = zorll
        self._zorli = zorli
        self._zorlw = zorlw
        self._zorlwav = zorlwav
        self._area = area
        #NOAHMP scalars
        self._tvxy = tvxy
        self._tgxy = tgxy
        self._tahxy = tahxy
        self._canicexy = canicexy
        self._canliqxy = canliqxy
        self._eahxy = eahxy
        self._cmxy = cmxy
        self._chxy = chxy
        self._fwetxy = fwetxy
        self._sneqvoxy = sneqvoxy
        self._alboldxy = alboldxy
        self._qsnowxy = qsnowxy
        self._wslakexy = wslakexy
        self._taussxy = taussxy
        self._waxy = waxy
        self._wtxy = wtxy
        self._zwtxy = zwtxy
        self._xlaixy = xlaixy
        self._xsaixy = xsaixy
        self._lfmassxy = lfmassxy
        self._stmassxy = stmassxy
        self._rtmassxy = rtmassxy
        self._woodxy = woodxy
        self._stblcpxy = stblcpxy
        self._fastcpxy = fastcpxy
        self._smcwtdxy = smcwtdxy
        self._deeprechxy = deeprechxy
        self._rechxy = rechxy
        self._snowxy = snowxy
        #RUC LSM scalars
        self._wetness = wetness
        self._clw_surf_land = clw_surf_land
        self._clw_surf_ice = clw_surf_ice
        self._qwv_surf_land = qwv_surf_land
        self._qwv_surf_ice = qwv_surf_ice
        self._tsnow_land = tsnow_land
        self._tsnow_ice = tsnow_ice
        self._snowfallac_land = snowfallac_land
        self._snowfallac_ice = snowfallac_ice
        self._sncovr_ice = sncovr_ice
        self._sfalb_lnd = sfalb_lnd
        self._sfalb_lnd_bck = sfalb_lnd_bck
        self._emis_ice = emis_ice
        self._lai = lai
        #orographic variables
        self._stddev = stddev
        self._convexity = convexity
        self._oa1 = oa1
        self._oa2 = oa2
        self._oa3 = oa3
        self._oa4 = oa4
        self._ol1 = ol1
        self._ol2 = ol2
        self._ol3 = ol3
        self._ol4 = ol4
        self._theta_oro = theta_oro
        self._gamma = gamma
        self._sigma = sigma
        self._elvmax = elvmax
        self._oro = oro
        self._oro_uf = oro_uf
        self._landfrac = landfrac
        self._lakefrac = lakefrac
        self._lakedepth = lakedepth
        #NSST vars
        self._tref = tref
        self._z_c = z_c
        self._c_0 = c_0
        self._c_d = c_d
        self._w_0 = w_0
        self._w_d = w_d
        self._xt = xt
        self._xs = xs
        self._xu = xu
        self._xv = xv
        self._xz = xz
        self._zm = zm
        self._xtts = xtts
        self._xzts = xzts
        self._d_conv = d_conv
        self._ifd = ifd
        self._dt_cool = dt_cool
        self._qrain = qrain
        #initial conditions (profiles for t0)
        self._height = height
        self._theta_il = theta_il
        self._t = t
        self._qt = qt
        self._ql = ql
        self._qi = qi
        self._u = u
        self._v = v
        self._tke = tke
        self._ozone = ozone
        self._stc = stc
        self._smc = smc
        self._slc = slc
        self._snicexy = snicexy
        self._snliqxy = snliqxy
        self._tsnoxy = tsnoxy
        self._smoiseq = smoiseq
        self._zsnsoxy = zsnsoxy
        self._tiice = tiice
        self._tslb = tslb
        self._smois = smois
        self._sh2o = sh2o
        self._smfr = smfr
        self._flfr = flfr
        #time series
        self._p_surf = p_surf
        self._T_surf = T_surf
        #2D forcing vars (vert, time)
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
        return f"""Case_Data
        Globals:
         name: {self._name}
         missing_value: {self._missing_value}
        Dimensions:
         time: {self._time} 
         levels: {self._levels} 
         soil_depth: {self._soil_depth} 
        Scalars:
         lat: {self._lat} 
         lon: {self._lon} 
         area: {self._area} 
         slmsk: {self._slmsk} 
         vegsrc: {self._vegsrc} 
         vegtyp: {self._vegtyp} 
         soiltyp: {self._soiltyp} 
         scolor: {self._scolor} 
         slopetyp: {self._slopetyp} 
         tsfco: {self._tsfco} 
         vegfrac: {self._vegfrac} 
         shdmin: {self._shdmin} 
         shdmax: {self._shdmax} 
         canopy: {self._canopy} 
         hice: {self._hice} 
         fice: {self._fice} 
         tisfc: {self._tisfc} 
         snowd: {self._snowd} 
         snoalb: {self._snoalb} 
         tg3: {self._tg3} 
         uustar: {self._uustar} 
         alvsf: {self._alvsf} 
         alnsf: {self._alnsf} 
         alvwf: {self._alvwf} 
         alnwf: {self._alnwf} 
         facsf: {self._facsf} 
         facwf: {self._facwf} 
         weasd: {self._weasd} 
         f10m: {self._f10m} 
         t2m: {self._t2m} 
         q2m: {self._q2m} 
         ffmm: {self._ffmm} 
         ffhh: {self._ffhh} 
         tprcp: {self._tprcp} 
         srflag: {self._srflag} 
         sncovr: {self._sncovr} 
         tsfcl: {self._tsfcl} 
         zorl: {self._zorl} 
         zorll: {self._zorll} 
         zorli: {self._zorli} 
         zorlw: {self._zorlw} 
         zorlwav: {self._zorlwav} 
         tvxy: {self._tvxy} 
         tgxy: {self._tgxy} 
         tahxy: {self._tahxy} 
         canicexy: {self._canicexy} 
         canliqxy: {self._canliqxy} 
         eahxy: {self._eahxy} 
         cmxy: {self._cmxy} 
         chxy: {self._chxy} 
         fwetxy: {self._fwetxy} 
         sneqvoxy: {self._sneqvoxy} 
         alboldxy: {self._alboldxy} 
         qsnowxy: {self._qsnowxy} 
         wslakexy: {self._wslakexy} 
         taussxy: {self._taussxy} 
         waxy: {self._waxy} 
         wtxy: {self._wtxy} 
         zwtxy: {self._zwtxy} 
         xlaixy: {self._xlaixy} 
         xsaixy: {self._xsaixy} 
         lfmassxy: {self._lfmassxy} 
         stmassxy: {self._stmassxy} 
         rtmassxy: {self._rtmassxy} 
         woodxy: {self._woodxy} 
         stblcpxy: {self._stblcpxy} 
         fastcpxy: {self._fastcpxy} 
         smcwtdxy: {self._smcwtdxy} 
         deeprechxy: {self._deeprechxy} 
         rechxy: {self._rechxy} 
         snowxy: {self._snowxy} 
         wetness: {self._wetness} 
         clw_surf_land: {self._clw_surf_land} 
         clw_surf_ice: {self._clw_surf_ice} 
         qwv_surf_land: {self._qwv_surf_land} 
         qwv_surf_ice: {self._qwv_surf_ice} 
         tsnow_land: {self._tsnow_land} 
         tsnow_ice: {self._tsnow_ice} 
         snowfallac_land: {self._snowfallac_land} 
         snowfallac_ice: {self._snowfallac_ice} 
         sncovr_ice: {self._sncovr_ice} 
         sfalb_lnd: {self._sfalb_lnd} 
         sfalb_lnd_bck: {self._sfalb_lnd_bck} 
         emis_ice: {self._emis_ice} 
         lai: {self._lai} 
         stddev: {self._stddev} 
         convexity: {self._convexity} 
         oa1: {self._oa1} 
         oa2: {self._oa2} 
         oa3: {self._oa3} 
         oa4: {self._oa4} 
         ol1: {self._ol1} 
         ol2: {self._ol2} 
         ol3: {self._ol3} 
         ol4: {self._ol4} 
         theta_oro: {self._theta_oro} 
         gamma: {self._gamma} 
         sigma: {self._sigma} 
         elvmax: {self._elvmax} 
         oro: {self._oro} 
         oro_uf: {self._oro_uf} 
         landfrac: {self._landfrac} 
         lakefrac: {self._lakefrac} 
         lakedepth: {self._lakedepth} 
         tref: {self._tref} 
         z_c: {self._z_c} 
         c_0: {self._c_0} 
         c_d: {self._c_d} 
         w_0: {self._w_0} 
         w_d: {self._w_d} 
         xt: {self._xt} 
         xs: {self._xs} 
         xu: {self._xu} 
         xv: {self._xv} 
         xz: {self._xz} 
         zm: {self._zm} 
         xtts: {self._xtts} 
         xzts: {self._xzts} 
         d_conv: {self._d_conv} 
         ifd: {self._ifd} 
         dt_cool: {self._dt_cool} 
         qrain: {self._qrain} 
        Initial: 
         height: {self._height} 
         theta_il: {self._theta_il} 
         t: {self._t} 
         qt: {self._qt} 
         ql: {self._ql} 
         qi: {self._qi} 
         u: {self._u} 
         v: {self._v} 
         tke: {self._tke} 
         ozone: {self._ozone} 
         stc: {self._stc} 
         smc: {self._smc} 
         slc: {self._slc} 
         snicexy: {self._snicexy} 
         snliqxy: {self._snliqxy} 
         tsnoxy: {self._tsnoxy} 
         smoiseq: {self._smoiseq} 
         zsnsoxy: {self._zsnsoxy} 
         tiice: {self._tiice} 
         tslb: {self._tslb} 
         smois: {self._smois} 
         sh2o: {self._sh2o} 
         smfr: {self._smfr} 
         flfr: {self._flfr} 
        Forcing: 
         p_surf: {self._p_surf} 
         T_surf: {self._T_surf} 
         w_ls (time avg): {np.mean(self._w_ls, axis=1)} 
         omega (time avg): {np.mean(self._omega, axis=1)} 
         u_g (time avg): {np.mean(self._u_g, axis=1)} 
         v_g (time avg): {np.mean(self._v_g, axis=1)} 
         u_nudge (time avg): {np.mean(self._u_nudge, axis=1)} 
         v_nudge (time avg): {np.mean(self._v_nudge, axis=1)} 
         T_nudge (time avg): {np.mean(self._T_nudge, axis=1)} 
         thil_nudge (time avg): {np.mean(self._thil_nudge, axis=1)} 
         qt_nudge (time avg): {np.mean(self._qt_nudge, axis=1)}
         dT_dt_rad (time avg): {np.mean(self._dT_dt_rad, axis=1)}
         h_advec_thil (time avg): {np.mean(self._h_advec_thil, axis=1)}
         v_advec_thil (time avg): {np.mean(self._v_advec_thil, axis=1)}
         h_advec_qt (time avg): {np.mean(self._h_advec_qt, axis=1)}
         v_advec_qt (time avg): {np.mean(self._v_advec_qt, axis=1)}
         sh_flux_sfc: {self._sh_flux_sfc}
         lh_flux_sfc: {self._lh_flux_sfc}"""

def get_case_nml(case_name):
    """Returns case configuration Fortran namelist"""
    
    filename = os.path.join(CASE_NML_DIR, case_name + '.nml')
    
    logging.debug(filename)
    
    error = False
    nml = ''
    if (os.path.exists(filename)):
        nml = f90nml.read(filename)
    else:
        error = True
    
    return (nml, error)

def get_case_data(case_name):
    """Returns proprietary CCPP SCM case data in NetCDF Dataset format"""
    
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
        try:
            soil_depth = nc_fid.variables['soil_depth'][:]
        except KeyError:
            soil_depth = [missing_value]
        
        #read variables from scalar group
        scalars_grp = nc_fid.groups['scalars']
        lat = scalars_grp.variables['lat'][:]
        lon = scalars_grp.variables['lon'][:]
        try:
            slmsk = scalars_grp.variables['slmsk'][:]
        except KeyError:
            slmsk = missing_value
        try:
            vegsrc = scalars_grp.variables['vegsrc'][:]
        except KeyError:
            vegsrc = missing_value
        try:
            vegtyp = scalars_grp.variables['vegtyp'][:]
        except KeyError:
            vegtyp = missing_value
        try:
            soiltyp = scalars_grp.variables['soiltyp'][:]
        except KeyError:
            soiltyp = missing_value
        try:
            scolor = scalars_grp.variables['scolor'][:]
        except KeyError:
            scolor = missing_value
        try:
            slopetyp = scalars_grp.variables['slopetyp'][:]
        except KeyError:
            slopetyp = missing_value
        try:
            tsfco = scalars_grp.variables['tsfco'][:]
        except KeyError:
            tsfco = missing_value
        try:
            vegfrac = scalars_grp.variables['vegfrac'][:]
        except KeyError:
            vegfrac = missing_value
        try:
            shdmin = scalars_grp.variables['shdmin'][:]
        except KeyError:
            shdmin = missing_value
        try:
            shdmax = scalars_grp.variables['shdmax'][:]
        except KeyError:
            shdmax = missing_value
        try:
            canopy = scalars_grp.variables['canopy'][:]
        except KeyError:
            canopy = missing_value
        try:
            hice = scalars_grp.variables['hice'][:]
        except KeyError:
            hice = missing_value
        try:
            fice = scalars_grp.variables['fice'][:]
        except KeyError:
            fice = missing_value
        try:
            tisfc = scalars_grp.variables['tisfc'][:]
        except KeyError:
            tisfc = missing_value
        try:
            snowd = scalars_grp.variables['snowd'][:]
        except KeyError:
            try:
                snowd = scalars_grp.variables['snwdph'][:]
            except KeyError:
                snowd = missing_value
        try:
            snoalb = scalars_grp.variables['snoalb'][:]
        except KeyError:
            snoalb = missing_value
        try:
            tg3 = scalars_grp.variables['tg3'][:]
        except KeyError:
            tg3 = missing_value
        try:
            uustar = scalars_grp.variables['uustar'][:]
        except KeyError:
            uustar = missing_value
        try:
            alvsf = scalars_grp.variables['alvsf'][:]
        except KeyError:
            alvsf = missing_value
        try:
            alnsf = scalars_grp.variables['alnsf'][:]
        except KeyError:
            alnsf = missing_value
        try:
            alvwf = scalars_grp.variables['alvwf'][:]
        except KeyError:
            alvwf = missing_value
        try:
            alnwf = scalars_grp.variables['alnwf'][:]
        except KeyError:
            alnwf = missing_value
        try:
            facsf = scalars_grp.variables['facsf'][:]
        except KeyError:
            facsf = missing_value
        try:
            facwf = scalars_grp.variables['facwf'][:]
        except KeyError:
            facwf = missing_value
        try:
            weasd = scalars_grp.variables['weasd'][:]
        except KeyError:
            weasd = missing_value
        try:
            f10m = scalars_grp.variables['f10m'][:]
        except KeyError:
            f10m = missing_value
        try:
            t2m = scalars_grp.variables['t2m'][:]
        except KeyError:
            t2m = missing_value
        try:
            q2m = scalars_grp.variables['q2m'][:]
        except KeyError:
            q2m = missing_value
        try:
            ffmm = scalars_grp.variables['ffmm'][:]
        except KeyError:
            ffmm = missing_value
        try:
            ffhh = scalars_grp.variables['ffhh'][:]
        except KeyError:
            ffhh = missing_value
        try:
            tprcp = scalars_grp.variables['tprcp'][:]
        except KeyError:
            tprcp = missing_value
        try:
            srflag = scalars_grp.variables['srflag'][:]
        except KeyError:
            srflag = missing_value
        try:
            sncovr = scalars_grp.variables['sncovr'][:]
        except KeyError:
            sncovr = missing_value
        try:
            tsfcl = scalars_grp.variables['tsfcl'][:]
        except KeyError:
            tsfcl = missing_value
        try:
            zorl = scalars_grp.variables['zorl'][:]
        except KeyError:
            zorl = missing_value
        try:
            zorll = scalars_grp.variables['zorll'][:]
        except KeyError:
            zorll = missing_value
        try:
            zorli = scalars_grp.variables['zorli'][:]
        except KeyError:
            zorli = missing_value
        try:
            zorlw = scalars_grp.variables['zorlw'][:]
        except KeyError:
            zorlw = missing_value
        try:
            zorlwav = scalars_grp.variables['zorlwav'][:]
        except KeyError:
            zorlwav = missing_value
        try:
            tvxy = scalars_grp.variables['tvxy'][:]
        except KeyError:
            tvxy = missing_value
        try:
            tgxy = scalars_grp.variables['tgxy'][:]
        except KeyError:
            tgxy = missing_value
        try:
            tahxy = scalars_grp.variables['tahxy'][:]
        except KeyError:
            tahxy = missing_value
        try:
            canicexy = scalars_grp.variables['canicexy'][:]
        except KeyError:
            canicexy = missing_value
        try:
            canliqxy = scalars_grp.variables['canliqxy'][:]
        except KeyError:
            canliqxy = missing_value
        try:
            eahxy = scalars_grp.variables['eahxy'][:]
        except KeyError:
            eahxy = missing_value
        try:
            cmxy = scalars_grp.variables['cmxy'][:]
        except KeyError:
            cmxy = missing_value
        try:
            chxy = scalars_grp.variables['chxy'][:]
        except KeyError:
            chxy = missing_value
        try:
            fwetxy = scalars_grp.variables['fwetxy'][:]
        except KeyError:
            fwetxy = missing_value
        try:
            sneqvoxy = scalars_grp.variables['sneqvoxy'][:]
        except KeyError:
            sneqvoxy = missing_value
        try:
            alboldxy = scalars_grp.variables['alboldxy'][:]
        except KeyError:
            alboldxy = missing_value
        try:
            qsnowxy = scalars_grp.variables['qsnowxy'][:]
        except KeyError:
            qsnowxy = missing_value
        try:
            wslakexy = scalars_grp.variables['wslakexy'][:]
        except KeyError:
            wslakexy = missing_value
        try:
            taussxy = scalars_grp.variables['taussxy'][:]
        except KeyError:
            taussxy = missing_value
        try:
            waxy = scalars_grp.variables['waxy'][:]
        except KeyError:
            waxy = missing_value
        try:
            wtxy = scalars_grp.variables['wtxy'][:]
        except KeyError:
            wtxy = missing_value
        try:
            zwtxy = scalars_grp.variables['zwtxy'][:]
        except KeyError:
            zwtxy = missing_value
        try:
            xlaixy = scalars_grp.variables['xlaixy'][:]
        except KeyError:
            xlaixy = missing_value
        try:
            xsaixy = scalars_grp.variables['xsaixy'][:]
        except KeyError:
            xsaixy = missing_value
        try:
            lfmassxy = scalars_grp.variables['lfmassxy'][:]
        except KeyError:
            lfmassxy = missing_value
        try:
            stmassxy = scalars_grp.variables['stmassxy'][:]
        except KeyError:
            stmassxy = missing_value
        try:
            rtmassxy = scalars_grp.variables['rtmassxy'][:]
        except KeyError:
            rtmassxy = missing_value
        try:
            woodxy = scalars_grp.variables['woodxy'][:]
        except KeyError:
            woodxy = missing_value
        try:
            stblcpxy = scalars_grp.variables['stblcpxy'][:]
        except KeyError:
            stblcpxy = missing_value
        try:
            fastcpxy = scalars_grp.variables['fastcpxy'][:]
        except KeyError:
            fastcpxy = missing_value
        try:
            smcwtdxy = scalars_grp.variables['smcwtdxy'][:]
        except KeyError:
            smcwtdxy = missing_value
        try:
            deeprechxy = scalars_grp.variables['deeprechxy'][:]
        except KeyError:
            deeprechxy = missing_value
        try:
            rechxy = scalars_grp.variables['rechxy'][:]
        except KeyError:
            rechxy = missing_value
        try:
            snowxy = scalars_grp.variables['snowxy'][:]
        except KeyError:
            snowxy = missing_value
        try:
            wetness = scalars_grp.variables['wetness'][:]
        except KeyError:
            wetness = missing_value
        try:
            clw_surf_land = scalars_grp.variables['clw_surf_land'][:]
        except KeyError:
            clw_surf_land = missing_value
        try:
            clw_surf_ice = scalars_grp.variables['clw_surf_ice'][:]
        except KeyError:
            clw_surf_ice = missing_value
        try:
            qwv_surf_land = scalars_grp.variables['qwv_surf_land'][:]
        except KeyError:
            qwv_surf_land = missing_value
        try:
            qwv_surf_ice = scalars_grp.variables['qwv_surf_ice'][:]
        except KeyError:
            qwv_surf_ice = missing_value
        try:
            tsnow_land = scalars_grp.variables['tsnow_land'][:]
        except KeyError:
            tsnow_land = missing_value
        try:
            tsnow_ice = scalars_grp.variables['tsnow_ice'][:]
        except KeyError:
            tsnow_ice = missing_value
        try:
            snowfallac_land = scalars_grp.variables['snowfallac_land'][:]
        except KeyError:
            snowfallac_land = missing_value
        try:
            snowfallac_ice = scalars_grp.variables['snowfallac_ice'][:]
        except KeyError:
            snowfallac_ice = missing_value
        try:
            sncovr_ice = scalars_grp.variables['sncovr_ice'][:]
        except KeyError:
            sncovr_ice = missing_value
        try:
            sfalb_lnd = scalars_grp.variables['sfalb_lnd'][:]
        except KeyError:
            sfalb_lnd = missing_value
        try:
            sfalb_lnd_bck = scalars_grp.variables['sfalb_lnd_bck'][:]
        except KeyError:
            sfalb_lnd_bck = missing_value
        try:
            emis_ice = scalars_grp.variables['emis_ice'][:]
        except KeyError:
            emis_ice = missing_value
        try:
            lai = scalars_grp.variables['lai'][:]
        except KeyError:
            lai = missing_value
        try:
            area = scalars_grp.variables['area'][:]
        except KeyError:
            area = missing_value
        try:
            stddev = scalars_grp.variables['stddev'][:]
        except KeyError:
            stddev = missing_value
        try:
            convexity = scalars_grp.variables['convexity'][:]
        except KeyError:
            convexity = missing_value
        try:
            oa1 = scalars_grp.variables['oa1'][:]
        except KeyError:
            oa1 = missing_value
        try:
            oa2 = scalars_grp.variables['oa2'][:]
        except KeyError:
            oa2 = missing_value
        try:
            oa3 = scalars_grp.variables['oa3'][:]
        except KeyError:
            oa3 = missing_value
        try:
            oa4 = scalars_grp.variables['oa4'][:]
        except KeyError:
            oa4 = missing_value
        try:
            ol1 = scalars_grp.variables['ol1'][:]
        except KeyError:
            ol1 = missing_value
        try:
            ol2 = scalars_grp.variables['ol2'][:]
        except KeyError:
            ol2 = missing_value
        try:
            ol3 = scalars_grp.variables['ol3'][:]
        except KeyError:
            ol3 = missing_value
        try:
            ol4 = scalars_grp.variables['ol4'][:]
        except KeyError:
            ol4 = missing_value
        try:
            theta_oro = scalars_grp.variables['theta_oro'][:]
        except KeyError:
            theta_oro = missing_value
        try:
            gamma = scalars_grp.variables['gamma'][:]
        except KeyError:
            gamma = missing_value
        try:
            sigma = scalars_grp.variables['sigma'][:]
        except KeyError:
            sigma = missing_value
        try:
            elvmax = scalars_grp.variables['elvmax'][:]
        except KeyError:
            elvmax = missing_value
        try:
            oro = scalars_grp.variables['oro'][:]
        except KeyError:
            oro = missing_value
        try:
            oro_uf = scalars_grp.variables['oro_uf'][:]
        except KeyError:
            oro_uf = missing_value
        try:
            landfrac = scalars_grp.variables['landfrac'][:]
        except KeyError:
            landfrac = missing_value
        try:
            lakefrac = scalars_grp.variables['lakefrac'][:]
        except KeyError:
            lakefrac = missing_value
        try:
            lakedepth = scalars_grp.variables['lakedepth'][:]
        except KeyError:
            lakedepth = missing_value
        try:
            tref = scalars_grp.variables['tref'][:]
        except KeyError:
            tref = missing_value
        try:
            z_c = scalars_grp.variables['z_c'][:]
        except KeyError:
            z_c = missing_value
        try:
            c_0 = scalars_grp.variables['c_0'][:]
        except KeyError:
            c_0 = missing_value
        try:
            c_d = scalars_grp.variables['c_d'][:]
        except KeyError:
            c_d = missing_value
        try:
            w_0 = scalars_grp.variables['w_0'][:]
        except KeyError:
            w_0 = missing_value
        try:
            w_d = scalars_grp.variables['w_d'][:]
        except KeyError:
            w_d = missing_value
        try:
            xt = scalars_grp.variables['xt'][:]
        except KeyError:
            xt = missing_value
        try:
            xs = scalars_grp.variables['xs'][:]
        except KeyError:
            xs = missing_value
        try:
            xu = scalars_grp.variables['xu'][:]
        except KeyError:
            xu = missing_value
        try:
            xv = scalars_grp.variables['xv'][:]
        except KeyError:
            xv = missing_value
        try:
            xz = scalars_grp.variables['xz'][:]
        except KeyError:
            xz = missing_value
        try:
            zm = scalars_grp.variables['zm'][:]
        except KeyError:
            zm = missing_value
        try:
            xtts = scalars_grp.variables['xtts'][:]
        except KeyError:
            xtts = missing_value
        try:
            xzts = scalars_grp.variables['xzts'][:]
        except KeyError:
            xzts = missing_value
        try:
            d_conv = scalars_grp.variables['d_conv'][:]
        except KeyError:
            d_conv = missing_value
        try:
            ifd = scalars_grp.variables['ifd'][:]
        except KeyError:
            ifd = missing_value
        try:
            dt_cool = scalars_grp.variables['dt_cool'][:]
        except KeyError:
            dt_cool = missing_value
        try:
            qrain = scalars_grp.variables['qrains'][:]
        except KeyError:
            qrain = missing_value        
          
        #read variables from initial group
        initial_grp = nc_fid.groups['initial']
        try:
            height = initial_grp.variables['height'][:]
        except KeyError:
            height = [missing_value]
        try:
            theta_il = initial_grp.variables['thetail'][:]
        except KeyError:
            theta_il = [missing_value]
        try:
            t = initial_grp.variables['temp'][:]
        except KeyError:
            t = [missing_value]
        qt = initial_grp.variables['qt'][:]
        ql = initial_grp.variables['ql'][:]
        qi = initial_grp.variables['qi'][:]
        u = initial_grp.variables['u'][:]
        v = initial_grp.variables['v'][:]
        tke = initial_grp.variables['tke'][:]
        ozone = initial_grp.variables['ozone'][:]
        try:
            stc = initial_grp.variables['stc'][:]
        except KeyError:
            stc = [missing_value]
        try:
            smc = initial_grp.variables['smc'][:]
        except KeyError:
            smc = [missing_value]
        try:
            slc = initial_grp.variables['slc'][:]
        except KeyError:
            slc = [missing_value]
        try:
            snicexy = initial_grp.variables['snicexy'][:]
        except KeyError:
            snicexy = [missing_value]
        try:
            snliqxy = initial_grp.variables['snliqxy'][:]
        except KeyError:
            snliqxy = [missing_value]
        try:
            tsnoxy = initial_grp.variables['tsnoxy'][:]
        except KeyError:
            tsnoxy = [missing_value]
        try:
            smoiseq = initial_grp.variables['smoiseq'][:]
        except KeyError:
            smoiseq = [missing_value]
        try:
            zsnsoxy = initial_grp.variables['zsnsoxy'][:]
        except KeyError:
            zsnsoxy = [missing_value]
        try:
            tiice = initial_grp.variables['tiice'][:]
        except KeyError:
            tiice = [missing_value]
        try:
            tslb = initial_grp.variables['tslb'][:]
        except KeyError:
            tslb = [missing_value]
        try:
            smois = initial_grp.variables['smois'][:]
        except KeyError:
            smois = [missing_value]
        try:
            sh2o = initial_grp.variables['sh2o'][:]
        except KeyError:
            sh2o = [missing_value]
        try:
            smfr = initial_grp.variables['smfr'][:]
        except KeyError:
            smfr = [missing_value]
        try:
            flfr = initial_grp.variables['flfr'][:]
        except KeyError:
            flfr = [missing_value]
        
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
    
    case_data = Case_Data(case_name, missing_value, time, levels, soil_depth, lat, lon, 
                          slmsk, vegsrc, vegtyp, soiltyp, scolor, slopetyp, tsfco, vegfrac, shdmin, shdmax,
                          canopy, hice, fice, tisfc, snowd, snoalb, tg3, uustar, 
                          alvsf, alnsf, alvwf, alnwf, facsf, facwf, weasd, f10m, t2m, q2m, ffmm, ffhh, 
                          tprcp, srflag, sncovr, tsfcl, zorl, zorll, zorli, zorlw, zorlwav, tvxy, tgxy, tahxy, canicexy, canliqxy, eahxy,
                          cmxy, chxy, fwetxy, sneqvoxy, alboldxy, qsnowxy, wslakexy, taussxy,
                          waxy, wtxy, zwtxy, xlaixy, xsaixy,
                          lfmassxy, stmassxy, rtmassxy, woodxy, stblcpxy, fastcpxy, smcwtdxy, deeprechxy, rechxy, snowxy,
                          wetness, clw_surf_land, clw_surf_ice, qwv_surf_land, qwv_surf_ice, tsnow_land, tsnow_ice,
                          snowfallac_land, snowfallac_ice, sncovr_ice, sfalb_lnd, sfalb_lnd_bck, emis_ice, lai, area, 
                          stddev, convexity, oa1, oa2, oa3, oa4, ol1, ol3, ol3, ol4, theta_oro, gamma, sigma, elvmax, oro, oro_uf, landfrac, lakefrac, lakedepth,
                          tref, z_c, c_0, c_d, w_0, w_d, xt, xs, xu, xv, xz, zm, xtts, xzts, d_conv, ifd, dt_cool, qrain, 
                          height, theta_il, t, qt, ql, qi, u, v, tke, ozone, stc, smc, slc, 
                          snicexy, snliqxy, tsnoxy, smoiseq, zsnsoxy, tiice, tslb, smois, sh2o, smfr, flfr, 
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
    
    nml_keys = case_nml['case_config'].todict().keys()
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
    
    if (case_nml['case_config']['sfc_type'] > 1.5):
        surface_string = 'ice'
    elif (case_nml['case_config']['sfc_type'] > 0.5):
        surface_string = 'land'
    else:
        surface_string = 'ocean'
    
    #override case nml with LSM/model data
    if ('lsm_ics' in nml_keys):
        if (case_nml['case_config']['lsm_ics']):
            if (case_data._slmsk > 1.5):
                surface_string = 'ice'
            elif (case_data._slmsk > 0.5):
                surface_string = 'land'
            else:
                surface_string = 'ocean'
    
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
    nc_file.surface_forcing_lsm      = 'none'
    
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
    
    if ('lsm_ics' in nml_keys or 'model_ics' in nml_keys):
        nc_file.surface_forcing_lsm      = 'lsm'
        if (('lsm_ics' in nml_keys and case_nml['case_config']['lsm_ics']) or ('model_ics' in nml_keys and case_nml['case_config']['model_ics'])):
            if (case_data._soil_depth[0] != case_data._missing_value):
                soil_dim   = nc_file.createDimension('nsoil', case_data._soil_depth.shape[0])
            else:
                message = 'LSM ICs are expected from the case_nml file, but no soil depth is provided.'
                logging.critical(message)
                raise Exception(message)
            
            if (case_data._snicexy[0] != case_data._missing_value):
                snow_dim   = nc_file.createDimension('nsnow', case_data._snicexy.shape[0])
                nslsnw_dim = nc_file.createDimension('nsoil_plus_nsnow',case_data._snicexy.shape[0] + case_data._soil_depth.shape[0])
            
            if (case_data._tiice[0] != case_data._missing_value):
                ice_dim    = nc_file.createDimension('nice',  case_data._tiice.shape[0])
    
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
    
    if (case_data._soil_depth[0] != case_data._missing_value):
        soil_depth_var               = nc_file.createVariable('soil_depth', wp, ('nsoil'))
        soil_depth_var.units         = 'm'
        soil_depth_var.standard_name = 'depth of bottom of soil layers'
        soil_depth_var[:]            = case_data._soil_depth[:]
    
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
    
    if (case_data._slmsk != case_data._missing_value):
        slmsk_var               = nc_file.createVariable('slmsk', wp)
        slmsk_var.units         = 'none'
        slmsk_var.standard_name = 'land_sea_ice_mask'
        slmsk_var[:]            = case_data._slmsk
    
    if (case_data._vegsrc != case_data._missing_value):
        vegsrc_var               = nc_file.createVariable('vegsrc', wp)
        vegsrc_var.units         = 'none'
        vegsrc_var.standard_name = 'vegetation source (1-2)'
        vegsrc_var[:]            = case_data._vegsrc
    
    if (case_data._vegtyp != case_data._missing_value):
        vegtyp_var               = nc_file.createVariable('vegtyp', wp)
        vegtyp_var.units         = 'none'
        vegtyp_var.standard_name = 'vegetation type (1-12)'
        vegtyp_var[:]            = case_data._vegtyp
    
    if (case_data._soiltyp != case_data._missing_value):
        soiltyp_var               = nc_file.createVariable('soiltyp', wp)
        soiltyp_var.units         = 'none'
        soiltyp_var.standard_name = 'soil type (1-12)'
        soiltyp_var[:]            = case_data._soiltyp
    
    if (case_data._scolor != case_data._missing_value):
        scolor_var               = nc_file.createVariable('scolor', wp)
        scolor_var.units         = 'none'
        scolor_var.standard_name = 'soil color'
        scolor_var[:]            = case_data._scolor
    
    if (case_data._slopetyp != case_data._missing_value):
        slopetyp_var               = nc_file.createVariable('slopetyp', wp)
        slopetyp_var.units         = 'none'
        slopetyp_var.standard_name = 'slope type (1-9)'
        slopetyp_var[:]            = case_data._slopetyp
    
    if (case_data._tsfco != case_data._missing_value):
        tsfco_var               = nc_file.createVariable('tsfco', wp)
        tsfco_var.units         = 'none'
        tsfco_var.standard_name = 'slope type (1-9)'
        tsfco_var[:]            = case_data._tsfco
    
    if (case_data._vegfrac != case_data._missing_value):
        vegfrac_var               = nc_file.createVariable('vegfrac', wp)
        vegfrac_var.units         = 'none'
        vegfrac_var.standard_name = 'slope type (1-9)'
        vegfrac_var[:]            = case_data._vegfrac
    
    if (case_data._shdmin != case_data._missing_value):
        shdmin_var               = nc_file.createVariable('shdmin', wp)
        shdmin_var.units         = 'none'
        shdmin_var.standard_name = 'slope type (1-9)'
        shdmin_var[:]            = case_data._shdmin
    
    if (case_data._shdmax != case_data._missing_value):
        shdmax_var               = nc_file.createVariable('shdmax', wp)
        shdmax_var.units         = 'none'
        shdmax_var.standard_name = 'slope type (1-9)'
        shdmax_var[:]            = case_data._shdmax
    
    if (case_data._canopy != case_data._missing_value):
        canopy_var               = nc_file.createVariable('canopy', wp)
        canopy_var.units         = 'kg m-2'
        canopy_var.standard_name = 'amount of water stored in canopy'
        canopy_var[:]            = case_data._canopy
    
    if (case_data._hice != case_data._missing_value):
        hice_var               = nc_file.createVariable('hice', wp)
        hice_var.units         = 'm'
        hice_var.standard_name = 'sea ice thickness'
        hice_var[:]            = case_data._hice
    
    if (case_data._fice != case_data._missing_value):
        fice_var               = nc_file.createVariable('fice', wp)
        fice_var.units         = 'none'
        fice_var.standard_name = 'ice fraction'
        fice_var[:]            = case_data._fice
    
    if (case_data._tisfc != case_data._missing_value):
        tisfc_var               = nc_file.createVariable('tisfc', wp)
        tisfc_var.units         = 'K'
        tisfc_var.standard_name = 'ice surface temperature'
        tisfc_var[:]            = case_data._tisfc
    
    if (case_data._snowd != case_data._missing_value):
        snowd_var               = nc_file.createVariable('snowd', wp)
        snowd_var.units         = 'mm'
        snowd_var.standard_name = 'water equivalent snow depth'
        snowd_var[:]            = case_data._snowd
    
    if (case_data._snoalb != case_data._missing_value):
        snoalb_var               = nc_file.createVariable('snoalb', wp)
        snoalb_var.units         = 'none'
        snoalb_var.standard_name = 'maximum snow albedo'
        snoalb_var[:]            = case_data._snoalb
    
    if (case_data._tg3 != case_data._missing_value):
        tg3_var               = nc_file.createVariable('tg3', wp)
        tg3_var.units         = 'K'
        tg3_var.standard_name = 'deep soil temperature'
        tg3_var[:]            = case_data._tg3
    
    if (case_data._uustar != case_data._missing_value):
        uustar_var               = nc_file.createVariable('uustar', wp)
        uustar_var.units         = 'm s-1'
        uustar_var.standard_name = 'surface_friction_velocity'
        uustar_var[:]            = case_data._uustar
    
    if (case_data._alvsf != case_data._missing_value):
        alvsf_var               = nc_file.createVariable('alvsf', wp)
        alvsf_var.units         = 'none'
        alvsf_var.standard_name = '60 degree vis albedo with strong cosz dependency'
        alvsf_var[:]            = case_data._alvsf
    
    if (case_data._alnsf != case_data._missing_value):
        alnsf_var               = nc_file.createVariable('alnsf', wp)
        alnsf_var.units         = 'none'
        alnsf_var.standard_name = '60 degree nir albedo with strong cosz dependency'
        alnsf_var[:]            = case_data._alnsf
    
    if (case_data._alvwf != case_data._missing_value):
        alvwf_var               = nc_file.createVariable('alvwf', wp)
        alvwf_var.units         = 'none'
        alvwf_var.standard_name = '60 degree vis albedo with weak cosz dependency'
        alvwf_var[:]            = case_data._alvwf
    
    if (case_data._alnwf != case_data._missing_value):
        alnwf_var               = nc_file.createVariable('alnwf', wp)
        alnwf_var.units         = 'none'
        alnwf_var.standard_name = '60 degree nir albedo with weak cosz dependency'
        alnwf_var[:]            = case_data._alnwf
    
    if (case_data._facsf != case_data._missing_value):
        facsf_var               = nc_file.createVariable('facsf', wp)
        facsf_var.units         = 'none'
        facsf_var.standard_name = 'fractional coverage with strong cosz dependency'
        facsf_var[:]            = case_data._facsf
    
    if (case_data._facwf != case_data._missing_value):
        facwf_var               = nc_file.createVariable('facwf', wp)
        facwf_var.units         = 'none'
        facwf_var.standard_name = 'fractional coverage with weak cosz dependency'
        facwf_var[:]            = case_data._facwf
    
    if (case_data._weasd != case_data._missing_value):
        weasd_var               = nc_file.createVariable('weasd', wp)
        weasd_var.units         = 'mm'
        weasd_var.standard_name = 'water equivalent accumulated snow depth'
        weasd_var[:]            = case_data._weasd
    
    if (case_data._f10m != case_data._missing_value):
        f10m_var               = nc_file.createVariable('f10m', wp)
        f10m_var.units         = 'none'
        f10m_var.standard_name = 'ratio of sigma level 1 wind and 10m wind'
        f10m_var[:]            = case_data._f10m
    
    if (case_data._t2m != case_data._missing_value):
        t2m_var               = nc_file.createVariable('t2m', wp)
        t2m_var.units         = 'K'
        t2m_var.standard_name = '2-meter absolute temperature'
        t2m_var[:]            = case_data._t2m
    
    if (case_data._q2m != case_data._missing_value):
        q2m_var               = nc_file.createVariable('q2m', wp)
        q2m_var.units         = 'kg kg-1'
        q2m_var.standard_name = '2-meter specific humidity'
        q2m_var[:]            = case_data._q2m
    
    if (case_data._ffmm != case_data._missing_value):
        ffmm_var               = nc_file.createVariable('ffmm', wp)
        ffmm_var.units         = 'none'
        ffmm_var.standard_name = 'Monin-Obukhov similarity function for momentum'
        ffmm_var[:]            = case_data._ffmm
    
    if (case_data._ffhh != case_data._missing_value):
        ffhh_var               = nc_file.createVariable('ffhh', wp)
        ffhh_var.units         = 'none'
        ffhh_var.standard_name = 'Monin-Obukhov similarity function for heat'
        ffhh_var[:]            = case_data._ffhh
    
    if (case_data._tprcp != case_data._missing_value):
        tprcp_var               = nc_file.createVariable('tprcp', wp)
        tprcp_var.units         = 'm'
        tprcp_var.standard_name = 'instantaneous total precipitation amount'
        tprcp_var[:]            = case_data._tprcp
    
    if (case_data._srflag != case_data._missing_value):
        srflag_var               = nc_file.createVariable('srflag', wp)
        srflag_var.units         = 'none'
        srflag_var.standard_name = 'snow/rain flag for precipitation'
        srflag_var[:]            = case_data._srflag
    
    if (case_data._sncovr != case_data._missing_value):
        sncovr_var               = nc_file.createVariable('sncovr', wp)
        sncovr_var.units         = 'none'
        sncovr_var.standard_name = 'surface snow area fraction'
        sncovr_var[:]            = case_data._sncovr
    
    if (case_data._tsfcl != case_data._missing_value):
        tsfcl_var               = nc_file.createVariable('tsfcl', wp)
        tsfcl_var.units         = 'K'
        tsfcl_var.standard_name = 'surface skin temperature over land'
        tsfcl_var[:]            = case_data._tsfcl
    
    if (case_data._zorl != case_data._missing_value):
        zorl_var               = nc_file.createVariable('zorl', wp)
        zorl_var.units         = 'cm'
        zorl_var.standard_name = 'surface roughness length'
        zorl_var[:]            = case_data._zorl
    
    if (case_data._zorll != case_data._missing_value):
        zorll_var               = nc_file.createVariable('zorll', wp)
        zorll_var.units         = 'cm'
        zorll_var.standard_name = 'surface roughness length over land'
        zorll_var[:]            = case_data._zorll
    
    if (case_data._zorli != case_data._missing_value):
        zorli_var               = nc_file.createVariable('zorli', wp)
        zorli_var.units         = 'cm'
        zorli_var.standard_name = 'surface roughness length over ice'
        zorli_var[:]            = case_data._zorli
    
    if (case_data._zorlw != case_data._missing_value):
        zorlw_var               = nc_file.createVariable('zorlw', wp)
        zorlw_var.units         = 'cm'
        zorlw_var.standard_name = 'surface roughness length over ocean'
        zorlw_var[:]            = case_data._zorlw
    
    if (case_data._zorlwav != case_data._missing_value):
        zorlwav_var               = nc_file.createVariable('zorlwav', wp)
        zorlwav_var.units         = 'cm'
        zorlwav_var.standard_name = 'surface_roughness_length_from_wave_model'
        zorlwav_var[:]            = case_data._zorlwav
    
    if (case_data._tvxy != case_data._missing_value):
        tvxy_var               = nc_file.createVariable('tvxy', wp)
        tvxy_var.units         = 'K'
        tvxy_var.standard_name = 'vegetation temperature for NoahMP'
        tvxy_var[:]            = case_data._tvxy
    
    if (case_data._tgxy != case_data._missing_value):
        tgxy_var               = nc_file.createVariable('tgxy', wp)
        tgxy_var.units         = 'K'
        tgxy_var.standard_name = 'ground temperature for NoahMP'
        tgxy_var[:]            = case_data._tgxy
    
    if (case_data._tahxy != case_data._missing_value):
        tahxy_var               = nc_file.createVariable('tahxy', wp)
        tahxy_var.units         = 'K'
        tahxy_var.standard_name = 'canopy air temperature for NoahMP'
        tahxy_var[:]            = case_data._tahxy
    
    if (case_data._canicexy != case_data._missing_value):
        canicexy_var               = nc_file.createVariable('canicexy', wp)
        canicexy_var.units         = 'mm'
        canicexy_var.standard_name = 'canopy intercepted ice mass for NoahMP'
        canicexy_var[:]            = case_data._canicexy
    
    if (case_data._canliqxy != case_data._missing_value):
        canliqxy_var               = nc_file.createVariable('canliqxy', wp)
        canliqxy_var.units         = 'mm'
        canliqxy_var.standard_name = 'canopy intercepted liquid water for NoahMP'
        canliqxy_var[:]            = case_data._canliqxy
    
    if (case_data._eahxy != case_data._missing_value):
        eahxy_var               = nc_file.createVariable('eahxy', wp)
        eahxy_var.units         = 'Pa'
        eahxy_var.standard_name = 'canopy air vapor pressure for NoahMP'
        eahxy_var[:]            = case_data._eahxy
    
    if (case_data._cmxy != case_data._missing_value):
        cmxy_var               = nc_file.createVariable('cmxy', wp)
        cmxy_var.units         = 'none'
        cmxy_var.standard_name = 'surface drag coefficient for momentum for NoahMP'
        cmxy_var[:]            = case_data._cmxy
    
    if (case_data._chxy != case_data._missing_value):
        chxy_var               = nc_file.createVariable('chxy', wp)
        chxy_var.units         = 'none'
        chxy_var.standard_name = 'surface exchange coeff heat & moisture for NoahMP'
        chxy_var[:]            = case_data._chxy
    
    if (case_data._fwetxy != case_data._missing_value):
        fwetxy_var               = nc_file.createVariable('fwetxy', wp)
        fwetxy_var.units         = 'none'
        fwetxy_var.standard_name = 'area fraction of canopy that is wetted/snowed for NoahMP'
        fwetxy_var[:]            = case_data._fwetxy
    
    if (case_data._sneqvoxy != case_data._missing_value):
        sneqvoxy_var               = nc_file.createVariable('sneqvoxy', wp)
        sneqvoxy_var.units         = 'mm'
        sneqvoxy_var.standard_name = 'snow mass at previous time step for NoahMP'
        sneqvoxy_var[:]            = case_data._sneqvoxy
    
    if (case_data._alboldxy != case_data._missing_value):
        alboldxy_var               = nc_file.createVariable('alboldxy', wp)
        alboldxy_var.units         = 'none'
        alboldxy_var.standard_name = 'snow albedo at previous time step for NoahMP'
        alboldxy_var[:]            = case_data._alboldxy
    
    if (case_data._qsnowxy != case_data._missing_value):
        qsnowxy_var               = nc_file.createVariable('qsnowxy', wp)
        qsnowxy_var.units         = 'mm s-1'
        qsnowxy_var.standard_name = 'snow precipitation rate at surface for NoahMP'
        qsnowxy_var[:]            = case_data._qsnowxy
    
    if (case_data._wslakexy != case_data._missing_value):
        wslakexy_var               = nc_file.createVariable('wslakexy', wp)
        wslakexy_var.units         = 'mm'
        wslakexy_var.standard_name = 'lake water storage for NoahMP'
        wslakexy_var[:]            = case_data._wslakexy
    
    if (case_data._taussxy != case_data._missing_value):
        taussxy_var               = nc_file.createVariable('taussxy', wp)
        taussxy_var.units         = 'none'
        taussxy_var.standard_name = 'non-dimensional snow age for NoahMP'
        taussxy_var[:]            = case_data._taussxy
    
    
    if (case_data._waxy != case_data._missing_value):
        waxy_var               = nc_file.createVariable('waxy', wp)
        waxy_var.units         = 'mm'
        waxy_var.standard_name = 'water storage in aquifer for NoahMP'
        waxy_var[:]            = case_data._waxy
    
    if (case_data._wtxy != case_data._missing_value):
        wtxy_var               = nc_file.createVariable('wtxy', wp)
        wtxy_var.units         = 'mm'
        wtxy_var.standard_name = 'ater storage in aquifer and saturated soil for NoahMP'
        wtxy_var[:]            = case_data._wtxy
    
    if (case_data._zwtxy != case_data._missing_value):
        zwtxy_var               = nc_file.createVariable('zwtxy', wp)
        zwtxy_var.units         = 'm'
        zwtxy_var.standard_name = 'water table depth for NoahMP'
        zwtxy_var[:]            = case_data._zwtxy
    
    if (case_data._xlaixy != case_data._missing_value):
        xlaixy_var               = nc_file.createVariable('xlaixy', wp)
        xlaixy_var.units         = 'none'
        xlaixy_var.standard_name = 'leaf area index for NoahMP'
        xlaixy_var[:]            = case_data._xlaixy
    
    if (case_data._xsaixy != case_data._missing_value):
        xsaixy_var               = nc_file.createVariable('xsaixy', wp)
        xsaixy_var.units         = 'none'
        xsaixy_var.standard_name = 'stem area index for NoahMP'
        xsaixy_var[:]            = case_data._xsaixy
    
    if (case_data._lfmassxy != case_data._missing_value):
        lfmassxy_var               = nc_file.createVariable('lfmassxy', wp)
        lfmassxy_var.units         = 'g m-2'
        lfmassxy_var.standard_name = 'leaf mass for NoahMP'
        lfmassxy_var[:]            = case_data._lfmassxy
    
    if (case_data._stmassxy != case_data._missing_value):
        stmassxy_var               = nc_file.createVariable('stmassxy', wp)
        stmassxy_var.units         = 'g m-2'
        stmassxy_var.standard_name = 'stem mass for NoahMP'
        stmassxy_var[:]            = case_data._stmassxy
    
    if (case_data._rtmassxy != case_data._missing_value):
        rtmassxy_var               = nc_file.createVariable('rtmassxy', wp)
        rtmassxy_var.units         = 'g m-2'
        rtmassxy_var.standard_name = 'fine root mass for NoahMP'
        rtmassxy_var[:]            = case_data._rtmassxy
    
    if (case_data._woodxy != case_data._missing_value):
        woodxy_var               = nc_file.createVariable('woodxy', wp)
        woodxy_var.units         = 'g m-2'
        woodxy_var.standard_name = 'wood mass including woody roots for NoahMP'
        woodxy_var[:]            = case_data._woodxy
    
    if (case_data._stblcpxy != case_data._missing_value):
        stblcpxy_var               = nc_file.createVariable('stblcpxy', wp)
        stblcpxy_var.units         = 'g m-2'
        stblcpxy_var.standard_name = 'stable carbon in deep soil for NoahMP'
        stblcpxy_var[:]            = case_data._stblcpxy
    
    if (case_data._fastcpxy != case_data._missing_value):
        fastcpxy_var               = nc_file.createVariable('fastcpxy', wp)
        fastcpxy_var.units         = 'g m-2'
        fastcpxy_var.standard_name = 'short-lived carbon in shallow soil for NoahMP'
        fastcpxy_var[:]            = case_data._fastcpxy
    
    if (case_data._smcwtdxy != case_data._missing_value):
        smcwtdxy_var               = nc_file.createVariable('smcwtdxy', wp)
        smcwtdxy_var.units         = 'm3 m-3'
        smcwtdxy_var.standard_name = 'oil water content between the bottom of the soil and the water table for NoahMP'
        smcwtdxy_var[:]            = case_data._smcwtdxy
    
    if (case_data._deeprechxy != case_data._missing_value):
        deeprechxy_var               = nc_file.createVariable('deeprechxy', wp)
        deeprechxy_var.units         = 'm'
        deeprechxy_var.standard_name = 'echarge to or from the water table when deep for NoahMP'
        deeprechxy_var[:]            = case_data._deeprechxy
    
    if (case_data._rechxy != case_data._missing_value):
        rechxy_var               = nc_file.createVariable('rechxy', wp)
        rechxy_var.units         = 'm'
        rechxy_var.standard_name = 'recharge to or from the water table when shallow for NoahMP'
        rechxy_var[:]            = case_data._rechxy
    
    if (case_data._snowxy != case_data._missing_value):
        snowxy_var               = nc_file.createVariable('snowxy', wp)
        snowxy_var.units         = 'none'
        snowxy_var.standard_name = 'number of snow layers for NoahMP'
        snowxy_var[:]            = case_data._snowxy
    
    if (case_data._wetness != case_data._missing_value):
        wetness_var               = nc_file.createVariable('wetness', wp)
        wetness_var.units         = 'none'
        wetness_var.standard_name = 'normalized soil wetness for RUC LSM'
        wetness_var[:]            = case_data._wetness
        
    if (case_data._clw_surf_land != case_data._missing_value):
        clw_surf_land_var               = nc_file.createVariable('clw_surf_land', wp)
        clw_surf_land_var.units         = 'kg kg-1'
        clw_surf_land_var.standard_name = 'cloud condensed water mixing ratio at surface over land for RUC LSM'
        clw_surf_land_var[:]            = case_data._clw_surf_land
        
    if (case_data._clw_surf_ice != case_data._missing_value):
        clw_surf_ice_var               = nc_file.createVariable('clw_surf_ice', wp)
        clw_surf_ice_var.units         = 'kg kg-1'
        clw_surf_ice_var.standard_name = 'cloud condensed water mixing ratio at surface over ice for RUC LSM'
        clw_surf_ice_var[:]            = case_data._clw_surf_ice
    
    if (case_data._qwv_surf_land != case_data._missing_value):
        qwv_surf_land_var               = nc_file.createVariable('qwv_surf_land', wp)
        qwv_surf_land_var.units         = 'kg kg-1'
        qwv_surf_land_var.standard_name = 'water vapor mixing ratio at surface over land for RUC LSM'
        qwv_surf_land_var[:]            = case_data._qwv_surf_land
    
    if (case_data._qwv_surf_ice != case_data._missing_value):
        qwv_surf_ice_var               = nc_file.createVariable('qwv_surf_ice', wp)
        qwv_surf_ice_var.units         = 'kg kg-1'
        qwv_surf_ice_var.standard_name = 'water vapor mixing ratio at surface over ice for RUC LSM'
        qwv_surf_ice_var[:]            = case_data._qwv_surf_ice
    
    if (case_data._tsnow_land != case_data._missing_value):
        tsnow_land_var               = nc_file.createVariable('tsnow_land', wp)
        tsnow_land_var.units         = 'K'
        tsnow_land_var.standard_name = 'snow temperature at the bottom of the first snow layer over land for RUC LSM'
        tsnow_land_var[:]            = case_data._tsnow_land
    
    if (case_data._tsnow_ice != case_data._missing_value):
        tsnow_ice_var               = nc_file.createVariable('tsnow_ice', wp)
        tsnow_ice_var.units         = 'K'
        tsnow_ice_var.standard_name = 'snow temperature at the bottom of the first snow layer over land for RUC LSM'
        tsnow_ice_var[:]            = case_data._tsnow_ice
    
    if (case_data._snowfallac_land != case_data._missing_value):
        snowfallac_land_var               = nc_file.createVariable('snowfallac_land', wp)
        snowfallac_land_var.units         = 'kg m-2'
        snowfallac_land_var.standard_name = 'run-total snow accumulation on the ground over land for RUC LSM'
        snowfallac_land_var[:]            = case_data._snowfallac_land
    
    if (case_data._snowfallac_ice != case_data._missing_value):
        snowfallac_ice_var               = nc_file.createVariable('snowfallac_ice', wp)
        snowfallac_ice_var.units         = 'kg m-2'
        snowfallac_ice_var.standard_name = 'run-total snow accumulation on the ground over land for RUC LSM'
        snowfallac_ice_var[:]            = case_data._snowfallac_ice
        
    if (case_data._sncovr_ice != case_data._missing_value):
        sncovr_ice_var               = nc_file.createVariable('sncovr_ice', wp)
        sncovr_ice_var.units         = 'none'
        sncovr_ice_var.standard_name = 'surface snow area fraction over ice'
        sncovr_ice_var[:]            = case_data._sncovr_ice
    
    if (case_data._sfalb_lnd != case_data._missing_value):
        sfalb_lnd_var               = nc_file.createVariable('sfalb_lnd', wp)
        sfalb_lnd_var.units         = 'none'
        sfalb_lnd_var.standard_name = 'surface albedo over land for RUC LSM'
        sfalb_lnd_var[:]            = case_data._sfalb_lnd
    
    if (case_data._sfalb_lnd_bck != case_data._missing_value):
        sfalb_lnd_bck_var               = nc_file.createVariable('sfalb_lnd_bckß', wp)
        sfalb_lnd_bck_var.units         = 'none'
        sfalb_lnd_bck_var.standard_name = 'surface snow-free albedo over land for RUC LSM'
        sfalb_lnd_bck_var[:]            = case_data._sfalb_lnd_bck
    
    if (case_data._emis_ice != case_data._missing_value):
        emis_ice_var               = nc_file.createVariable('emis_ice', wp)
        emis_ice_var.units         = 'none'
        emis_ice_var.standard_name = 'surface emissivity over ice for RUC LSM'
        emis_ice_var[:]            = case_data._emis_ice
    
    if (case_data._lai != case_data._missing_value):
        lai_var               = nc_file.createVariable('lai', wp)
        lai_var.units         = 'none'
        lai_var.standard_name = 'leaf area index for RUC LSM'
        lai_var[:]            = case_data._lai
    
    area_var               = nc_file.createVariable('area', wp)
    area_var.units         = 'm2'
    area_var.standard_name = 'grid cell area'
    if ('column_area' in nml_keys and case_nml['case_config']['column_area']):
        area_var[:]            = case_nml['case_config']['column_area']
        message = 'Since column_area was supplied in the case namelist, it will be used instead of the data from the case data file (if it exists).'
        logging.info(message)
    elif (case_data._area != case_data._missing_value):
        area_var[:]            = case_data._area
    
    if (case_data._stddev != case_data._missing_value):
        stddev_var               = nc_file.createVariable('stddev', wp)
        stddev_var.units         = 'm'
        stddev_var.standard_name = 'standard deviation of subgrid orography'
        stddev_var[:]            = case_data._stddev
        
    if (case_data._convexity != case_data._missing_value):
        convexity_var               = nc_file.createVariable('convexity', wp)
        convexity_var.units         = 'none'
        convexity_var.standard_name = 'convexity of subgrid orography'
        convexity_var[:]            = case_data._convexity
    
    if (case_data._oa1 != case_data._missing_value):
        oa1_var               = nc_file.createVariable('oa1', wp)
        oa1_var.units         = 'none'
        oa1_var.standard_name = 'assymetry of subgrid orography 1'
        oa1_var[:]            = case_data._oa1
    
    if (case_data._oa2 != case_data._missing_value):
        oa2_var               = nc_file.createVariable('oa2', wp)
        oa2_var.units         = 'none'
        oa2_var.standard_name = 'assymetry of subgrid orography 2'
        oa2_var[:]            = case_data._oa2
    
    if (case_data._oa3 != case_data._missing_value):
        oa3_var               = nc_file.createVariable('oa3', wp)
        oa3_var.units         = 'none'
        oa3_var.standard_name = 'assymetry of subgrid orography 3'
        oa3_var[:]            = case_data._oa3
    
    if (case_data._oa4 != case_data._missing_value):
        oa4_var               = nc_file.createVariable('oa4', wp)
        oa4_var.units         = 'none'
        oa4_var.standard_name = 'assymetry of subgrid orography 4'
        oa4_var[:]            = case_data._oa4
        
    if (case_data._ol1 != case_data._missing_value):
        ol1_var               = nc_file.createVariable('ol1', wp)
        ol1_var.units         = 'none'
        ol1_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 1'
        ol1_var[:]            = case_data._ol1
        
    if (case_data._ol2 != case_data._missing_value):
        ol2_var               = nc_file.createVariable('ol2', wp)
        ol2_var.units         = 'none'
        ol2_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 2'
        ol2_var[:]            = case_data._ol2
        
    if (case_data._ol3 != case_data._missing_value):
        ol3_var               = nc_file.createVariable('ol3', wp)
        ol3_var.units         = 'none'
        ol3_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 3'
        ol3_var[:]            = case_data._ol3
        
    if (case_data._ol4 != case_data._missing_value):
        ol4_var               = nc_file.createVariable('ol4', wp)
        ol4_var.units         = 'none'
        ol4_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 4'
        ol4_var[:]            = case_data._ol4
        
    if (case_data._theta_oro != case_data._missing_value):
        theta_oro_var               = nc_file.createVariable('theta_oro', wp)
        theta_oro_var.units         = 'deg'
        theta_oro_var.standard_name = 'angle with respect to east of maximum subgrid orographic variations'
        theta_oro_var[:]            = case_data._theta_oro
        
    if (case_data._gamma != case_data._missing_value):
        gamma_var               = nc_file.createVariable('gamma', wp)
        gamma_var.units         = 'none'
        gamma_var.standard_name = 'anisotropy of subgrid orography'
        gamma_var[:]            = case_data._gamma
        
    if (case_data._sigma != case_data._missing_value):
        sigma_var               = nc_file.createVariable('sigma', wp)
        sigma_var.units         = 'none'
        sigma_var.standard_name = 'slope of subgrid orography'
        sigma_var[:]            = case_data._sigma
    
    if (case_data._elvmax != case_data._missing_value):
        elvmax_var               = nc_file.createVariable('elvmax', wp)
        elvmax_var.units         = 'm'
        elvmax_var.standard_name = 'maximum of subgrid orography'
        elvmax_var[:]            = case_data._elvmax
        
    if (case_data._oro != case_data._missing_value):
        oro_var               = nc_file.createVariable('oro', wp)
        oro_var.units         = 'm'
        oro_var.standard_name = 'orography'
        oro_var[:]            = case_data._oro
        
    if (case_data._oro_uf != case_data._missing_value):
        oro_uf_var               = nc_file.createVariable('oro_uf', wp)
        oro_uf_var.units         = 'm'
        oro_uf_var.standard_name = 'unfiltered orography'
        oro_uf_var[:]            = case_data._oro_uf
        
    if (case_data._landfrac != case_data._missing_value):
        landfrac_var               = nc_file.createVariable('landfrac', wp)
        landfrac_var.units         = 'none'
        landfrac_var.standard_name = 'fraction of horizontal grid area occupied by land'
        landfrac_var[:]            = case_data._landfrac
        
    if (case_data._lakefrac != case_data._missing_value):
        lakefrac_var               = nc_file.createVariable('lakefrac', wp)
        lakefrac_var.units         = 'none'
        lakefrac_var.standard_name = 'fraction of horizontal grid area occupied by lake'
        lakefrac_var[:]            = case_data._lakefrac
        
    if (case_data._lakedepth != case_data._missing_value):
        lakedepth_var               = nc_file.createVariable('lakedepth', wp)
        lakedepth_var.units         = 'none'
        lakedepth_var.standard_name = 'lake depth'
        lakedepth_var[:]            = case_data._lakedepth
    
    if (case_data._tref != case_data._missing_value):
        tref_var               = nc_file.createVariable('tref', wp)
        tref_var.units         = 'K'
        tref_var.standard_name = 'sea surface reference temperature for NSST'
        tref_var[:]            = case_data._tref
    
    if (case_data._z_c != case_data._missing_value):
        z_c_var               = nc_file.createVariable('z_c', wp)
        z_c_var.units         = 'm'
        z_c_var.standard_name = 'sub-layer cooling thickness for NSST'
        z_c_var[:]            = case_data._z_c
        
    if (case_data._c_0 != case_data._missing_value):
        c_0_var               = nc_file.createVariable('c_0', wp)
        c_0_var.units         = 'none'
        c_0_var.standard_name = 'coefficient 1 to calculate d(Tz)/d(Ts) for NSST'
        c_0_var[:]            = case_data._c_0
        
    if (case_data._c_d != case_data._missing_value):
        c_d_var               = nc_file.createVariable('c_d', wp)
        c_d_var.units         = 'none'
        c_d_var.standard_name = 'coefficient 2 to calculate d(Tz)/d(Ts) for NSST'
        c_d_var[:]            = case_data._c_d
        
    if (case_data._w_0!= case_data._missing_value):
        w_0_var               = nc_file.createVariable('w_0', wp)
        w_0_var.units         = 'none'
        w_0_var.standard_name = 'coefficient 3 to calculate d(Tz)/d(Ts) for NSST'
        w_0_var[:]            = case_data._w_0
        
    if (case_data._w_d != case_data._missing_value):
        w_d_var               = nc_file.createVariable('w_d', wp)
        w_d_var.units         = 'none'
        w_d_var.standard_name = 'coefficient 4 to calculate d(Tz)/d(Ts) for NSST'
        w_d_var[:]            = case_data._w_d
        
    if (case_data._xt != case_data._missing_value):
        xt_var               = nc_file.createVariable('xt', wp)
        xt_var.units         = 'K m'
        xt_var.standard_name = 'heat content in diurnal thermocline layer for NSST'
        xt_var[:]            = case_data._xt
        
    if (case_data._xs != case_data._missing_value):
        xs_var               = nc_file.createVariable('xs', wp)
        xs_var.units         = 'ppt m'
        xs_var.standard_name = 'salinity content in diurnal thermocline layer for NSST'
        xs_var[:]            = case_data._xs
        
    if (case_data._xu != case_data._missing_value):
        xu_var               = nc_file.createVariable('xu', wp)
        xu_var.units         = 'm2 s-1'
        xu_var.standard_name = 'u-current in diurnal thermocline layer for NSST'
        xu_var[:]            = case_data._xu
        
    if (case_data._xv != case_data._missing_value):
        xv_var               = nc_file.createVariable('xv', wp)
        xv_var.units         = 'm2 s-1'
        xv_var.standard_name = 'v-current in diurnal thermocline layer for NSST'
        xv_var[:]            = case_data._xv
        
    if (case_data._xz != case_data._missing_value):
        xz_var               = nc_file.createVariable('xz', wp)
        xz_var.units         = 'm'
        xz_var.standard_name = 'thickness of diurnal thermocline layer for NSST'
        xz_var[:]            = case_data._xz
        
    if (case_data._zm != case_data._missing_value):
        zm_var               = nc_file.createVariable('zm', wp)
        zm_var.units         = 'm'
        zm_var.standard_name = 'thickness of ocean mixed layer for NSST'
        zm_var[:]            = case_data._zm
        
    if (case_data._xtts != case_data._missing_value):
        xtts_var               = nc_file.createVariable('xtts', wp)
        xtts_var.units         = 'm'
        xtts_var.standard_name = 'sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST'
        xtts_var[:]            = case_data._xtts
    
    if (case_data._xzts != case_data._missing_value):
        xzts_var               = nc_file.createVariable('xzts', wp)
        xzts_var.units         = 'm K-1'
        xzts_var.standard_name = 'sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST'
        xzts_var[:]            = case_data._xzts
        
    if (case_data._d_conv != case_data._missing_value):
        d_conv_var               = nc_file.createVariable('d_conv', wp)
        d_conv_var.units         = 'm'
        d_conv_var.standard_name = 'thickness of free convection layer for NSST'
        d_conv_var[:]            = case_data._d_conv
        
    if (case_data._ifd != case_data._missing_value):
        ifd_var               = nc_file.createVariable('ifd', wp)
        ifd_var.units         = 'none'
        ifd_var.standard_name = 'index to start DTM run for NSST'
        ifd_var[:]            = case_data._ifd
        
    if (case_data._dt_cool != case_data._missing_value):
        dt_cool_var               = nc_file.createVariable('dt_cool', wp)
        dt_cool_var.units         = 'K'
        dt_cool_var.standard_name = 'sub-layer cooling amount for NSST'
        dt_cool_var[:]            = case_data._dt_cool
        
    if (case_data._qrain != case_data._missing_value):
        qrain_var               = nc_file.createVariable('qrain', wp)
        qrain_var.units         = 'W m-2'
        qrain_var.standard_name = 'sensible heat due to rainfall for NSST'
        qrain_var[:]            = case_data._qrain
              
    if (case_data._theta_il[0] != case_data._missing_value):
        thetal_var                   = nc_file.createVariable('thetal', wp, ('t0','lev'))
        thetal_var.units             = 'K'
        thetal_var.standard_name     = 'air_liquid_potential_temperature'
        thetal_var[:]                = case_data._theta_il[:]
    
    if (case_data._t[0] != case_data._missing_value):
        t_var                   = nc_file.createVariable('t', wp, ('t0','lev'))
        t_var.units             = 'K'
        t_var.standard_name     = 'absolute temperature'
        t_var[:]                = case_data._t[:]
    
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
    
    if (case_data._stc[0] != case_data._missing_value):
        stc_var               = nc_file.createVariable('stc', wp, ('t0','nsoil'))
        stc_var.units         = 'K'
        stc_var.standard_name = 'initial profile of soil temperature'
        stc_var[:]            = case_data._stc[:]
    
    if (case_data._smc[0] != case_data._missing_value):
        smc_var               = nc_file.createVariable('smc', wp, ('t0','nsoil'))
        smc_var.units         = 'kg'
        smc_var.standard_name = 'initial profile of soil moisture'
        smc_var[:]            = case_data._smc[:]
        
    if (case_data._slc[0] != case_data._missing_value):
        slc_var               = nc_file.createVariable('slc', wp, ('t0','nsoil'))
        slc_var.units         = 'kg'
        slc_var.standard_name = 'initial profile of soil liquid moisture'
        slc_var[:]            = case_data._slc[:]
    
    if (case_data._snicexy[0] != case_data._missing_value):
        snicexy_var               = nc_file.createVariable('snicexy', wp, ('t0','nsnow'))
        snicexy_var.units         = 'mm'
        snicexy_var.standard_name = 'initial profile of snow layer ice'
        snicexy_var[:]            = case_data._snicexy[:]
    
    if (case_data._snliqxy[0] != case_data._missing_value):
        snliqxy_var               = nc_file.createVariable('snliqxy', wp, ('t0','nsnow'))
        snliqxy_var.units         = 'mm'
        snliqxy_var.standard_name = 'initial profile of snow layer liquid'
        snliqxy_var[:]            = case_data._snliqxy[:]

    if (case_data._tsnoxy[0] != case_data._missing_value):
        tsnoxy_var               = nc_file.createVariable('tsnoxy', wp, ('t0','nsnow'))
        tsnoxy_var.units         = 'K'
        tsnoxy_var.standard_name = 'initial profile of snow layer temperature'
        tsnoxy_var[:]            = case_data._tsnoxy[:]
    
    if (case_data._smoiseq[0] != case_data._missing_value):
        smoiseq_var               = nc_file.createVariable('smoiseq', wp, ('t0','nsoil'))
        smoiseq_var.units         = 'm3 m-3'
        smoiseq_var.standard_name = 'initial profile of equilibrium soil water content'
        smoiseq_var[:]            = case_data._smoiseq[:]
    
    if (case_data._zsnsoxy[0] != case_data._missing_value):
        zsnsoxy_var               = nc_file.createVariable('zsnxosy', wp, ('t0','nsoil_plus_nsnow'))
        zsnsoxy_var.units         = 'm'
        zsnsoxy_var.standard_name = 'layer bottom depth from snow surface'
        zsnsoxy_var[:]            = case_data._zsnsoxy[:]
    
    if (case_data._tiice[0] != case_data._missing_value):
        tiice_var               = nc_file.createVariable('tiice', wp, ('t0','nice'))
        tiice_var.units         = 'K'
        tiice_var.standard_name = 'sea ice internal temperature'
        tiice_var[:]            = case_data._tiice[:]
    
    if (case_data._tslb[0] != case_data._missing_value):
        tslb_var               = nc_file.createVariable('tslb', wp, ('t0','nsoil'))
        tslb_var.units         = 'K'
        tslb_var.standard_name = 'soil temperature for RUC LSM'
        tslb_var[:]            = case_data._tslb[:]
    
    if (case_data._smois[0] != case_data._missing_value):
        smois_var               = nc_file.createVariable('smois', wp, ('t0','nsoil'))
        smois_var.units         = 'none'
        smois_var.standard_name = 'volume fraction of soil moisture for RUC LSM'
        smois_var[:]            = case_data._smois[:]
    
    if (case_data._sh2o[0] != case_data._missing_value):
        sh2o_var               = nc_file.createVariable('sh2o', wp, ('t0','nsoil'))
        sh2o_var.units         = 'none'
        sh2o_var.standard_name = 'volume fraction of unfrozen soil moisture for RUC LSM'
        sh2o_var[:]            = case_data._sh2o[:]
    
    if (case_data._smfr[0] != case_data._missing_value):
        smfr_var               = nc_file.createVariable('smfr', wp, ('t0','nsoil'))
        smfr_var.units         = 'none'
        smfr_var.standard_name = 'volume fraction of frozen soil moisture for RUC LSM'
        smfr_var[:]            = case_data._smfr[:]
    
    if (case_data._flfr[0] != case_data._missing_value):
        flfr_var               = nc_file.createVariable('flfr', wp, ('t0','nsoil'))
        flfr_var.units         = 'none'
        flfr_var.standard_name = 'flag for frozen soil physics for RUC LSM'
        flfr_var[:]            = case_data._flfr[:]
        
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
        logging.critical(message)
        raise Exception(message)
        # tnta_adv_var                    = nc_file.createVariable('tnta_adv', wp, ('time','lev'))
        # tnta_adv_var.units              = 'K s-1'
        # tnta_adv_var.standard_name      = 'tendency_of_air_temperature_due_to_advection'
        # tnta_adv_var[:]                 = np.swapaxes(case_data._tnta_adv[:],0,1)
        
    if (nc_file.adv_qv == forcing_on):
        message = 'adv_qv is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnqv_adv_var                    = nc_file.createVariable('tnqv_adv', wp, ('time','lev'))
        # tnqv_adv_var.units              = 'kg kg-1 s-1'
        # tnqv_adv_var.standard_name      = 'tendency_of_specific_humidity_due_to_advection'
        # tnqv_adv_var[:]                 = np.swapaxes(case_data._tnqv_adv[:],0,1)
    
    if (nc_file.adv_ua == forcing_on):
        message = 'adv_ua is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnua_adv_var                    = nc_file.createVariable('tnua_adv', wp, ('time','lev'))
        # tnua_adv_var.units              = 'm s-2'
        # tnua_adv_var.standard_name      = 'tendency_of_eastward_wind_due_to_advection'
        # tnua_adv_var[:]                 = np.swapaxes(case_data._tnua_adv[:],0,1)
    
    if (nc_file.adv_va == forcing_on):
        message = 'adv_va is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnva_adv_var                    = nc_file.createVariable('tnva_adv', wp, ('time','lev'))
        # tnva_adv_var.units              = 'm s-2'
        # tnva_adv_var.standard_name      = 'tendency_of_northward_wind_due_to_advection'
        # tnva_adv_var[:]                 = np.swapaxes(case_data._tnva_adv[:],0,1)
    
    if (nc_file.adv_theta == forcing_on):
        message = 'adv_theta is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
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
        logging.critical(message)
        raise Exception(message)
        # tnrv_adv_var                    = nc_file.createVariable('tnrv_adv', wp, ('time','lev'))
        # tnrv_adv_var.units              = 'kg kg-1 s-1'
        # tnrv_adv_var.standard_name      = 'tendency_of_humidity_mixing_ratio_due_to_advection'
        # tnrv_adv_var[:]                 = np.swapaxes(case_data._tnrv_adv[:],0,1)
    
    if (nc_file.adv_rt == forcing_on):
        message = 'adv_rt is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
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
def main():
    (case_name, use_area, debug) = parse_arguments()
    
    setup_logging(debug)
    
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
    
if __name__ == '__main__':
    main()
