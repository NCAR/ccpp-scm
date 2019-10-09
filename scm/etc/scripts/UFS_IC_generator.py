#!/usr/bin/env python

import argparse
import logging
import os
import fnmatch
import logging
from netCDF4 import Dataset
import numpy as np
from scipy.spatial.distance import euclidean
from shapely.geometry import Point, Polygon
import copy
import math
import f90nml
import re

###############################################################################
# Global settings                                                             #
###############################################################################

#Physical constants
earth_radius = 6371000.0 #m
rdgas  = 287.05
rvgas  = 461.50
zvir = rvgas/rdgas - 1.
grav=9.80665

missing_value = 9.99e20

# Path to the directory containing processed case input files ()
PROCESSED_CASE_DIR = '../../data/processed_case_input'

NOAHMP_TABLES_DIR = '../../data/raw_case_input/NoahMP_tables'

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-l', '--location',   help='longitude and latitude in degress E and N, respectively, separated by a space', nargs=2, type=float)
group1.add_argument('-ij','--index',      help='i,j indices within the tile (if known - bypasses search for closest model point to lon/lat location)', nargs=2, type=int)
parser.add_argument('-d', '--date',       help='date corresponding to initial conditions in YYYYMMDDHHMM format', required=True)
parser.add_argument('-i', '--in_dir',     help='input directory path containing FV3 input files', required=True)
parser.add_argument('-g', '--grid_dir',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-t', '--tile',       help='tile of desired point (if known - bypasses tile search if present)', type=int, choices=range(1,7))
parser.add_argument('-a', '--area',       help='area of grid cell in m^2', type=float)
parser.add_argument('-n', '--case_name',  help='name of case', required=True)

###############################################################################
# Functions and subroutines                                                   #
###############################################################################

def parse_arguments():
    """Parse command line arguments"""
    args = parser.parse_args()
    location = args.location
    index = args.index
    date = args.date
    in_dir = args.in_dir
    grid_dir = args.grid_dir
    tile = args.tile
    area = args.area
    case_name = args.case_name
    
    #validate args
    if not os.path.exists(in_dir):
        message = 'The directory {0} does not exist'.format(in_dir)
        logging.critical(message)
        raise Exception(message)
    
    if not index:
        if not 0 <= location[0] <= 360 :
            message = 'The longitude {0} is outside of the range {1}'.format(location[0], '[0,360]')
            logging.critical(message)
            raise Exception(message)
        
        if not -90 <= location[1] <= 90:
            message = 'The latitude {0} is outside of the range {1}'.format(location[1], '[-90,90]')
            logging.critical(message)
            raise Exception(message)
    
    date_dict = {}
    if len(date) != 12:
        message = 'The entered date {0} does not have the 12 characters expected in the format YYYYMMDDHHMM'.format(date)
        logging.critical(message)
        raise Exception(message)
    else:
        date_dict["year"] = np.int(date[0:4])
        date_dict["month"] = np.int(date[4:6])
        date_dict["day"] = np.int(date[6:8])
        date_dict["hour"] = np.int(date[8:10])
        date_dict["minute"] = np.int(date[10:])
        
    return (location, index, date_dict, in_dir, grid_dir, tile, area, case_name)

def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)
    
def find_tile(loc, dir):
    """Find the FV3 tile with the given lon/lat"""
    
    filename_pattern = '*grid.tile*.nc'
    
    #find all grid files in the directory
    grid_fnames = []
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          grid_fnames.append(f_name)
    if not grid_fnames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    #non-polar tiles can use traditional 2D point-in-polygon methods; if a point is not in a non-polar tile,
    #it is in one of the polar tiles, and the tile can be distinguished by the sign of latitude of the point
    polar_tile_filenames = []
    found_tile = False
    for f_name in grid_fnames:
        if not found_tile:
            nc_file = Dataset('{0}/{1}'.format(dir,f_name))
            longitude = np.array(nc_file['x']).swapaxes(0,1)
            #longitude[longitude > 180] -= 360
            latitude = np.array(nc_file['y']).swapaxes(0,1)
            nc_file.close()
            
            adj_long = False        
            #look for reversal of longitude; if found, adjust longitude so that 0-360 transition doesn't exist
            for row in longitude:
                if not (np.all(np.diff(row) >= 0) or np.all(np.diff(row) <= 0)):
                    adj_long = True
            if adj_long:
                longitude[longitude < 180] += 360
            
            #print f_name
            edge_1_lon = longitude[0,:]
            edge_1_lat = latitude[0,:]
            edge_1 = zip(edge_1_lon, edge_1_lat)
            #print edge_1[0], edge_1[-1]
            #print len(edge_1)
            #print edge_1
            
            edge_2_lon = longitude[:,-1]
            edge_2_lat = latitude[:,-1]
            edge_2 = zip(edge_2_lon, edge_2_lat)
            #print edge_2[0], edge_2[-1]
            #print len(edge_2)
            #print edge_2
            
            edge_3_lon = longitude[-1,:]
            edge_3_lat = latitude[-1,:]
            edge_3 = zip(edge_3_lon, edge_3_lat)
            edge_3.reverse()
            #print edge_3[0], edge_3[-1]
            #print len(edge_3)
            #print edge_3
            
            edge_4_lon = longitude[:,0]
            edge_4_lat = latitude[:,0]
            edge_4 = zip(edge_4_lon, edge_4_lat)
            edge_4.reverse()
            #print edge_4[0], edge_4[-1]
            #print len(edge_4)
            #print edge_4
            
            polygon_points = edge_1 + edge_2 + edge_3 + edge_4
            #print len(polygon_points)
            
            tile_polygon = Polygon(polygon_points)
            tile_polygon = tile_polygon.simplify(0)
            
            if tile_polygon.is_valid:
                temp_loc = copy.deepcopy(loc)
                if adj_long:
                    if loc[0] < 180:
                        temp_loc[0] += 360
                loc_point = Point(temp_loc)
                if tile_polygon.contains(loc_point):
                    found_tile = True
                    return f_name.split('tile')[1].split('.nc')[0] 
            else:
                polar_tile_filenames.append(f_name)
                
    for f_name in polar_tile_filenames:
        nc_file = Dataset('{0}/{1}'.format(dir,f_name))
        latitude = np.array(nc_file['y']).swapaxes(0,1)
        nc_file.close()
        
        if np.sign(np.mean(latitude)) == np.sign(loc[1]):
            found_tile = True
            return f_name.split('tile')[1].split('.nc')[0]        
    return -1

def find_loc_indices(loc, dir, tile):
    filename_pattern = 'sfc_data.tile{0}.nc'.format(tile)
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    longitude = np.array(nc_file['geolon'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    latitude = np.array(nc_file['geolat'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    nc_file.close()
    
    adj_long = False        
    #look for reversal of longitude; if found, adjust longitude so that 0-360 transition doesn't exist
    temp_loc = copy.deepcopy(loc)
    for row in longitude:
        if not (np.all(np.diff(row) >= 0) or np.all(np.diff(row) <= 0)):
            adj_long = True
    if adj_long:
        longitude[longitude < 180] += 360
        if loc[0] < 180:
            temp_loc[0] += 360
    
    eucl_dist = np.zeros((longitude.shape[0],longitude.shape[1]))
    cart_loc = np.array(sph2cart(math.radians(temp_loc[0]), math.radians(temp_loc[1]), earth_radius))
    for i in range(len(longitude)):
        for j in range(len(longitude[i])):
            cart_cell = np.array(sph2cart(math.radians(longitude[i,j]), math.radians(latitude[i,j]), earth_radius))
            eucl_dist[i,j] = np.linalg.norm(cart_loc - cart_cell)
    i,j = np.unravel_index(eucl_dist.argmin(), eucl_dist.shape)
    
    # filename_pattern = '*grid.tile{0}.nc'.format(tile)
    # for f_name in os.listdir(dir):
    #    if fnmatch.fnmatch(f_name, filename_pattern):
    #       filename = f_name
    # nc_file = Dataset('{0}/{1}'.format(dir,filename))
    # super_lon = np.array(nc_file['x'])#.swapaxes(0,1)
    # super_lat = np.array(nc_file['y'])#.swapaxes(0,1)
    # nc_file.close()
    
    #print super_lon[i*2+1,j*2+1], super_lat[i*2+1,j*2+1] => lon/lat in supergrid (C96_grid.tile*.nc) are same values as in geolon/geolat 
    
    return (i,j,longitude[i,j]%360.0, latitude[i,j], eucl_dist[i,j])

def find_lon_lat_of_indices(indices, dir, tile):
    filename_pattern = 'sfc_data.tile{0}.nc'.format(tile)
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    longitude = np.array(nc_file['geolon'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    latitude = np.array(nc_file['geolat'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    nc_file.close()
    
    return (longitude[indices[1],indices[0]], latitude[indices[1],indices[0]])
    
    

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return (x, y, z)    

def get_UFS_IC_data(dir, tile, i, j):
    
    state_data = get_UFS_state_data(dir, tile, i, j)
    surface_data = get_UFS_surface_data(dir, tile, i, j)
    oro_data = get_UFS_oro_data(dir, tile, i, j)
    vgrid_data = get_UFS_vgrid_data(dir)
    
    #calculate derived quantities
    
    #temperature
    nlevs = state_data["nlevs"]
    gz=state_data["z"]*grav
    pn1=np.zeros([nlevs+1])
    temp=np.zeros([nlevs])
    for k in range(nlevs+1):
      pn1[k]=np.log(vgrid_data["ak"][k]+state_data["p_surf"]*vgrid_data["bk"][k])
    for k in range(nlevs):
      temp[k] = (gz[k]-gz[k+1])/( rdgas*(pn1[k+1]-pn1[k])*(1.+zvir*state_data["qv"][k]) )
    state_data["T"] = temp
    state_data["pres"] = np.exp(pn1[0:nlevs])
    
    return (state_data, surface_data, oro_data)
    
def get_UFS_state_data(dir, tile, i, j):
    
    nc_file = Dataset('{0}/{1}'.format(dir,'gfs_data.tile{0}.nc'.format(tile)))
    
    # assume model contains one less level than the cold start spectral GFS initial conditions
    nlevs=len(nc_file.dimensions['lev'])-1
    
    # upper air fields from initial conditions
    zh=nc_file['zh'][::-1,j,i]
    uw1=nc_file['u_w'][::-1,j,i]
    uw2=nc_file['u_w'][::-1,j,i+1]
    us1=nc_file['u_s'][::-1,j,i]
    us2=nc_file['u_s'][::-1,j+1,i]
    vw1=nc_file['v_w'][::-1,j,i]
    vw2=nc_file['v_w'][::-1,j,i+1]
    vs1=nc_file['v_s'][::-1,j,i]
    vs2=nc_file['v_s'][::-1,j+1,i]
    ucomp=0.25*(uw1+uw2+us1+us2)  # estimate u winds on the a grid
    vcomp=0.25*(vw1+vw2+vs1+vs2)  # estimate v winds on the a grid
    sphum=nc_file['sphum'][::-1,j,i]
    # o3 and qv are taken from ics. 
    o3=nc_file['o3mr'][::-1,j,i]
    liqwat=nc_file['liq_wat'][:-1,j,i]

    # surface pressure and skin temperature
    ps=nc_file['ps'][j,i]
    
    nc_file.close()
    
    #put data in a dictionary
    state = {
        "nlevs": nlevs,
        "z": zh,
        "u": ucomp,
        "v": vcomp,
        "qv": sphum,
        "o3": o3,
        "ql": liqwat,
        "p_surf": ps
    }
    return state

def get_UFS_surface_data(dir, tile, i, j):
    
    nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.tile{0}.nc'.format(tile)))
    
    ts_in=nc_file['tsea'][j,i]

    # land state
    stc_in=nc_file['stc'][:,j,i]
    smc_in=nc_file['smc'][:,j,i]
    slc_in=nc_file['slc'][:,j,i]
    tg3_in=nc_file['tg3'][j,i]

    # surface properties
    uustar_in=nc_file['uustar'][j,i]
    alvsf_in=nc_file['alvsf'][j,i]
    alvwf_in=nc_file['alvwf'][j,i]
    alnsf_in=nc_file['alnsf'][j,i]
    alnwf_in=nc_file['alnwf'][j,i]
    facsf_in=nc_file['facsf'][j,i]
    facwf_in=nc_file['facwf'][j,i]
    styp_in=nc_file['stype'][j,i]
    slope_in=nc_file['slope'][j,i]
    vtyp_in=nc_file['vtype'][j,i]
    vfrac_in=nc_file['vfrac'][j,i]
    shdmin_in=nc_file['shdmin'][j,i]
    shdmax_in=nc_file['shdmax'][j,i]
    zorl_in=nc_file['zorl'][j,i]
    slmsk_in=nc_file['slmsk'][j,i]
    canopy_in=nc_file['canopy'][j,i]
    hice_in=nc_file['hice'][j,i]
    fice_in=nc_file['fice'][j,i]
    tisfc_in=nc_file['tisfc'][j,i]
    snwdph_in=nc_file['snwdph'][j,i]
    snoalb_in=nc_file['snoalb'][j,i]
    sheleg_in=nc_file['sheleg'][j,i]
    
    nc_file.close()
    
    #put data in a dictionary
    surface = {
        "T_surf": ts_in,
        "stc": stc_in,
        "smc": smc_in,
        "slc": slc_in,
        "tg3": tg3_in,
        "uustar": uustar_in,
        "alvsf": alvsf_in,
        "alvwf": alvwf_in,
        "alnsf": alnsf_in,
        "alnwf": alnwf_in,
        "facsf": facsf_in,
        "facwf": facwf_in,
        "styp": styp_in,
        "slope": slope_in,
        "vtyp": vtyp_in,
        "vfrac": vfrac_in,
        "shdmin": shdmin_in,
        "shdmax": shdmax_in,
        "zorl": zorl_in,
        "slmsk": slmsk_in,
        "canopy": canopy_in,
        "hice": hice_in,
        "fice": fice_in,
        "tisfc": tisfc_in,
        "snwdph": snwdph_in,
        "snoalb": snoalb_in,
        "sheleg": sheleg_in
    }
    return surface

def get_UFS_oro_data(dir, tile, i, j):
    
    filename_pattern = 'oro_data.tile{0}.nc'.format(tile)
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    # orographyic properties
    stddev_in=nc_file['stddev'][j,i]
    convexity_in=nc_file['convexity'][j,i]
    oa1_in=nc_file['oa1'][j,i]
    oa2_in=nc_file['oa2'][j,i]
    oa3_in=nc_file['oa3'][j,i]
    oa4_in=nc_file['oa4'][j,i]
    ol1_in=nc_file['ol1'][j,i]
    ol2_in=nc_file['ol2'][j,i]
    ol3_in=nc_file['ol3'][j,i]
    ol4_in=nc_file['ol4'][j,i]
    theta_in=nc_file['theta'][j,i]
    gamma_in=nc_file['gamma'][j,i]
    sigma_in=nc_file['sigma'][j,i]
    elvmax_in=nc_file['elvmax'][j,i]
    
    nc_file.close()
    
    #put data in a dictionary
    oro = {
        "stddev": stddev_in,
        "convexity": convexity_in,
        "oa1": oa1_in,
        "oa2": oa2_in,
        "oa3": oa3_in,
        "oa4": oa4_in,
        "ol1": ol1_in,
        "ol2": ol2_in,
        "ol3": ol3_in,
        "ol4": ol4_in,
        "theta": theta_in,
        "gamma": gamma_in,
        "sigma": sigma_in,
        "elvmax": elvmax_in
    }
    return oro

def get_UFS_vgrid_data(dir):
    
    nc_file = Dataset('{0}/{1}'.format(dir,'gfs_ctrl.nc'))
    
    # vertical coordinate definition
    ak=nc_file['vcoord'][0,::-1]
    bk=nc_file['vcoord'][1,::-1]
    
    nc_file.close()
    
    vgrid = {
        "ak": ak,
        "bk": bk
    }
    
    return vgrid    

def get_UFS_grid_area(dir, tile, i, j):
    filename_pattern = '*grid.tile{0}.nc'.format(tile)
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          file_name = f_name
    if not file_name:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,file_name))
    
    # extract out area of grid cell
    jpt2 = j*2+1
    ipt2 = i*2+1
    area_in=nc_file['area'][jpt2-1:jpt2+1,ipt2-1:ipt2+1]
    
    return area_in.sum()

def get_UFS_forcing_data(nlevs):
    
    ntimes = 1
        
    time = np.zeros(ntimes)
    w_ls = np.zeros((nlevs,ntimes),dtype=float)
    omega = np.zeros((nlevs,ntimes),dtype=float)
    u_g = np.zeros((nlevs,ntimes),dtype=float)
    v_g = np.zeros((nlevs,ntimes),dtype=float)
    u_nudge = np.zeros((nlevs,ntimes),dtype=float)
    v_nudge = np.zeros((nlevs,ntimes),dtype=float)
    T_nudge = np.zeros((nlevs,ntimes),dtype=float)
    thil_nudge = np.zeros((nlevs,ntimes),dtype=float)
    qt_nudge = np.zeros((nlevs,ntimes),dtype=float)
    rad_heating = np.zeros((nlevs,ntimes),dtype=float)
    h_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
    v_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
    h_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
    v_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
    
    forcing = {
        "time": time,
        "w_ls": w_ls,
        "omega": omega,
        "u_g": u_g,
        "v_g": v_g,
        "u_nudge": u_nudge,
        "v_nudge": v_nudge,
        "T_nudge": T_nudge,
        "thil_nudge": thil_nudge,
        "qt_nudge": qt_nudge,
        "rad_heating": rad_heating,
        "h_advec_thil": h_advec_thil,
        "v_advec_thil": v_advec_thil,
        "h_advec_qt": h_advec_qt,
        "v_advec_qt": v_advec_qt
    }
    
    return forcing

def add_noahmp_coldstart(surface, date):
    #use cold start section of FV3/io/FV3GFS_io.F90 to initialize NoahMP-specific variables
    print surface["slmsk"]
    
    mptable_nml_all = f90nml.read(os.path.join(NOAHMP_TABLES_DIR, 'MPTABLE.TBL'))
    mptable_nml_active = mptable_nml_all['noah_mp_modis_parameters']
    
    n_snow_layers = 3
    n_soil_layers = 4
    
    dzs = np.array([0.1,0.3,0.6,1.0])
    zsoil = np.array([-0.1,-0.4,-1.0,-2.0])
    
    surface["tvxy"]     = missing_value
    surface["tgxy"]     = missing_value
    surface["tahxy"]    = missing_value
    surface["canicexy"] = missing_value
    surface["canliqxy"] = missing_value
    surface["eahxy"]    = missing_value
    surface["cmxy"]     = missing_value
    surface["chxy"]     = missing_value
    surface["fwetxy"]   = missing_value
    surface["sneqvoxy"] = missing_value
    surface["alboldxy"] = missing_value
    surface["qsnowxy"]  = missing_value
    surface["wslakexy"] = missing_value
    surface["taussxy"]  = missing_value
    surface["waxy"]     = missing_value
    surface["wtxy"]     = missing_value
    surface["zwtxy"]    = missing_value
    surface["xlaixy"]   = missing_value
    surface["xsaixy"]   = missing_value

    surface["lfmassxy"] = missing_value
    surface["stmassxy"] = missing_value
    surface["rtmassxy"] = missing_value
    surface["woodxy"]   = missing_value
    surface["stblcpxy"] = missing_value
    surface["fastcpxy"] = missing_value
    surface["smcwtdxy"] = missing_value
    surface["deeprechxy"] = missing_value
    surface["rechxy"]   = missing_value

    surface["snowxy"]   = missing_value
    surface["snicexy"]  = np.ones(n_snow_layers)*missing_value
    surface["snliqxy"]  = np.ones(n_snow_layers)*missing_value
    surface["tsnoxy"]   = np.ones(n_snow_layers)*missing_value
    surface["smoiseq"]  = np.ones(n_soil_layers)*missing_value
    surface["zsnsoxy"]  = np.ones(n_snow_layers + n_soil_layers)*missing_value
    
    if surface["slmsk"] > 0.01:
        surface["tvxy"] = surface["T_surf"]
        surface["tgxy"] = surface["T_surf"]
        surface["tahxy"] = surface["T_surf"]
        
        if (surface["snwdph"] > 0.01 and surface["T_surf"] > 273.15 ):
            surface["tvxy"] = 273.15
            surface["tgxy"] = 273.15
            surface["tahxy"]= 273.15
            
        surface["canicexy"] = 0.0
        surface["canliqxy"] = surface["canopy"]
        surface["eahxy"]    = 2000.0
        
        #      eahxy = psfc*qv/(0.622+qv); qv is mixing ratio, converted from sepcific
        #      humidity specific humidity /(1.0 - specific humidity)
        
        surface["cmxy"]     = 0.0
        surface["chxy"]     = 0.0
        surface["fwetxy"]   = 0.0
        surface["sneqvoxy"] = surface["sheleg"]     # mm
        surface["alboldxy"] = 0.65
        surface["qsnowxy"]  = 0.0
        
        surface["wslakexy"] = 0.0
        surface["taussxy"]  = 0.0
        surface["waxy"]     = 4900.0
        surface["wtxy"]     = surface["waxy"]
        surface["zwtxy"]    = (25.0 + 2.0) - surface["waxy"] / 1000.0 /0.2
        
        vegtyp = np.int(surface['vtyp'])
        if (vegtyp == 0):
            vegtyp = 7
        if ((vegtyp == mptable_nml_active['ISBARREN']) or (vegtyp == mptable_nml_active['ISSNOW']) or  (vegtyp == mptable_nml_active['ISURBAN']) or (vegtyp == mptable_nml_active['ISWATER'])) :
            surface["xlaixy"] = 0.0
            surface["xsaixy"] = 0.0

            surface["lfmassxy"] = 0.0
            surface["stmassxy"] = 0.0
            surface["rtmassxy"] = 0.0

            surface["woodxy"] = 0.0       
            surface["stblcpxy"] = 0.0      
            surface["fastcpxy"] = 0.0
        else:
            laim = np.array(mptable_nml_active['LAIM']).reshape(12,20)
            surface["xlaixy"] = np.amax([laim[date["month"]-1,vegtyp-1],0.05])
            surface["xsaixy"] = np.amax([surface["xlaixy"]*0.1,0.05])
            
            sla = np.array(mptable_nml_active['SLA'])
            masslai = 1000.0 / np.amax([sla[vegtyp],1.0])
            surface["lfmassxy"] = surface["xlaixy"]*masslai
            masssai = 1000.0 / 3.0
            surface["stmassxy"] = surface["xsaixy"]*masssai
            
            surface["rtmassxy"] = 500.0      

            surface["woodxy"] = 500.0       
            surface["stblcpxy"] = 1000.0      
            surface["fastcpxy"] = 1000.0
            
        if ( vegtyp == mptable_nml_active['ISSNOW'] ):
            for k in range(n_soil_layers):
                surface["stc"][k] = np.amin([surface["stc"][k],np.amin([surface["tg3"],263.15])])
                surface["smc"][k] = 1
                surface["slc"][k] = 0
        
        snd = surface["snwdph"]/1000.0  # go to m from snwdph
        
        if (surface["sheleg"] != 0.0 and snd == 0.0 ):
            snd = surface["sheleg"]/1000.0
            
        if (vegtyp == 15):                       # land ice in MODIS/IGBP
            if ( surface["sheleg"] < 0.1):
                surface["sheleg"] = 0.1
                snd = 0.01
        
        dzsno = np.zeros(n_snow_layers)
        if (snd < 0.025 ):
            surface["snowxy"] = 0.0
            dzsno[:]          = 0.0
        elif (snd >= 0.025 and snd <= 0.05 ):
            surface["snowxy"] = -1.0
            dzsno[-1]         = snd
        elif (snd > 0.05 and snd <= 0.10 ):
            surface["snowxy"] = -2.0
            dzsno[-2] = 0.5*snd
            dzsno[-1] = 0.5*snd
        elif (snd > 0.10 and snd <= 0.25 ):
            surface["snowxy"] = -2.0
            dzsno[-2] = 0.05
            dzsno[-1] = snd - 0.05
        elif (snd > 0.25 and snd <= 0.45 ):
            surface["snowxy"] = -3.0
            dzsno[-3] = 0.05
            dzsno[-2] = 0.5*(snd-0.05)
            dzsno[-1] = 0.5*(snd-0.05)
        elif (snd > 0.45): 
            surface["snowxy"] = -3.0
            dzsno[-3] = 0.05
            dzsno[-2] = 0.20
            dzsno[-1] = snd - 0.05 - 0.20
        else:
            message = 'problem with the logic assigning snow layers.'
            logging.critical(message)
            raise Exception(message)
        
        surface["tsnoxy"][:]  = 0.0
        surface["snicexy"][:] = 0.0
        surface["snliqxy"][:] = 0.0
        surface["zsnsoxy"][:] = 0.0
        
        isnow = np.int(surface["snowxy"] + n_snow_layers)
        dzsnso = np.zeros(n_snow_layers + n_soil_layers)
        for k in range(isnow, n_snow_layers):
            surface["tsnoxy"][k]  = surface["tgxy"]
            surface["snliqxy"][k] = 0.0
            surface["snicexy"][k] = 1.00 * dzsno[k] * surface["sheleg"]/snd  #this line causes a warning
            
            dzsnso[k] = -dzsno[k]
        
        for k in range(n_snow_layers, n_snow_layers + n_soil_layers):
            dzsnso[k] = -dzs[k - n_snow_layers]
        
        surface["zsnsoxy"][isnow] = dzsnso[isnow]
        for k in range(isnow+1,n_snow_layers + n_soil_layers):
            surface["zsnsoxy"][k] = surface["zsnsoxy"][k-1] + dzsnso[k]
        
        soilparm = read_noahmp_soil_table()
        
        soiltyp  = int(surface["styp"])
        if (soiltyp != 0):
            index = soilparm["index"].index(soiltyp)
            bexp   = soilparm["bb"][index]
            smcmax = soilparm["maxsmc"][index]
            smcwlt = soilparm["wltsmc"][index]
            dwsat  = soilparm["satdw"][index]
            dksat  = soilparm["satdk"][index]
            psisat = -soilparm["satpsi"][index]
        
        if (vegtyp == mptable_nml_active['ISURBAN']):
            smcmax = 0.45
            smcwlt = 0.40
        
        if ((bexp > 0.0) and (smcmax > 0.0) and (-psisat > 0.0 )):
            for k in range(n_soil_layers):
                if ( k == 0 ):
                    ddz = -zsoil[k+1] * 0.5
                elif ( k < n_soil_layers-1 ):
                    ddz = (zsoil[k-1] - zsoil[k+1] ) * 0.5
                else:
                    ddz = zsoil[k-1] - zsoil[k]
# !
# ! Use newton-raphson method to find eq soil moisture
# !
                expon = bexp +1.
                aa = dwsat/ddz
                bb = dksat / smcmax ** expon

                smc = 0.5 * smcmax
                for iter in range(100):
                    func = (smc - smcmax) * aa +  bb * smc ** expon
                    dfunc = aa + bb * expon * smc ** bexp
                    dx  = func/dfunc
                    smc = smc - dx
                    if ( abs (dx) < 1.e-6):
                        break

                surface["smoiseq"][k] = np.amin([np.amax([smc,1.e-4]),smcmax*0.99])
        else:
            surface["smoiseq"][:] = smcmax

        surface["smcwtdxy"]   = smcmax
        surface["deeprechxy"] = 0.0
        surface["rechxy"]     = 0.0
        
    return surface

def read_noahmp_soil_table():
    
    choices = ["STAS","STAS-RUC"]
    
    with open(os.path.join(NOAHMP_TABLES_DIR, 'SOILPARM.TBL'), 'r') as f:
        lineList = f.readlines()
    f.close()
    
    line_index = 0
    for line in lineList:
        line_index += 1
        m = re.match(choices[0],line) and not re.match(choices[1],line)
        if m:
            start_index = line_index
            break
    
    n_soil_types = int(lineList[start_index].split()[0].split(',')[0])
    soil_index = []
    bb = []
    drysmc = []
    f11 = []
    maxsmc = []
    refsmc = []
    satpsi = []
    satdk = []
    satdw = []
    wltsmc = []
    qtz = []
    name = []
    for line in lineList[start_index+1:start_index+n_soil_types+1]:
        values = line.strip().split(',')
        soil_index.append(int(values[0]))
        bb.append(float(values[1]))
        drysmc.append(float(values[2]))
        f11.append(float(values[3]))
        maxsmc.append(float(values[4]))
        refsmc.append(float(values[5]))
        satpsi.append(float(values[6]))
        satdk.append(float(values[7]))
        satdw.append(float(values[8]))
        wltsmc.append(float(values[9]))
        qtz.append(float(values[10]))
        name.append(values[11].strip())
    
    soilparm = {
        "index": soil_index,
        "bb": bb,
        "drysmc": drysmc,
        "f11": f11,
        "maxsmc": maxsmc,
        "refsmc": refsmc,
        "satpsi": satpsi,
        "satdk": satdk,
        "satdw": satdw,
        "wltsmc": wltsmc,
        "qtz": qtz,
        "name": name        
    }
    
    return soilparm

def write_SCM_case_file(state, surface, oro, forcing, case, date):
    nlevs = state["nlevs"]
    nsoil = len(surface["stc"])
    
    noahmp = False
    if "snicexy" in surface.keys():
        noahmp = True
        nsnow = len(surface["snicexy"])
    
    nc_file = Dataset(os.path.join(PROCESSED_CASE_DIR, case + '.nc'), 'w', format='NETCDF4')
    nc_file.description = "FV3GFS model profile input (no forcing)"
    
    #create groups for scalars, intitialization, and forcing

    scalar_grp = nc_file.createGroup("scalars")
    initial_grp = nc_file.createGroup("initial")
    forcing_grp = nc_file.createGroup("forcing")
    
    #create dimensions and write them out

    time_dim = nc_file.createDimension('time', None)
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var[:] = forcing["time"]
    time_var.units = 's'
    time_var.description = 'elapsed time since the beginning of the simulation'

    levels_dim = nc_file.createDimension('levels', None)
    levels_var = nc_file.createVariable('levels', 'f4', ('levels',))
    levels_var[:] = state["pres"]
    levels_var.units = 'Pa'
    levels_var.description = 'pressure levels'
    
    soil_dim  = nc_file.createDimension('nsoil',None)
    soil_depth_var = nc_file.createVariable('soil_depth', 'f4', ('nsoil',))
    soil_depth_var[:] = [0.1,0.4,1.0,2.0]
    soil_depth_var.units = 'm'
    soil_depth_var.description = 'depth of bottom of soil layers'
    
    if noahmp:
        snow_dim = nc_file.createDimension('nsnow',None)
        soil_plus_snow_dim = nc_file.createDimension('nsoil_plus_nsnow',None)
        
    #initial group

    temperature_var = initial_grp.createVariable('temp', 'f4', ('levels',))
    temperature_var[:] = state["T"][0:nlevs]
    temperature_var.units = 'K'
    temperature_var.description = 'initial profile of absolute temperature'

    qt_var = initial_grp.createVariable('qt', 'f4', ('levels',))
    qt_var[:] = state["qv"][0:nlevs]
    qt_var.units = 'kg kg^-1'
    qt_var.description = 'initial profile of total water specific humidity'

    ql_var = initial_grp.createVariable('ql', 'f4', ('levels',))
    ql_var[:] = state["ql"][0:nlevs]
    ql_var.units = 'kg kg^-1'
    ql_var.description = 'initial profile of liquid water specific humidity'

    qi_var = initial_grp.createVariable('qi', 'f4', ('levels',))
    qi_var[:] = 0.0
    qi_var.units = 'kg kg^-1'
    qi_var.description = 'initial profile of ice water specific humidity'

    u_var = initial_grp.createVariable('u', 'f4', ('levels',))
    u_var[:] = state["u"][0:nlevs]
    u_var.units = 'm s^-1'
    u_var.description = 'initial profile of E-W horizontal wind'

    v_var = initial_grp.createVariable('v', 'f4', ('levels',))
    v_var[:] = state["v"][0:nlevs]
    v_var.units = 'm s^-1'
    v_var.description = 'initial profile of N-S horizontal wind'

    tke_var = initial_grp.createVariable('tke', 'f4', ('levels',))
    tke_var[:] = 0.0
    tke_var.units = 'm^2 s^-2'
    tke_var.description = 'initial profile of turbulence kinetic energy'

    ozone_var = initial_grp.createVariable('ozone', 'f4', ('levels',))
    ozone_var[:] = state["o3"][0:nlevs]
    ozone_var.units = 'kg kg^-1'
    ozone_var.description = 'initial profile of ozone mass mixing ratio'
    
    stc_var = initial_grp.createVariable('stc',np.float64,('nsoil',))
    stc_var[:] = surface['stc'][0:nsoil]
    stc_var.units = "K"
    stc_var.description = "initial profile of soil temperature"
    
    smc_var = initial_grp.createVariable('smc',np.float64,('nsoil',))
    smc_var[:] = surface['smc'][0:nsoil]
    smc_var.units = "kg"
    smc_var.description = "initial profile of soil moisture"
    
    slc_var = initial_grp.createVariable('slc',np.float64,('nsoil',))
    slc_var[:] = surface['slc'][0:nsoil]
    slc_var.units = "kg"
    slc_var.description = "initial profile of soil liquid moisture"
    
    if noahmp:
        snicexy_var = initial_grp.createVariable('snicexy',np.float64,('nsnow',))
        snicexy_var[:] = surface['snicexy'][0:nsnow]
        snicexy_var.units = ""
        snicexy_var.description = ""
        
        snliqxy_var = initial_grp.createVariable('snliqxy',np.float64,('nsnow',))
        snliqxy_var[:] = surface['snliqxy'][0:nsnow]
        snliqxy_var.units = ""
        snliqxy_var.description = ""
        
        tsnoxy_var = initial_grp.createVariable('tsnoxy',np.float64,('nsnow',))
        tsnoxy_var[:] = surface['tsnoxy'][0:nsnow]
        tsnoxy_var.units = ""
        tsnoxy_var.description = ""
        
        smoiseq_var = initial_grp.createVariable('smoiseq',np.float64,('nsoil',))
        smoiseq_var[:] = surface['smoiseq'][0:nsoil]
        smoiseq_var.units = ""
        smoiseq_var.description = ""
        
        zsnsoxy_var = initial_grp.createVariable('zsnsoxy',np.float64,('nsoil_plus_nsnow',))
        zsnsoxy_var[:] = surface['zsnsoxy'][0:nsoil + nsnow]
        zsnsoxy_var.units = ""
        zsnsoxy_var.description = ""
    
    #forcing group

    p_surf_var = forcing_grp.createVariable('p_surf', 'f4', ('time',))
    p_surf_var[:] = state["p_surf"]
    p_surf_var.units = 'Pa'
    p_surf_var.description = 'surface pressure'

    T_surf_var = forcing_grp.createVariable('T_surf', 'f4', ('time',))
    T_surf_var[:] = surface["T_surf"]
    T_surf_var.units = 'K'
    T_surf_var.description = 'surface absolute temperature'

    w_ls_var = forcing_grp.createVariable('w_ls', 'f4', ('levels','time',))
    w_ls_var[:] = forcing["w_ls"]
    w_ls_var.units = 'm s^-1'
    w_ls_var.description = 'large scale vertical velocity'
    
    omega_var = forcing_grp.createVariable('omega', 'f4', ('levels','time',))
    omega_var[:] = forcing["omega"]
    omega_var.units = 'Pa s^-1'
    omega_var.description = 'large scale pressure vertical velocity'
    
    u_g_var = forcing_grp.createVariable('u_g', 'f4', ('levels','time',))
    u_g_var[:] = forcing["u_g"]
    u_g_var.units = 'm s^-1'
    u_g_var.description = 'large scale geostrophic E-W wind'
    
    v_g_var = forcing_grp.createVariable('v_g', 'f4', ('levels','time',))
    v_g_var[:] = forcing["v_g"]
    v_g_var.units = 'm s^-1'
    v_g_var.description = 'large scale geostrophic N-S wind'
    
    u_nudge_var = forcing_grp.createVariable('u_nudge', 'f4', ('levels','time',))
    u_nudge_var[:] = forcing["u_nudge"]
    u_nudge_var.units = 'm s^-1'
    u_nudge_var.description = 'E-W wind to nudge toward'
    
    v_nudge_var = forcing_grp.createVariable('v_nudge', 'f4', ('levels','time',))
    v_nudge_var[:] = forcing["v_nudge"]
    v_nudge_var.units = 'm s^-1'
    v_nudge_var.description = 'N-S wind to nudge toward'
    
    T_nudge_var = forcing_grp.createVariable('T_nudge', 'f4', ('levels','time',))
    T_nudge_var[:] = forcing["T_nudge"]
    T_nudge_var.units = 'K'
    T_nudge_var.description = 'absolute temperature to nudge toward'
     
    thil_nudge_var = forcing_grp.createVariable('thil_nudge', 'f4', ('levels','time',))
    thil_nudge_var[:] = forcing["thil_nudge"]
    thil_nudge_var.units = 'K'
    thil_nudge_var.description = 'potential temperature to nudge toward'
    
    qt_nudge_var = forcing_grp.createVariable('qt_nudge', 'f4', ('levels','time',))
    qt_nudge_var[:] = forcing["qt_nudge"]
    qt_nudge_var.units = 'kg kg^-1'
    qt_nudge_var.description = 'q_t to nudge toward'
    
    rad_heating_var = forcing_grp.createVariable('dT_dt_rad', 'f4', ('levels','time',))
    rad_heating_var[:] = forcing["rad_heating"]
    rad_heating_var.units = 'K s^-1'
    rad_heating_var.description = 'prescribed radiative heating rate'
    
    h_advec_thil_var = forcing_grp.createVariable('h_advec_thetail', 'f4', ('levels','time',))
    h_advec_thil_var[:] = forcing["h_advec_thil"]
    h_advec_thil_var.units = 'K s^-1'
    h_advec_thil_var.description = 'prescribed theta_il tendency due to horizontal advection'
    
    v_advec_thil_var = forcing_grp.createVariable('v_advec_thetail', 'f4', ('levels','time',))
    v_advec_thil_var[:] = forcing["v_advec_thil"]
    v_advec_thil_var.units = 'K s^-1'
    v_advec_thil_var.description = 'prescribed theta_il tendency due to vertical advection'
    
    h_advec_qt_var = forcing_grp.createVariable('h_advec_qt', 'f4', ('levels','time',))
    h_advec_qt_var[:] = forcing["h_advec_qt"]
    h_advec_qt_var.units = 'kg kg^-1 s^-1'
    h_advec_qt_var.description = 'prescribed q_t tendency due to horizontal advection'
    
    v_advec_qt_var = forcing_grp.createVariable('v_advec_qt', 'f4', ('levels','time',))
    v_advec_qt_var[:] = forcing["v_advec_qt"]
    v_advec_qt_var.units = 'kg kg^-1 s^-1'
    v_advec_qt_var.description = 'prescribed q_t tendency due to vertical advection'
    
    #scalar group
    year_var = scalar_grp.createVariable('init_year',np.int32)
    year_var[:] = date["year"]
    year_var.units = "years"
    year_var.description = "year at time of initial values"
    
    month_var = scalar_grp.createVariable('init_month',np.int32)
    month_var[:] = date["month"]
    month_var.units = "months"
    month_var.description = "month at time of initial values"
    
    day_var = scalar_grp.createVariable('init_day',np.int32)
    day_var[:] = date["day"]
    day_var.units = "days"
    day_var.description = "day at time of initial values"
    
    hour_var = scalar_grp.createVariable('init_hour',np.int32)
    hour_var[:] = date["hour"]
    hour_var.units = "hours"
    hour_var.description = "hour at time of initial values"
    
    minute_var = scalar_grp.createVariable('init_minute',np.int32)
    minute_var[:] = date["minute"]
    minute_var.units = "minutes"
    minute_var.description = "minute at time of initial values"
    
    second_var = scalar_grp.createVariable('init_second',np.int32)
    second_var[:] = 0.0
    second_var.units = "seconds"
    second_var.description = "second at time of initial values"
    
    lat_var = scalar_grp.createVariable('lat', 'f4')
    lat_var[:] = surface["lat"]
    lat_var.units = 'degrees N'
    lat_var.description = 'latitude of column'

    lon_var = scalar_grp.createVariable('lon', 'f4')
    lon_var[:] = surface["lon"]
    lon_var.units = 'degrees E'
    lon_var.description = 'longitude of column'
    
    area = scalar_grp.createVariable('area',np.float64)
    area[:] = surface["area"]
    area.units = "m^2" 
    area.description = "grid cell area"
    
    #Noah initial parameters
    
    vegsrc  = scalar_grp.createVariable('vegsrc',np.int32)
    vegsrc[:] = 1 #when would this be 2?
    vegsrc.description = "vegetation soure (1-2)"
    
    vegtyp  = scalar_grp.createVariable('vegtyp',np.int32)
    vegtyp[:] = surface["vtyp"]
    vegtyp.description = "vegetation type (1-12)"

    soiltyp = scalar_grp.createVariable('soiltyp',np.int32)
    soiltyp[:] = surface["styp"]
    soiltyp.description = "soil type (1-12)"
    
    slopetyp = scalar_grp.createVariable('slopetyp',np.int32)
    slopetyp[:] = surface["slope"]
    slopetyp.description = "slope type (1-9)"
    
    vegfrac = scalar_grp.createVariable('vegfrac',np.float)
    vegfrac[:] = surface["vfrac"]
    vegfrac.description = "vegetation fraction"
    
    shdmin = scalar_grp.createVariable('shdmin',np.float)
    shdmin[:] = surface["shdmin"]
    shdmin.description = "minimum vegetation fraction"
    
    shdmax = scalar_grp.createVariable('shdmax',np.float)
    shdmax[:] = surface["shdmax"]
    shdmax.description = "maximum vegetation fraction"
    
    zorl = scalar_grp.createVariable('zorl',np.float)
    zorl[:] = surface["zorl"]
    zorl.units = ""
    zorl.description = "surface roughness length"
    
    islmsk = scalar_grp.createVariable('slmsk',np.float)
    islmsk[:] = surface["slmsk"]
    islmsk.description = "land-sea-ice mask"
    
    canopy = scalar_grp.createVariable('canopy',np.float)
    canopy[:] = surface["canopy"]
    canopy.units = ""
    canopy.description = "canopy moisture"
    
    hice = scalar_grp.createVariable('hice',np.float)
    hice[:] = surface["hice"]
    hice.units = ""
    hice.description = "ice thickness"
    
    fice = scalar_grp.createVariable('fice',np.float)
    fice[:] = surface["fice"]
    fice.units = ""
    fice.description = "ice fraction"
    
    tisfc = scalar_grp.createVariable('tisfc',np.float)
    tisfc[:] = surface["tisfc"]
    tisfc.units = ""
    tisfc.description = "ice temperature"
    
    snwdph = scalar_grp.createVariable('snwdph',np.float)
    snwdph[:] = surface["snwdph"]
    snwdph.units = ""
    snwdph.description = "snow depth"
    
    snoalb = scalar_grp.createVariable('snoalb',np.float)
    snoalb[:] = surface["snoalb"]
    snoalb.units = ""
    snoalb.description = "snow albedo"
    
    sncovr = scalar_grp.createVariable('sncovr',np.float)
    sncovr[:] = 0.0
    sncovr.units = ""
    sncovr.description = "snow cover"
        
    tg3 = scalar_grp.createVariable('tg3',np.float)
    tg3[:] = surface["tg3"]
    tg3.units = "K"  
    tg3.description = "deep soil temperature"
    
    uustar = scalar_grp.createVariable('uustar',np.float)
    uustar[:] = surface["uustar"]
    uustar.units = "m2s-2?"  
    uustar.description = "friction velocity"
    
    alvsf = scalar_grp.createVariable('alvsf',np.float64)
    alvsf[:] = surface["alvsf"]
    alvsf.units = "None" 
    alvsf.description = "uv+visible black sky albedo (z=60 degree)"
    
    alnsf = scalar_grp.createVariable('alnsf',np.float64)
    alnsf[:] = surface["alnsf"]
    alnsf.units = "None"
    alnsf.description = "near IR black sky albedo (z=60 degree)"
    
    alvwf = scalar_grp.createVariable('alvwf',np.float64)
    alvwf[:] = surface["alvwf"]
    alvwf.units = "None"
    alvwf.description = "uv+visible white sky albedo"
    
    alnwf = scalar_grp.createVariable('alnwf',np.float64)
    alnwf[:] = surface["alnwf"]
    alnwf.units = "None"
    alnwf.description = "near IR white sky albedo"
    
    facsf = scalar_grp.createVariable('facsf',np.float64)
    facsf[:] = surface["facsf"]
    facsf.units = "None" 
    facsf.description = "fraction of grid cell with strong sun angle albedo dependence"
    
    facwf = scalar_grp.createVariable('facwf',np.float64)
    facwf[:] = surface["facwf"]
    facwf.units = "None" 
    facwf.description = "fraction of grid cell with weak sun angle albedo dependence"
    
    #Orography initial parameters
    
    stddev = scalar_grp.createVariable('stddev',np.float64)
    stddev[:] = oro["stddev"]
    stddev.units = ""
    stddev.description = "surface orography standard deviation"
    
    convexity = scalar_grp.createVariable('convexity',np.float64)
    convexity[:] = oro["convexity"]
    convexity.units = ""
    convexity.description = ""
    
    oa1 = scalar_grp.createVariable('oa1',np.float64)
    oa1[:] = oro["oa1"]
    oa1.units = ""
    oa1.description = ""
    
    oa2 = scalar_grp.createVariable('oa2',np.float64)
    oa2[:] = oro["oa2"]
    oa2.units = ""
    oa2.description = ""
    
    oa3 = scalar_grp.createVariable('oa3',np.float64)
    oa3[:] = oro["oa3"]
    oa3.units = ""
    oa3.description = ""
    
    oa4 = scalar_grp.createVariable('oa4',np.float64)
    oa4[:] = oro["oa4"]
    oa4.units = ""
    oa4.description = ""
    
    ol1 = scalar_grp.createVariable('ol1',np.float64)
    ol1[:] = oro["ol1"]
    ol1.units = ""
    ol1.description = ""
    
    ol2 = scalar_grp.createVariable('ol2',np.float64)
    ol2[:] = oro["ol2"]
    ol2.units = ""
    ol2.description = ""
    
    ol3 = scalar_grp.createVariable('ol3',np.float64)
    ol3[:] = oro["ol3"]
    ol3.units = ""
    ol3.description = ""
    
    ol4 = scalar_grp.createVariable('ol4',np.float64)
    ol4[:] = oro["ol4"]
    ol4.units = ""
    ol4.description = ""
    
    theta = scalar_grp.createVariable('theta',np.float64)
    theta[:] = oro["theta"]
    theta.units = ""
    theta.description = ""
    
    gamma = scalar_grp.createVariable('gamma',np.float64)
    gamma[:] = oro["gamma"]
    gamma.units = ""
    gamma.description = ""
    
    sigma = scalar_grp.createVariable('sigma',np.float64)
    sigma[:] = oro["sigma"]
    sigma.units = ""
    sigma.description = ""
    
    elvmax = scalar_grp.createVariable('elvmax',np.float64)
    elvmax[:] = oro["elvmax"]
    elvmax.units = ""
    elvmax.description = ""
    
    #NoahMP initial scalar parameters
    if noahmp:
        tvxy = scalar_grp.createVariable('tvxy',np.float64)
        tvxy[:] = surface["tvxy"]
        tvxy.units = ""
        tvxy.description = ""
        
        tgxy = scalar_grp.createVariable('tgxy',np.float64)
        tgxy[:] = surface["tgxy"]
        tgxy.units = ""
        tgxy.description = ""
        
        tahxy = scalar_grp.createVariable('tahxy',np.float64)
        tahxy[:] = surface["tahxy"]
        tahxy.units = ""
        tahxy.description = ""
        
        canicexy = scalar_grp.createVariable('canicexy',np.float64)
        canicexy[:] = surface["canicexy"]
        canicexy.units = ""
        canicexy.description = ""
        
        canliqxy = scalar_grp.createVariable('canliqxy',np.float64)
        canliqxy[:] = surface["canliqxy"]
        canliqxy.units = ""
        canliqxy.description = ""
        
        eahxy = scalar_grp.createVariable('eahxy',np.float64)
        eahxy[:] = surface["eahxy"]
        eahxy.units = ""
        eahxy.description = ""
        
        cmxy = scalar_grp.createVariable('cmxy',np.float64)
        cmxy[:] = surface["cmxy"]
        cmxy.units = ""
        cmxy.description = ""
        
        chxy = scalar_grp.createVariable('chxy',np.float64)
        chxy[:] = surface["chxy"]
        chxy.units = ""
        chxy.description = ""
    
        fwetxy = scalar_grp.createVariable('fwetxy',np.float64)
        fwetxy[:] = surface["fwetxy"]
        fwetxy.units = ""
        fwetxy.description = ""
        
        sneqvoxy = scalar_grp.createVariable('sneqvoxy',np.float64)
        sneqvoxy[:] = surface["sneqvoxy"]
        sneqvoxy.units = ""
        sneqvoxy.description = ""
        
        alboldxy = scalar_grp.createVariable('alboldxy',np.float64)
        alboldxy[:] = surface["alboldxy"]
        alboldxy.units = ""
        alboldxy.description = ""
        
        qsnowxy = scalar_grp.createVariable('qsnowxy',np.float64)
        qsnowxy[:] = surface["qsnowxy"]
        qsnowxy.units = ""
        qsnowxy.description = ""
        
        wslakexy = scalar_grp.createVariable('wslakexy',np.float64)
        wslakexy[:] = surface["wslakexy"]
        wslakexy.units = ""
        wslakexy.description = ""
        
        taussxy = scalar_grp.createVariable('taussxy',np.float64)
        taussxy[:] = surface["taussxy"]
        taussxy.units = ""
        taussxy.description = ""
        
        waxy = scalar_grp.createVariable('waxy',np.float64)
        waxy[:] = surface["waxy"]
        waxy.units = ""
        waxy.description = ""
        
        wtxy = scalar_grp.createVariable('wtxy',np.float64)
        wtxy[:] = surface["wtxy"]
        wtxy.units = ""
        wtxy.description = ""
    
        zwtxy = scalar_grp.createVariable('zwtxy',np.float64)
        zwtxy[:] = surface["zwtxy"]
        zwtxy.units = ""
        zwtxy.description = ""
        
        xlaixy = scalar_grp.createVariable('xlaixy',np.float64)
        xlaixy[:] = surface["xlaixy"]
        xlaixy.units = ""
        xlaixy.description = ""
        
        xsaixy = scalar_grp.createVariable('xsaixy',np.float64)
        xsaixy[:] = surface["xsaixy"]
        xsaixy.units = ""
        xsaixy.description = ""

        lfmassxy = scalar_grp.createVariable('lfmassxy',np.float64)
        lfmassxy[:] = surface["lfmassxy"]
        lfmassxy.units = ""
        lfmassxy.description = ""
        
        stmassxy = scalar_grp.createVariable('stmassxy',np.float64)
        stmassxy[:] = surface["stmassxy"]
        stmassxy.units = ""
        stmassxy.description = ""
        
        rtmassxy = scalar_grp.createVariable('rtmassxy',np.float64)
        rtmassxy[:] = surface["rtmassxy"]
        rtmassxy.units = ""
        rtmassxy.description = ""
        
        woodxy = scalar_grp.createVariable('woodxy',np.float64)
        woodxy[:] = surface["woodxy"]
        woodxy.units = ""
        woodxy.description = ""
        
        stblcpxy = scalar_grp.createVariable('stblcpxy',np.float64)
        stblcpxy[:] = surface["stblcpxy"]
        stblcpxy.units = ""
        stblcpxy.description = ""
        
        fastcpxy = scalar_grp.createVariable('fastcpxy',np.float64)
        fastcpxy[:] = surface["fastcpxy"]
        fastcpxy.units = ""
        fastcpxy.description = ""
    
        smcwtdxy = scalar_grp.createVariable('smcwtdxy',np.float64)
        smcwtdxy[:] = surface["smcwtdxy"]
        smcwtdxy.units = ""
        smcwtdxy.description = ""
        
        deeprechxy = scalar_grp.createVariable('deeprechxy',np.float64)
        deeprechxy[:] = surface["deeprechxy"]
        deeprechxy.units = ""
        deeprechxy.description = ""
        
        rechxy = scalar_grp.createVariable('rechxy',np.float64)
        rechxy[:] = surface["rechxy"]
        rechxy.units = ""
        rechxy.description = ""
        
        snowxy = scalar_grp.createVariable('snowxy',np.float64)
        snowxy[:] = surface["snowxy"]
        snowxy.units = ""
        snowxy.description = ""
    
    nc_file.close()

def main():
    setup_logging()
    
    #read in arguments
    (location, indices, date, in_dir, grid_dir, tile, area, case_name) = parse_arguments()
        
    #find tile containing the point using the supergrid if no tile is specified 
    if not tile:
        tile = find_tile(location, grid_dir)
        if tile < 0:
            message = 'No tile was found for location {0}'.format(location)
            logging.critical(message)
            raise Exception(message)
        print 'Tile found: {0}'.format(tile)
    
    #find index of closest point in the tile if indices are not specified
    if not indices:
        (tile_j, tile_i, point_lon, point_lat, dist_min) = find_loc_indices(location, in_dir, tile)
        print 'The closest point in tile {0} has indices [{1},{2}]'.format(tile,tile_i,tile_j)
        print 'This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat)
        print 'This grid cell is approximately {0} km away from the desired location of {1} {2}'.format(dist_min/1.0E3,location[0],location[1])
    else:
        tile_i = indices[0]
        tile_j = indices[1]
        (point_lon, point_lat) = find_lon_lat_of_indices(indices, in_dir, tile)
        
        print 'This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat)
    
    #get UFS IC data (or restart data!)
    (state_data, surface_data, oro_data) = get_UFS_IC_data(in_dir, tile, tile_i, tile_j)
    
    #cold start NoahMP variables
    noahmp = False
    if (noahmp):
        surface_data = add_noahmp_coldstart(surface_data, date)
    
    #get grid cell area if not given
    if not area:
        area = get_UFS_grid_area(grid_dir, tile, tile_i, tile_j)
    surface_data["area"] = area
    
    surface_data["lon"] = point_lon
    surface_data["lat"] = point_lat
        
    #get UFS forcing data (zeros for now)
    forcing_data = get_UFS_forcing_data(state_data["nlevs"])
    
    #write SCM case file
    write_SCM_case_file(state_data, surface_data, oro_data, forcing_data, case_name, date)
    
    
if __name__ == '__main__':
    main()
