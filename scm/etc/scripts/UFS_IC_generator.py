#!/usr/bin/env python

import argparse
import os
import fnmatch
import logging
from netCDF4 import Dataset
import numpy as np
from shapely.geometry import Point, Polygon
import copy
import math
import f90nml
import re
import fv3_remap
import xesmf
from datetime import datetime, timedelta

###############################################################################
# Global settings                                                             #
###############################################################################

#Physical constants
earth_radius = 6371000.0 #m
rdgas        = 287.05
rvgas        = 461.50
cp           = 1004.6
zvir         = rvgas/rdgas - 1.
rocp         = rdgas/cp
grav         = 9.80665
deg_to_rad   = math.pi/180.0
kappa        = rdgas/cp
p0           = 100000.0

missing_value = -9999.0 #9.99e20

n_lam_halo_points = 3

missing_variable_snow_layers = 3
missing_variable_soil_layers = 4
missing_variable_ice_layers = 2

# Path to the directory containing processed case input files
PROCESSED_CASE_DIR = '../../data/processed_case_input'

# Path to the directory containing comparison data files
COMPARISON_DATA_DIR = '../../data/comparison_data'

# Path to the directory containing NoahMP table files (need MPTABLE.TBL and SOILPARM.TBL)
NOAHMP_TABLES_DIR = '../../data/raw_case_input/NoahMP_tables'

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-l',   '--location',   help='longitude and latitude in degress E and N, respectively, separated by a space', nargs=2, type=float)
group1.add_argument('-ij',  '--index',      help='i,j indices within the tile (if known - bypasses search for closest model point to lon/lat location)', nargs=2, type=int)
parser.add_argument('-d',   '--date',       help='date corresponding to initial conditions in YYYYMMDDHHMMSS format', required=False)
parser.add_argument('-i',   '--in_dir',     help='input directory path containing FV3 input files', required=True)
parser.add_argument('-g',   '--grid_dir',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-f',   '--forcing_dir',help='directory path containing physics diag files')
parser.add_argument('-t',   '--tile',       help='tile of desired point (if known - bypasses tile search if present)', type=int, choices=range(1,8))
parser.add_argument('-a',   '--area',       help='area of grid cell in m^2', type=float)
parser.add_argument('-mp',  '--noahmp',     help='flag to generate cold-start ICs for NoahMP LSM from Noah LSM ICs', action='store_true')
parser.add_argument('-lam', '--lam',        help='flag to signal that the ICs and forcing is from a limited-area model run', action='store_true')
parser.add_argument('-lami','--lam_ic',     help='flag to signal that the ICs are from a limited-area model run', action='store_true')
parser.add_argument('-lamf','--lam_frc',    help='flag to signal that the forcings are from a limited-area model run', action='store_true')
parser.add_argument('-n',   '--case_name',  help='name of case', required=True)
parser.add_argument('-oc',  '--old_chgres', help='flag to denote that the initial conditions use an older data format (pre-chgres_cube)', action='store_true')
parser.add_argument('-sc',  '--save_comp',       help='flag to save and write out a file with UFS output data to compare SCM simulations with', action='store_true')
parser.add_argument('-lsm', '--add_UFS_NOAH_lsm', help='flag to include UFS NOAH LSM surface forcing', action='store_true')
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
    forcing_dir = args.forcing_dir
    tile = args.tile
    area = args.area
    case_name = args.case_name
    noahmp = args.noahmp
    old_chgres = args.old_chgres
    lam = args.lam
    lami = args.lam_ic
    lamf = args.lam_frc
    save_comp = args.save_comp
    add_UFS_NOAH_lsm = args.add_UFS_NOAH_lsm

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
    if date:
        if len(date) != 14:
            message = 'The entered date {0} does not have the 14 characters expected in the format YYYYMMDDHHMMSS'.format(date)
            logging.critical(message)
            raise Exception(message)
        else:
            date_dict["year"] = int(date[0:4])
            date_dict["month"] = int(date[4:6])
            date_dict["day"] = int(date[6:8])
            date_dict["hour"] = int(date[8:10])
            date_dict["minute"] = int(date[10:12])
            date_dict["second"] = int(date[12:])
    
    if tile:
        if (not lam and tile > 6):
            message = 'The entered tile {0} is not compatibile with the global cubed-sphere grid'.format(date)
            logging.critical(message)
            raise Exception(message)
    
    return (location, index, date_dict, in_dir, grid_dir, forcing_dir, tile, area, noahmp, case_name, old_chgres, lam, lami, lamf, save_comp, add_UFS_NOAH_lsm)

def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)
    
def find_tile(loc, dir, lam):
    """Find the FV3 tile with the given lon/lat"""
    #returns the integer tile number
    
    # should be looking in the directory with supergrid data (probably "fix" directory)
    filename_pattern = '*grid.tile*.nc'
    
    #find all supergrid files in the directory
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
            longitude = np.asarray(nc_file['x']).swapaxes(0,1)
            latitude = np.asarray(nc_file['y']).swapaxes(0,1)
            nc_file.close()
            
            adj_long = False        
            #look for reversal of longitude; if found, adjust longitude so that 0-360 transition doesn't exist
            for row in longitude:
                if not (np.all(np.diff(row) >= 0) or np.all(np.diff(row) <= 0)):
                    adj_long = True
            if adj_long:
                longitude[longitude < 180] += 360
            
            #get lon/lat pairs for all edges of the tiles
            
            edge_1_lon = longitude[0,:]
            edge_1_lat = latitude[0,:]
            edge_1 = list(zip(edge_1_lon, edge_1_lat))
                        
            edge_2_lon = longitude[:,-1]
            edge_2_lat = latitude[:,-1]
            edge_2 = list(zip(edge_2_lon, edge_2_lat))
                        
            edge_3_lon = longitude[-1,:]
            edge_3_lat = latitude[-1,:]
            edge_3 = list(zip(edge_3_lon, edge_3_lat))
            edge_3.reverse() #need to reverse the direction of this edge to form a regular polygon
            
            edge_4_lon = longitude[:,0]
            edge_4_lat = latitude[:,0]
            edge_4 = list(zip(edge_4_lon, edge_4_lat))
            edge_4.reverse() #need to reverse the direction of this edge to form a regular polygon
                        
            polygon_points = edge_1 + edge_2 + edge_3 + edge_4
                        
            tile_polygon = Polygon(polygon_points)
            tile_polygon = tile_polygon.simplify(0)
            
            if tile_polygon.is_valid:  #this will be True unless the tile is a polar tile, which will not form a regular polygon in Cartesian space using lon/lat data
                temp_loc = copy.deepcopy(loc)
                if adj_long:
                    if loc[0] < 180:
                        temp_loc[0] += 360
                loc_point = Point(temp_loc)
                if tile_polygon.contains(loc_point):
                    found_tile = True
                    if (lam):
                        return f_name.split('tile')[1].split('.halo')[0] 
                    else:
                        return f_name.split('tile')[1].split('.nc')[0] 
            else:
                polar_tile_filenames.append(f_name)
                
    #if the tile hasn't been found by this point, it must be contained within a polar tile
    for f_name in polar_tile_filenames:
        nc_file = Dataset('{0}/{1}'.format(dir,f_name))
        latitude = np.asarray(nc_file['y']).swapaxes(0,1)
        nc_file.close()
        
        #if the sign of the mean latitude of the tile is the same as that of the point, the tile has been found
        if np.sign(np.mean(latitude)) == np.sign(loc[1]):
            found_tile = True
            if (lam):
                return f_name.split('tile')[1].split('.halo')[0] 
            else:
                return f_name.split('tile')[1].split('.nc')[0]        
    return -1

def find_loc_indices(loc, dir, tile, lam):
    """Find the nearest neighbor FV3 grid point given a lon/lat pair and the tile number"""
    #returns the indices of the nearest neighbor point in the given tile, the lon/lat of the nearest neighbor, 
    #and the distance (m) from the given point to the nearest neighbor grid cell
    
    if (lam):
        filename_pattern = '*grid.tile7.halo{}.nc'.format(n_lam_halo_points)
    else:
        filename_pattern = '*grid.tile{0}.nc'.format(tile)
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    #read in supergrid longitude and latitude
    lon_super = np.asarray(nc_file['x'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    lat_super = np.asarray(nc_file['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    if (lam):
        #strip ghost/halo points and return central (A-grid) points
        #assuming n_lam_halo_points
        lon_super_no_halo = lon_super[2*n_lam_halo_points:lon_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lon_super.shape[1]-2*n_lam_halo_points]
        lat_super_no_halo = lat_super[2*n_lam_halo_points:lat_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lat_super.shape[1]-2*n_lam_halo_points]
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super_no_halo[1::2,1::2]
        latitude = lat_super_no_halo[1::2,1::2]
    else:
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super[1::2,1::2]
        latitude = lat_super[1::2,1::2]
    
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
    
    #set up an array to hold the euclidean distance between the given point and every grid cell
    eucl_dist = np.zeros((longitude.shape[0],longitude.shape[1]))
    
    #get the Cartesian location of the given point
    cart_loc = np.asarray(sph2cart(math.radians(temp_loc[0]), math.radians(temp_loc[1]), earth_radius))
    
    for i in range(len(longitude)):
        for j in range(len(longitude[i])):
            #get the Cartesian location of all grid points
            cart_cell = np.asarray(sph2cart(math.radians(longitude[i,j]), math.radians(latitude[i,j]), earth_radius))
            
            #calculate the euclidean distance from the given point to the current grid cell
            eucl_dist[i,j] = np.linalg.norm(cart_loc - cart_cell)
    
    #get the indices of the grid point with the minimum euclidean distance to the given point
    i,j = np.unravel_index(eucl_dist.argmin(), eucl_dist.shape)
    
    return (i,j,longitude[i,j]%360.0, latitude[i,j], eucl_dist[i,j])

def find_lon_lat_of_indices(indices, dir, tile, lam):
    """Find the longitude and latitude of the given indices within the given tile."""
    
    if (lam):
        filename_pattern = '*grid.tile{0}.halo{1}.nc'.format(tile, n_lam_halo_points)
    else:
        filename_pattern = '*grid.tile{0}.nc'.format(tile)
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    #read in supergrid longitude and latitude
    lon_super = np.asarray(nc_file['x'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    lat_super = np.asarray(nc_file['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    if (lam):
        #strip ghost/halo points and return central (A-grid) points
        #assuming n_lam_halo_points
        lon_super_no_halo = lon_super[2*n_lam_halo_points:lon_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lon_super.shape[1]-2*n_lam_halo_points]
        lat_super_no_halo = lat_super[2*n_lam_halo_points:lat_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lat_super.shape[1]-2*n_lam_halo_points]
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super_no_halo[1::2,1::2]
        latitude = lat_super_no_halo[1::2,1::2]
    else:
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super[1::2,1::2]
        latitude = lat_super[1::2,1::2]
    
    nc_file.close()
    
    return (longitude[indices[1],indices[0]], latitude[indices[1],indices[0]])

def get_initial_lon_lat_grid(dir, tile, lam):
    if (lam):
        filename_pattern = '*grid.tile{0}.halo{1}.nc'.format(tile, n_lam_halo_points)
    else:
        filename_pattern = '*grid.tile{0}.nc'.format(tile)
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    #read in supergrid longitude and latitude
    lon_super = np.asarray(nc_file['x'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    lat_super = np.asarray(nc_file['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    if (lam):
        #strip ghost/halo points and return central (A-grid) points
        #assuming n_lam_halo_points
        lon_super_no_halo = lon_super[2*n_lam_halo_points:lon_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lon_super.shape[1]-2*n_lam_halo_points]
        lat_super_no_halo = lat_super[2*n_lam_halo_points:lat_super.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lat_super.shape[1]-2*n_lam_halo_points]
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super_no_halo[1::2,1::2]
        latitude = lat_super_no_halo[1::2,1::2]
    else:
        #get the longitude and latitude data for the grid centers by slicing the supergrid 
        #and taking only odd-indexed values
        longitude = lon_super[1::2,1::2]
        latitude = lat_super[1::2,1::2]
    
    nc_file.close()
    
    return (longitude, latitude)

def sph2cart(az, el, r):
    """Calculate the Cartesian coordiates from spherical coordinates"""
    
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    
    return (x, y, z)    

def read_NetCDF_var(nc_file, var_name, i, j):
    try:
        var = nc_file[var_name][j,i]
    except (KeyError, IndexError):
        message = "Variable {0} is not found in {1}. Filling with missing value.".format(var_name,nc_file.filepath())
        logging.debug(message)
        var = missing_value
    return var

def read_NetCDF_surface_var(nc_file, var_name, i, j, old_chgres, vert_dim):
    if old_chgres:
        if vert_dim > 0:
            try:
                var = nc_file[var_name][:,j,i]
            except (KeyError, IndexError):
                message = "Variable {0} is not found in {1}. Filling with array of size {2} with missing value.".format(var_name,nc_file.filepath(),vert_dim)
                logging.debug(message)
                var = missing_value*np.ma.ones(vert_dim)
        else:
            try:
                var = nc_file[var_name][j,i]
            except (KeyError, IndexError):
                message = "Variable {0} is not found in {1}. Filling with missing value.".format(var_name,nc_file.filepath())
                logging.debug(message)
                var = missing_value
    else:
        #the sfc_data.tileX.nc files created from chgres_cube have an extra time dimension in front compared to those created from global_chgres
        if vert_dim > 0:
            try:
                var = nc_file[var_name][0,:,j,i]
            except (KeyError, IndexError):
                message = "Variable {0} is not found in {1}. Filling with array of size {2} with missing value.".format(var_name,nc_file.filepath(),vert_dim)
                logging.debug(message)
                var = missing_value*np.ma.ones(vert_dim)
        else:
            try:
                var = nc_file[var_name][0,j,i]
            except (KeyError, IndexError):
                message = "Variable {0} is not found in {1}. Filling with missing value.".format(var_name,nc_file.filepath())
                logging.debug(message)
                var = missing_value
    return var

def get_UFS_IC_data(dir, grid_dir, forcing_dir, tile, i, j, old_chgres, lam):
    """Get the state, surface, and orographic data for the given tile and indices"""
    #returns dictionaries with the data
    
    vgrid_data = get_UFS_vgrid_data(grid_dir) #only needed for ak, bk to calculate pressure
    state_data = get_UFS_state_data(vgrid_data, dir, tile, i, j, old_chgres, lam)
    surface_data = get_UFS_surface_data(dir, tile, i, j, old_chgres, lam)
    oro_data = get_UFS_oro_data(dir, tile, i, j, lam)
    
    return (state_data, surface_data, oro_data)
    
def get_UFS_state_data(vgrid, dir, tile, i, j, old_chgres, lam):
    """Get the state data for the given tile and indices"""
    
    if lam:
        nc_file_data = Dataset('{0}/{1}'.format(dir,'gfs_data.nc'))
    else:
        nc_file_data = Dataset('{0}/{1}'.format(dir,'gfs_data.tile{0}.nc'.format(tile)))
    
    # get nlevs from the gfs_ctrl.nc data
    nlevs_model=vgrid["nlevs"]
    
    # upper air fields from initial conditions (all data are top-first)
    zh_rev=nc_file_data['zh'][:,j,i]
    sphum_rev=nc_file_data['sphum'][:,j,i]
    # o3 and qv are taken from ics. 
    o3_rev=nc_file_data['o3mr'][:,j,i]
    liqwat_rev=nc_file_data['liq_wat'][:,j,i]
    ps_data = nc_file_data['ps'][j,i]
    
    #The 3D fields above are apparently on grid vertical interfaces. In the file external_ic.F90/get_nggps_ic subroutine in FV3, these fields
    #are further processed to get to the vertical grid centers/means.
    
    # following remap_scalar_nggps in external_ic.F90
    levp_data = len(sphum_rev)
        
    ak_rev = vgrid["ak"][::-1]
    bk_rev = vgrid["bk"][::-1]
    ak_rev[0] = np.max([1.0E-9, ak_rev[0]])
    
    ptop_data = ak_rev[1]
            
    pressure_from_data_rev = ak_rev + bk_rev*ps_data
    log_pressure_from_data_rev = np.log(pressure_from_data_rev)
    
    gz_rev = np.zeros(2*levp_data +1)
    pn_rev = np.zeros(2*levp_data +1)
        
    for k in range(0,levp_data+1):
        gz_rev[k] = zh_rev[k]*grav
        pn_rev[k] = log_pressure_from_data_rev[k]
    k2 = int(np.max([10, levp_data/2]))
    for k in range(levp_data+1,levp_data+k2):
        #do k=km+2, km+k2
        l = 2*(levp_data) - k
        gz_rev[k] = 2.*gz_rev[levp_data] - gz_rev[l]
        pn_rev[k] = 2.*pn_rev[levp_data] - pn_rev[l]
    
    phis = zh_rev[-1]*grav
    
    for k in range(levp_data+k2-2,0,-1):
        #do k=km+k2-1, 2, -1
        if (phis <= gz_rev[k] and phis >= gz_rev[k+1]):
            log_ps_calc = pn_rev[k] + (pn_rev[k+1]-pn_rev[k])*(gz_rev[k]-phis)/(gz_rev[k]-gz_rev[k+1])
            break
    
    ps_calc = np.exp(log_ps_calc)
    
    pressure_model_interfaces_rev = np.zeros(nlevs_model+1)
    log_pressure_model_interfaces_rev = np.zeros(nlevs_model+1)
    pressure_model_interfaces_rev[0] = ak_rev[1]
    log_pressure_model_interfaces_rev[0] = np.log(pressure_model_interfaces_rev[0])
    for k in range(1,nlevs_model+1):
        pressure_model_interfaces_rev[k] = ak_rev[k+1] + bk_rev[k+1]*ps_calc
        log_pressure_model_interfaces_rev[k] = np.log(pressure_model_interfaces_rev[k])
    
    pressure_thickness_model_rev = np.zeros(nlevs_model)
    for k in range(0,nlevs_model):
        pressure_thickness_model_rev[k] = pressure_model_interfaces_rev[k+1] - pressure_model_interfaces_rev[k]
    
    sphum_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], sphum_rev[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, 0, 8, ptop_data)
    sphum_model_rev_3d = fv3_remap.fillq(1, nlevs_model, 1, np.expand_dims(sphum_model_rev, axis=2), pressure_thickness_model_rev[np.newaxis, :])
    sphum_model_rev = sphum_model_rev_3d[:,:,0]
    
    o3_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], o3_rev[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, 0, 8, ptop_data)
    o3_model_rev_3d = fv3_remap.fillz(1, nlevs_model, 1, np.expand_dims(o3_model_rev, axis=2), pressure_thickness_model_rev[np.newaxis, :])
    o3_model_rev = o3_model_rev_3d[:,:,0]
    
    liqwat_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], liqwat_rev[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, 0, 8, ptop_data)
    liqwat_model_rev_3d = fv3_remap.fillz(1, nlevs_model, 1, np.expand_dims(liqwat_model_rev, axis=2), pressure_thickness_model_rev[np.newaxis, :])
    liqwat_model_rev = liqwat_model_rev_3d[:,:,0]
    
    if old_chgres:
        gz_fv = np.zeros(nlevs_model+1)
        gz_fv[-1] = phis
        m = 0
        for k in range(0,nlevs_model):
            for l in range(m, levp_data+k2-1):
                if ( (log_pressure_model_interfaces_rev[k] <= pn_rev[l+1]) and (log_pressure_model_interfaces_rev[k] >= pn_rev[l]) ):
                    gz_fv[k] = gz_rev[l] + (gz_rev[l+1]-gz_rev[l])*(log_pressure_model_interfaces_rev[k]-pn_rev[l])/(pn_rev[l+1]-pn_rev[l])
                    break
            m = l
        
        temp_model_rev = np.zeros((1,nlevs_model))
        for k in range(0, nlevs_model):
            temp_model_rev[0,k] = (gz_fv[k]-gz_fv[k+1])/(rdgas*(log_pressure_model_interfaces_rev[k+1]-log_pressure_model_interfaces_rev[k])*(1.+zvir*sphum_model_rev[0,k]) )
    else:
        temp_rev = nc_file_data['t'][:,j,i]
        
        temp_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], temp_rev[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, 2, 4, ptop_data)
        
    
    icewat_model_rev = np.zeros(nlevs_model)
    all_liquid_threshold = 273.16
    all_ice_threshold = 233.16
    intermediate_threshold = 258.16
    cloud_ice_mixing_ratio_threshold = 1.0E-5
    for k in range(0, nlevs_model):
        cloud_water = liqwat_model_rev[0,k]
        if (temp_model_rev[0,k] > all_liquid_threshold):
            liqwat_model_rev[0,k] = cloud_water
            icewat_model_rev[k] = 0.0
        elif (temp_model_rev[0,k] < all_ice_threshold):
            liqwat_model_rev[0,k] = 0.0
            icewat_model_rev[k] = cloud_water
        else:
            if k == 0:
                liqwat_model_rev[0,k] = cloud_water*(temp_model_rev[0,k]-all_ice_threshold)/(all_liquid_threshold - all_ice_threshold)
                icewat_model_rev[k] = cloud_water - liqwat_model_rev[0,k]
            else:
                if (temp_model_rev[0,k] < intermediate_threshold and icewat_model_rev[k-1] > cloud_ice_mixing_ratio_threshold):
                    liqwat_model_rev[0,k] = 0.0
                    icewat_model_rev[k] = cloud_water
                else:
                    liqwat_model_rev[0,k] = cloud_water*(temp_model_rev[0,k]-all_ice_threshold)/(all_liquid_threshold - all_ice_threshold)
                    icewat_model_rev[k] = cloud_water - liqwat_model_rev[0,k]
        (liqwat_model_rev[0,k], dummy_rain, icewat_model_rev[k], dummy_snow) = fv3_remap.mp_auto_conversion(liqwat_model_rev[0,k], icewat_model_rev[k])
    
    [u_s, u_n, v_w, v_e] = get_zonal_and_meridional_winds_on_cd_grid(tile, dir, i, j, nc_file_data, lam)
    
    #put C/D grid zonal/meridional winds on model pressure levels
    u_s_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], u_s[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, -1, 8, ptop_data)
    u_n_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], u_n[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, -1, 8, ptop_data)
    v_w_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], v_w[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, -1, 8, ptop_data)
    v_e_model_rev = fv3_remap.mappm(levp_data, pressure_from_data_rev[np.newaxis, :], v_e[np.newaxis, :], nlevs_model, pressure_model_interfaces_rev[np.newaxis, :], 1, 1, -1, 8, ptop_data)
    
    #put C/D grid zonal/meridional winds on A grid (simple averaging for now, but FV3 has more complex methods that should be implemented)
    
    u_model_rev = np.zeros(nlevs_model)
    v_model_rev = np.zeros(nlevs_model)
    u_model_rev = 0.5*(u_s_model_rev + u_n_model_rev)
    v_model_rev = 0.5*(v_w_model_rev + v_e_model_rev)
    
    nc_file_data.close()
    
    pressure_model_interfaces = pressure_model_interfaces_rev[::-1]
    pressure_model = np.zeros(nlevs_model)
    for k in range(0,nlevs_model):
        #from gmtb_scm_vgrid
        pressure_model[k] = ((1.0/(rocp+1.0))*(pressure_model_interfaces[k]**(rocp+1.0) - pressure_model_interfaces[k+1]**(rocp+1.0))/(pressure_model_interfaces[k] - pressure_model_interfaces[k+1]))**(1.0/rocp)
    
    #put data in a dictionary
    state = {
        "time": 1,
        "nlevs": nlevs_model,
        "z": zh_rev[::-1],
        "u": u_model_rev[0][::-1],
        "v": v_model_rev[0][::-1],
        "qv": sphum_model_rev[0][::-1],
        "o3": o3_model_rev[0][::-1],
        "ql": liqwat_model_rev[0][::-1],
        "qi": icewat_model_rev[::-1],
        "p_surf": ps_calc,
        "T_surf": temp_model_rev[0,0],
        "T": temp_model_rev[0,::-1],
        "pres": pressure_model,
        "pres_i": pressure_model_interfaces
    }
        
    return state

def get_zonal_and_meridional_winds_on_cd_grid(tile, dir, i, j, nc_file_data, lam):
    if lam:
        filename_pattern = '*grid.tile{0}.halo{1}.nc'.format(tile, n_lam_halo_points)
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, filename_pattern):
              filename = f_name
        if not filename:
            message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
            logging.critical(message)
            raise Exception(message)
    else:
        filename_pattern = '*grid.tile{0}.nc'.format(tile)
        
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, filename_pattern):
              filename = f_name
        if not filename:
            message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
            logging.critical(message)
            raise Exception(message)
    
    nc_file_grid = Dataset('{0}/{1}'.format(dir,filename))
    
    if (lam):
        #strip ghost/halo points and return supergrid
        lon_super_data = np.asarray(nc_file_grid['x'])
        lat_super_data = np.asarray(nc_file_grid['y'])
        #assuming n_lam_halo_points
        lon_super = lon_super_data[2*n_lam_halo_points:lon_super_data.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lon_super_data.shape[1]-2*n_lam_halo_points]
        lat_super = lat_super_data[2*n_lam_halo_points:lat_super_data.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:lat_super_data.shape[1]-2*n_lam_halo_points]
    else:
        lon_super = np.asarray(nc_file_grid['x'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
        lat_super = np.asarray(nc_file_grid['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    
    num_agrid_x = int(0.5*(lon_super.shape[1]-1))
    num_agrid_y = int(0.5*(lon_super.shape[0]-1))
    
    #find orientation
    #A-grid point
    agrid_super_i_index = 2*i + 1
    agrid_super_j_index = 2*j + 1
    point_on_agrid = np.asarray((lon_super[agrid_super_j_index,agrid_super_i_index],lat_super[agrid_super_j_index,agrid_super_i_index]))
    
    test_dgrid_points = [(lon_super[agrid_super_j_index,agrid_super_i_index+1],lat_super[agrid_super_j_index,agrid_super_i_index+1]),\
                         (lon_super[agrid_super_j_index,agrid_super_i_index-1],lat_super[agrid_super_j_index,agrid_super_i_index-1]),\
                         (lon_super[agrid_super_j_index+1,agrid_super_i_index],lat_super[agrid_super_j_index+1,agrid_super_i_index]),\
                         (lon_super[agrid_super_j_index-1,agrid_super_i_index],lat_super[agrid_super_j_index-1,agrid_super_i_index])]
    
    test_lon_diff = [p[0] - point_on_agrid[0] for p in test_dgrid_points]
    test_lat_diff = [p[1] - point_on_agrid[1] for p in test_dgrid_points]
    
    east_test_point = np.argmax(test_lon_diff)
    north_test_point = np.argmax(test_lat_diff)

    gridNotKnown = False
    if east_test_point == 0:
        #longitude increases most along the positive i axis
        if north_test_point == 2:
            #latitude increases most along the positive j axis
            #     ---> j+ north
            #     |
            #     V
            #     i+ east
            
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i,]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i+1)],lat_super[2*(j+1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j+1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j+1,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i+1)],lat_super[2*(j+1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j,i+1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i+1]*fv3_remap.inner_prod(e1, ey)
        elif north_test_point == 3:
            #latitude increases most along the negative j axis
            # <--- j- north
            #    |
            #    V
            #    i+ east
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i+1)],lat_super[2*(j-1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j-1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j-1,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i+1)],lat_super[2*(j-1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j,i+1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i+1]*fv3_remap.inner_prod(e1, ey)
        else:
            gridNotKnown = True
            print('unknown grid orientation A')
    elif east_test_point == 1:
        #longitude increases most along the negative i axis
        if north_test_point == 2:
            #latitude increases most along the positive j axis
            #     i- east
            #     ^
            #     |
            #     ---> j+ north
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i-1)],lat_super[2*(j+1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j+1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j+1,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i-1)],lat_super[2*(j+1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j,i-1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i-1]*fv3_remap.inner_prod(e1, ey)
        elif north_test_point == 3:
            #latitude increases most along the negative j axis
            #     i- east
            #     ^
            #     |
            # <--- j- north
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i-1)],lat_super[2*(j-1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j-1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j-1,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i-1)],lat_super[2*(j-1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j,i-1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i-1]*fv3_remap.inner_prod(e1, ey)
        else:
            gridNotKnown = True
            print('unknown grid orientation B')
    elif east_test_point == 2:
        #longitude increases most along the positive j axis
        if north_test_point == 0:
            #latitude increases most along the positive i axis
            #     ---> j+ east
            #     |
            #     V
            #     i+ north
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i+1)],lat_super[2*(j+1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j,i+1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i+1]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i+1)],lat_super[2*(j+1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ey)
        elif north_test_point == 1:
            #latitude increases most along the negative i axis
            #     i- north
            #     ^
            #     |
            #     ---> j+ east
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i-1)],lat_super[2*(j+1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j,i-1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i-1]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*(j+1),2*i],lat_super[2*(j+1),2*i]))
            p2 = np.asarray((lon_super[2*(j+1),2*(i-1)],lat_super[2*(j+1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ey)
        else:
            gridNotKnown = True
            print('unknown grid orientation C')
    elif east_test_point == 3:
        #longitude increases most along the negative j axis
        if north_test_point == 0:
            #latitude increases most along the positive i axis
            # <--- j- east
            #    |
            #    V
            #    i+ north
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i+1)],lat_super[2*(j-1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j,i+1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i+1]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i+1)],lat_super[2*j,2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i+1)],lat_super[2*(j-1),2*(i+1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j-1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j-1,i]*fv3_remap.inner_prod(e1, ey)
        elif north_test_point == 1:
            #latitude increases most along the negative i axis
            #     i- north
            #     ^
            #     |
            # <--- j- east
            #calculation of zonal wind on first (south) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_s = nc_file_data['u_s'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of zonal wind on second (north) D-grid point
            p1 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i-1)],lat_super[2*(j-1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            u_n = nc_file_data['u_s'][:,j,i-1]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_s'][:,j,i-1]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on first (west) D-grid point
            p1 = np.asarray((lon_super[2*j,2*i],lat_super[2*j,2*i]))
            p2 = np.asarray((lon_super[2*j,2*(i-1)],lat_super[2*j,2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_w = nc_file_data['u_w'][:,j,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j,i]*fv3_remap.inner_prod(e1, ey)
            
            #calculation of meridionial wind on second (east) D-grid point
            p1 = np.asarray((lon_super[2*(j-1),2*i],lat_super[2*(j-1),2*i]))
            p2 = np.asarray((lon_super[2*(j-1),2*(i-1)],lat_super[2*(j-1),2*(i-1)]))
            p3 = fv3_remap.mid_pt_sphere(p1*deg_to_rad, p2*deg_to_rad)
            e1 = fv3_remap.get_unit_vect2(p1*deg_to_rad, p2*deg_to_rad)
            (ex, ey) = fv3_remap.get_latlon_vector(p3)
            v_e = nc_file_data['u_w'][:,j-1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j-1,i]*fv3_remap.inner_prod(e1, ey)
        else:
            gridNotKnown = True
            print('unknown grid orientation D')    
    
    nc_file_grid.close()

    if gridNotKnown: exit()
    
    return [u_s, u_n, v_w, v_e]

def get_UFS_surface_data(dir, tile, i, j, old_chgres, lam):
    """Get the surface data for the given tile and indices"""
    
    if lam:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.nc'))
    else:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.tile{0}.nc'.format(tile)))

    vars_in_sfc_data = nc_file.variables.keys()

    #FV3/io/FV3GFS_io.F90/sfc_prop_restart_read was used as reference for variables that can be read in    
    
    #read in scalars (would be 2D variables in a 3D model)
    #{"name": "", "dim":}
    surface_vars = [{"name": "tsea",              "dim":0}, {"name": "tg3",              "dim":0}, \
                    {"name": "uustar",            "dim":0}, {"name": "alvsf",            "dim":0}, \
                    {"name": "alvwf",             "dim":0}, {"name": "alnsf",            "dim":0}, \
                    {"name": "alnwf",             "dim":0}, {"name": "facsf",            "dim":0}, \
                    {"name": "facwf",             "dim":0}, {"name": "stype",            "dim":0}, \
                    {"name": "slope",             "dim":0}, {"name": "vtype" ,           "dim":0}, \
                    {"name": "vfrac",             "dim":0}, {"name": "shdmin",           "dim":0}, \
                    {"name": "shdmax",            "dim":0}, {"name": "zorl",             "dim":0}, \
                    {"name": "slmsk",             "dim":0}, {"name": "canopy",           "dim":0}, \
                    {"name": "hice",              "dim":0}, {"name": "fice",             "dim":0}, \
                    {"name": "tisfc",             "dim":0}, {"name": "snwdph",           "dim":0}, \
                    {"name": "snoalb",            "dim":0}, {"name": "sheleg",           "dim":0}, \
                    {"name": "f10m",              "dim":0}, {"name": "t2m",              "dim":0}, \
                    {"name": "q2m",               "dim":0}, {"name": "ffmm",             "dim":0}, \
                    {"name": "ffhh",              "dim":0}, {"name": "tprcp",            "dim":0}, \
                    {"name": "srflag",            "dim":0}, {"name": "sncovr",           "dim":0}, \
                    {"name": "tsfcl",             "dim":0}, {"name": "zorlwav",          "dim":0}, \
                    #{"name": "zorlo",            "dim":0}, {"name": "zorll" ,            "dim":0}, \
                    #{"name": "zorli",            "dim":0}, {"name": "zorlw" ,            "dim":0}, \
                    {"name": "tref",              "dim":0}, {"name": "z_c",              "dim":0}, \
                    {"name": "c_0",               "dim":0}, {"name": "c_d",              "dim":0}, \
                    {"name": "w_0",               "dim":0}, {"name": "w_d",              "dim":0}, \
                    {"name": "xt",                "dim":0}, {"name": "xs",               "dim":0}, \
                    {"name": "xu",                "dim":0}, {"name": "xv",               "dim":0}, \
                    {"name": "zm",                "dim":0}, {"name": "xtts",             "dim":0}, \
                    {"name": "xzts",              "dim":0}, {"name": "d_conv",           "dim":0}, \
                    {"name": "xz",                "dim":0}, \
                    {"name": "ifd",               "dim":0}, {"name": "dt_cool",          "dim":0}, \
                    {"name": "qrain",             "dim":0}, {"name": "snowxy",           "dim":0}, \
                    {"name": "tvxy",              "dim":0}, {"name": "tgxy",             "dim":0}, \
                    {"name": "canicexy",          "dim":0}, {"name": "canliqxy",         "dim":0}, \
                    {"name": "eahxy",             "dim":0}, {"name": "tahxy",            "dim":0}, \
                    {"name": "cmxy",              "dim":0}, {"name": "chxy",             "dim":0}, \
                    {"name": "fwetxy",            "dim":0}, {"name": "sneqvoxy",         "dim":0}, \
                    {"name": "alboldxy",          "dim":0}, {"name": "qsnowxy",          "dim":0}, \
                    {"name": "wslakexy",          "dim":0}, {"name": "zwtxy",            "dim":0}, \
                    {"name": "waxy",              "dim":0}, {"name": "wtxy",             "dim":0}, \
                    {"name": "lfmassxy",          "dim":0}, {"name": "rtmassxy",         "dim":0}, \
                    {"name": "stmassxy",          "dim":0}, {"name": "woodxy",           "dim":0}, \
                    {"name": "stblcpxy",          "dim":0}, {"name": "fastcpxy",         "dim":0}, \
                    {"name": "xsaixy",            "dim":0}, {"name": "xlaixy",           "dim":0}, \
                    {"name": "taussxy",           "dim":0}, {"name": "smcwtdxy",         "dim":0}, \
                    {"name": "deeprechxy",        "dim":0}, {"name": "rechxy",           "dim":0}, \
                    {"name": "albdvis",           "dim":0}, {"name": "albdnir",          "dim":0}, \
                    {"name": "albivis",           "dim":0}, {"name": "albinir",          "dim":0}, \
                    {"name": "emiss",             "dim":0}, {"name": "wetness",          "dim":0}, \
                    {"name": "clw_surf_land",     "dim":0}, {"name": "clw_surf_ice",     "dim":0}, \
                    {"name": "qwv_surf_land",     "dim":0}, {"name": "qwv_surf_ice",     "dim":0}, \
                    {"name": "tsnow_land",        "dim":0}, {"name": "tsnow_ice",        "dim":0}, \
                    {"name": "snowfall_acc_land", "dim":0}, {"name": "snowfall_acc_ice", "dim":0}, \
                    {"name": "sncovr_ice",        "dim":0}, {"name": "lai",              "dim":0}, \
                    {"name": "swe_snowfall_acc",  "dim":0}, \
                    {"name": "stc",               "dim": missing_variable_soil_layers}, \
                    {"name": "smc",               "dim": missing_variable_soil_layers}, \
                    {"name": "slc",               "dim": missing_variable_soil_layers}, \
                    {"name": "snicexy",           "dim": missing_variable_snow_layers}, \
                    {"name": "snliqxy",           "dim": missing_variable_snow_layers}, \
                    {"name": "tsnoxy",            "dim": missing_variable_snow_layers}, \
                    {"name": "smoiseq",           "dim": missing_variable_soil_layers}, \
                    {"name": "zsnsoxy",           "dim": missing_variable_soil_layers + missing_variable_snow_layers}, \
                    {"name": "tslb",              "dim": missing_variable_soil_layers}, \
                    {"name": "smois",             "dim": missing_variable_soil_layers}, \
                    {"name": "sh2o",              "dim": missing_variable_soil_layers}, \
                    {"name": "smfr",              "dim": missing_variable_soil_layers}, \
                    {"name": "flfr",              "dim": missing_variable_soil_layers}, \
                    {"name": "tiice",             "dim": missing_variable_ice_layers}]

    # Using surface_vars dictionary, read-in netcdf data and store into "surface" dictionary.
    surface = {}
    for var in surface_vars:
        if var["name"] in vars_in_sfc_data:
            surface[var["name"]] = read_NetCDF_surface_var(nc_file, var["name"], i, j, old_chgres, var["dim"])
        else:
            if var["dim"] == 0:
                surface[var["name"]] = 0.0
            else:
                surface[var["name"]] = np.full(var["dim"],0.0)
    nc_file.close()

    return surface

def get_UFS_oro_data(dir, tile, i, j, lam):
    """Get the orographic data for the given tile and indices"""
    
    if lam:
        nc_file = Dataset('{0}/{1}'.format(dir,'oro_data.nc'))
    else:
        filename_pattern = 'oro_data.tile{0}.nc'.format(tile)
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, filename_pattern):
              filename = f_name
              
        nc_file = Dataset('{0}/{1}'.format(dir,filename))

    # orographic properties
    stddev_in     = read_NetCDF_var(nc_file, "stddev",     i, j)
    convexity_in  = read_NetCDF_var(nc_file, "convexity",  i, j)
    oa1_in        = read_NetCDF_var(nc_file, "oa1",        i, j)
    oa2_in        = read_NetCDF_var(nc_file, "oa2",        i, j)
    oa3_in        = read_NetCDF_var(nc_file, "oa3",        i, j)
    oa4_in        = read_NetCDF_var(nc_file, "oa4",        i, j)
    ol1_in        = read_NetCDF_var(nc_file, "ol1",        i, j)
    ol2_in        = read_NetCDF_var(nc_file, "ol2",        i, j)
    ol3_in        = read_NetCDF_var(nc_file, "ol3",        i, j)
    ol4_in        = read_NetCDF_var(nc_file, "ol4",        i, j)
    theta_in      = read_NetCDF_var(nc_file, "theta",      i, j)
    gamma_in      = read_NetCDF_var(nc_file, "gamma",      i, j)
    sigma_in      = read_NetCDF_var(nc_file, "sigma",      i, j)
    elvmax_in     = read_NetCDF_var(nc_file, "elvmax",     i, j)
    orog_filt_in  = read_NetCDF_var(nc_file, "orog_filt",  i, j)
    orog_raw_in   = read_NetCDF_var(nc_file, "orog_raw",   i, j)
    land_frac_in  = read_NetCDF_var(nc_file, "land_frac",  i, j)
    lake_frac_in  = read_NetCDF_var(nc_file, "lake_frac",  i, j)
    lake_depth_in = read_NetCDF_var(nc_file, "lake_depth", i, j)
    nc_file.close()
    
    # Store data in dictionary.
    oro = {"stddev":     stddev_in,
           "convexity":  convexity_in,
           "oa1":        oa1_in,
           "oa2":        oa2_in,
           "oa3":        oa3_in,
           "oa4":        oa4_in,
           "ol1":        ol1_in,
           "ol2":        ol2_in,
           "ol3":        ol3_in,
           "ol4":        ol4_in,
           "theta_oro":  theta_in,
           "gamma":      gamma_in,
           "sigma":      sigma_in,
           "elvmax":     elvmax_in,
           "orog_filt":  orog_filt_in,
           "orog_raw":   orog_raw_in,
           "land_frac":  land_frac_in,
           "lake_frac":  lake_frac_in,
           "lake_depth": lake_depth_in}

    return oro

def get_UFS_vgrid_data(dir):
    """Get the vertical grid data for resolution of the data within the IC directory"""
    
    nc_file = Dataset('{0}/{1}'.format(dir,'gfs_ctrl.nc'))
    
    # vertical coordinate definition
    # GJF: it looks like there is an extra level on top that represents 0 Pa, otherwise these values are for vertical grid interfaces
    ak=nc_file['vcoord'][0,::-1]
    bk=nc_file['vcoord'][1,::-1]
    
    #GJF: in external_ic.F90, when external_eta is true (which it apparently is for FV3GFS runs), the top value is ignored
    #ak = ak[0:len(ak)-1]
    #bk = bk[0:len(bk)-1]
    
    nc_file.close()
    
    vgrid = {
        "ak": ak,
        "bk": bk,
        "nlevs": len(ak)-2  #full grid levels are interfaces - 1 and there is an extra level on top (subtract 2)
    }
    
    return vgrid    

def get_UFS_grid_area(dir, tile, i, j, lam):
    """Get the horizontal grid cell area for the given tile and indices"""
    #this information is in the supergrid files
    
    if lam:
        filename_pattern = '*grid.tile{0}.halo{1}.nc'.format(tile, n_lam_halo_points)
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, filename_pattern):
              filename = f_name
        if not filename:
            message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
            logging.critical(message)
            raise Exception(message)
    else:
        filename_pattern = '*grid.tile{0}.nc'.format(tile)
        
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, filename_pattern):
              filename = f_name
        if not filename:
            message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
            logging.critical(message)
            raise Exception(message)
        
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    
    # extract out area of grid cell
    
    #calculate supergrid indices from regular grid indices
    jpt2 = j*2+1
    ipt2 = i*2+1
    
    #from Phil Pegion: the area is calculated by adding up the 4 components of the contained supergrid cells
    if lam:
        area_data = nc_file['area'][:,:]
        area_data_no_halo = area_data[2*n_lam_halo_points:area_data.shape[0]-2*n_lam_halo_points,2*n_lam_halo_points:area_data.shape[1]-2*n_lam_halo_points]
        area_in = area_data_no_halo[jpt2-1:jpt2+1,ipt2-1:ipt2+1]
    else:
        area_in=nc_file['area'][jpt2-1:jpt2+1,ipt2-1:ipt2+1]
    
    return area_in.sum()

def get_UFS_forcing_data2(nlevs, ic_state, forcing_dir, grid_dir, tile, i, j, lam, lami, lamf, save_comp_data):
    """Get the horizontal and vertical advective tendencies for the given tile and indices"""

    regrid_output = 'point'
    #regrid_output = 'all'

    if lamf:
        atm_filename_pattern = 'atmf*.tile{0}.nc'.format(tile)
        sfc_filename_pattern = 'sfcf*.tile{0}.nc'.format(tile)
    else:
        atm_filename_pattern = 'atmf*.nc'
        sfc_filename_pattern = 'sfcf*.nc'

    if lami:
        print("Initial Conditions are on limited area model grid")
    else:
        print("Initial Conditions are NOT on limited area model grid")
    if lamf:
        print("Forcing data are on limited area model grid")
    else:
        print("Forcing data are NOT on limited area model grid")
    if lami == lamf:
        print("Both the initial conditions and forcings are on the same grids")

    # Create list of input forcing files
    atm_filenames = []
    sfc_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, atm_filename_pattern):
          atm_filenames.append(f_name)
       if fnmatch.fnmatch(f_name, sfc_filename_pattern):
          sfc_filenames.append(f_name)
    if not atm_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(atm_filename_pattern,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    if not sfc_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(sfc_filename_pattern,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    atm_filenames = sorted(atm_filenames)
    sfc_filenames = sorted(sfc_filenames)
    if (len(atm_filenames) != len(sfc_filenames)):
        message = 'The number of atm files and sfc files in {0} matching the patterns does not match.'.format(forcing_dir)
        logging.critical(message)
        raise Exception(message)

    #atm_filenames = atm_filenames[0:2]
    #sfc_filenames = sfc_filenames[0:2]

    kord_tm = -9
    kord_mt = 9
    kord_tr = 9
    t_min = 184.0
    q_min = 0.0
    secinhr = 3600.0

    ##################################################################################################################
    # Read in UFS-state (u,v,T,qv), stored in atmf*.nc history files.
    ##################################################################################################################
    #
    natmf = len(atm_filenames)

    #
    ps             = []
    p_interfaces   = []
    p_layers       = []
    t_layers       = []
    qv_layers      = []
    u_layers       = []
    v_layers       = []
    time_atm_hours = []
    (ic_grid_lon, ic_grid_lat) = get_initial_lon_lat_grid(grid_dir, tile, lam)
    for count, filename in enumerate(atm_filenames, start=1):
        # Open file
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)

        #check if output grid is different than initial (native) grid
        if lam:
            data_grid_lon = nc_file['grid_xt'][:,:]
            data_grid_lat = nc_file['grid_yt'][:,:]
        else:
            data_grid_lon = nc_file['lon'][:,:]
            data_grid_lat = nc_file['lat'][:,:]

        equal_grids = False
        if (ic_grid_lon.shape == data_grid_lon.shape and ic_grid_lat.shape == ic_grid_lat.shape):
            if (np.equal(ic_grid_lon,data_grid_lon).all() and np.equal(ic_grid_lat,data_grid_lat).all()):
                equal_grids = True
        
        if (not equal_grids):
            grid_in = {'lon': data_grid_lon, 'lat': data_grid_lat}
            if regrid_output == 'all':
                grid_out = {'lon': ic_grid_lon, 'lat': ic_grid_lat}
                i_get = i
                j_get = j
            elif regrid_output == 'point':
                grid_out = {'lon': np.reshape(ic_grid_lon[j,i],(-1,1)), 'lat': np.reshape(ic_grid_lat[j,i],(-1,1))}
                i_get = 0
                j_get = 0
            else:
                print('Unrecognized regrid_output variable. Exiting...')
                exit()
            
            print('Regridding {} onto native grid: regridding progress = {}%'.format(filename, 100.0*count/(len(atm_filenames) + len(sfc_filenames))))
            regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
            ps_data   = regridder(nc_file['pressfc'][0,:,:])
            t_data    = regridder(nc_file['tmp'][0,::-1,:,:])
            qv_data   = regridder(nc_file['spfh'][0,::-1,:,:])
            u_data    = regridder(nc_file['ugrd'][0,::-1,:,:])
            v_data    = regridder(nc_file['vgrd'][0,::-1,:,:])
        else:
            ps_data = nc_file['pressfc'][0,:,:]
            t_data  = nc_file['tmp'][0,::-1,:,:]
            qv_data = nc_file['spfh'][0,::-1,:,:]
            u_data  = nc_file['ugrd'][0,::-1,:,:]
            v_data  = nc_file['vgrd'][0,::-1,:,:]
            i_get = i
            j_get = j

        # Scalars
        ak    = getattr(nc_file, "ak")[::-1]
        bk    = getattr(nc_file, "bk")[::-1]

        # Extract fields for current time
        ps.append(       ps_data[   j_get, i_get])
        t_layers.append( t_data[ :, j_get, i_get])
        qv_layers.append(qv_data[:, j_get, i_get])
        u_layers.append( u_data[ :, j_get, i_get])
        v_layers.append( v_data[ :, j_get, i_get])

        time_atm_hours.append(nc_file['time'][0])

        # Compute pressure at model-level (edges)
        p_interface = np.zeros(nlevs+1)
        for k in range(nlevs+1):
            p_interface[k]=ak[k]+ps[-1]*bk[k]        
        p_interfaces.append(p_interface)

        # Compute pressure at model-layers (centers)
        p_layer = np.zeros(nlevs)
        for k in range(nlevs):
            p_layer[k] = ((1.0/(rocp+1.0))*(p_interface[k]**(rocp+1.0) - p_interface[k+1]**(rocp+1.0))/(p_interface[k] - p_interface[k+1]))**(1.0/rocp)
        p_layers.append(p_layer)
        #
        nc_file.close()

    #
    ufs_state = {"time": np.asarray(time_atm_hours)*secinhr, \
                 "qv":   np.asarray(qv_layers),      \
                 "T":    np.asarray(t_layers),       \
                 "u":    np.asarray(u_layers),       \
                 "v":    np.asarray(v_layers),       \
                 "pres": np.asarray(p_layers),       \
                 "presi":np.asarray(p_interfaces),   \
                 "ps":   np.asarray(ps),             \
                 "tv":   np.asarray(t_layers)*(1.0 + zvir*np.asarray(qv_layers))}

    ##################################################################################################################
    # Regrid UFS-state (u,v,T,qv), for all times, using the same vertical-grid from ICs
    ##################################################################################################################
    #
    pres_init_rev          = np.zeros([1,nlevs+1])
    pres_init_rev[0,:]     = ic_state["pres_i"][::-1]
    log_pres_init_rev      = np.zeros([1,nlevs+1])
    log_pres_init_rev[0,:] = np.log(pres_init_rev)
    #
    dp2             = np.zeros([1,nlevs])
    tv_layers_remap = np.zeros([natmf,nlevs])
    qv_layers_remap = np.zeros([natmf,nlevs])
    u_layers_remap  = np.zeros([natmf,nlevs])
    v_layers_remap  = np.zeros([natmf,nlevs])
    p_layers_remap  = np.zeros([natmf,nlevs])
    for t in range(natmf):
        for k in range(0,nlevs): dp2[0,k] = pres_init_rev[0,k+1] - pres_init_rev[0,k]
        qv_rev_new = fv3_remap.map1_q2(   nlevs, ufs_state["presi"][t:t+1,::-1],  ufs_state["qv"][t:t+1,::-1],              \
                                          nlevs, pres_init_rev, dp2, 0, 0,  0, kord_tr,         q_min)
        u_rev_new  = fv3_remap.map1_ppm(  nlevs, ufs_state["presi"][t:t+1,::-1],  ufs_state["u"][ t:t+1,::-1], 0.0,         \
                                          nlevs, pres_init_rev,      0, 0, -1, kord_tm               )
        v_rev_new  = fv3_remap.map1_ppm(  nlevs, ufs_state["presi"][t:t+1,::-1],  ufs_state["v"][ t:t+1,::-1], 0.0,         \
                                          nlevs, pres_init_rev,      0, 0, -1, kord_tm               )
        tv_rev_new = fv3_remap.map_scalar(nlevs, np.log(ufs_state["presi"][t:t+1,::-1]), ufs_state["tv"][t:t+1,::-1], np.zeros(1), \
                                          nlevs, log_pres_init_rev,  0, 0,  1, np.abs(kord_tm), t_min)
        
        tv_layers_remap[t,:] = tv_rev_new[0,::-1]
        qv_layers_remap[t,:] = qv_rev_new[0,::-1]
        u_layers_remap[t,:]  = u_rev_new[0,::-1]
        v_layers_remap[t,:]  = v_rev_new[0,::-1]
        p_layers_remap[t,:]  = ic_state["pres"]
    t_layers_remap           = tv_layers_remap/(1.0 + zvir*qv_layers_remap)

    #
    ufs_state_remap = {"time": ufs_state["time"],\
                       "qv":   qv_layers_remap, \
                       "T":    t_layers_remap,  \
                       "u":    u_layers_remap,  \
                       "v":    v_layers_remap,  \
                       "tv":   tv_layers_remap, \
                       "pres": p_layers_remap,  \
                       "ps":   ufs_state["ps"]}

    ##################################################################################################################
    # Read in diagnostic tendencies, stored in sfcf*.nc files
    ##################################################################################################################
    #
    nsfcf = len(sfc_filenames)
    #
    time_dyn_hours = []

    #
    # Dynamic tendency dictionary
    #
    vars_comp_nophys = [{"name":"dtend_temp_nophys"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_qv_nophys"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_nophys"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_nophys"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_cld_amt_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_graupel_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_ice_wat_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_liq_wat_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_nophys"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_rainwat_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_sgs_tke_nophys"   , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_snowwat_nophys"   , "values":[], "dims":[nsfcf,nlevs]}]
    #
    # Physics tendency dictionaries (Saved to comparision file, alogng with "ufs_state")
    #
    vars_comp_pbl    = [{"name":"dtend_qv_pbl"           , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_pbl"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_pbl"            , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_pbl"            , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_cld_amt_pbl"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_graupel_pbl"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_liq_wat_pbl"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_pbl"           , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_rainwat_pbl"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_sgs_tke_pbl"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_snowwat_pbl"      , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_dpcnv  = [{"name":"dtend_qv_deepcnv"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_deepcnv"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_deepcnv"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_deepcnv"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_sgs_tke_deepcnv"  , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_shlcnv = [{"name":"dtend_qv_shalcnv"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_shalcnv"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_shalcnv"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_shalcnv"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_sgs_tke_shalcnv"  , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_rad    = [{"name":"dtend_temp_lw"          , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_sw"          , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_mp     = [{"name":"dtend_qv_mp"            , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_mp"          , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_cld_amt_mp"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_graupel_mp"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_ice_wat_mp"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_liq_wat_mp"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_rainwat_mp"       , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_snowwat_mp"       , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_gwd    = [{"name":"dtend_temp_cnvgwd"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_cnvgwd"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_cnvgwd"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_temp_orogwd"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_orogwd"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_orogwd"         , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_totcld = [{"name":"dtend_cld_amt_cnvtrans" , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_o3     = [{"name":"dtend_o3_o3column"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_o3mix"         , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_photochem"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_prodloss"      , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_temp"          , "values":[], "dims":[nsfcf,nlevs]}]
    vars_comp_phys   = [{"name":"dtend_temp_phys"        , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_u_phys"           , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_v_phys"           , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_qv_phys"          , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_snowwat_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_sgs_tke_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_rainwat_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_o3_phys"          , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_liq_wat_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_ice_wat_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_graupel_phys"     , "values":[], "dims":[nsfcf,nlevs]},\
                        {"name":"dtend_cld_amt_phys"     , "values":[], "dims":[nsfcf,nlevs]}]

    vars_comp = []
    vars_comp.extend(vars_comp_pbl)
    #vars_comp.extend(vars_comp_dpcnv)
    #vars_comp.extend(vars_comp_shlcnv)
    #vars_comp.extend(vars_comp_rad)
    #vars_comp.extend(vars_comp_mp)
    #vars_comp.extend(vars_comp_gwd)
    #vars_comp.extend(vars_comp_totcld)
    #vars_comp.extend(vars_comp_o3)
    #vars_comp.extend(vars_comp_phys)

    #
    # Loop through all sfc_filenames...
    #
    init = True
    time_dyn_hours = []
    for count, filename in enumerate(sfc_filenames, start=1):
        # Open file
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)

        # Initialize dictionaries
        if (init):
            # Get list of variables stored in file.
            var_av = nc_file.variables.keys()

            # Loop through dictionary, if requested field present, initialize field.
            for var in vars_comp_nophys:
                if var["name"] in var_av:
                    var["values"]    = np.zeros(var["dims"])
                    var["units"]     = nc_file[var["name"]].getncattr(name="units")
                    var["long_name"] = nc_file[var["name"]].getncattr(name="long_name")
                    var["active"]    = True
                else:
                    var["active"]    = False
            # If requested, initialize comparison data dictionary...
            if (save_comp_data):
                for var in vars_comp:
                    if var["name"] in var_av:
                        var["values"]    = np.zeros(var["dims"])
                        var["units"]     = nc_file[var["name"]].getncattr(name="units")
                        var["long_name"] = nc_file[var["name"]].getncattr(name="long_name")
                        var["active"]    = True
                    else:
                        var["active"]    = False
            init=False

        #
        # Read in output grid
        #
        if lam:
            data_grid_lon = nc_file['grid_xt'][:,:]
            data_grid_lat = nc_file['grid_yt'][:,:]
        else:
            data_grid_lon = nc_file['lon'][:,:]
            data_grid_lat = nc_file['lat'][:,:]

        # Check if output grid is different than initial (native) grid.
        # If grids are different, regrid output grid to IC grid
        equal_grids = False
        if (ic_grid_lon.shape == data_grid_lon.shape and ic_grid_lat.shape == ic_grid_lat.shape):
            if (np.equal(ic_grid_lon,data_grid_lon).all() and np.equal(ic_grid_lat,data_grid_lat).all()):
                equal_grids = True
        
        if (not equal_grids):
            # Save point indices
            if regrid_output == 'all':
                grid_out = {'lon': ic_grid_lon, 'lat': ic_grid_lat}
                i_get = i
                j_get = j
            elif regrid_output == 'point':
                grid_out = {'lon': np.reshape(ic_grid_lon[j,i],(-1,1)), 'lat': np.reshape(ic_grid_lat[j,i],(-1,1))}
                i_get = 0
                j_get = 0
            else:
                print('Unrecognized regrid_output variable. Exiting...')
                exit()
            
            #
            # Regrid to requested UFS grid-point..
            #
            print('Regridding {} onto native grid: regridding progress = {}%'.format(filename, 100.0*(count + len(atm_filenames))/(len(atm_filenames) + len(sfc_filenames))))
            grid_in = {'lon': data_grid_lon, 'lat': data_grid_lat}
            regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
            
            #
            # Regrid and store in tendency dictionary...
            #
            for var in vars_comp_nophys:
                if var["active"]:
                    var["values"][count-1,:] = regridder(nc_file[var["name"]][0,::-1,:,:])[:,0,0]

            #
            # If requested, regrid and store comparision dictionary...
            #
            if (save_comp_data):
                for var in vars_comp:
                    if var["active"]:
                        var["values"][count-1,:] = regridder(nc_file[var["name"]][0,::-1,:,:])[:,0,0]
        else:
            # Save point indices
            i_get = i
            j_get = j

            #
            # Store in tendency dictionary...
            #
            for var in vars_comp_nophys:
                if var["active"]: var["values"][count-1,:] = nc_file[var["name"]][0,::-1,:,:]
            #
            # If requested, store in comparision dictionary...
            #
            if (save_comp_data):
                for var in vars_comp:
                    if var["active"]: var["values"][count-1,:] = nc_file[var["name"]][0,::-1,:,:]

        # Save time
        time_dyn_hours.append(nc_file['time'][0])

        # Close file
        nc_file.close()

    ##################################################################################################################
    # Compute advective terms.
    #
    # Use IC state at time = 0, and atmf* and sfcf* for subsequent timesteps.
    # Compute tendencies at time = 0 using the IC state @ t=0 and the state within atmf001.nc.
    # Need to omit state within sfcf000.nc and atmf000.nc, they are on a non-standard interval, and
    # we need to compute the tendency over the whole interval.
    ##################################################################################################################
    #
    # Initialize
    #
    dqvdt_adv = np.zeros([natmf, nlevs])
    dtdt_adv  = np.zeros([natmf, nlevs])
    dudt_adv  = np.zeros([natmf, nlevs])
    dvdt_adv  = np.zeros([natmf, nlevs])

    #
    # @ time = 0
    #
    dqvdt_adv[0,:] = (ufs_state_remap["qv"][1,:] - ic_state["qv"][:]) / (secinhr*ufs_state_remap["time"][1])
    dtdt_adv[0,:]  = (ufs_state_remap["T"][1,:]  - ic_state["T"][:])  / (secinhr*ufs_state_remap["time"][1])
    dudt_adv[0,:]  = (ufs_state_remap["u"][1,:]  - ic_state["u"][:])  / (secinhr*ufs_state_remap["time"][1])
    dvdt_adv[0,:]  = (ufs_state_remap["v"][1,:]  - ic_state["v"][:])  / (secinhr*ufs_state_remap["time"][1])

    #
    # @ time > 0
    #
    for t in range(natmf-1):
        #dtdt_adv[t+1,:]  = vars_comp_nophys[0]["values"][t+1,:]
        #dqvdt_adv[t+1,:] = vars_comp_nophys[1]["values"][t+1,:]
        #dudt_adv[t+1,:]  = vars_comp_nophys[2]["values"][t+1,:]
        #dvdt_adv[t+1,:]  = vars_comp_nophys[3]["values"][t+1,:]
        dt = secinhr*(ufs_state_remap["time"][t+1] - ufs_state_remap["time"][t])
        dqvdt_adv[t+1,:] = (ufs_state_remap["qv"][t+1,:] - ufs_state_remap["qv"][t,:]) / dt
        dtdt_adv[t+1,:]  = (ufs_state_remap["T"][t+1,:]  - ufs_state_remap["T"][t,:])  / dt
        dudt_adv[t+1,:]  = (ufs_state_remap["u"][t+1,:]  - ufs_state_remap["u"][t,:])  / dt
        dvdt_adv[t+1,:]  = (ufs_state_remap["v"][t+1,:]  - ufs_state_remap["v"][t,:])  / dt

    ##################################################################################################################
    #
    # Create comparision state data 
    #
    ##################################################################################################################
    #
    # Initialize
    #
    comp_time = np.zeros(natmf)
    comp_p  = np.zeros([natmf,nlevs])
    comp_qv = np.zeros([natmf,nlevs])
    comp_T  = np.zeros([natmf,nlevs])
    comp_u  = np.zeros([natmf,nlevs])
    comp_v  = np.zeros([natmf,nlevs])

    #
    # @ time = 0
    #
    comp_p[0,:]  = ic_state["pres"][:]
    comp_qv[0,:] = ic_state["qv"][:]
    comp_T[0,:]  = ic_state["T"][:]
    comp_u[0,:]  = ic_state["u"][:]
    comp_v[0,:]  = ic_state["v"][:]

    #
    # @ time > 0
    #
    for t in range(1,natmf):
        comp_time[t] = ufs_state["time"][t]
        comp_p[t,:]  = ufs_state["pres"][t,:]
        comp_qv[t,:] = ufs_state["qv"][t,:]
        comp_T[t,:]  = ufs_state["T"][t,:]
        comp_u[t,:]  = ufs_state["u"][t,:]
        comp_v[t,:]  = ufs_state["v"][t,:]
    
    #
    # Store in dictionary
    #
    state_comp = {"pres": comp_p,    "qv":comp_qv, "T": comp_T,\
                  "time": comp_time, "u": comp_u,  "v": comp_v}

    #if we had dynf,phyf files at every timestep (and the SCM timestep is made to match the UFS), then dqvdt_adv should be
    #applied uninterpolated for each time step. If dynf and phyf files represent time averages over the previous diagnostic period,
    #and if forcing terms are interpolatd in time in the SCM, then dqvdt_adv should represent the forcing values in the 
    #middle of time[t] and time[t+1] from dynf/phyf. That way, the time-averaged applied forcing from time[t] to time[t+1] in the SCM will 
    #be equal to what is derived from dynf/phyf. (preference should be to have option to remove time-interpolation of forcing such
    #that the constant forcing applied converged to time-step values as the diag interval approaches the time step)    
    
    time_method = 'constant_simple' #this is not implemented in the SCM code yet
    #time_method = 'constant_interp'
    #time_method = 'gradient' #this produced wonky results in the SCM; avoid until investigated more
    
    ##################################################################################################################
    # Interpolate tendencies in time
    ##################################################################################################################
    
    if (time_method == 'constant_simple'):
        print('Forcing should not be interpolated in time. Rather, forcing should held constant at their current values until the next forcing interval is reached.')
        ntimes = natmf

        #
        # Initialize
        #
        time          = np.zeros(natmf)
        p_s           = np.zeros((natmf),dtype=float)
        pressure_forc = np.zeros((nlevs,natmf),dtype=float)
        tot_advec_T   = np.zeros((nlevs,natmf),dtype=float)
        tot_advec_qv  = np.zeros((nlevs,natmf),dtype=float)
        tot_advec_u   = np.zeros((nlevs,natmf),dtype=float)
        tot_advec_v   = np.zeros((nlevs,natmf),dtype=float)

        # 
        # @ time = 0
        #
        p_s[0]             = ic_state["p_surf"]
        pressure_forc[:,0] = ic_state["pres"][:]
        tot_advec_T[:,0]   = dtdt_adv[0,:]
        tot_advec_qv[:,0]  = dqvdt_adv[0,:]
        tot_advec_u[:,0]   = dudt_adv[0,:]
        tot_advec_v[:,0]   = dvdt_adv[0,:]

        #
        # @ time > 0
        #
        for t in range(1,natmf):
            time[t]            = secinhr*time_dyn_hours[t]
            p_s[t]             = ufs_state["ps"][t]
            pressure_forc[:,t] = ufs_state["pres"][t,:]
            tot_advec_T[:,t]   = dtdt_adv[t,:]
            tot_advec_qv[:,t]  = dqvdt_adv[t,:]
            tot_advec_u[:,t]   = dudt_adv[t,:]
            tot_advec_v[:,t]   = dvdt_adv[t,:]

    elif (time_method == 'constant_interp'):
        print('Forcing can be interpolated in time, but the time values are chosen such that forcing will effectively be held consant during a diagnostic time interval.')
        ntimes       = 2*natmf
        time_setback = 1.0

        #
        # Initialize
        #
        time          = np.zeros(ntimes,        dtype=float)
        p_s           = np.zeros((ntimes),      dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T   = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv  = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u   = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v   = np.zeros((nlevs,ntimes),dtype=float)

        #
        # @ time = 0
        #
        p_s[0] = ic_state["p_surf"]
        p_s[1] = p_s[0]
        pressure_forc[:,0] = ic_state["pres"][:]
        pressure_forc[:,1] = pressure_forc[:,0]
        tot_advec_T[:,0]   = dtdt_adv[0,:]
        tot_advec_T[:,1]   = tot_advec_T[:,0]
        tot_advec_qv[:,0]  = dqvdt_adv[0,:]
        tot_advec_qv[:,1]  = tot_advec_qv[:,0]
        tot_advec_u[:,0]   = dudt_adv[0,:]
        tot_advec_u[:,1]   = tot_advec_u[:,0]
        tot_advec_v[:,0]   = dvdt_adv[0,:]
        tot_advec_v[:,1]   = tot_advec_v[:,0]

        #
        # @ time > 0
        #
        time[1] = 3600.0*time_dyn_hours[1] - time_setback
        for t in range(1,natmf):
            time[2*t]   = secinhr*time_dyn_hours[t]
            time[2*t+1] = secinhr*(time_dyn_hours[t]+1) - time_setback
            p_s[2*t]    = ps[t]
            p_s[2*t+1]  = p_s[2*t]
            pressure_forc[:,2*t]   = ufs_state["pres"][t,:]
            pressure_forc[:,2*t+1] = pressure_forc[:,2*t]
            tot_advec_T[:,2*t]     = dtdt_adv[t,:]
            tot_advec_T[:,2*t+1]   = tot_advec_T[:,2*t]
            tot_advec_qv[:,2*t]    = dqvdt_adv[t,:]
            tot_advec_qv[:,2*t+1]  = tot_advec_qv[:,2*t]
            tot_advec_u[:,2*t]     = dudt_adv[t,:]
            tot_advec_u[:,2*t+1]   = tot_advec_u[:,2*t]
            tot_advec_v[:,2*t]     = dvdt_adv[t,:]
            tot_advec_v[:,2*t+1]   = tot_advec_v[:,2*t]

    elif (time_method == 'gradient'): #this produced wonky results in the SCM; avoid until investigated more
        print('Forcing can be interpolated in time since the forcing terms are assumed to follow a constant time-gradient.')
        
        ntimes = 2*natmf + 1
        # Initialize
        time          = np.zeros(ntimes)
        p_s           = np.zeros((ntimes),dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T   = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv  = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u   = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v   = np.zeros((nlevs,ntimes),dtype=float)
        
        p_s[0]             = ic_state['p_surf']
        pressure_forc[:,0] = ic_state['pres']
        tot_advec_T[:,0]   = 0.0
        tot_advec_qv[:,0]  = 0.0
        tot_advec_u[:,0]   = 0.0
        tot_advec_v[:,0]   = 0.0
        for t in range(natmf):
            time[2*t + 1] = time[2*t] + 0.5*(secinhr*time_dyn_hours[t] - time[2*t])
            time[2*t + 2] = secinhr*time_dyn_hours[t]
            p_s[2*t+1]    = ufs_state_remap["ps"][t]
            pressure_forc[:,2*t+1] = valid_pres_adv[t,:]
            tot_advec_T[:,2*t+1]   = dtdt_adv[t,:]
            tot_advec_qv[:,2*t+1]  = dqvdt_adv[t,:]
            tot_advec_u[:,2*t+1]   = dudt_adv[t,:]
            tot_advec_v[:,2*t+1]   = dvdt_adv[t,:]
            
            #calculate gradient in time and extrapolate for time (2t + 2)
            grad = (p_s[2*t + 1] - p_s[2*t])/(time[2*t + 1] - time[2*t])
            p_s[2*t + 2] = p_s[2*t + 1] + grad*(time[2*t + 2] - time[2*t + 1])
            
            for k in range(nlevs):
                grad = (pressure_forc[k,2*t + 1] - pressure_forc[k, 2*t])/(time[2*t + 1] - time[2*t])
                pressure_forc[k,2*t + 2] = pressure_forc[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                    
                grad = (tot_advec_T[k,2*t + 1] - tot_advec_T[k, 2*t])/(time[2*t + 1] - time[2*t])
                tot_advec_T[k,2*t + 2] = tot_advec_T[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                
                grad = (tot_advec_qv[k,2*t + 1] - tot_advec_qv[k, 2*t])/(time[2*t + 1] - time[2*t])
                tot_advec_qv[k,2*t + 2] = tot_advec_qv[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                
                grad = (tot_advec_u[k,2*t + 1] - tot_advec_u[k, 2*t])/(time[2*t + 1] - time[2*t])
                tot_advec_u[k,2*t + 2] = tot_advec_u[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                
                grad = (tot_advec_v[k,2*t + 1] - tot_advec_v[k, 2*t])/(time[2*t + 1] - time[2*t])
                tot_advec_v[k,2*t + 2] = tot_advec_v[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
    else:
        print('Unrecognized forcing time method. Exiting.')
        exit()

    forcing = {"time":          time,
               "w_ls":          np.zeros((nlevs,ntimes),dtype=float),
               "omega":         np.zeros((nlevs,ntimes),dtype=float),
               "u_g":           np.zeros((nlevs,ntimes),dtype=float),
               "v_g":           np.zeros((nlevs,ntimes),dtype=float),
               "u_nudge":       np.zeros((nlevs,ntimes),dtype=float),
               "v_nudge":       np.zeros((nlevs,ntimes),dtype=float),
               "T_nudge":       np.zeros((nlevs,ntimes),dtype=float),
               "thil_nudge":    np.zeros((nlevs,ntimes),dtype=float),
               "qt_nudge":      np.zeros((nlevs,ntimes),dtype=float),
               "rad_heating":   np.zeros((nlevs,ntimes),dtype=float),
               "h_advec_thil":  np.zeros((nlevs,ntimes),dtype=float),
               "v_advec_thil":  np.zeros((nlevs,ntimes),dtype=float),
               "h_advec_qt":    np.zeros((nlevs,ntimes),dtype=float),
               "v_advec_qt":    np.zeros((nlevs,ntimes),dtype=float),
               "h_advec_u":     np.zeros((nlevs,ntimes),dtype=float),
               "h_advec_v":     np.zeros((nlevs,ntimes),dtype=float),
               #"ps_forc": p_s, #uncomment when SCM switches to semi-Lagrangian vertical coordinate
               "ps_forc":       np.ones(ntimes)*ufs_state_remap["ps"][0],
               "tot_advec_T":   tot_advec_T,
               "tot_advec_qv":  tot_advec_qv,
               "tot_advec_u":   tot_advec_u,
               "tot_advec_v":   tot_advec_v,
               "pressure_forc": pressure_forc}

    return (forcing, vars_comp, state_comp)

def add_noahmp_coldstart(surface, date):
    """Add cold-start ICs for the NoahMP LSM from Noah LSM variables"""
    
    #use cold start section of FV3/io/FV3GFS_io.F90 to initialize NoahMP-specific variables (this is a python port of the Fortran code in that file)
    
    #MPTABLE.TBL uses a namelist format, so can use f90nml to read it in
    mptable_nml_all = f90nml.read(os.path.join(NOAHMP_TABLES_DIR, 'MPTABLE.TBL'))
    #MPTABLE.TBL contains data (with distinct namelists) for USGS and MODIS data; looks like MODIS is the operational
    mptable_nml_active = mptable_nml_all['noah_mp_modis_parameters'] #alternative is mptable_nml_all['noah_mp_usgs_parameters']
    
    #operational values; change if necessary (or read from somewhere?)
    n_snow_layers = 3
    n_soil_layers = 4
    
    #thickness of each soil level
    dzs = np.asarray([0.1,0.3,0.6,1.0])
    
    #bottom depth of each soil level
    zsoil = np.asarray([-0.1,-0.4,-1.0,-2.0])
    
    #initialize all NoahMP vars as missing
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
        surface["tvxy"] = surface["tsea"]
        surface["tgxy"] = surface["tsea"]
        surface["tahxy"] = surface["tsea"]
        
        if (surface["snwdph"] > 0.01 and surface["tsea"] > 273.15 ):
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
        
        vegtyp = int(surface['vtype'])
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
            #laim gives monthly values for each of the vegetation types
            laim = np.asarray(mptable_nml_active['LAIM']).reshape(12,20)
            
            #be sure to use month-1, vegtyp-1 since python is 0-indexed
            surface["xlaixy"] = np.amax([laim[date["month"]-1,vegtyp-1],0.05])
            surface["xsaixy"] = np.amax([surface["xlaixy"]*0.1,0.05])
            
            sla = np.asarray(mptable_nml_active['SLA'])
            masslai = 1000.0 / np.amax([sla[vegtyp-1],1.0])
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
        
        isnow = int(surface["snowxy"] + n_snow_layers)
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
        
        soiltyp  = int(surface["stype"])
        if (soiltyp != 0):
            #find the index of the soiltype from the "index" field
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
    """Read values from SOILPARM.TBL for NoahMP LSM ICs"""
    #returns a dictionary with data
    
    #two different datasets are included in the table
    choices = ["STAS","STAS-RUC"]
    
    #get all lines of the file
    with open(os.path.join(NOAHMP_TABLES_DIR, 'SOILPARM.TBL'), 'r') as f:
        lineList = f.readlines()
    f.close()
    
    #find the line where the desired data starts 
    line_index = 0
    for line in lineList:
        line_index += 1
        #hardcoded to look for choices[0]; swap choices[0] for choices[1] to use choices[1] below
        m = re.match(choices[0],line) and not re.match(choices[1],line)
        if m:
            start_index = line_index
            break
    
    #get the data for each variable from the lines    
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

def write_comparison_file(vars_comp, state_comp, case_name, date, surface, add_UFS_dyn_tend, add_UFS_NOAH_lsm):
    """Write UFS history file data to netCDF file for comparison"""
    
    #
    real_type = np.float64
    int_type  = np.int32
    #
    start_date        = datetime(date["year"],date["month"],date["day"],date["hour"],date["minute"],date["second"])
    start_date_string = start_date.strftime("%Y%m%d%H%M%S")
    #
    loc_string  = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string
    #
    fileOUT           = os.path.join(COMPARISON_DATA_DIR, case_name + '_comp_data.nc')
    nc_file           = Dataset(fileOUT, 'w', format='NETCDF4')
    nc_file.case      = case_string
    nc_file.title     = 'UFS history file data for ' + case_string
    nc_file.reference = ''
    nc_file.author    = 'Grant J. Firl and Dustin Swales'
    nc_file.version   = 'Created on ' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    nc_file.script    = os.path.basename(__file__)
    nc_file.startDate = start_date_string

    #
    # Set file dimensions
    #
    lon_dim  = nc_file.createDimension('lon', 1)
    lat_dim  = nc_file.createDimension('lat', 1)
    if (add_UFS_dyn_tend):
        lev_dim  = nc_file.createDimension('lev', state_comp["pres"].shape[1])
        time_dim = nc_file.createDimension('t0',  state_comp["time"].shape[0])
    else:
        lev_dim  = nc_file.createDimension('lev', state_comp["pres"].shape[0])
        time_dim = nc_file.createDimension('t0',  1)

    # Variable dictionary to create netcdf file (DEPHY format).
    var_dict = [{"dict": date,       "name": "year",   "type":int_type,  "dimd": ( ),                         "units": "year",          "description":"year at time of initial values", "alias": "init_year"}, \
                {"dict": date,       "name": "month",  "type":int_type,  "dimd": ( ),                         "units": "month",         "description": "month at time of initial values", "alias": "init_month"}, \
                {"dict": date,       "name": "day",    "type":int_type,  "dimd": ( ),                         "units": "day",           "description": "day at time of initial values", "alias": "init_day"}, \
                {"dict": date,       "name": "hour",   "type":int_type,  "dimd": ( ),                         "units": "hour",          "description": "hour at time of initial values", "alias": "init_hour"}, \
                {"dict": date,       "name": "minute", "type":int_type,  "dimd": ( ),                         "units": "minute",        "description": "minute at time of initial values", "alias": "init_minute"}, \
                {"dict": date,       "name": "second", "type":int_type,  "dimd": ( ),                         "units": "second",        "description": "second at time of initial values", "alias": "init_second"}, \
                {"dict": state_comp, "name": "time",   "type":real_type, "dimd": ('t0'),                      "units": 'seconds since ' + str(start_date), "description": "history file time"}, 
                {"dict": state_comp, "name": "pres",   "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": 'Pa', "description": "pressure"}, \
                {"dict": state_comp, "name": "v",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m s-1",         "description": "meridional wind"},\
                {"dict": state_comp, "name": "u",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m s-1",         "description": "zonal wind"},\
                {"dict": state_comp, "name": "qv",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "specific humidity"},\
                {"dict": state_comp, "name": "T",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "K",             "description": "Temperature"}]

    # Write dictionaries to output...

    #
    # var_dict (local)
    #
    for var in var_dict:
        var_name = var["name"]
        if "alias" in var: var_name = var["alias"]
        var_temp             = nc_file.createVariable(var_name, var["type"], var["dimd"])
        var_temp.units       = var["units"]
        var_temp.description = var["description"]
        var_temp[:]          = var["dict"][var["name"]]

    #
    # vars_comp (input; created in get_UFS_forcing_data2)
    #
    for var in vars_comp:
        if var["active"]:
            var_temp             = nc_file.createVariable(var["name"], real_type, ('t0', 'lev', 'lat', 'lon'))
            var_temp.units       = var["units"]
            var_temp.description = var["long_name"]
            var_temp[:]          = var["values"]

    #
    # Close file
    #
    nc_file.close()
    
    return (fileOUT)

def write_SCM_case_file(state, surface, oro, forcing, case, date, add_UFS_dyn_tend, add_UFS_NOAH_lsm):
    """Write all data to a netCDF file in the DEPHY-SCM format"""

    # Working types
    real_type = np.float64
    int_type  = np.int32
    
    # Dimensions
    nlevs = state["nlevs"]
    nsoil = len(surface["stc"])
    nsnow = len(surface["snicexy"])
    nice  = len(surface["tiice"])

    # Local switches
    forcing_on  = 1
    forcing_off = 0

    # Output file
    fileOUT = os.path.join(PROCESSED_CASE_DIR, case + '.nc')
    nc_file = Dataset(fileOUT, 'w', format='NETCDF4')
    if (add_UFS_dyn_tend):
        nc_file.description = "FV3GFS model profile input (UFS dynamic tendencies, SCM-UFS replay mode.)"
    elif (add_UFS_NOAH_lsm):
        nc_file.description = "FV3GFS model profile input (With NOAH Land surface moodel surface forcings)"
    else:
        nc_file.description = "FV3GFS model profile input (no forcings)"

    print("SCM case file created: ",fileOUT)

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
    if (add_UFS_dyn_tend): 
        runtime = timedelta(seconds=forcing['time'][-1])
    else:
        runtime = timedelta(seconds=0) 
    end_date   = start_date + runtime
    start_date_string = start_date.strftime("%Y%m%d%H%M%S")
    #
    loc_string  = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string
    #
    nc_file.case           = case_string
    nc_file.title          = 'Forcing and Initial Conditions for ' + case_string
    nc_file.reference      = ''
    nc_file.author         = 'Grant J. Firl and Dustin Swales'
    nc_file.version        = 'Created on ' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    nc_file.format_version = '1.0'
    #nc_file.modifications  = 'contains initial conditions for Noah LSM'
    nc_file.script         = os.path.basename(__file__)
    nc_file.comment        = ''
    nc_file.startDate      = start_date_string
    nc_file.endDate        = end_date.strftime("%Y%m%d%H%M%S")

    # Start with all forcing switches OFF
    nc_file.adv_theta        = forcing_off
    nc_file.adv_thetal       = forcing_off
    nc_file.rad_temp         = forcing_off
    nc_file.rad_theta        = forcing_off
    nc_file.rad_thetal       = forcing_off
    nc_file.adv_qt           = forcing_off
    nc_file.adv_rv           = forcing_off
    nc_file.adv_rt           = forcing_off
    nc_file.forc_w           = forcing_off
    nc_file.forc_omega       = forcing_off
    nc_file.forc_geo         = forcing_off
    nc_file.nudging_u        = forcing_off
    nc_file.nudging_v        = forcing_off
    nc_file.nudging_temp     = forcing_off
    nc_file.nudging_theta    = forcing_off
    nc_file.nudging_thetal   = forcing_off
    nc_file.nudging_qv       = forcing_off
    nc_file.nudging_qt       = forcing_off
    nc_file.nudging_rv       = forcing_off
    nc_file.nudging_rt       = forcing_off
    nc_file.z_nudging_temp   = forcing_off
    nc_file.z_nudging_theta  = forcing_off
    nc_file.z_nudging_thetal = forcing_off
    nc_file.z_nudging_qv     = forcing_off
    nc_file.z_nudging_qt     = forcing_off
    nc_file.z_nudging_rv     = forcing_off
    nc_file.z_nudging_rt     = forcing_off
    nc_file.z_nudging_u      = forcing_off
    nc_file.z_nudging_v      = forcing_off
    nc_file.p_nudging_temp   = forcing_off
    nc_file.p_nudging_theta  = forcing_off
    nc_file.p_nudging_thetal = forcing_off
    nc_file.p_nudging_qv     = forcing_off
    nc_file.p_nudging_qt     = forcing_off
    nc_file.p_nudging_rv     = forcing_off
    nc_file.p_nudging_rt     = forcing_off
    nc_file.p_nudging_u      = forcing_off
    nc_file.p_nudging_v      = forcing_off
    nc_file.adv_temp         = forcing_off
    nc_file.adv_qv           = forcing_off
    nc_file.adv_u            = forcing_off
    nc_file.adv_v            = forcing_off
    #
    nc_file.zorog              = oro['orog_filt']
    nc_file.z0                 = surface['zorl']
    nc_file.surfaceType        = surface_string

    if (add_UFS_dyn_tend):
        nc_file.adv_temp = forcing_on
        nc_file.adv_qv   = forcing_on
        nc_file.adv_u    = forcing_on
        nc_file.adv_v    = forcing_on
        #
        time_dim                   = nc_file.createDimension('time', len(forcing['time']))
        time_var                   = nc_file.createVariable('time', real_type, ('time',))
        time_var.units             = 'seconds since ' + str(start_date)
        time_var.long_name         = 'Forcing time'
        time_var[:]                = forcing['time']
    else:
        time_dim                   = nc_file.createDimension('time', 1)
        time_var                   = nc_file.createVariable('time', real_type, ('time',))
        time_var[:]                = 0.
    if (add_UFS_NOAH_lsm):
        nc_file.surfaceForcing     = 'lsm'
        nc_file.surfaceForcingWind = 'lsm'
    else:
        nc_file.surfaceForcing     = 'none'
        nc_file.surfaceForcingWind = 'none'

    # Set file dimension
    initial_time_dim = nc_file.createDimension('t0',    1)
    lev_dim          = nc_file.createDimension('lev',   nlevs)
    lon_dim          = nc_file.createDimension('lon',   1)
    lat_dim          = nc_file.createDimension('lat',   1)
    soil_dim         = nc_file.createDimension('nsoil', nsoil)
    snow_dim         = nc_file.createDimension('nsnow', nsnow)
    nslsnw_dim       = nc_file.createDimension('nsoil_plus_nsnow',nsnow + nsoil)
    ice_dim          = nc_file.createDimension('nice',  nice)
    
    #
    init_time_var              = nc_file.createVariable('t0', real_type, ('t0',))
    init_time_var.units        = 'seconds since ' + str(start_date)
    init_time_var.description  = 'Initial time'
    init_time_var.calendar     = 'gregorian'
    init_time_var[:]           = 0.0
    #
    lev_var                    = nc_file.createVariable('lev', real_type, ('lev',))
    lev_var.units              = 'Pa'
    lev_var.description        = 'pressure'
    lev_var[:]                 = state["pres"]
    #
    soil_depth_var             = nc_file.createVariable('soil_depth', real_type, ('nsoil',))
    soil_depth_var.units       = 'm'
    soil_depth_var.description = 'depth of bottom of soil layers'
    soil_depth_var[:]          = [0.1,0.4,1.0,2.0]
    #    
    lat_var                    = nc_file.createVariable('lat', real_type, ('lat',))
    lat_var.units              = 'degrees_north'
    lat_var.description        = "Latitude"
    lat_var[:]                 = surface["lat"]
    #
    lon_var                    = nc_file.createVariable('lon', real_type, ('lon',))
    lon_var.units              = 'degrees_east'
    lon_var.description        = "Longitude"
    lon_var[:]                 = surface["lon"]
    #
    zorlw_var                  = nc_file.createVariable('zorlw', real_type, ('t0', 'lat', 'lon'))
    zorlw_var.units            = "cm"
    zorlw_var.description      = "surface roughness length over ocean"
    zorlw_var[:]               = missing_value
    #
    zorll_var                  = nc_file.createVariable('zorll', real_type, ('t0', 'lat', 'lon'))
    zorll_var.units            = "cm"
    zorll_var.description      = "surface roughness length over land"
    zorll_var[:]               = missing_value
    #
    zorli_var                  = nc_file.createVariable('zorli', real_type, ('t0', 'lat', 'lon'))
    zorli_var.units            = "cm"
    zorli_var.description      = "surface roughness length over ice"
    zorli_var[:]               = missing_value
    
    #
    if (surface_string == "ice"):
        zorli_var[:] = surface["zorl"]
    elif (surface_string == "land"):
        zorll_var[:] = surface["zorl"]
    else:
        zorlw_var[:] = surface["zorl"]

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
    # {"name": "", "type", "dimd": (), "units": "", "description": "", "alias": ""}
    #
    ######################################################################################## 
    var_dict = [\
                #
                # Required SCM inputs
                #
                {"name": "p_surf", "type":real_type, "dimd": ('t0', 'lat', 'lon'),        "units": "Pa",            "description": "inital surface pressure", "alias":"ps"}, \
                {"name": "T",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "K",             "description": "Temperature", "alias":"temp"}, \
                {"name": "u",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m s-1",         "description": "Zonal wind"}, \
                {"name": "v",      "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m s-1",         "description": "Meridional wind"}, \
                {"name": "pres",   "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "Pa",            "description": "Pressure", "alias": "pressure"}, \
                {"name": "qt",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Total water content"}, \
                {"name": "qv",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Total pecific humidity"}, \
                {"name": "ql",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Liquid water specific humidity"}, \
                {"name": "qi",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Ice water specific humidity", "default_value": 0.0}, \
                {"name": "rv",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Water vapor mixing ratio"}, \
                {"name": "rt",     "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "kg kg-1",       "description": "Total water mixing ratio"}, \
                {"name": "tke",    "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m2 s-2",        "description": "Turbulen kinetic energy", "default_value": 0.0}, \
                {"name": "height", "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "m",             "description": "Height above ground"},\
                {"name": "theta",  "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "K",             "description": "Potential temperature"},\
                {"name": "thetal", "type":real_type, "dimd": ('t0', 'lev', 'lat', 'lon'), "units": "K",             "description": "Liquid potential temperature"}, \
                {"name": "o3",         "type":real_type, "dimd": ('t0', 'lev', 'lat','lon'), "units": "kg kg-1", "description": "initial profile of ozone mass mixing ratio", "alias": "ozone"}, \
                {"name": "area",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m^2",     "description": "grid cell area"}, \
                {"name": "stddev",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "standard deviation of subgrid orography"}, \
                {"name": "convexity",  "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "convexity of subgrid orography"}, \
                {"name": "oa1",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "assymetry of subgrid orography 1"}, \
                {"name": "oa2",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "assymetry of subgrid orography 2"}, \
                {"name": "oa3",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "assymetry of subgrid orography 3"}, \
                {"name": "oa4",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "assymetry of subgrid orography 4"}, \
                {"name": "ol1",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of grid box with subgrid orography higher than critical height 1"}, \
                {"name": "ol2",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of grid box with subgrid orography higher than critical height 2"}, \
                {"name": "ol3",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of grid box with subgrid orography higher than critical height 3"}, \
                {"name": "ol4",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of grid box with subgrid orography higher than critical height 4"}, \
                {"name": "sigma",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "slope of subgrid orography"}, \
                {"name": "theta_oro",  "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "deg",     "description": "angle with respect to east of maximum subgrid orographic variations"}, \
                {"name": "gamma",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "anisotropy of subgrid orography"}, \
                {"name": "elvmax",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "maximum of subgrid orography"}, \
                {"name": "orog_filt",  "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "orography", "alias": "oro"}, \
                {"name": "orog_raw",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "unfiltered orography", "alias": "oro_uf"}, \
                {"name": "land_frac",  "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of horizontal grid area occupied by land", "alias": "landfrac"}, \
                {"name": "lake_frac" , "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "fraction of horizontal grid area occupied by lake", "alias": "lakefrac", "default_value":0}, \
                {"name": "lake_depth", "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "lake depth", "alias": "lakedepth", "default_value":0}, \
                {"name": "tref",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "K",       "description": "sea surface reference temperature for NSST"}, \
                {"name": "z_c",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "sub-layer cooling thickness for NSST"}, \
                {"name": "c_0",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "coefficient 1 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "c_d",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "nonw",    "description": "coefficient 2 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "w_0",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "coefficient 3 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "w_d",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "coefficient 4 to calculate d(Tz)/d(Ts) for NSST"}, \
                {"name": "xt",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "K m",     "description": "heat content in diurnal thermocline layer for NSST"}, \
                {"name": "xs",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "ppt m",   "description": "salinity content in diurnal thermocline layer for NSST"}, \
                {"name": "xu",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m2 s-1",  "description": "u-current in diurnal thermocline layer for NSST"}, \
                {"name": "xv",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m2 s-1",  "description": "v-current in diurnal thermocline layer for NSST"}, \
                {"name": "xz",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "thickness of diurnal thermocline layer for NSST"}, \
                {"name": "zm"   ,      "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "thickness of ocean mixed layer for NSST"}, \
                {"name": "xtts",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST"}, \
                {"name": "xzts",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m K-1",   "description": "sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST"}, \
                {"name": "d_conv",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "m",       "description": "thickness of free convection layer for NSST"}, \
                {"name": "ifd",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "none",    "description": "index to start DTM run for NSST"}, \
                {"name": "dt_cool",    "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "K",       "description": "sub-layer cooling amount for NSST"}, \
                {"name": "qrain",         "type":real_type, "dimd": ('t0', 'lat', 'lon'),       "units": "W m-2",   "description": "sensible heat due to rainfall for NSST"},\
                {"name": "ps_forc",       "type":real_type, "dimd": ('time',        'lat', 'lon'),       "units": 'Pa',      "description": 'Surface pressure for forcing'},\
                {"name": "pressure_forc", "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'),       "units": 'Pa',      "description": 'Pressure for forcing'},\
                {"name": "height_forc",   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'),       "units": 'm',       "description": 'Height above the ground for forcing',"default_value": 1.}]
    var_lsm_ics = [\
                   #
                   # Fields needed if model_ics = True AND input_surfaceForcing == "lsm"
                   #
                   {"name": "stc",        "type":real_type, "dimd": ('t0','nsoil','lat','lon'),            "units": "K",       "description": "initial profile of soil liquid moisture"}, \
                   {"name": "smc",        "type":real_type, "dimd": ('t0','nsoil','lat','lon'),            "units": "kg",      "description": "initial profile of soil moisture"}, \
                   {"name": "slc",        "type":real_type, "dimd": ('t0','nsoil','lat','lon'),            "units": "kg",      "description": "initial profile of soil temperature"}, \
                   {"name": "snicexy",    "type":real_type, "dimd": ('t0','nsnow','lat','lon'),            "units": "mm",      "description": "initial profile of snow layer ice"}, \
                   {"name": "snliqxy",    "type":real_type, "dimd": ('t0','nsnow','lat','lon'),            "units": "mm",      "description": "initial profile of snow layer liquid"}, \
                   {"name": "tsnoxy",     "type":real_type, "dimd": ('t0','nsnow','lat','lon'),            "units": "K",       "description": "initial profile of snow layer temperature"}, \
                   {"name": "smoiseq",    "type":real_type, "dimd": ('t0','nsoil','lat','lon'),            "units": "m3 m-3",  "description": "initial profile of equilibrium soil water content"}, \
                   {"name": "zsnsoxy",    "type":real_type, "dimd": ('t0','nsoil_plus_nsnow','lat','lon'), "units": "m",       "description": "layer bottom depth from snow surface"},\
                   {"name": "tiice",      "type":real_type, "dimd": ('t0','nice', 'lat','lon'),            "units": "K",       "description": "sea ice internal temperature"}, \
                   {"name": "tslb",       "type":real_type, "dimd": ('t0','nsoil','lat','lon',),           "units": "K",       "description": "soil temperature for RUC LSM"}, \
                   {"name": "smois",      "type":real_type, "dimd": ('t0','nsoil','lat','lon',),           "units": "none",    "description": "volume fraction of soil moisture for RUC LSM"}, \
                   {"name": "sh2o",       "type":real_type, "dimd": ('t0','nsoil','lat','lon',),           "units": "none",    "description": "volume fraction of unfrozen soil moisture for RUC LSM"}, \
                   {"name": "smfr",       "type":real_type, "dimd": ('t0','nsoil','lat','lon',),           "units": "none",    "description": "volume fraction of frozen soil moisture for RUC LSM"}, \
                   {"name": "flfr",       "type":real_type, "dimd": ('t0','nsoil','lat','lon',),           "units": "none",    "description": "flag for frozen soil physics for RUC LSM"}, \
                   {"name": "tsea",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "sea/skin/ice surface temperature", "alias":"tsfco"}, \
                   {"name": "vegsrc",     "type":int_type,  "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "vegetation soure (1-2)", "default_value": 1}, \
                   {"name": "vtype",      "type":int_type,  "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "vegetation type (1-12)", "alias": "vegtyp"}, \
                   {"name": "stype",      "type":int_type,  "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "soil type (1-12)", "alias":"soiltyp"}, \
                   {"name": "slope",      "type":int_type,  "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "slope type (1-9)", "alias":"slopetyp"}, \
                   {"name": "vfrac",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "vegetation fraction", "alias":"vegfrac"}, \
                   {"name": "shdmin",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "minimum vegetation fraction"}, \
                   {"name": "shdmax",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "maximum vegetation fraction"}, \
                   {"name": "slmsk",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "land-sea-ice mask"}, \
                   {"name": "canopy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "kg m-2",  "description": "amount of water stored in canopy"}, \
                   {"name": "hice",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m",       "description": "sea ice thickness"}, \
                   {"name": "fice",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "ice fraction"}, \
                   {"name": "tisfc",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "ice surface temperature"}, \
                   {"name": "snwdph",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "water equivalent snow depth"}, \
                   {"name": "snoalb",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "maximum snow albedo"}, \
                   {"name": "sncovr",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "surface snow area fraction"}, \
                   {"name": "tg3",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "deep soil temperature"}, \
                   {"name": "uustar",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m s-1",   "description": "friction velocity"}, \
                   {"name": "alvsf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "60 degree vis albedo with strong cosz dependency"}, \
                   {"name": "alnsf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "60 degree nir albedo with strong cosz dependency"}, \
                   {"name": "alvwf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "60 degree vis albedo with weak cosz dependency"}, \
                   {"name": "alnwf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "60 degree nir albedo with weak cosz dependency"}, \
                   {"name": "facsf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "fractional coverage with strong cosz dependency"}, \
                   {"name": "facwf",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "fractional coverage with weak cosz dependency"}, \
                   {"name": "sheleg",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "water equivalent accumulated snow depth", "alias": "weasd"}, \
                   {"name": "f10m",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "ratio of sigma level 1 wind and 10m wind"}, \
                   {"name": "t2m",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "2-meter absolute temperature"}, \
                   {"name": "q2m",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "kg kg-1", "description": "2-meter specific humidity"}, \
                   {"name": "ffmm",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "Monin-Obukhov similarity function for momentum"}, \
                   {"name": "ffhh",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "Monin-Obukhov similarity function for heat"}, \
                   {"name": "tprcp",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m",       "description": "instantaneous total precipitation amount"}, \
                   {"name": "srflag",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "snow/rain flag for precipitation"}, \
                   {"name": "tsfcl",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "surface skin temperature over land"}, \
                   # NoahMP
                   {"name": "tvxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "vegetation temperature for NoahMP"}, \
                   {"name": "tgxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "ground temperature for NoahMP"}, \
                   {"name": "tahxy",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "K",       "description": "canopy air temperature for NoahMP"}, \
                   {"name": "canicexy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "canopy intercepted ice mass for NoahMP"}, \
                   {"name": "canliqxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "canopy intercepted liquid water for NoahMP"}, \
                   {"name": "eahxy",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "Pa",      "description": "canopy air vapor pressure for NoahMP"}, \
                   {"name": "cmxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "surface drag coefficient for momentum for NoahMP"}, \
                   {"name": "chxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "surface exchange coeff heat & moisture for NoahMP"}, \
                   {"name": "fwetxy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "area fraction of canopy that is wetted/snowed for NoahMP"}, \
                   {"name": "sneqvoxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "snow mass at previous time step for NoahMP"}, \
                   {"name": "alboldxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "snow albedo at previous time step for NoahMP"}, \
                   {"name": "qsnowxy",    "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm s-1",  "description": "snow precipitation rate at surface for NoahMP"}, \
                   {"name": "wslakexy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "lake water storage for NoahMP"}, \
                   {"name": "taussxy",    "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "non-dimensional snow age for NoahMP"}, \
                   {"name": "waxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "water storage in aquifer for NoahMP"}, \
                   {"name": "wtxy",       "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "mm",      "description": "water storage in aquifer and saturated soil for NoahMP"}, \
                   {"name": "zwtxy",      "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m",       "description": "water table depth for NoahMP"}, \
                   {"name": "xlaixy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "leaf area index for NoahMP"}, \
                   {"name": "xsaixy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "stem area index for NoahMP"}, \
                   {"name": "lfmassxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "leaf mass for NoahMP"}, \
                   {"name": "stmassxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "stem mass for NoahMP"}, \
                   {"name": "rtmassxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "fine root mass for NoahMP"}, \
                   {"name": "woodxy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "wood mass including woody roots for NoahMP"}, \
                   {"name": "stblcpxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "stable carbon in deep soil for NoahMP"}, \
                   {"name": "fastcpxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "g m-2",   "description": "short-lived carbon in shallow soil for NoahMP"}, \
                   {"name": "smcwtdxy",   "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m3 m-3",  "description": "soil water content between the bottom of the soil and the water table for NoahMP"}, \
                   {"name": "deeprechxy", "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m",       "description": "recharge to or from the water table when deep for NoahMP"}, \
                   {"name": "rechxy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "m",       "description": "recharge to or from the water table when shallow for NoahMP"}, \
                   {"name": "snowxy",     "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "number of snow layers for NoahMP"}, \
                   # RUC LSM
                   {"name": "wetness",    "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "normalized soil wetness for RUC LSM"}, \
                   {"name": "lai",        "type":real_type, "dimd": ('t0', 'lat', 'lon'),                  "units": "none",    "description": "leaf area index for RUC LSM"}]

    var_forcing = [{"name": "w_ls",         "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "m s-1",       "description": "large scale vertical velocity"}, \
                   {"name": "omega",        "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "Pa s-1",      "description": "large scale pressure vertical velocity"}, \
                   {"name": "u_g",          "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "m s-1",       "description": "large scale geostrophic E-W wind"}, \
                   {"name": "v_g",          "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "m s-1",       "description": "large scale geostrophic N-S wind"}, \
                   {"name": "u_nudge",      "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "m s-1",       "description": "E-W wind to nudge toward"}, \
                   {"name": "v_nudge",      "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "m s-1",       "description": "N-S wind to nudge toward"}, \
                   {"name": "T_nudge",      "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "K",           "description": "absolute temperature to nudge toward"}, \
                   {"name": "thil_nudge",   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "K",           "description": "potential temperature to nudge toward"}, \
                   {"name": "qt_nudge",     "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "kg kg-1",     "description": "total water content to nudge toward"}, \
                   {"name": "rad_heating",  "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "K s-1",       "description": "prescribed radiative heating rate", "alias": "dT_dt_rad"}, \
                   {"name": "h_advec_thil", "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "K s-1",       "description": "prescribed theta_il tendency due to horizontal advection", "alias": "h_advec_thetail"}, \
                   {"name": "v_advec_thil", "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "K s-1",       "description": "prescribed theta_il tendency due to vertical advection",   "alias": "v_advec_thetail"}, \
                   {"name": "h_advec_qt",   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "kg kg-1 s-1", "description": "prescribed q_t tendency due to horizontal advection"}, \
                   {"name": "v_advec_qt",   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": "kg kg-1 s-1", "description": "prescribed q_t tendency due to vertical advection"}, \
                   {"name":'tot_advec_T',   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": 'K s-1',       "description": 'Temperature large-scale advection', "alias": "temp_adv"},\
                   {"name":'tot_advec_qv',  "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": 'kg kg-1 s-1', "description": 'Specific humidity large-scale advection', "alias": "qv_adv"},\
                   {"name":'tot_advec_u',   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": 'm s-2',       "description": 'Zonal wind large-scale advection', "alias": "u_adv"},\
                   {"name":'tot_advec_v',   "type":real_type, "dimd": ('time', 'lev', 'lat', 'lon'), "units": 'm s-2',       "description": 'Meridional wind large-scale advection', "alias": "v_adv"}, \
                   {"name": "T_surf",       "type":real_type, "dimd": ('time',        'lat', 'lon'), "units": "K",           "description": "surface temperature"}]

    #
    # Include dynamic forcing tendencies?
    #
    if (add_UFS_dyn_tend):
        var_dict.extend(var_forcing)

    #
    # Include surface forcing from NOAH LSM?
    #
    if (add_UFS_NOAH_lsm):
        var_dict.extend(var_lsm_ics)

    #
    # Write all fields in "var_dict" to SCM input file.
    #
    for var in var_dict:
        #
        # Do we need to rename this variable? (e.g does it have an alias in variable dictionary)
        #
        var_name = var["name"]
        if "alias" in var: var_name = var["alias"]

        #
        # Create variable
        #
        var_temp = nc_file.createVariable(var_name, var["type"], var["dimd"])

        #
        # Set variable attributes (copy from dictionary)
        #
        var_temp.units       = var["units"]
        var_temp.description = var["description"]

        #
        # Set variable (copy from input dictionary? Set to default_value? Only if provided in dictionary. Otherwise, set to missing value)
        #
        if (var_name in dict) or ("alias" in var):
            var_temp[:] = dict[var["name"]]
        else:
            # If field is not present, set to default_value if provided in dictionary, otherwise set to missing_value.
            var_temp[:]     = missing_value 
        if "default_value" in var: var_temp[:] = var["default_value"]

    #
    # Close file
    #
    nc_file.close()

    return

def find_date(forcing_dir, lam):
    
    dyn_filename_pattern = 'atmf*.nc'
    
    dyn_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, dyn_filename_pattern):
          dyn_filenames.append(f_name)
    if not dyn_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(dyn_filename_pattern,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    dyn_filenames = sorted(dyn_filenames)
    
    nc_file = Dataset('{0}/{1}'.format(forcing_dir,dyn_filenames[0]))
    
    #starting date is in the units attribute of time
    
    date_string = nc_file['time'].getncattr('units').split('since ')[1] #should be in format YYYY-MM-DD HH:MM:SS
    
    nc_file.close()
    
    date_dict = {}
    date_dict["year"] = int(date_string[0:4])
    date_dict["month"] = int(date_string[5:7])
    date_dict["day"] = int(date_string[8:10])
    date_dict["hour"] = int(date_string[11:13])
    date_dict["minute"] = int(date_string[14:16])
    date_dict["second"] = int(date_string[17:])
    
    return date_dict
    
def main():
    setup_logging()
    
    #
    # Read in arguments
    #
    (location, indices, date, in_dir, grid_dir, forcing_dir, tile, area, noahmp, case_name, old_chgres, lam, lami, lamf, save_comp, add_UFS_NOAH_lsm) = parse_arguments()
    
    #
    # Are we provided forcing data?
    #
    add_UFS_dyn_tend  = False
    forcing_data = {}
    if forcing_dir: 
        print("Forcing data directory provided. Adding dyanmic tendencies to SCM input file.")
        add_UFS_dyn_tend = True
    
    #
    # Find tile containing the point using the supergrid if no tile is specified 
    #
    if not tile:
        tile = int(find_tile(location, grid_dir, lam))
        if tile < 0:
            message = 'No tile was found for location {0}'.format(location)
            logging.critical(message)
            raise Exception(message)
        print('Tile found: {0}'.format(tile))
    
    #
    # Find index of closest point in the tile if indices are not specified
    #
    if not indices:
        (tile_j, tile_i, point_lon, point_lat, dist_min) = find_loc_indices(location, grid_dir, tile, lam)
        print('The closest point in tile {0} has indices [{1},{2}]'.format(tile,tile_i,tile_j))
        print('This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat))
        print('This grid cell is approximately {0} km away from the desired location of {1} {2}'.format(dist_min/1.0E3,location[0],location[1]))
    else:
        tile_i = indices[0]
        tile_j = indices[1]
        #still need to grab the lon/lat if the tile and indices are supplied
        (point_lon, point_lat) = find_lon_lat_of_indices(indices, grid_dir, tile, lam)
        
        print('This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat))
    
    #
    # Get UFS IC data (TODO: flag to read in RESTART data rather than IC data and implement different file reads)
    #
    (state_data, surface_data, oro_data) = get_UFS_IC_data(in_dir, grid_dir, forcing_dir, tile, tile_i, tile_j, old_chgres, lam)
    
    if not date:
        #date was not included on command line; look in atmf* file for initial date
        date = find_date(forcing_dir, lam)
    
    #
    # Cold start NoahMP variables (shouldn't be necessary to use this anymore, since same capability exists in SCM code if given Noah ICs only)
    #
    if (noahmp):
        surface_data = add_noahmp_coldstart(surface_data, date)

    #
    # Get grid cell area if not given
    #
    if not area:
        area = get_UFS_grid_area(grid_dir, tile, tile_i, tile_j, lam)
    surface_data["area"] = area
    surface_data["lon"]  = point_lon
    surface_data["lat"]  = point_lat
    
    #
    # Get UFS forcing data
    #
    if (add_UFS_dyn_tend):
        (forcing_data, vars_comp, state_comp) = get_UFS_forcing_data2(state_data["nlevs"], state_data, forcing_dir, grid_dir, tile, tile_i, tile_j, lam, lami, lamf, save_comp)
    else:
        forcing_data = {}
        forcing_data["ps_forc"]       = state_data["p_surf"]
        forcing_data["pressure_forc"] = state_data["pres"][:]
        forcing_data["height_forc"]   = 1
        vars_comp    = {}
        state_comp   = state_data

    #
    # Write SCM case file
    #
    write_SCM_case_file(state_data, surface_data, oro_data, forcing_data, case_name, date, add_UFS_dyn_tend, add_UFS_NOAH_lsm)

    #
    # Create comparison file for UFS/SCM diagnostic scripts.
    #
    if (save_comp):
        fileOUT = write_comparison_file(vars_comp, state_comp, case_name, date, surface_data, add_UFS_dyn_tend, add_UFS_NOAH_lsm)

        # Write name of output file to text file (used by ensemble wrapper)
        fileID  = open("nameout.txt", 'w')
        fileID.write(fileOUT)
        fileID.write('\n')
        fileID.close()

if __name__ == '__main__':
    main()
