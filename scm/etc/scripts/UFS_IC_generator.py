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
group1.add_argument('-l', '--location',   help='longitude and latitude in degress E and N, respectively, separated by a space', nargs=2, type=float)
group1.add_argument('-ij','--index',      help='i,j indices within the tile (if known - bypasses search for closest model point to lon/lat location)', nargs=2, type=int)
parser.add_argument('-d', '--date',       help='date corresponding to initial conditions in YYYYMMDDHHMMSS format', required=False)
parser.add_argument('-i', '--in_dir',     help='input directory path containing FV3 input files', required=True)
parser.add_argument('-g', '--grid_dir',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-f', '--forcing_dir',help='directory path containing physics diag files', required=True)
parser.add_argument('-t', '--tile',       help='tile of desired point (if known - bypasses tile search if present)', type=int, choices=range(1,8))
parser.add_argument('-a', '--area',       help='area of grid cell in m^2', type=float)
parser.add_argument('-mp','--noahmp',     help='flag to generate cold-start ICs for NoahMP LSM from Noah LSM ICs', action='store_true')
parser.add_argument('-lam','--lam',       help='flag to signal that the ICs and forcing is from a limited-area model run', action='store_true')
parser.add_argument('-n', '--case_name',  help='name of case', required=True)
parser.add_argument('-oc','--old_chgres', help='flag to denote that the initial conditions use an older data format (pre-chgres_cube)', action='store_true')
parser.add_argument('-nf','--no_force',   help='flag to additionally write out a case data file with forcing turned off', action='store_true')
parser.add_argument('-sc','--save_comp',  help='flag to save and write out a file with UFS output data to compare SCM simulations with', action='store_true')
parser.add_argument('-lsm', '--add_UFS_NOAH_lsm', help='flag to include UFS NOAH LSM surface forcing', action='store_true')
parser.add_argument('-ufsf','--add_UFS_dyn_tend', help='flag to include UFS dynamic tendencies for SCM forcing', action='store_true')

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
    no_force = args.no_force
    save_comp = args.save_comp
    add_UFS_NOAH_lsm = args.add_UFS_NOAH_lsm
    add_UFS_dyn_tend = args.add_UFS_dyn_tend
    
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
    
    return (location, index, date_dict, in_dir, grid_dir, forcing_dir, tile, \
            area, noahmp, case_name, old_chgres, lam, no_force, save_comp,   \
            add_UFS_NOAH_lsm, add_UFS_dyn_tend)

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
        "nlevs": nlevs_model,
        "zh": zh_rev[nlevs_model:0:-1],
        "ua": u_model_rev[0][::-1],
        "va": v_model_rev[0][::-1],
        "qv": sphum_model_rev[0][::-1],
        "o3": o3_model_rev[0][::-1],
        "ql": liqwat_model_rev[0][::-1],
        "qi": icewat_model_rev[::-1],
        "ps": ps_calc,
        "ta": temp_model_rev[0,::-1],
        "pa": pressure_model,
        "pa_i": pressure_model_interfaces
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
            print('unknown grid orientation')
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
            print('unknown grid orientation')
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
            print('unknown grid orientation')
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
            print('unknown grid orientation')
    
    
    nc_file_grid.close()
    
    return [u_s, u_n, v_w, v_e]

def get_UFS_surface_data(dir, tile, i, j, old_chgres, lam):
    """Get the surface data for the given tile and indices"""
    
    if lam:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.nc'))
    else:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.tile{0}.nc'.format(tile)))
    
    #FV3/io/FV3GFS_io.F90/sfc_prop_restart_read was used as reference for variables that can be read in    
    
    #read in scalars (would be 2D variables in a 3D model)
    
    # surface properties (assuming Noah LSM; may contain variables needed for fractional land fraction)
    tsfco_in = read_NetCDF_surface_var(nc_file, 'tsea', i, j, old_chgres, 0)
    tg3_in = read_NetCDF_surface_var(nc_file, 'tg3', i, j, old_chgres, 0)
    uustar_in = read_NetCDF_surface_var(nc_file, 'uustar', i, j, old_chgres, 0)
    alvsf_in = read_NetCDF_surface_var(nc_file, 'alvsf', i, j, old_chgres, 0)
    alvwf_in = read_NetCDF_surface_var(nc_file, 'alvwf', i, j, old_chgres, 0)
    alnsf_in = read_NetCDF_surface_var(nc_file, 'alnsf', i, j, old_chgres, 0)
    alnwf_in = read_NetCDF_surface_var(nc_file, 'alnwf', i, j, old_chgres, 0)
    facsf_in = read_NetCDF_surface_var(nc_file, 'facsf', i, j, old_chgres, 0)
    facwf_in = read_NetCDF_surface_var(nc_file, 'facwf', i, j, old_chgres, 0)
    styp_in = read_NetCDF_surface_var(nc_file, 'stype', i, j, old_chgres, 0)
    slope_in = read_NetCDF_surface_var(nc_file, 'slope', i, j, old_chgres, 0)
    vtyp_in = read_NetCDF_surface_var(nc_file, 'vtype', i, j, old_chgres, 0)
    vfrac_in = read_NetCDF_surface_var(nc_file, 'vfrac', i, j, old_chgres, 0)
    shdmin_in = read_NetCDF_surface_var(nc_file, 'shdmin', i, j, old_chgres, 0)
    shdmax_in = read_NetCDF_surface_var(nc_file, 'shdmax', i, j, old_chgres, 0)
    zorlo_in = read_NetCDF_surface_var(nc_file, 'zorl', i, j, old_chgres, 0)
    slmsk_in = read_NetCDF_surface_var(nc_file, 'slmsk', i, j, old_chgres, 0)
    canopy_in = read_NetCDF_surface_var(nc_file, 'canopy', i, j, old_chgres, 0)
    hice_in = read_NetCDF_surface_var(nc_file, 'hice', i, j, old_chgres, 0)
    fice_in = read_NetCDF_surface_var(nc_file, 'fice', i, j, old_chgres, 0)
    tisfc_in = read_NetCDF_surface_var(nc_file, 'tisfc', i, j, old_chgres, 0)
    snwdph_in = read_NetCDF_surface_var(nc_file, 'snwdph', i, j, old_chgres, 0)
    snoalb_in = read_NetCDF_surface_var(nc_file, 'snoalb', i, j, old_chgres, 0)
    sheleg_in = read_NetCDF_surface_var(nc_file, 'sheleg', i, j, old_chgres, 0)
    f10m_in = read_NetCDF_surface_var(nc_file, 'f10m', i, j, old_chgres, 0)
    t2m_in = read_NetCDF_surface_var(nc_file, 't2m', i, j, old_chgres, 0)
    q2m_in = read_NetCDF_surface_var(nc_file, 'q2m', i, j, old_chgres, 0)
    ffmm_in = read_NetCDF_surface_var(nc_file, 'ffmm', i, j, old_chgres, 0)
    ffhh_in = read_NetCDF_surface_var(nc_file, 'ffhh', i, j, old_chgres, 0)
    tprcp_in = read_NetCDF_surface_var(nc_file, 'tprcp', i, j, old_chgres, 0)
    srflag_in = read_NetCDF_surface_var(nc_file, 'srflag', i, j, old_chgres, 0)
    sncovr_in = read_NetCDF_surface_var(nc_file, 'sncovr', i, j, old_chgres, 0)
    tsfcl_in = read_NetCDF_surface_var(nc_file, 'tsfcl', i, j, old_chgres, 0)
    zorll_in = read_NetCDF_surface_var(nc_file, 'zorll', i, j, old_chgres, 0)   
    zorli_in = read_NetCDF_surface_var(nc_file, 'zorli', i, j, old_chgres, 0)
    
    #present when cplwav = T
    zorlw_in = read_NetCDF_surface_var(nc_file, 'zorlw', i, j, old_chgres, 0)
    
    #NSST variables that may be in the surface file
    tref_in = read_NetCDF_surface_var(nc_file, 'tref', i, j, old_chgres, 0)
    z_c_in = read_NetCDF_surface_var(nc_file, 'z_c', i, j, old_chgres, 0)
    c_0_in = read_NetCDF_surface_var(nc_file, 'c_0', i, j, old_chgres, 0)
    c_d_in = read_NetCDF_surface_var(nc_file, 'c_d', i, j, old_chgres, 0)
    w_0_in = read_NetCDF_surface_var(nc_file, 'w_0', i, j, old_chgres, 0)
    w_d_in = read_NetCDF_surface_var(nc_file, 'w_d', i, j, old_chgres, 0)
    xt_in = read_NetCDF_surface_var(nc_file, 'xt', i, j, old_chgres, 0)
    xs_in = read_NetCDF_surface_var(nc_file, 'xs', i, j, old_chgres, 0)
    xu_in = read_NetCDF_surface_var(nc_file, 'xu', i, j, old_chgres, 0)
    xv_in = read_NetCDF_surface_var(nc_file, 'xv', i, j, old_chgres, 0)
    xz_in = read_NetCDF_surface_var(nc_file, 'xz', i, j, old_chgres, 0)
    zm_in = read_NetCDF_surface_var(nc_file, 'zm', i, j, old_chgres, 0)
    xtts_in = read_NetCDF_surface_var(nc_file, 'xtts', i, j, old_chgres, 0)
    xzts_in = read_NetCDF_surface_var(nc_file, 'xzts', i, j, old_chgres, 0)
    d_conv_in = read_NetCDF_surface_var(nc_file, 'd_conv', i, j, old_chgres, 0)
    ifd_in = read_NetCDF_surface_var(nc_file, 'ifd', i, j, old_chgres, 0)
    dt_cool_in = read_NetCDF_surface_var(nc_file, 'dt_cool', i, j, old_chgres, 0)
    qrain_in = read_NetCDF_surface_var(nc_file, 'qrain', i, j, old_chgres, 0)

    #NoahMP variables that may be in the surface file
    snowxy_in = read_NetCDF_surface_var(nc_file, 'snowxy', i, j, old_chgres, 0)
    tvxy_in = read_NetCDF_surface_var(nc_file, 'tvxy', i, j, old_chgres, 0)
    tgxy_in = read_NetCDF_surface_var(nc_file, 'tgxy', i, j, old_chgres, 0)
    canicexy_in = read_NetCDF_surface_var(nc_file, 'canicexy', i, j, old_chgres, 0)
    canliqxy_in = read_NetCDF_surface_var(nc_file, 'canliqxy', i, j, old_chgres, 0)
    eahxy_in = read_NetCDF_surface_var(nc_file, 'eahxy', i, j, old_chgres, 0)
    tahxy_in = read_NetCDF_surface_var(nc_file, 'tahxy', i, j, old_chgres, 0)
    cmxy_in = read_NetCDF_surface_var(nc_file, 'cmxy', i, j, old_chgres, 0)
    chxy_in = read_NetCDF_surface_var(nc_file, 'chxy', i, j, old_chgres, 0)
    fwetxy_in = read_NetCDF_surface_var(nc_file, 'fwetxy', i, j, old_chgres, 0)
    sneqvoxy_in = read_NetCDF_surface_var(nc_file, 'sneqvoxy', i, j, old_chgres, 0)
    alboldxy_in = read_NetCDF_surface_var(nc_file, 'alboldxy', i, j, old_chgres, 0)
    qsnowxy_in = read_NetCDF_surface_var(nc_file, 'qsnowxy', i, j, old_chgres, 0)
    wslakexy_in = read_NetCDF_surface_var(nc_file, 'wslakexy', i, j, old_chgres, 0)
    zwtxy_in = read_NetCDF_surface_var(nc_file, 'zwtxy', i, j, old_chgres, 0)
    waxy_in = read_NetCDF_surface_var(nc_file, 'waxy', i, j, old_chgres, 0)
    wtxy_in = read_NetCDF_surface_var(nc_file, 'wtxy', i, j, old_chgres, 0)
    lfmassxy_in = read_NetCDF_surface_var(nc_file, 'lfmassxy', i, j, old_chgres, 0)
    rtmassxy_in = read_NetCDF_surface_var(nc_file, 'rtmassxy', i, j, old_chgres, 0)
    stmassxy_in = read_NetCDF_surface_var(nc_file, 'stmassxy', i, j, old_chgres, 0)
    woodxy_in = read_NetCDF_surface_var(nc_file, 'woodxy', i, j, old_chgres, 0)
    stblcpxy_in = read_NetCDF_surface_var(nc_file, 'stblcpxy', i, j, old_chgres, 0)
    fastcpxy_in = read_NetCDF_surface_var(nc_file, 'fastcpxy', i, j, old_chgres, 0)
    xsaixy_in = read_NetCDF_surface_var(nc_file, 'xsaixy', i, j, old_chgres, 0)
    xlaixy_in = read_NetCDF_surface_var(nc_file, 'xlaixy', i, j, old_chgres, 0)
    taussxy_in = read_NetCDF_surface_var(nc_file, 'taussxy', i, j, old_chgres, 0)
    smcwtdxy_in = read_NetCDF_surface_var(nc_file, 'smcwtdxy', i, j, old_chgres, 0)
    deeprechxy_in = read_NetCDF_surface_var(nc_file, 'deeprechxy', i, j, old_chgres, 0)
    rechxy_in = read_NetCDF_surface_var(nc_file, 'rechxy', i, j, old_chgres, 0)
    
    # RUC LSM variables
    wetness_in           = read_NetCDF_surface_var(nc_file, 'wetness',           i, j, old_chgres, 0)
    clw_surf_land_in     = read_NetCDF_surface_var(nc_file, 'clw_surf_land',     i, j, old_chgres, 0)
    clw_surf_ice_in      = read_NetCDF_surface_var(nc_file, 'clw_surf_ice',      i, j, old_chgres, 0)
    qwv_surf_land_in     = read_NetCDF_surface_var(nc_file, 'qwv_surf_land',     i, j, old_chgres, 0)
    qwv_surf_ice_in      = read_NetCDF_surface_var(nc_file, 'qwv_surf_ice',      i, j, old_chgres, 0)
    tsnow_land_in        = read_NetCDF_surface_var(nc_file, 'tsnow_land',        i, j, old_chgres, 0)
    tsnow_ice_in         = read_NetCDF_surface_var(nc_file, 'tsnow_ice',         i, j, old_chgres, 0)
    snowfall_acc_land_in = read_NetCDF_surface_var(nc_file, 'snowfall_acc_land', i, j, old_chgres, 0)
    snowfall_acc_ice_in  = read_NetCDF_surface_var(nc_file, 'snowfall_acc_ice',  i, j, old_chgres, 0)
    sfalb_lnd_in         = read_NetCDF_surface_var(nc_file, 'sfalb_lnd',         i, j, old_chgres, 0)
    sfalb_lnd_bck_in     = read_NetCDF_surface_var(nc_file, 'sfalb_lnd_bck',     i, j, old_chgres, 0)
    sfalb_ice_in         = read_NetCDF_surface_var(nc_file, 'sfalb_ice',         i, j, old_chgres, 0)
    lai_in               = read_NetCDF_surface_var(nc_file, 'lai',               i, j, old_chgres, 0)
    albdirvis_ice_in     = read_NetCDF_surface_var(nc_file, 'albdirvis_ice',     i, j, old_chgres, 0)
    albdirnir_ice_in     = read_NetCDF_surface_var(nc_file, 'albdirnir_ice',     i, j, old_chgres, 0)
    albdifvis_ice_in     = read_NetCDF_surface_var(nc_file, 'albdifvis_ice',     i, j, old_chgres, 0)
    albdifnir_ice_in     = read_NetCDF_surface_var(nc_file, 'albdifnir_ice',     i, j, old_chgres, 0)
    emis_lnd_in          = read_NetCDF_surface_var(nc_file, 'emis_lnd',          i, j, old_chgres, 0)
    emis_ice_in          = read_NetCDF_surface_var(nc_file, 'emis_ice',          i, j, old_chgres, 0)
    
    #read in profiles (would be 3D variables in a 3D model)
    
    #land_state
    stc_in = read_NetCDF_surface_var(nc_file, 'stc', i, j, old_chgres, missing_variable_soil_layers)
    smc_in = read_NetCDF_surface_var(nc_file, 'smc', i, j, old_chgres, missing_variable_soil_layers)
    slc_in = read_NetCDF_surface_var(nc_file, 'slc', i, j, old_chgres, missing_variable_soil_layers)
    
    #NoahMP 3D variables
    snicexy_in = read_NetCDF_surface_var(nc_file, 'snicexy', i, j, old_chgres, missing_variable_snow_layers)
    snliqxy_in = read_NetCDF_surface_var(nc_file, 'snliqxy', i, j, old_chgres, missing_variable_snow_layers)
    tsnoxy_in = read_NetCDF_surface_var(nc_file, 'tsnoxy', i, j, old_chgres, missing_variable_snow_layers)
    smoiseq_in = read_NetCDF_surface_var(nc_file, 'smoiseq', i, j, old_chgres, missing_variable_soil_layers)
    zsnsoxy_in = read_NetCDF_surface_var(nc_file, 'zsnsoxy', i, j, old_chgres, missing_variable_soil_layers + missing_variable_snow_layers)
     
    #RUC LSM 3D variables
    tslb_in = read_NetCDF_surface_var(nc_file, 'tslb', i, j, old_chgres, missing_variable_soil_layers)
    smois_in = read_NetCDF_surface_var(nc_file, 'smois', i, j, old_chgres, missing_variable_soil_layers)
    sh2o_in = read_NetCDF_surface_var(nc_file, 'sh2o', i, j, old_chgres, missing_variable_soil_layers)
    smfr_in = read_NetCDF_surface_var(nc_file, 'smfr', i, j, old_chgres, missing_variable_soil_layers)
    flfr_in = read_NetCDF_surface_var(nc_file, 'flfr', i, j, old_chgres, missing_variable_soil_layers)
    
    #fractional grid 3D variables
    tiice_in = read_NetCDF_surface_var(nc_file, 'tiice', i, j, old_chgres, missing_variable_ice_layers)

    #print("zorlw_in = {}".format(zorlw_in))
    
    nc_file.close()
    
    #put data in a dictionary
    surface = {
        #Noah LSM
        "tsfco": tsfco_in,  
        "tg3": tg3_in,
        "ustar": uustar_in,
        "alvsf": alvsf_in,
        "alvwf": alvwf_in,
        "alnsf": alnsf_in,
        "alnwf": alnwf_in,
        "facsf": facsf_in,
        "facwf": facwf_in,
        "soiltyp": styp_in,
        "slopetyp": slope_in,
        "vegtyp": vtyp_in,
        "vegfrac": vfrac_in,
        "shdmin": shdmin_in,
        "shdmax": shdmax_in,
        "z0": zorlo_in,
        "slmsk": slmsk_in,
        "canopy": canopy_in,
        "hice": hice_in,
        "fice": fice_in,
        "tisfc": tisfc_in,
        "snwdph": snwdph_in,
        "snoalb": snoalb_in,
        "weasd": sheleg_in,
        "f10m": f10m_in,
        "t2m": t2m_in,
        "q2m": q2m_in,
        "ffmm": ffmm_in,
        "ffhh": ffhh_in,
        "tprcp": tprcp_in,
        "srflag": srflag_in,
        "sncovr": sncovr_in,
        "tsfcl": tsfcl_in,
        "zorll": zorll_in,
        "zorli": zorli_in,
        #cplwav
        "zorlw": zorlw_in,
        #NSST
        "tref": tref_in,
        "z_c": z_c_in,
        "c_0": c_0_in,
        "c_d": c_d_in,
        "w_0": w_0_in,
        "w_d": w_d_in,
        "xt": xt_in,
        "xs": xs_in,
        "xu": xu_in,
        "xv": xv_in,
        "xz": xz_in,
        "zm": zm_in,
        "xtts": xtts_in,
        "xzts": xzts_in,
        "d_conv": d_conv_in,
        "ifd": ifd_in,
        "dt_cool": dt_cool_in,
        "qrain": qrain_in,
        #NoahMP
        "snowxy": snowxy_in,
        "tvxy": tvxy_in,
        "tgxy": tgxy_in,
        "canicexy": canicexy_in,
        "canliqxy": canliqxy_in,
        "eahxy": eahxy_in,
        "tahxy": tahxy_in,
        "cmxy": cmxy_in,
        "chxy": chxy_in,
        "fwetxy": fwetxy_in,
        "sneqvoxy": sneqvoxy_in,
        "alboldxy": alboldxy_in,
        "qsnowxy": qsnowxy_in,
        "wslakexy": wslakexy_in,
        "zwtxy": zwtxy_in,
        "waxy": waxy_in,
        "wtxy": wtxy_in,
        "lfmassxy": lfmassxy_in,
        "rtmassxy": rtmassxy_in,
        "stmassxy": stmassxy_in,
        "woodxy": woodxy_in,
        "stblcpxy": stblcpxy_in,
        "fastcpxy": fastcpxy_in,
        "xsaixy": xsaixy_in,
        "xlaixy": xlaixy_in,
        "taussxy": taussxy_in,
        "smcwtdxy": smcwtdxy_in,
        "deeprechxy": deeprechxy_in,
        "rechxy": rechxy_in,
        # RUC LSM 2D
        "wetness": wetness_in,
        "clw_surf_land": clw_surf_land_in,
        "clw_surf_ice": clw_surf_ice_in,
        "qwv_surf_land": qwv_surf_land_in,
        "qwv_surf_ice": qwv_surf_ice_in,
        "tsnow_land": tsnow_land_in,
        "tsnow_ice": tsnow_ice_in,
        "snowfall_acc_land": snowfall_acc_land_in,
        "snowfall_acc_ice": snowfall_acc_ice_in,
        "sfalb_lnd": sfalb_lnd_in,
        "sfalb_lnd_bck": sfalb_lnd_bck_in,
        "sfalb_ice": sfalb_ice_in,
        "lai": lai_in,
        "albdirvis_ice": albdirvis_ice_in,
        "albdirnir_ice": albdirnir_ice_in,
        "albdifvis_ice": albdifvis_ice_in,
        "albdifnir_ice": albdifnir_ice_in,
        "emis_lnd": emis_lnd_in,
        "emis_ice": emis_ice_in,
        #Noah LSM 3D
        "stc": stc_in,
        "smc": smc_in,
        "slc": slc_in,
        #NoahMP LSM 3D
        "snicexy": snicexy_in,
        "snliqxy": snliqxy_in,
        "tsnoxy": tsnoxy_in,
        "smoiseq": smoiseq_in,
        "zsnsoxy": zsnsoxy_in,         
        #RUC LSM 3D variables
        "tslb": tslb_in,
        "smois": smois_in,
        "sh2o": sh2o_in,
        "smfr": smfr_in,
        "flfr": flfr_in,        
        #fractional grid 3D variables
        "tiice": tiice_in,
    }
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
    
    
    # orographyic properties
    stddev_in = read_NetCDF_var(nc_file, "stddev", i, j)
    convexity_in = read_NetCDF_var(nc_file, "convexity", i, j)
    oa1_in = read_NetCDF_var(nc_file, "oa1", i, j)
    oa2_in = read_NetCDF_var(nc_file, "oa2", i, j)
    oa3_in = read_NetCDF_var(nc_file, "oa3", i, j)
    oa4_in = read_NetCDF_var(nc_file, "oa4", i, j)
    ol1_in = read_NetCDF_var(nc_file, "ol1", i, j)
    ol2_in = read_NetCDF_var(nc_file, "ol2", i, j)
    ol3_in = read_NetCDF_var(nc_file, "ol3", i, j)
    ol4_in = read_NetCDF_var(nc_file, "ol4", i, j)
    theta_in = read_NetCDF_var(nc_file, "theta", i, j)
    gamma_in = read_NetCDF_var(nc_file, "gamma", i, j)
    sigma_in = read_NetCDF_var(nc_file, "sigma", i, j)
    elvmax_in = read_NetCDF_var(nc_file, "elvmax", i, j)
    orog_filt_in = read_NetCDF_var(nc_file, "orog_filt", i, j)
    orog_raw_in = read_NetCDF_var(nc_file, "orog_raw", i, j)
    #fractional landmask variables
    land_frac_in = read_NetCDF_var(nc_file, "land_frac", i, j)
    #lake variables (optional)
    lake_frac_in = read_NetCDF_var(nc_file, "lake_frac", i, j)
    lake_depth_in = read_NetCDF_var(nc_file, "lake_depth", i, j)
    
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
        "elvmax": elvmax_in,
        "oro": orog_filt_in,
        "oro_uf": orog_raw_in,
        "landfrac": land_frac_in,
        "lakefrac": lake_frac_in,
        "lakedepth": lake_depth_in
    }
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

def search_in_dict(listin,name):
    for count, dictionary in enumerate(listin):
        if dictionary["name"] == name:
            return count

def get_UFS_forcing_data(nlevs, state_IC, forcing_dir, grid_dir, tile, i, j, lam, save_comp_data):
    """Get the horizontal and vertical advective tendencies for the given tile and indices"""
    
    regrid_output = 'point'
    #regrid_output = 'all'
    
    orig_SCM_format = False
    
    #if lam:
    dyn_filename_pattern = 'atmf*.nc'
    phy_filename_pattern = 'sfcf*.nc'
    #else:
    #    dyn_filename_pattern = 'atmf*.tile{0}.nc'.format(tile)
    #    phy_filename_pattern = 'sfcf*.tile{0}.nc'.format(tile)
    
    dyn_filenames = []
    phy_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, dyn_filename_pattern):
          dyn_filenames.append(f_name)
       if fnmatch.fnmatch(f_name, phy_filename_pattern):
          phy_filenames.append(f_name)
    if not dyn_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(dyn_filename_pattern,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    if not phy_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(phy_filename_pattern,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    dyn_filenames = sorted(dyn_filenames)
    phy_filenames = sorted(phy_filenames)
    
    if (len(dyn_filenames) != len(phy_filenames)):
        message = 'The number of dyn files and phy files in {0} matching the patterns does not match.'.format(forcing_dir)
        logging.critical(message)
        raise Exception(message)
    
    n_files = len(dyn_filenames)
    
    kord_tm = -9
    kord_mt = 9
    kord_tr = 9
    t_min = 184.0
    q_min = 0.0
    
    ####################################################################################
    #
    # Read in atmospheric state_UFS (atmf*.nc)
    #
    ####################################################################################
    ps     = []
    p_lev  = []
    p_lay  = []
    t_lay  = []
    qv_lay = []
    u_lay  = []
    v_lay  = []
    time_dyn_hours = []
    (ic_grid_lon, ic_grid_lat) = get_initial_lon_lat_grid(grid_dir, tile, lam)
    for count, filename in enumerate(dyn_filenames, start=1):
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)
        
        #check if output grid is different than initial (native) grid
        try:
            data_grid_lon = nc_file['lon'][:,:]
            data_grid_lat = nc_file['lat'][:,:]
        except:
            data_grid_lon = nc_file['grid_xt'][:,:]
            data_grid_lat = nc_file['grid_yt'][:,:]
        
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
            
            print('Regridding {} onto native grid: regridding progress = {}%'.format(filename, 100.0*count/(len(dyn_filenames) + len(phy_filenames))))
            regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
            #print(regridder.weights)            
            
            ps_data = regridder(nc_file['pressfc'][0,:,:])
            
            t_data = regridder(nc_file['tmp'][0,::-1,:,:])
            qv_data = regridder(nc_file['spfh'][0,::-1,:,:])
            u_data = regridder(nc_file['ugrd'][0,::-1,:,:])
            v_data = regridder(nc_file['vgrd'][0,::-1,:,:])
        else:
            ps_data = nc_file['pressfc'][0,:,:]
            t_data = nc_file['tmp'][0,::-1,:,:]
            qv_data = nc_file['spfh'][0,::-1,:,:]
            u_data = nc_file['ugrd'][0,::-1,:,:]
            v_data = nc_file['vgrd'][0,::-1,:,:]
            i_get = i
            j_get = j
        
        nlevs=len(nc_file.dimensions['pfull'])
        
        ak = getattr(nc_file, "ak")[::-1]
        bk = getattr(nc_file, "bk")[::-1]
    
        ps.append(ps_data[j_get,i_get])
        
        p_interface = np.zeros(nlevs+1)
        for k in range(nlevs+1):
            p_interface[k]=ak[k]+ps[-1]*bk[k]
        
        p_lev.append(p_interface)
        
        p_layer = np.zeros(nlevs)
        for k in range(nlevs):
            p_layer[k] = ((1.0/(rocp+1.0))*(p_interface[k]**(rocp+1.0) - p_interface[k+1]**(rocp+1.0))/(p_interface[k] - p_interface[k+1]))**(1.0/rocp)
        
        p_lay.append(p_layer)
        
        t_lay.append(t_data[:,j_get,i_get])
        qv_lay.append(qv_data[:,j_get,i_get])
        u_lay.append(u_data[:,j_get,i_get])
        v_lay.append(v_data[:,j_get,i_get])
            
        time_dyn_hours.append(nc_file['time'][0])
        
        nc_file.close()
    ps = np.asarray(ps)
    p_lev = np.asarray(p_lev)
    p_lay = np.asarray(p_lay)
    t_lay = np.asarray(t_lay)
    qv_lay = np.asarray(qv_lay)
    u_lay = np.asarray(u_lay)
    v_lay = np.asarray(v_lay)
    time_dyn_hours = np.asarray(time_dyn_hours)
    tv_lay = t_lay*(1.0 + zvir*qv_lay)

    ####################################################################################
    #
    # Read in tendencies (sfcf*.nc)
    #
    ####################################################################################

    # Dictionary containing dynamic tendency names.
    vars_dyn  =[{"name":"dtend_qv_nophys"},  {"name":"dtend_temp_nophys"}, \
                {"name":"dtend_u_nophys"},   {"name":"dtend_v_nophys"}]

    # Dictionary of diagnostics to save for compariaion to SCM output.
    vars_comp= [{"name":"dtend_qv_pbl"},     {"name":"dtend_temp_pbl"},    \
                {"name":"dtend_u_pbl"},      {"name":"dtend_v_pbl"},       \
                {"name":"dtend_qv_deepcnv"}, {"name":"dtend_temp_deepcnv"},\
                {"name":"dtend_u_deepcnv"},  {"name":"dtend_v_deepcnv"},   \
                {"name":"dtend_qv_shalcnv"}, {"name":"dtend_temp_shalcnv"},\
                {"name":"dtend_u_shalcnv"},  {"name":"dtend_v_shalcnv"},   \
                {"name":"dtend_temp_lw"},    {"name":"dtend_temp_sw"},     \
                {"name":"dtend_qv_mp"},      {"name":"dtend_temp_mp"},     \
                {"name":"dtend_qv_cnvgwd"},  {"name":"dtend_temp_cnvgwd"}, \
                {"name":"dtend_u_cnvgwd"},   {"name":"dtend_v_cnvgwd"},    \
                {"name":"dtend_qv_phys"},    {"name":"dtend_temp_phys"},   \
                {"name":"dtend_u_phys"},     {"name":"dtend_v_phys"}]

    # Initialize dictionary "values" to empty list.
    for vars in vars_dyn: vars["values"] = []
    for vars in vars_comp: vars["values"] = []

    # Read in each sfcf***.nc
    time_phys_hours = []
    for count, filename in enumerate(phy_filenames, start=1):
        # Open file
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)
        
        # check if output grid is different than initial (native) grid
        try:
            data_grid_lon = nc_file['lon'][:,:]
            data_grid_lat = nc_file['lat'][:,:]
        except:
            data_grid_lon = nc_file['grid_xt'][:,:]
            data_grid_lat = nc_file['grid_yt'][:,:]
        
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
            
            print('Regridding {} onto native grid: regridding progress = {}%'.\
                  format(filename, 100.0*(count + len(dyn_filenames))/(len(dyn_filenames) + len(phy_filenames))))
            regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
            
            # Read in each available fields in dictionaries. Copy variable attributes to dictionaries (use
            # for output).
            for dtend in vars_dyn:
                try:
                    data = regridder(nc_file[dtend["name"]][0,::-1,:,:])
                    dtend["values"].append(data[:,j_get,i_get])
                    dtend["units"]     = nc_file[dtend["name"]].getncattr(name="units")
                    dtend["long_name"] = nc_file[dtend["name"]].getncattr(name="long_name")
                except:
                    logging.info(dtend["name"] + ' not found in ' + filename)
                    exit()

            if (save_comp_data):
                for dtend in vars_comp:
                    try:
                        data = regridder(nc_file[dtend["name"]][0,::-1,:,:])
                        dtend["values"].append(data[:,j_get,i_get])
                        dtend["units"]     = nc_file[dtend["name"]].getncattr(name="units")
                        dtend["long_name"] = nc_file[dtend["name"]].getncattr(name="long_name")
                    except:
                        logging.debug(dtend["name"] + ' not found in ' + filename)
        #
        # Regridding not necessary
        #
        else:
            i_get = i
            j_get = j
            # Read in fields. Copy variable attributes
            for dtend in vars_dyn:
                try:
                    data = nc_file[dtend["name"]][0,::-1,:,:]
                    dtend["values"].append(data[:,j_get,i_get])
                    dtend["units"]     = nc_file[dtend["name"]].getncattr(name="units")
                    dtend["long_name"] = nc_file[dtend["name"]].getncattr(name="long_name")
                except:
                    logging.debug(dtend["name"] + ' not found in ' + filename)
            #
            if (save_comp_data):
                for dtend in vars_comp:
                    try:
                        data = nc_file[dtend["name"]][0,::-1,:,:]
                        dtend["values"].append(data[:,j_get,i_get])
                        dtend["units"]     = nc_file[dtend["name"]].getncattr(name="units")
                        dtend["long_name"] = nc_file[dtend["name"]].getncattr(name="long_name")
                    except:
                        logging.debug(dtend["name"] + ' not available in ' + filename)
        
        # Save 
        time_phys_hours.append(nc_file['time'][0])

        # Close file
        nc_file.close()

    # Convert to numpy arrays
    for dtend in vars_dyn:
        try:    dtend["values"] = np.asarray(dtend["values"])
        except: logging.debug(dtend["name"] + ' not available')
    if (save_comp_data):
        for dtend in vars_comp:
            try:    dtend["values"] = np.asarray(dtend["values"])
            except: logging.debug(dtend["name"] + ' not available')
    time_phys_hours = np.asarray(time_phys_hours)

    # Indices for dynamic tendencies in dictionary (used when referencing vars_dyn dict)
    itemp = search_in_dict(vars_dyn,"dtend_temp_nophys")
    iqv   = search_in_dict(vars_dyn,"dtend_qv_nophys")
    iu    = search_in_dict(vars_dyn,"dtend_u_nophys")
    iv    = search_in_dict(vars_dyn,"dtend_v_nophys")

    ####################################################################################
    #
    # Handle forcing from initialization to the first history file
    #
    ####################################################################################
    nlevs           = len(p_lay[0,:])
    dummy           = np.zeros(1)
    tv_temp         = np.zeros([1,nlevs])
    qv_temp         = np.zeros([1,nlevs])
    u_temp          = np.zeros([1,nlevs])
    v_temp          = np.zeros([1,nlevs])
    dtdt_temp       = np.zeros([1,nlevs])
    dqvdt_temp      = np.zeros([1,nlevs])
    dudt_temp       = np.zeros([1,nlevs])
    dvdt_temp       = np.zeros([1,nlevs])
    from_p          = np.zeros([1,nlevs+1])
    to_p            = np.zeros([1,nlevs+1])
    log_from_p      = np.zeros([1,nlevs+1])
    log_to_p        = np.zeros([1,nlevs+1])
    tv_rev          = np.zeros([1,nlevs])
    qv_rev          = np.zeros([1,nlevs])
    dp2             = np.zeros([1,nlevs])
    u_rev           = np.zeros([1,nlevs])
    v_rev           = np.zeros([1,nlevs])
    tv_layers_remap = np.zeros([n_files,nlevs])
    qv_layers_remap = np.zeros([n_files,nlevs])
    u_layers_remap  = np.zeros([n_files,nlevs])
    v_layers_remap  = np.zeros([n_files,nlevs])
    dtdt_adv        = np.zeros([n_files,nlevs])
    dqvdt_adv       = np.zeros([n_files,nlevs])
    dudt_adv        = np.zeros([n_files,nlevs])
    dvdt_adv        = np.zeros([n_files,nlevs])
    valid_pres_adv  = np.zeros([n_files,nlevs])
    valid_pres_i_adv= np.zeros([n_files,nlevs+1])

    # Interpolation range
    from_p[0,:]     = state_IC["pa_i"][::-1]
    to_p[0,:]       = p_lev[0,::-1]
    log_from_p[0,:] = np.log(state_IC["pa_i"][::-1])
    log_to_p[0,:]   = np.log(p_lev[0,::-1])

    # Virtual Temperature @ time = 0
    tv_init      = state_IC["ta"]*(1.0 + zvir*state_IC["qv"])
    tv_init_rev  = tv_init[::-1]
    tv_rev_new   = fv3_remap.map_scalar(nlevs, log_from_p, tv_init_rev[np.newaxis, :], dummy, \
                                        nlevs, log_to_p, 0, 0, 1, np.abs(kord_tm), t_min)
    tv_temp[0,:] = tv_rev_new[0,::-1]

    # Specific humidity @ time = 0
    qv_init_rev  = state_IC["qv"][::-1]
    for k in range(0,nlevs): dp2[0,k] = from_p[0,k+1] - from_p[0,k]
    qv_rev_new   = fv3_remap.map1_q2(nlevs, from_p, qv_init_rev[np.newaxis, :],               \
                                     nlevs, to_p, dp2, 0, 0, 0, kord_tr, q_min)
    qv_temp[0,:] = qv_rev_new[0,::-1]

    # Temperature @ time = 0
    t_temp = tv_temp/(1.0 + zvir*qv_temp)

    # Zonal wind @ time = 0
    u_init_rev  = state_IC["ua"][::-1]
    u_rev_new   = fv3_remap.map1_ppm(nlevs, from_p, u_init_rev[np.newaxis, :], 0.0,           \
                                     nlevs, to_p, 0, 0, -1, kord_tm )
    u_temp[0,:] = u_rev_new[0,::-1]

    # Meridional wind @ time = 0
    v_init_rev  = state_IC["va"][::-1]
    v_rev_new   = fv3_remap.map1_ppm(nlevs, from_p, v_init_rev[np.newaxis, :], 0.0,           \
                                     nlevs, to_p, 0, 0, -1, kord_tm )
    v_temp[0,:] = v_rev_new[0,::-1]

    # Pressure advection @ time = 0
    valid_pres_adv[0,:]   = state_IC["pa"]
    valid_pres_i_adv[0,:] = state_IC["pa_i"]

    # Moisture advection @ time = 0
    dqvdt_temp[0,:] = (qv_temp[0,:] - state_IC["qv"][:])/(3600.0*(time_dyn_hours[0]))
    dqvdt_adv[0,:]  = vars_dyn[iqv]["values"][0,:] - dqvdt_temp[0,:]

    # Temperature advection @ time = 0
    dtdt_temp[0,:]  = (t_temp[0,:] - state_IC["ta"][:])/(3600.0*(time_dyn_hours[0]))
    dtdt_adv[0,:]   = vars_dyn[itemp]["values"][0,:] - dtdt_temp[0,:]

    # u-momentum advection @ time = 0
    dudt_temp[0,:]  = (u_temp[0,:] - state_IC["ua"][:])/(3600.0*(time_dyn_hours[0]))
    dudt_adv[0,:]   = vars_dyn[iu]["values"][0,:] - dudt_temp[0,:]

    # v-momentum advection @ time = 0
    dvdt_temp[0,:]  = (v_temp[0,:] - state_IC["va"][:])/(3600.0*(time_dyn_hours[0]))
    dvdt_adv[0,:]   = vars_dyn[iv]["values"][0,:] - dvdt_temp[0,:]

    # @ time > 0
    for t in range(n_files-1):
        #
        from_p[0,:]     = p_lev[t,::-1]
        to_p[0,:]       = p_lev[t+1,::-1]
        log_from_p[0,:] = np.log(p_lev[t,::-1])
        log_to_p[0,:]   = np.log(p_lev[t+1,::-1])

        # Virtual Temperature @ time > 0
        tv_rev      = np.zeros([1,nlevs])
        tv_rev[0,:] = tv_lay[t,::-1]
        tv_rev_new  = fv3_remap.map_scalar(nlevs, log_from_p, tv_rev, dummy, nlevs,     \
                                           log_to_p, 0, 0, 1, np.abs(kord_tm), t_min)

        # Specific humidity @ time > 0
        qv_rev[0,:] = qv_lay[t,::-1]
        for k in range(0,nlevs): dp2[0,k] = to_p[0,k+1] - to_p[0,k]
        qv_rev_new = fv3_remap.map1_q2(nlevs, from_p, qv_rev, nlevs, to_p, dp2,         \
                                       0, 0, 0, kord_tr, q_min)

        # Zonal wind  @ time > 0
        u_rev[0,:] = u_lay[t,::-1]
        u_rev_new  = fv3_remap.map1_ppm(nlevs, from_p, u_rev, 0.0, nlevs, to_p,         \
                                        0, 0, -1, kord_tm )

        # Meridional wind @ time > 0
        v_rev[0,:] = v_lay[t,::-1]
        v_rev_new  = fv3_remap.map1_ppm(nlevs, from_p, v_rev, 0.0, nlevs, to_p,         \
                                        0, 0, -1, kord_tm )

        # Store
        tv_layers_remap[t+1,:] = tv_rev_new[0,::-1]
        qv_layers_remap[t+1,:] = qv_rev_new[0,::-1]
        u_layers_remap[t+1,:]  = u_rev_new[0,::-1]
        v_layers_remap[t+1,:]  = v_rev_new[0,::-1]

    # Temperature at all times (using tv/qv at all times)
    t_layers_remap = tv_layers_remap/(1.0 + zvir*qv_layers_remap)

    # @ time > 0
    for t in range(n_files-1):
        dtime = 3600.0*(time_dyn_hours[t+1] - time_dyn_hours[t])

        # Pressure advection @ time > 0
        valid_pres_adv[t+1,:]   = p_lay[t]
        valid_pres_i_adv[t+1,:] = p_lev[t]
        
        # Moisture advection @ time > 0
        dqvdt_temp[0,:] = (qv_layers_remap[t+1,:] - qv_lay[t,:])/dtime
        dqvdt_adv[t+1,:] = vars_dyn[iqv]["values"][t+1,:] - dqvdt_temp[0,:]
        
        # Temperature advection @ time > 0
        dtdt_temp[0,:] = (t_layers_remap[t+1,:] - t_lay[t,:])/dtime
        dtdt_adv[t+1,:] = vars_dyn[itemp]["values"][t+1,:] - dtdt_temp[0,:]

        # u-momentum advection @ time > 0
        dudt_temp[0,:] = (u_layers_remap[t+1,:] - u_lay[t,:])/dtime
        dudt_adv[t+1,:] = vars_dyn[iu]["values"][t+1,:] - dudt_temp[0,:]
        
        # v-momentum advection @ time > 0
        dvdt_temp[0,:] = (v_layers_remap[t+1,:] - v_lay[t,:])/dtime
        dvdt_adv[t+1,:] = vars_dyn[iv]["values"][t+1,:] - dvdt_temp[0,:]

    ####################################################################################
    # for the original SCM forcing input file, all forcing terms should be valid on the
    # initial pressure levels; one should be able to use fv3_remap.map_scalar to remap
    # these forcing terms to the initial pressure profile, rather than resort to linear
    # interpolation or something else. (this interpolation can be removed when using
    # DEPHY, that allows for varying pressure levels for forcing)
    ####################################################################################
    if orig_SCM_format:
        dqvdt_adv_at_init_pres_rev = np.zeros([n_files,1,nlevs])
        dtdt_adv_at_init_pres_rev  = np.zeros([n_files,1,nlevs])
        dudt_adv_at_init_pres_rev  = np.zeros([n_files,1,nlevs])
        dvdt_adv_at_init_pres_rev  = np.zeros([n_files,1,nlevs])
        from_log_p = np.zeros([1,nlevs+1])
        to_log_p   = np.zeros([1,nlevs+1])

        # @ time = 0
        # don't need to remap the first time interval because it is already valid at the
        # initial pressure levels.
        dqvdt_adv_at_init_pres_rev[0,0,:] = dqvdt_adv[0,::-1]
        dtdt_adv_at_init_pres_rev[0,0,:]  = dtdt_adv[0,::-1]
        dudt_adv_at_init_pres_rev[0,0,:]  = dudt_adv[0,::-1]
        dvdt_adv_at_init_pres_rev[0,0,:]  = dvdt_adv[0,::-1]

        # @ time > 0
        for t in range(1,n_files):
            from_log_p[0,:] = np.log(valid_pres_i_adv[t,::-1])
            to_log_p[0,:] = np.log(valid_pres_i_adv[0,::-1])
            dqvdt_adv_at_init_pres_rev[t,:,:] = fv3_remap.map_scalar(nlevs, from_log_p, \
                dqvdt_adv[t,np.newaxis,::-1], dummy, nlevs, to_log_p, 0, 0, 1, 1, 0.0)
            dtdt_adv_at_init_pres_rev[t,:,:] = fv3_remap.map_scalar( nlevs, from_log_p, \
                dtdt_adv[t,np.newaxis,::-1], dummy, nlevs, to_log_p, 0, 0, 1, 1, 0.0)
            dudt_adv_at_init_pres_rev[t,:,:] = fv3_remap.map_scalar( nlevs, from_log_p, \
                dudt_adv[t,np.newaxis,::-1], dummy, nlevs, to_log_p, 0, 0, 1, 1, 0.0)
            dvdt_adv_at_init_pres_rev[t,:,:] = fv3_remap.map_scalar( nlevs, from_log_p, \
                dvdt_adv[t,np.newaxis,::-1], dummy, nlevs, to_log_p, 0, 0, 1, 1, 0.0)
        # Store
        dtdt_adv_at_init_pres  = dtdt_adv_at_init_pres_rev[:,0,::-1]
        dqvdt_adv_at_init_pres = dqvdt_adv_at_init_pres_rev[:,0,::-1]
        dudt_adv_at_init_pres  = dudt_adv_at_init_pres_rev[:,0,::-1]
        dvdt_adv_at_init_pres  = dvdt_adv_at_init_pres_rev[:,0,::-1]

    if save_comp_data:
        #
        # State variables
        #
        from_log_p = np.zeros([1,nlevs+1])
        to_log_p   = np.zeros([1,nlevs+1])
        from_p     = np.zeros([1,nlevs+1])
        to_p       = np.zeros([1,nlevs+1])
        t_layr     = np.zeros([n_files,1,nlevs])
        qv_layr    = np.zeros([n_files,1,nlevs])
        u_layr     = np.zeros([n_files,1,nlevs])
        v_layr     = np.zeros([n_files,1,nlevs])

        # @ time = 0
        t_layr[0,0,:]  = t_lay[0,::-1]
        qv_layr[0,0,:] = qv_lay[0,::-1]
        u_layr[0,0,:]  = u_lay[0,::-1]
        v_layr[0,0,:]  = v_lay[0,::-1]

        # @ time > 0
        for t in range(1,n_files):
            from_log_p[0,:] = np.log(p_lev[t,::-1])
            to_log_p[0,:]   = np.log(p_lev[0,::-1])
            from_p[0,:]     = p_lev[t,::-1]
            to_p[0,:]       = p_lev[0,::-1]
            for k in range(0,nlevs): dp2[0,k]    = to_p[0,k+1] - to_p[0,k]
            t_layr[t,:,:]   = fv3_remap.map_scalar(nlevs, from_log_p, t_lay[t,np.newaxis,::-1],  \
                                                dummy, nlevs, to_log_p, 0, 0, 1, np.abs(kord_tm), t_min)
            qv_layr[t,:,:]  = fv3_remap.map1_q2(nlevs,    from_p,     qv_lay[t,np.newaxis,::-1], \
                                                nlevs, to_p, dp2, 0, 0, 0, kord_tr, q_min)
            u_layr[t,:,:]   = fv3_remap.map1_ppm(nlevs,   from_p,     u_lay[t,np.newaxis,::-1],  \
                                                 0.0, nlevs, to_p, 0, 0, -1, kord_tm)
            v_layr[t,:,:]   = fv3_remap.map1_ppm(nlevs,   from_p,     v_lay[t,np.newaxis,::-1],  \
                                                 0.0, nlevs, to_p, 0, 0, -1, kord_tm)
        #
        # Physics tendencies
        #
        to_tend   = np.zeros([n_files,1,nlevs])
        from_tend = np.zeros([n_files,1,nlevs]) 
        for diag in vars_comp:
            try:
                # @ time = 0
                to_tend[0,0,:] = diag["values"][0,::-1]
                # @ time > 0
                for itime in range(1,n_files):
                    from_log_p[0,:]  = np.log(p_lev[itime,::-1])
                    to_log_p[0,:]    = np.log(p_lev[0,::-1])
                    from_p[0,:]      = p_lev[itime,::-1]
                    to_p[0,:]        = p_lev[0,::-1]
                    from_tend[0,0,:] = diag["values"][itime,np.newaxis,::-1]
                    to_tend[t,:,:]   = fv3_remap.map_scalar(nlevs, from_log_p, from_tend,\
                                            dummy, nlevs, to_log_p, 0, 0, 1, 1, -9999.99)
                # Store
                diag["values_regrid"] = to_tend[:,0,::-1]
            except:
                logging.debug(dtend["name"] + ' not available for interpolation')

    ####################################################################################
    # if we had atmf,sfcf files at every timestep (and the SCM timestep is made to match
    # the UFS), then dqvdt_adv should be applied uninterpolated for each time step. If
    # atmf and sfcf files represent time averages over the previous diagnostic period,
    # and if forcing terms are interpolatd in time in the SCM, then dqvdt_adv should
    # represent the forcing values in the middle of time[t] and time[t+1] from atmf/sfcf.
    # That way, the time-averaged applied forcing from time[t] to time[t+1] in the SCM
    # will be equal to what is derived from atmf/sfcf. (preference should be to have
    # option to remove time-interpolation of forcing such that the constant forcing
    # applied converged to time-step values as the diag interval approaches the time
    # step)
    ####################################################################################
    #time_method = 'constant_simple' #this is not implemented in the SCM code yet
    time_method = 'constant_interp'
    #time_method = 'gradient' #this produced wonky results in the SCM; avoid until investigated more
    
    if (time_method == 'constant_simple'):
        print('Forcing should not be interpolated in time. Rather, forcing should held constant at their current values until the next forcing interval is reached.')
        ntimes = n_files
        time = np.zeros(ntimes)
        
        if orig_SCM_format:
            h_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_v = np.zeros((nlevs,ntimes),dtype=float)
        else:
            p_s = np.zeros((ntimes),dtype=float)
            pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)        
        
        if orig_SCM_format:
            h_advec_qt[:,0] = dqvdt_adv_at_init_pres[0,:]
            for k in range(nlevs):
                h_advec_thil[k,0] = (p0/state_IC["pa"][k])**kappa*dtdt_adv_at_init_pres[0,k]
            h_advec_u[:,0] = dudt_adv_at_init_pres[0,:]
            h_advec_v[:,0] = dvdt_adv_at_init_pres[0,:]
        else:
            p_s[0] = ps[0]
            pressure_forc[:,0] = valid_pres_adv[0,:]
            tot_advec_T[:,0] = dtdt_adv[0,:]
            tot_advec_qv[:,0] = dqvdt_adv[0,:]
            tot_advec_u[:,0] = dudt_adv[0,:]
            tot_advec_v[:,0] = dvdt_adv[0,:]
        
        for t in range(1,n_files):
            time[t] = 3600.0*time_dyn_hours[t-1]
            
            if orig_SCM_format:
                h_advec_qt[:,t] = dqvdt_adv_at_init_pres[t,:]
                for k in range(nlevs):
                    h_advec_thil[k,t] = (p0/state_IC["pa"][k])**kappa*dtdt_adv_at_init_pres[t,k]
                h_advec_u[:,t] = dudt_adv_at_init_pres[t,:]
                h_advec_v[:,t] = dvdt_adv_at_init_pres[t,:]
            else:
                p_s[t] = ps[t]
                pressure_forc[:,t] = valid_pres_adv[t,:]
                tot_advec_T[:,t] = dtdt_adv[t,:]
                tot_advec_qv[:,t] = dqvdt_adv[t,:]
                tot_advec_u[:,t] = dudt_adv[t,:]
                tot_advec_v[:,t] = dvdt_adv[t,:]
    elif (time_method == 'constant_interp'):
        print('Forcing can be interpolated in time, but the time values are chosen such that forcing will effectively be held consant during a diagnostic time interval.')
        ntimes = 2*n_files
        
        time_setback = 1.0 #s
        
        time = np.zeros(ntimes)
        
        if orig_SCM_format:
            h_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_v = np.zeros((nlevs,ntimes),dtype=float)
        else:
            p_s = np.zeros((ntimes),dtype=float)
            pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)
            
        time[0] = 0.0
        time[1] = 3600.0*time_dyn_hours[0] - time_setback #forcing period should extend from beginning of diagnostic period to right BEFORE the next one
        
        if orig_SCM_format:
            h_advec_qt[:,0] = dqvdt_adv_at_init_pres[0,:]
            h_advec_qt[:,1] = h_advec_qt[:,0]
            for k in range(nlevs):
                h_advec_thil[k,0] = (p0/state_IC["pa"][k])**kappa*dtdt_adv_at_init_pres[0,k]
            h_advec_thil[:,1] = h_advec_thil[:,0]
            h_advec_u[:,0] = dudt_adv_at_init_pres[0,:]
            h_advec_u[:,1] = h_advec_u[:,0]
            h_advec_v[:,0] = dvdt_adv_at_init_pres[0,:]
            h_advec_v[:,1] = h_advec_v[:,0]
        else:
            p_s[0] = ps[0]
            p_s[1] = p_s[0]
            pressure_forc[:,0] = valid_pres_adv[0,:]
            pressure_forc[:,1] = pressure_forc[:,0]
            tot_advec_T[:,0] = dtdt_adv[0,:]
            tot_advec_T[:,1] = tot_advec_T[:,0]
            tot_advec_qv[:,0] = dqvdt_adv[0,:]
            tot_advec_qv[:,1] = tot_advec_qv[:,0]
            tot_advec_u[:,0] = dudt_adv[0,:]
            tot_advec_u[:,1] = tot_advec_u[:,0]
            tot_advec_v[:,0] = dvdt_adv[0,:]
            tot_advec_v[:,1] = tot_advec_v[:,0]
        
        for t in range(1,n_files):
            time[2*t] = 3600.0*time_dyn_hours[t-1]
            time[2*t+1] = 3600*time_dyn_hours[t] - time_setback
            
            if orig_SCM_format:
                h_advec_qt[:,2*t] = dqvdt_adv_at_init_pres[t,:]
                h_advec_qt[:,2*t+1] = h_advec_qt[:,2*t]
                for k in range(nlevs):
                    h_advec_thil[k,2*t] = (p0/state_IC["pa"][k])**kappa*dtdt_adv_at_init_pres[t,k]
                h_advec_thil[:,2*t+1] = h_advec_thil[:,2*t]
                h_advec_u[:,2*t] = dudt_adv_at_init_pres[t,:]
                h_advec_u[:,2*t+1] = h_advec_u[:,2*t]
                h_advec_v[:,2*t] = dvdt_adv_at_init_pres[t,:]
                h_advec_v[:,2*t+1] = h_advec_v[:,2*t]
            else:
                p_s[2*t] = ps[t]
                p_s[2*t+1] = p_s[2*t]
                pressure_forc[:,2*t] = valid_pres_adv[t,:]
                pressure_forc[:,2*t+1] = pressure_forc[:,2*t]
                tot_advec_T[:,2*t] = dtdt_adv[t,:]
                tot_advec_T[:,2*t+1] = tot_advec_T[:,2*t]
                tot_advec_qv[:,2*t] = dqvdt_adv[t,:]
                tot_advec_qv[:,2*t+1] = tot_advec_qv[:,2*t]
                tot_advec_u[:,2*t] = dudt_adv[t,:]
                tot_advec_u[:,2*t+1] = tot_advec_u[:,2*t]
                tot_advec_v[:,2*t] = dvdt_adv[t,:]
                tot_advec_v[:,2*t+1] = tot_advec_v[:,2*t]
        
    elif (time_method == 'gradient'): #this produced wonky results in the SCM; avoid until investigated more
        print('Forcing can be interpolated in time since the forcing terms are assumed to follow a constant time-gradient.')
        
        ntimes = 2*n_files + 1
        time = np.zeros(ntimes)
        
        if orig_SCM_format:
            h_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            h_advec_v = np.zeros((nlevs,ntimes),dtype=float)
        else:
            p_s = np.zeros((ntimes),dtype=float)
            pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
            tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)
        
        if orig_SCM_format:
            h_advec_qt[:,0] = 0.0
            h_advec_thil[:,0] = 0.0
            h_advec_u[:,0] = 0.0
            h_advec_v[:,0] = 0.0
        else:
            p_s[0] = state_IC['ps']
            pressure_forc[:,0] = state_IC['pa']
            tot_advec_T[:,0] = 0.0
            tot_advec_qv[:,0] = 0.0
            tot_advec_u[:,0] = 0.0
            tot_advec_v[:,0] = 0.0
        
        for t in range(n_files):
            time[2*t + 1] = time[2*t] + 0.5*(3600*time_dyn_hours[t] - time[2*t])
            time[2*t + 2] = 3600.0*time_dyn_hours[t]
            
            if orig_SCM_format:
                h_advec_qt[:,2*t + 1] = dqvdt_adv_at_init_pres[t,:]
                for k in range(nlevs):
                    h_advec_thil[k,2*t + 1] = (p0/state_IC["pa"][k])**kappa*dtdt_adv_at_init_pres[t,k]
                h_advec_u[:,2*t + 1] = dudt_adv_at_init_pres[t,:]
                h_advec_v[:,2*t + 1] = dvdt_adv_at_init_pres[t,:]
            else:
                p_s[2*t+1] = ps[t]
                pressure_forc[:,2*t+1] = valid_pres_adv[t,:]
                tot_advec_T[:,2*t+1] = dtdt_adv[t,:]
                tot_advec_qv[:,2*t+1] = dqvdt_adv[t,:]
                tot_advec_u[:,2*t+1] = dudt_adv[t,:]
                tot_advec_v[:,2*t+1] = dvdt_adv[t,:]
            
            if orig_SCM_format:
                #calculate gradient in time and extrapolate for time (2t + 2)
                for k in range(nlevs):
                    grad = (h_advec_qt[k,2*t + 1] - h_advec_qt[k, 2*t])/(time[2*t + 1] - time[2*t])
                    h_advec_qt[k,2*t + 2] = h_advec_qt[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                    
                    grad = (h_advec_thil[k,2*t + 1] - h_advec_thil[k, 2*t])/(time[2*t + 1] - time[2*t])
                    h_advec_thil[k,2*t + 2] = h_advec_thil[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                    
                    grad = (h_advec_u[k,2*t + 1] - h_advec_u[k, 2*t])/(time[2*t + 1] - time[2*t])
                    h_advec_u[k,2*t + 2] = h_advec_u[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
                    
                    grad = (h_advec_v[k,2*t + 1] - h_advec_v[k, 2*t])/(time[2*t + 1] - time[2*t])
                    h_advec_v[k,2*t + 2] = h_advec_v[k,2*t+1] + grad*(time[2*t + 2] - time[2*t + 1])
            else:
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
    v_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
    v_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
    if orig_SCM_format:
        p_s = np.zeros((ntimes),dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)
    else:
        h_advec_thil = np.zeros((nlevs,ntimes),dtype=float)
        h_advec_qt = np.zeros((nlevs,ntimes),dtype=float)
        h_advec_u = np.zeros((nlevs,ntimes),dtype=float)
        h_advec_v = np.zeros((nlevs,ntimes),dtype=float)

    forcing = {
        "time":     time,
        "wa":       w_ls,
        "wap":      omega,
        "tnta_rad": rad_heating,
        "ps_forc":  np.ones(ntimes)*ps[0],
        "pa_forc":  pressure_forc.swapaxes(0,1),
        "tnta_adv": tot_advec_T.swapaxes(0,1),
        "tnqv_adv": tot_advec_qv.swapaxes(0,1),
        "tnua_adv": tot_advec_u.swapaxes(0,1),
        "tnva_adv": tot_advec_v.swapaxes(0,1)
    }
    
    if (save_comp_data):
        comp_data = {
            "time": time_dyn_hours*3600.0,
            "pa": p_lay[0,:],
            "ta" : t_layr[:,0,::-1], #(time,nlevs)
            "qv" : qv_layr[:,0,::-1],
            "ua" : u_layr[:,0,::-1],
            "va" : v_layr[:,0,::-1],
            "vars_comp": vars_comp
        }
    else:
        comp_data = {}

    return (forcing, comp_data)
    
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
        surface["tvxy"] = surface["tsfco"]
        surface["tgxy"] = surface["tsfco"]
        surface["tahxy"] = surface["tsfco"]
        
        if (surface["snwdph"] > 0.01 and surface["tsfco"] > 273.15 ):
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
        surface["sneqvoxy"] = surface["weasd"]     # mm
        surface["alboldxy"] = 0.65
        surface["qsnowxy"]  = 0.0
        
        surface["wslakexy"] = 0.0
        surface["taussxy"]  = 0.0
        surface["waxy"]     = 4900.0
        surface["wtxy"]     = surface["waxy"]
        surface["zwtxy"]    = (25.0 + 2.0) - surface["waxy"] / 1000.0 /0.2
        
        vegtyp = int(surface['vegtyp'])
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
        
        if (surface["weasd"] != 0.0 and snd == 0.0 ):
            snd = surface["weasd"]/1000.0
            
        if (vegtyp == 15):                       # land ice in MODIS/IGBP
            if ( surface["weasd"] < 0.1):
                surface["weasd"] = 0.1
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
            surface["snicexy"][k] = 1.00 * dzsno[k] * surface["weasd"]/snd  #this line causes a warning
            
            dzsnso[k] = -dzsno[k]
        
        for k in range(n_snow_layers, n_snow_layers + n_soil_layers):
            dzsnso[k] = -dzs[k - n_snow_layers]
        
        surface["zsnsoxy"][isnow] = dzsnso[isnow]
        for k in range(isnow+1,n_snow_layers + n_soil_layers):
            surface["zsnsoxy"][k] = surface["zsnsoxy"][k-1] + dzsnso[k]
        
        soilparm = read_noahmp_soil_table()
        
        soiltyp  = int(surface["soiltyp"])
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

def write_SCM_case_file(state, surface, oro, forcing, case, date, add_UFS_dyn_tend, add_UFS_NOAH_lsm):
    """Write all data to a netCDF file in the DEPHY-SCM format"""

    # Working types
    wp = np.float64
    wi = np.int32
    
    # Dimensions
    nlevs = state["nlevs"]
    nsoil = len(surface["stc"])
    nsnow = len(surface["snicexy"])
    nice  = len(surface["tiice"])

    # Local switches
    forcing_on  = 1
    forcing_off = 0

    # Output file
    if (add_UFS_dyn_tend): 
        fileOUT = os.path.join(PROCESSED_CASE_DIR, case + '.nc')
    else:
        fileOUT = os.path.join(PROCESSED_CASE_DIR, case + '_noforce.nc')

    nc_file = Dataset(fileOUT, 'w', format='NETCDF4')
    if (add_UFS_dyn_tend):
        nc_file.description = "FV3GFS model profile input (UFS dynamic tendencies, SCM-UFS replay mode.)"
        nc_file.modifications = 'contains dynamic forcing from UFS'
    elif (add_UFS_NOAH_lsm):
        nc_file.description   = "FV3GFS model profile input (With NOAH Land surface moodel surface forcings.)"
        nc_file.modifications = 'contains initial conditions for Noah LSM' 
    elif (add_UFS_dyn_tend and add_UFS_NOAH_lsm):
        nc_file.description = "FV3GFS model profile input (UFS dynamic tendencies and NOAH Land surface moodel surface forcings)"
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
    runtime           = timedelta(seconds=forcing['time'][-1])
    end_date          = start_date + runtime
    end_date_string   = end_date.strftime("%Y-%m-%d %H:%M:%S")
    start_date_string = start_date.strftime("%Y-%m-%d %H:%M:%S")
    #
    loc_string  = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string
    #
    nc_file.case           = case_string
    nc_file.title          = 'Forcing and Initial Conditions for ' + case_string
    nc_file.reference      = ''
    nc_file.author         = 'Grant J. Firl and Dustin Swales'
    nc_file.version        = 'Created on ' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    nc_file.format_version = '1.0'
    nc_file.modifications  = ''
    nc_file.script         = os.path.basename(__file__)
    nc_file.comment        = ''
    nc_file.start_date     = start_date_string
    nc_file.end_date       = end_date_string
    nc_file.radiation      = "off"
    nc_file.adv_ta            = forcing_off
    nc_file.adv_qv            = forcing_off
    nc_file.adv_ua            = forcing_off
    nc_file.adv_va            = forcing_off
    nc_file.adv_theta         = forcing_off
    nc_file.adv_thetal        = forcing_off
    nc_file.adv_qt            = forcing_off
    nc_file.adv_rv            = forcing_off
    nc_file.adv_rt            = forcing_off
    nc_file.forc_wa           = forcing_off
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
    if (add_UFS_dyn_tend):
        nc_file.adv_ta = forcing_on
        nc_file.adv_qv = forcing_on
        nc_file.adv_ua = forcing_on
        nc_file.adv_va = forcing_on
    #
    nc_file.surface_forcing_temp     = 'none'
    nc_file.surface_forcing_moisture = 'none'
    nc_file.surface_forcing_wind     = 'none'
    nc_file.surface_forcing_lsm      = 'noah' #'noahmp' #'ruc'

    # Set file dimension
    time_dim   = nc_file.createDimension('time', len(forcing['time']))
    timei_dim  = nc_file.createDimension('t0',    1)
    lev_dim    = nc_file.createDimension('lev',   nlevs)
    soil_dim   = nc_file.createDimension('nsoil', nsoil)
    snow_dim   = nc_file.createDimension('nsnow', nsnow)
    nslsnw_dim = nc_file.createDimension('nsoil_plus_nsnow',nsnow + nsoil)
    ice_dim    = nc_file.createDimension('nice',  nice)
    
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
    soil_depth_var               = nc_file.createVariable('soil_depth', wp, ('nsoil'))
    soil_depth_var.units         = 'm'
    soil_depth_var.standard_name = 'depth of bottom of soil layers'
    soil_depth_var[:]            = [0.1,0.4,1.0,2.0]
    #
    z0_var                       = nc_file.createVariable('z0', wp ('time'))
    z0_var.units                 =  "m"
    z0_var.standard_name         = 'surface_roughness_length_for_momentum_in_air'
    z0_var[:]                    = surface["z0"]
    #
    theta_oro                    = nc_file.createVariable('theta_oro',wp, ('t0'))
    theta_oro.units              = "deg"
    theta_oro.standard_name      = "angle with respect to east of maximum subgrid orographic variations"
    theta_oro[:]                 = oro["theta"]

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
    var_dict = [{"name": "lat",          "type":wp, "dimd": ('t0'         ),    "units": "degrees_north", "desc": "latitude"},\
                {"name": "lon",          "type":wp, "dimd": ('t0'         ),    "units": "degrees_east",  "desc": "longitude"},\
                {"name": "orog",         "type":wp, "dimd": ('t0'         ),    "units": "m",             "desc": "surface_altitude"},\
                {"name": "zh",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m",             "desc": "height"},\
                {"name": "pa",           "type":wp, "dimd": ('t0',   'lev'),    "units": "Pa",            "desc": "air_ressure"}, \
                {"name": "zh_forc",      "type":wp, "dimd": ('time', 'lev'),    "units": "m",             "desc": "height_forcing","default_value": 1.},\
                {"name": "pa_forc",      "type":wp, "dimd": ('time', 'lev'),    "units": "Pa",            "desc": "air_pressure_forcing"}, \
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
                {"name": "hfss",         "type":wp, "dimd": ('time',      ),    "units": "W m-2",         "desc": "surface_upward_sensible_heat_flux"}, \
                {"name": "hfls",         "type":wp, "dimd": ('time',      ),    "units": "W m-2",         "desc": "surface_upward_latent_heat_flux"}, \
                {"name": "wpthetap_s",   "type":wp, "dimd": ('time',      ),    "units": "K m s-1",       "desc": "surface_upward_potential_temperature_flux"}, \
                {"name": "wpqvp_s",      "type":wp, "dimd": ('time',      ),    "units": "m s-1",         "desc": "surface_upward_specific_humidity_flux"}, \
                {"name": "wpqtp_s",      "type":wp, "dimd": ('time',      ),    "units": "m s-1",         "desc": "surface_upward_water_mass_fraction_flux"}, \
                {"name": "wprvp_s",      "type":wp, "dimd": ('time',      ),    "units": "m s-1",         "desc": "surface_upward_humidity_mixing_ratio_flux"}, \
                {"name": "wprtp_s",      "type":wp, "dimd": ('time',      ),    "units": "m s-1",         "desc": "surface_upward_water_mixing_ratio_flux"}, \
                {"name": "ts",           "type":wp, "dimd": ('t0',        ),    "units": "K",             "desc": "surface_temperature"},\
                {"name": "ts_forc",      "type":wp, "dimd": ('time',      ),    "units": "K",             "desc": "forcing_surface_temperature"},\
                {"name": "tskin",        "type":wp, "dimd": ('t0',        ),    "units": "K",             "desc": "surface_skin_pressure"}, \
                {"name": "ps",           "type":wp, "dimd": ('t0',        ),    "units": "Pa",            "desc": "surface_air_pressure"}, \
                {"name": "ps_forc",      "type":wp, "dimd": ('time',      ),    "units": "Pa",            "desc": "forcing_surface_air_pressure"},\
                {"name": "ustar",        "type":wp, "dimd": ('time',      ),    "units": "m s-1",         "desc": "surface_friction_velocity"}, \
                #{"name": "z0",           "type":wp, "dimd": ('time',      ),    "units": "m",             "desc": "surface_roughness_length_for_momentum_in_air"}, \
                {"name": "z0h",          "type":wp, "dimd": ('time',      ),    "units": "m",             "desc": "surface_roughness_length_for_heat_in_air"}, \
                {"name": "z0q",          "type":wp, "dimd": ('time',      ),    "units": "m",             "desc": "surface_roughness_length_for_humidity_in_air"}, \
                {"name": "beta",         "type":wp, "dimd": ('t0',        ),    "units": "m",             "desc": "soil_water_stress_factor"}, \
                {"name": "mrsos",        "type":wp, "dimd": ('t0',        ),    "units": "kg m-2",        "desc": "mass_content_of_water_in_soil_layer"}, \
                {"name": "mrsos_forc",   "type":wp, "dimd": ('time',      ),    "units": "kg m-2",        "desc": "forcing_mass_content_of_water_in_soil_layer"}, \
                {"name": "o3",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mole_fraction_of_ozone_in_air"}, \
                {"name": "sza",          "type":wp, "dimd": ('t0',        ),    "units": "degree",        "desc": "solar_zenith_angle"}, \
                {"name": "io",           "type":wp, "dimd": ('t0',        ),    "units": "W m-2",         "desc": "solar_irradiance"}, \
                {"name": "alb",          "type":wp, "dimd": ('t0',        ),    "units": "1",             "desc": "surface_albedo"}, \
                {"name": "emis",         "type":wp, "dimd": ('t0',        ),    "units": "1",             "desc": "surface_longwave_emissivity"}]

    var_oro  = [{"name": "area",         "type":wp, "dimd": ('t0'),             "units": "m 2-1",   "desc": "grid_cell_area"},\
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
    var_lsm =  [{"name": "tiice",        "type":wp, "dimd": ('t0','nice'),      "units": "K",       "desc": "sea ice internal temperature"}]
    var_noah = [{"name": "vegsrc",       "type":wi,  "dimd": ('t0'),             "units": "none",    "desc": "vegetation source (1-2)", "default_value": 1}, \
                {"name": "slmsk",        "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "land-sea-ice mask"}, \
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
                {"name": "vegtyp",       "type":wi,  "dimd": ('t0'),             "units": "none",    "desc": "vegetation type (1-12)"}, \
                {"name": "soiltyp",      "type":wi,  "dimd": ('t0'),             "units": "none",    "desc": "soil type (1-12)"}, \
                {"name": "ffmm",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "Monin-Obukhov similarity function for momentum"}, \
                {"name": "ffhh",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "Monin-Obukhov similarity function for heat"}, \
                {"name": "hice",         "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "sea ice thickness"}, \
                {"name": "fice",         "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "ice fraction"}, \
                {"name": "tisfc",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "ice surface temperature"}, \
                {"name": "tprcp",        "type":wp, "dimd": ('t0'),             "units": "m",       "desc": "instantaneous total precipitation amount"}, \
                {"name": "srflag",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "snow/rain flag for precipitation"}, \
                {"name": "snwdph",       "type":wp, "dimd": ('t0'),             "units": "mm",      "desc": "water equivalent snow depth"}, \
                {"name": "shdmin",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "minimum vegetation fraction"}, \
                {"name": "shdmax",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "maximum vegetation fraction"}, \
                {"name": "slopetyp",     "type":wi,  "dimd": ('t0'),             "units": "none",    "desc": "slope type (1-9)"}, \
                {"name": "snoalb",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "maximum snow albedo"}, \
                {"name": "sncovr",       "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "surface snow area fraction"}, \
                {"name": "tsfcl",        "type":wp, "dimd": ('t0'),             "units": "K",       "desc": "surface skin temperature over land"}, \
                {"name": "stc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "K",       "desc": "initial profile of soil liquid moisture"}, \
                {"name": "smc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "kg",      "desc": "initial profile of soil moisture"}, \
                {"name": "slc",          "type":wp, "dimd": ('t0','nsoil'),     "units": "kg",      "desc": "initial profile of soil temperature"}]
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
    var_dict.extend(var_oro)
    var_dict.extend(var_nsst)
    var_dict.extend(var_lsm)
    var_dict.extend(var_ruc)
    var_dict.extend(var_noah)
    var_dict.extend(var_noahmp)

    #
    # Write all fields in "var_dict" to SCM input file.
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
    #
    # Close file
    #
    nc_file.close()

    return

def write_comparison_file(comp_data, case_name, date, surface):
    """Write UFS history file data to netCDF file for comparison"""

    wp = np.float64
    wi = np.int32

    nlevs = comp_data["pa"].shape[0]
    ntime = comp_data["time"].shape[0]

    start_date = datetime(date["year"],date["month"],date["day"],date["hour"],date["minute"],date["second"])
    start_date_string = start_date.strftime("%Y%m%d%H%M%S")

    loc_string = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string

    nc_file = Dataset(os.path.join(COMPARISON_DATA_DIR, case_name + '_comp_data.nc'), 'w', format='NETCDF3_CLASSIC')
    nc_file.case = case_string
    nc_file.title = 'UFS history file data for ' + case_string
    nc_file.reference = ''
    nc_file.author = 'Grant J. Firl'
    nc_file.version = 'Created on ' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    nc_file.script = os.path.basename(__file__)
    nc_file.startDate = start_date_string

    lev_dim = nc_file.createDimension('lev', size=nlevs)
    lev_var = nc_file.createVariable('lev', wp, ('lev',))
    lev_var.units = 'Pa'
    lev_var.long_name = 'pressure'
    lev_var[:] = comp_data["pa"]

    time_dim = nc_file.createDimension('time', size=ntime)
    time_var = nc_file.createVariable('time', wp, ('time',))
    time_var.units = 'second'# since ' + str(start_date)
    time_var.long_name = 'history file time'
    time_var[:] = comp_data['time']

    init_year_var = nc_file.createVariable('init_year', wi)
    init_year_var.units = 'year'
    init_year_var.description = 'UFS initialization year'
    init_year_var[:] = date["year"]

    init_month_var = nc_file.createVariable('init_month', wi)
    init_month_var.units = 'month'
    init_month_var.description = 'UFS initialization month'
    init_month_var[:] = date["month"]

    init_day_var = nc_file.createVariable('init_day', wi)
    init_day_var.units = 'day'
    init_day_var.description = 'UFS initialization day'
    init_day_var[:] = date["day"]

    init_hour_var = nc_file.createVariable('init_hour', wi)
    init_hour_var.units = 'hour'
    init_hour_var.description = 'UFS initialization hour'
    init_hour_var[:] = date["hour"]

    init_minute_var = nc_file.createVariable('init_minute', wi)
    init_minute_var.units = 'minute'
    init_minute_var.description = 'UFS initialization minute'
    init_minute_var[:] = date["minute"]

    init_second_var = nc_file.createVariable('init_second', wi)
    init_second_var.units = 'second'
    init_second_var.description = 'UFS initialization second'
    init_second_var[:] = date["second"]

    temperature_var = nc_file.createVariable('temp', wp, ('time', 'lev',))
    temperature_var.units = 'K'
    temperature_var.long_name = 'Temperature'
    temperature_var[:] = comp_data["ta"]

    qv_var = nc_file.createVariable('qv', wp, ('time', 'lev',))
    qv_var.units = 'kg kg-1'
    qv_var.long_name = 'specific humidity'
    qv_var[:] = comp_data["qv"]

    u_var = nc_file.createVariable('u', wp, ('time', 'lev',))
    u_var.units = 'm s-1'
    u_var.long_name = 'zonal wind'
    u_var[:] = comp_data["ua"]

    v_var = nc_file.createVariable('v', wp, ('time', 'lev',))
    v_var.units = 'm s-1'
    v_var.long_name = 'meridional wind'
    v_var[:] = comp_data["va"]

    for dtend in comp_data["vars_comp"]:
        try:
            tempVar           = nc_file.createVariable(dtend["name"], wp, ('time', 'lev',))
            tempVar.units     = dtend["units"]
            tempVar.long_name = dtend["long_name"]
            tempVar[:]        = dtend["values"]
        except:
            logging.debug(dtend["name"] + ' not available for output')
    nc_file.close()

    return

def find_date(forcing_dir):
    
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
    
    #read in arguments
    (location, indices, date, in_dir, grid_dir, forcing_dir, tile, area, noahmp, \
     case_name, old_chgres, lam, no_force, save_comp, add_UFS_NOAH_lsm,          \
     add_UFS_dyn_tend) = parse_arguments()
        
    #find tile containing the point using the supergrid if no tile is specified 
    #if not tile and not lam:
    if not tile:
        tile = int(find_tile(location, grid_dir, lam))
        if tile < 0:
            message = 'No tile was found for location {0}'.format(location)
            logging.critical(message)
            raise Exception(message)
        print('Tile found: {0}'.format(tile))
    
    #find index of closest point in the tile if indices are not specified
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
    
    # get UFS IC data (TODO: flag to read in RESTART data rather than IC data and implement
    # different file reads)
    (state_data, surface_data, oro_data) = get_UFS_IC_data(in_dir, grid_dir, forcing_dir, \
                                                    tile, tile_i, tile_j, old_chgres, lam)
    
    if not date:
        # date was not included on command line; look in atmf* file for initial date
        date = find_date(forcing_dir)
    
    # Cold start NoahMP variables (shouldn't be necessary to use this anymore, since same 
    # capability exists in SCM code if given Noah ICs only)
    if (noahmp):
        surface_data = add_noahmp_coldstart(surface_data, date)
    
    #get grid cell area if not given
    if not area:
        area = get_UFS_grid_area(grid_dir, tile, tile_i, tile_j, lam)
    surface_data["area"] = area
    
    surface_data["lon"] = point_lon
    surface_data["lat"] = point_lat
    
    # Get UFS forcing data
    (forcing_data, comp_data) = get_UFS_forcing_data(state_data["nlevs"], state_data,      \
                                                     forcing_dir, grid_dir, tile, tile_i,  \
                                                     tile_j, lam, save_comp)
    
    # Write SCM case file
    write_SCM_case_file(state_data, surface_data, oro_data, forcing_data, case_name, date, \
                        add_UFS_dyn_tend, add_UFS_NOAH_lsm)

    # Create file w/o forcings for testing (only if forcings are being used)
    if (add_UFS_dyn_tend):
        write_SCM_case_file(state_data, surface_data, oro_data, forcing_data, case_name,   \
                            date, False, add_UFS_NOAH_lsm)

    # read in and remap the state variables to the first history file pressure profile and 
    # write them out to compare SCM output to (atmf for state variables and sfcf for physics 
    # tendencies)
    if (save_comp):
        write_comparison_file(comp_data, case_name, date, surface_data)
    
if __name__ == '__main__':
    main()
