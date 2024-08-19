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
Omega        = 7.292E-5
rdgas        = 287.05
rvgas        = 461.50
cp           = 1004.6
zvir         = rvgas/rdgas - 1.
rocp         = rdgas/cp
grav         = 9.80665
deg_to_rad   = math.pi/180.0
kappa        = rdgas/cp
p0           = 100000.0
one_twelfth  = 1.0/12.0
two_thirds   = 2.0/3.0

use_oro_for_geos_wind = False  #use surface geopotential in the calculation of geostrophic winds (if True, the gradient of geopotential will include the surface geopotential height gradient)
use_actual_geos_wind = False   #use geostrophic winds calculated from the large-scale UFS-modeled geopotential (the scale is determined by the modeled winds and the choice of critical Rossby number for which geostrophic balance is assumed to hold)
critical_Rossby = 0.1          #the Rossby number under which geostrophic balance is assumed to be valid (a lower number means a larger horizontal area is required for a given wind speed); geostrophic balance can be assumed when Ro << 1.

PBL_top_for_geos = 85000.0 # (Pa) when using the mean UFS-modeled winds as geostrophic winds, below this pressure, geostrophic winds return linearly to 0 at the surface

missing_value = -9999.0 #9.99e20

top_pres_for_forcing = 20000.0 # (Pa) pressure above which to vertically smooth wind profiles to handle noisy model output at upper UFS levels

const_nudging_time = 3600.0

n_lam_halo_points = 3

n_forcing_halo_points = 2 #number of points in each direction used to retrieve forcing (1 = 3x3, 2 = 5x5, 3 = 7x7); must be >=1

#dist_method = 'euclidean'
dist_method = 'haversine' #faster

missing_variable_snow_layers = 3
missing_variable_soil_layers = 4
missing_variable_ice_layers = 2
missing_variable_vegetation_categories = 20
missing_variable_soil_categories = 16

# Path to the directory containing processed case input files
PROCESSED_CASE_DIR = '../../data/processed_case_input'

# Path to the directory containing comparison data files
COMPARISON_DATA_DIR = '../../data/comparison_data'

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('-l',   '--location',    help='longitude and latitude in degress E and N, respectively, separated by a space', nargs=2, type=float)
group1.add_argument('-ij',  '--index',       help='i,j indices within the tile (if known - bypasses search for closest model point to lon/lat location)', nargs=2, type=int)
parser.add_argument('-d',   '--date',        help='date corresponding to initial conditions in YYYYMMDDHHMMSS format', required=False)
parser.add_argument('-i',   '--in_dir',      help='path to the directory containing the UFS initial conditions',                    required=True)
parser.add_argument('-g',   '--grid_dir',    help='path to the directory containing the UFS supergrid files (AKA "fix" directory)', required=True)
parser.add_argument('-f',   '--forcing_dir', help='path to the directory containing the UFS history files',                         required=True)
parser.add_argument('-n',   '--case_name',   help='name of case',                                                                   required=True)
parser.add_argument('-t',   '--tile',        help='tile of desired point (if known - bypasses tile search if present)', type=int, choices=range(1,8))
parser.add_argument('-a',   '--area',        help='area of grid cell in m^2', type=float)
parser.add_argument('-oc',  '--old_chgres',  help='flag to denote that the initial conditions use an older data format (pre-chgres_cube)', action='store_true')
parser.add_argument('-lam', '--lam',         help='flag to signal that the ICs and forcing is from a limited-area model run',              action='store_true')
parser.add_argument('-sc',  '--save_comp',   help='flag to create UFS reference file for comparison',                                      action='store_true')
parser.add_argument('-near','--use_nearest', help='flag to indicate using the nearest UFS history file gridpoint',                         action='store_true')
parser.add_argument('-fm',  '--forcing_method',  help='method used to calculate forcing (1=total tendencies from UFS dycore, 2=advective terms calculated from UFS history files, 3=total time tendency terms calculated)', type=int, choices=range(1,4), default=2)
parser.add_argument('-vm',  '--vertical_method', help='method used to calculate vertical advective forcing (1=vertical advective terms calculated from UFS history files and added to total, 2=smoothed vertical velocity provided)', type=int, choices=range(1,3), default=2)
parser.add_argument('-wn',  '--wind_nudge',  help='flag to turn on wind nudging to UFS profiles', action='store_true')
parser.add_argument('-geos','--geostrophic', help='flag to turn on geostrophic wind forcing', action='store_true')

########################################################################################
#
########################################################################################
def parse_arguments():
    """Parse command line arguments"""
    args        = parser.parse_args()
    location    = args.location
    index       = args.index
    date        = args.date
    in_dir      = args.in_dir
    grid_dir    = args.grid_dir
    forcing_dir = args.forcing_dir
    tile        = args.tile
    area        = args.area
    case_name   = args.case_name
    old_chgres  = args.old_chgres
    lam         = args.lam
    save_comp   = args.save_comp
    use_nearest = args.use_nearest
    forcing_method = args.forcing_method
    vertical_method = args.vertical_method
    geos_wind_forcing = args.geostrophic
    wind_nudge  = args.wind_nudge
    
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
            area, case_name, old_chgres, lam, save_comp, use_nearest, forcing_method, \
            vertical_method, geos_wind_forcing, wind_nudge)

########################################################################################
#
########################################################################################
def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)

########################################################################################
#
########################################################################################
def haversine_distance(lat1, lon1, lat2, lon2):
    #lats and lons in degrees
    p = 0.017453292519943295
    hav = 0.5 - math.cos((lat2-lat1)*p)/2 + math.cos(lat1*p)*math.cos(lat2*p) * (1-math.cos((lon2-lon1)*p)) / 2
    return 12742 * math.asin(math.sqrt(hav)) #km

########################################################################################
#
########################################################################################
def sph2cart(az, el, r):
    """Calculate the Cartesian coordiates from spherical coordinates"""
    
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    
    return (x, y, z)
    
########################################################################################
#
########################################################################################
def moving_average(interval, window_size):
    convolve_mode = 'same'
    window= np.ones(int(window_size))/float(window_size)
    
    output = np.zeros((interval.shape[0]))
    if (convolve_mode == 'same'):
        convolution = np.convolve(interval, window, mode='same')
        output = convolution
    elif (convolve_mode == 'valid'):
        convolution = np.convolve(interval, window, mode='valid')
        num_not_valid_points = interval.shape[0] - (max(interval.shape[0], window_size) - min(interval.shape[0], window_size) + 1)
        output[num_not_valid_points:] = convolution
        for k in range(num_not_valid_points):
            output[k] = interval[k]
            output[-1-k] = interval[-1-k]
        
    return output

########################################################################################
#
########################################################################################
def centered_diff_derivative(y, dx):
    #y is the variable of interest, dx is the space between contiguous points (starting from the first point)
    cent_diff_coef_three = [-0.5, 0.0, 0.5] #second-order accuracy in space (assuming constant grid)
    cent_diff_coef_five = [one_twelfth, -two_thirds, 0.0, two_thirds, -one_twelfth] #fourth-order accuracy in space (assuming constant grid)
    
    if y.shape[0] % 2:
        center = int((y.shape[0]-1)/2)
        if y.shape[0] == 3:
            deriv = np.sum(cent_diff_coef_three*y)/np.mean(dx)
        elif y.shape[0] == 5:
            deriv = np.sum(cent_diff_coef_five*y)/np.mean(dx)
        else:
            message = 'centered_diff_derivative can only use n_forcing_halo_points = 1 or 2'
            logging.exception(message)
            raise Exception(message)
    else:
        message = 'An even number of points was passed into centered_diff_derivative'
        logging.exception(message)
        raise Exception(message)
    
    return deriv

########################################################################################
#
########################################################################################
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

########################################################################################
#
########################################################################################
def find_loc_indices(loc, dir, tile, lam):
    """Find the nearest neighbor FV3 grid point given a lon/lat pair and the tile number"""
    #returns the indices of the nearest neighbor point in the given tile, the lon/lat of the nearest neighbor, 
    #and the distance (m) from the given point to the nearest neighbor grid cell
    
    if (tile >= 0):
        if (lam):
            filename_pattern = '*grid.tile7.halo{}.nc'.format(n_lam_halo_points)
        else:
            filename_pattern = '*grid.tile{0}.nc'.format(tile)
    else:
        filename_pattern = 'atmf000.nc'
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    if (tile >= 0):
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
    else:
        longitude = np.asarray(nc_file['lon'])
        latitude  = np.asarray(nc_file['lat'])

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
    
    #set up an array to hold the distance between the given point and every grid cell
    dist = np.zeros((longitude.shape[0],longitude.shape[1]))
    if (dist_method == 'euclidean'):
        
        #get the Cartesian location of the given point
        cart_loc = np.asarray(sph2cart(math.radians(loc[0]), math.radians(loc[1]), earth_radius))
        
        for i in range(len(longitude)):
            for j in range(len(longitude[i])):
                #get the Cartesian location of all grid points
                cart_cell = np.asarray(sph2cart(math.radians(longitude[i,j]), math.radians(latitude[i,j]), earth_radius))
                
                #calculate the euclidean distance from the given point to the current grid cell
                dist[i,j] = np.linalg.norm(cart_loc - cart_cell)
    else:
        for i in range(len(longitude)):
            for j in range(len(longitude[i])):
                dist[i,j] = haversine_distance(latitude[i,j],longitude[i,j],loc[1],loc[0])
    #get the indices of the grid point with the minimum euclidean distance to the given point
    i,j = np.unravel_index(dist.argmin(), dist.shape)
    
    #get the direction of the closest point from the given point
    theta = math.atan2(math.sin(math.radians(longitude[i,j] - loc[0]))*math.cos(math.radians(latitude[i,j])), math.cos(math.radians(loc[1]))*math.sin(math.radians(latitude[i,j])) - math.sin(math.radians(loc[1]))*math.cos(math.radians(latitude[i,j]))*math.cos(math.radians(longitude[i,j] - loc[0])))
    theta_deg = math.fmod((math.degrees(theta) + 360.0), 360.0)
    
    return (i,j,longitude[i,j]%360.0, latitude[i,j], dist[i,j], theta_deg)

########################################################################################
#
########################################################################################
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

########################################################################################
#
########################################################################################
def find_loc_indices_UFS_history(loc, dir, lam):
    """Find the nearest neighbor UFS history file grid point given a lon/lat pair"""
    #returns the indices of the nearest neighbor point in the given tile, the lon/lat of the nearest neighbor, 
    #and the distance (m) from the given point to the nearest neighbor grid cell

    if lam:
        filename_pattern = 'dynf000.nc'
    else: 
        filename_pattern = 'atmf000.nc'
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    longitude = np.asarray(nc_file['lon'])
    latitude  = np.asarray(nc_file['lat'])
    
    nc_file.close()
    
    #set up an array to hold the distance between the given point and every grid cell
    dist = np.zeros((longitude.shape[0],longitude.shape[1]))
    if (dist_method == 'euclidean'):
        
        #get the Cartesian location of the given point
        cart_loc = np.asarray(sph2cart(math.radians(loc[0]), math.radians(loc[1]), earth_radius))
        
        for i in range(len(longitude)):
            for j in range(len(longitude[i])):
                #get the Cartesian location of all grid points
                cart_cell = np.asarray(sph2cart(math.radians(longitude[i,j]), math.radians(latitude[i,j]), earth_radius))
                
                #calculate the euclidean distance from the given point to the current grid cell
                dist[i,j] = np.linalg.norm(cart_loc - cart_cell)
    else:
        for i in range(len(longitude)):
            for j in range(len(longitude[i])):
                dist[i,j] = haversine_distance(latitude[i,j],longitude[i,j],loc[1],loc[0])
    #get the indices of the grid point with the minimum euclidean distance to the given point
    i,j = np.unravel_index(dist.argmin(), dist.shape)
    
    #get the direction of the closest point from the given point
    theta = math.atan2(math.sin(math.radians(longitude[i,j] - loc[0]))*math.cos(math.radians(latitude[i,j])), math.cos(math.radians(loc[1]))*math.sin(math.radians(latitude[i,j])) - math.sin(math.radians(loc[1]))*math.cos(math.radians(latitude[i,j]))*math.cos(math.radians(longitude[i,j] - loc[0])))
    theta_deg = math.fmod((math.degrees(theta) + 360.0), 360.0)
    
    if i < n_forcing_halo_points or i > latitude.shape[0]-1-n_forcing_halo_points:
        message = 'Chosen point too close to the poles for this algorithm'
        logging.critical(message)
        raise Exception(message)
    else:
        #longitude = (n rows of latitude, 2*n columns of longitude), i index increases toward EAST; rows are circles of latitude, columns are meridians
        #latitude  = (n rows of latitude, 2*n columns of longitude), j index increases toward NORTH; rows are circles of latitude, columns are meridians
        neighbors = np.mgrid[-n_forcing_halo_points:n_forcing_halo_points+1,-n_forcing_halo_points:n_forcing_halo_points+1]
        neighbors[0,:,:] = neighbors[0,:,:] + i
        neighbors[1,:,:] = neighbors[1,:,:] + j
        
        
        # the closest point is near the western edge of the prime meridian, but longitude should be cyclic, so get neighbors on the eastern/western side of the prime meridian as necessary
        neighbors[1,:,:] = neighbors[1,:,:]%longitude.shape[1]
        
        # if j > n_forcing_halo_points:
        #     if j < latitude.shape[1]-1-n_forcing_halo_points:
        #         neighbors[1,:,:] = neighbors[1,:,:] + j
        #     else:
        #         
        #         n_points_on_eastern_side = (j + n_forcing_halo_points)%latitude.shape[1]
                 
                
    
    message = 'nearest history file indices (i,j): ({},{})'.format(i, j)
    message += '\n lon/lat : {},{}'.format(longitude[i,j], latitude[i,j])
    message += '\n distance between history file point and desired point: {} km'.format(dist[i,j])
    message += '\n bearing from desired point to history file point: {} degrees'.format(theta_deg)
    logging.debug(message)
    
    message = 'neighboring points'
    for ii in neighbors[0,:,0]:
        for jj in neighbors[1,0,:]:
            message += '\n (i,j,lon,lat,dist): ({:>},{:>},{:f},{:f},{:f})'.format(ii,jj,longitude[ii,jj], latitude[ii,jj], dist[ii,jj])
    logging.debug(message)
    
    #calc dx, dy (n_forcing_halo_points*2 x n_forcing_halo_points*2)
    dx = np.zeros((n_forcing_halo_points*2+1,n_forcing_halo_points*2))
    dy = np.zeros((n_forcing_halo_points*2,n_forcing_halo_points*2+1))
    for ii in range(2*n_forcing_halo_points+1):
         for jj in range(2*n_forcing_halo_points):
             current_ii_index = neighbors[0,ii,jj]
             current_jj_index = neighbors[1,ii,jj]
             eastern_jj_index = neighbors[1,ii,jj+1]
             dx[ii,jj] = haversine_distance(latitude[current_ii_index,current_jj_index],longitude[current_ii_index,current_jj_index],latitude[current_ii_index,eastern_jj_index],longitude[current_ii_index,eastern_jj_index])
        
    for ii in range(2*n_forcing_halo_points):
         for jj in range(2*n_forcing_halo_points+1):
             current_ii_index  = neighbors[0,ii,jj]
             northern_ii_index = neighbors[0,ii+1,jj]
             current_jj_index  = neighbors[1,ii,jj]
             dy[ii,jj] = haversine_distance(latitude[current_ii_index,current_jj_index],longitude[current_ii_index,current_jj_index],latitude[northern_ii_index,current_jj_index],longitude[northern_ii_index,current_jj_index])
    
    return (i,j,longitude[i,j], latitude[i,j], dist[i,j], theta_deg, neighbors, dx, dy)

########################################################################################
#
########################################################################################
def find_loc_indices_UFS_IC(loc, dir, lam, tile, indices):
    
    #find tile containing the point using the supergrid if no tile is specified 
    #if not tile and not lam:
    if not tile:
        tile = int(find_tile(loc, dir, lam))
        if tile < 0:
            message = 'No tile was found for location {0}'.format(location)
            logging.critical(message)
            raise Exception(message)
        message = 'Tile found: {0}'.format(tile)
        logging.debug(message)
    
    #find index of closest point in the tile if indices are not specified
    theta = 0.0
    if not indices:
        (tile_j, tile_i, point_lon, point_lat, dist_min, theta) = find_loc_indices(loc, dir, tile, lam)
        message = 'nearest IC file indices in tile {} - (i,j): ({},{})'.format(tile, tile_i, tile_j)
        message += '\n lon/lat : {},{}'.format(point_lon, point_lat)
        message += '\n distance between history file point and desired point: {} km'.format(dist_min)
        message += '\n bearing from desired point to history file point: {} degrees'.format(theta)
        logging.debug(message)    
    else:
        tile_i = indices[0]
        tile_j = indices[1]
        #still need to grab the lon/lat if the tile and indices are supplied
        (point_lon, point_lat) = find_lon_lat_of_indices(indices, grid_dir, tile, lam)
        
        message = 'This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat)
        logging.debug(message)
        
    return (tile_i, tile_j, tile, point_lon, point_lat, dist_min, theta)

########################################################################################
#
########################################################################################
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

def compare_hist_and_IC_points(location, hist_lon, hist_lat, IC_lon, IC_lat, hist_dist, angle_to_hist_point, IC_dist, angle_to_IC_point):
    #determine distance and angle from IC point to hist point
    dist = haversine_distance(IC_lat,IC_lon,hist_lat,hist_lon)
    theta = math.atan2(math.sin(math.radians(hist_lon - IC_lon))*math.cos(math.radians(hist_lat)), math.cos(math.radians(IC_lat))*math.sin(math.radians(hist_lat)) - math.sin(math.radians(IC_lat))*math.cos(math.radians(hist_lat))*math.cos(math.radians(hist_lon - IC_lon)))
    theta_deg = math.fmod((math.degrees(theta) + 360.0), 360.0)
    
    message = 'Location summary'
    message += '\n The location as entered is {} deg E {} deg N'.format(location[0],location[1])
    message += '\n The closest point in the UFS history file is {} deg E {} deg N, or {} km away at a bearing of {} deg'.format(hist_lon, hist_lat, hist_dist, angle_to_hist_point)
    message += '\n The closest point in the UFS IC files (native grid) is {} deg E {} deg N, or {} km away at a bearing of {} deg'.format(IC_lon, IC_lat, IC_dist, angle_to_IC_point)
    message += '\n Therefore, the history point (used to define the SCM grid column) is {} km away from the closest UFS native grid poing at a bearing of {} deg'.format(dist, theta_deg)
    logging.info(message)
     
    return

def check_IC_hist_surface_compatibility(dir, i, j, surface_data, lam, old_chgres, tile):
    
    # Determine UFS history file format (tiled/quilted)
    if lam:
        filename_pattern = '*phyf000.nc'
    else:
        filename_pattern = '*sfcf000.nc'
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    slmsk_hist = read_NetCDF_surface_var(nc_file, "land",     j, i, old_chgres, 0)
    
    if (slmsk_hist != surface_data["slmsk"]):
        message = 'There is a mismatch between the UFS IC file and history file land masks: IC={}, history={}'.format(surface_data["slmsk"],slmsk_hist)
        logging.critical(message)
        raise Exception(message)
    else:
        message = 'UFS IC file and history file land masks match: IC={}, history={}'.format(surface_data["slmsk"],slmsk_hist)
        logging.info(message)
            
    nc_file.close()
    
    return

########################################################################################
#
########################################################################################
def read_NetCDF_var(nc_file, var_name, i, j):
    try:
        var = nc_file[var_name][j,i]
    except (KeyError, IndexError):
        message = "Variable {0} is not found in {1}. Filling with missing value.".format(var_name,nc_file.filepath())
        logging.debug(message)
        var = missing_value
    return var

########################################################################################
#
########################################################################################
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

def get_IC_data_from_UFS_history(dir, i, j, lam, tile):
    
    # Determine UFS history file format (tiled/quilted)
    if lam:
        filename_pattern = '*dynf000.nc'
    else:
        filename_pattern = '*atmf000.nc'
    
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(dir,filename))
    
    pfull   = nc_file['pfull'][::-1]*100.0
    phalf   = nc_file['phalf'][::-1]*100.0
    nlevs   = len(pfull)
    time    = nc_file['time'][0] #model hours since being initialized
    lon     = nc_file['lon'][i,j]
    lat     = nc_file['lat'][i,j]
    logging.debug('Grabbing data from history file point with lon/lat: {}/{}'.format(lon,lat))
    delz    = -1*nc_file['delz'][0,::-1,i,j] #derive zh
    ugrd    = nc_file['ugrd'][0,::-1,i,j]
    vgrd    = nc_file['vgrd'][0,::-1,i,j]
    spfh    = nc_file['spfh'][0,::-1,i,j]
    o3mr    = nc_file['o3mr'][0,::-1,i,j]
    clwmr   = nc_file['clwmr'][0,::-1,i,j]
    icmr    = nc_file['icmr'][0,::-1,i,j]
    pressfc = nc_file['pressfc'][0,i,j]
    tmp     = nc_file['tmp'][0,::-1,i,j]
    
    delz = np.asarray(delz)
    zh = np.zeros(nlevs)
    zh[0] = 0.5*delz[0]
    for k in range(1,nlevs):
        zh[k] = zh[k-1] + 0.5*delz[k-1] + 0.5*delz[k]
    
    dz = np.zeros(nlevs-1)
    for k in range(nlevs-1):
        dz[k] = zh[k+1]-zh[k]
    
    state = {
        "nlevs": nlevs,
        "zh": zh,
        "dz": dz,
        "ua": np.asarray(ugrd),
        "va": np.asarray(vgrd),
        "qv": np.asarray(spfh),
        "o3": np.asarray(o3mr),
        "ql": np.asarray(clwmr),
        "qi": np.asarray(icmr),
        "ps": pressfc,
        "ta": np.asarray(tmp),
        "pa": np.asarray(pfull),
        "pa_i": np.asarray(phalf)
    }

    nc_file.close()
    
    return state

########################################################################################
#
########################################################################################
def get_IC_data_from_UFS_ICs(dir, grid_dir, tile, i, j, old_chgres, lam):
    """Get the state, surface, and orographic data for the given tile and indices"""
    #returns dictionaries with the data
    
    vgrid_data = get_UFS_vgrid_data(grid_dir) #only needed for ak, bk to calculate pressure
    (state_data, error_msg) = get_UFS_state_data(vgrid_data, dir, tile, i, j, old_chgres, lam)
    
    return (state_data, error_msg)

########################################################################################
#
########################################################################################
def get_UFS_state_data(vgrid, dir, tile, i, j, old_chgres, lam):
    """Get the state data for the given tile and indices"""

    state = {}
    error_msg=None
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
    
    [u_s, u_n, v_w, v_e, unknown_grid] = get_zonal_and_meridional_winds_on_cd_grid(tile, dir, i, j, nc_file_data, lam)
    if unknown_grid:
        error_msg='unknown grid orientation'
        return(state,error_msg)
    
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
        "ps": np.asarray(ps_calc),
        "ta": temp_model_rev[0,::-1],
        "pa": pressure_model,
        "pa_i": pressure_model_interfaces
    }
        
    return (state,error_msg)

########################################################################################
#
########################################################################################
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

    # Get grid dimension
    nz,nx,ny = np.shape(nc_file_data['u_w'])

    
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
    unknown=False
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
            u_s = np.zeros(nz)
            u_n = np.zeros(nz)
            v_w = np.zeros(nz)
            v_e = np.zeros(nz)
            unknown=True
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
            u_s = np.zeros(nz)
            u_n = np.zeros(nz)
            v_w = np.zeros(nz)
            v_e = np.zeros(nz)
            unknown=True
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
            if (j < nx -1):
                v_e = nc_file_data['u_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ex) + nc_file_data['v_w'][:,j+1,i]*fv3_remap.inner_prod(e1, ey)
            else:
                v_e = nc_file_data['v_w'][:,j,i]
        else:
            u_s = np.zeros(nz)
            u_n = np.zeros(nz)
            v_w = np.zeros(nz)
            v_e = np.zeros(nz)
            unknown=True
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
            u_s = np.zeros(nz)
            u_n = np.zeros(nz)
            v_w = np.zeros(nz)
            v_e = np.zeros(nz)
            unknown=True
    
    nc_file_grid.close()
    
    return [u_s, u_n, v_w, v_e, unknown]

########################################################################################
#
########################################################################################
def get_UFS_surface_data(dir, tile, i, j, old_chgres, lam):
    """Get the surface data for the given tile and indices"""
    
    if lam:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.nc'))
    else:
        nc_file = Dataset('{0}/{1}'.format(dir,'sfc_data.tile{0}.nc'.format(tile)))
    
    #FV3/io/FV3GFS_io.F90/sfc_prop_restart_read was used as reference for variables that can be read in    
    
    #
    # read in scalars (2D)
    #

    # Noah LSM ( may contain variables needed for fractional land fraction)
    tsfco_in  = read_NetCDF_surface_var(nc_file, 'tsea',   i, j, old_chgres, 0)
    tg3_in    = read_NetCDF_surface_var(nc_file, 'tg3',    i, j, old_chgres, 0)
    uustar_in = read_NetCDF_surface_var(nc_file, 'uustar', i, j, old_chgres, 0)
    alvsf_in  = read_NetCDF_surface_var(nc_file, 'alvsf',  i, j, old_chgres, 0)
    alvwf_in  = read_NetCDF_surface_var(nc_file, 'alvwf',  i, j, old_chgres, 0)
    alnsf_in  = read_NetCDF_surface_var(nc_file, 'alnsf',  i, j, old_chgres, 0)
    alnwf_in  = read_NetCDF_surface_var(nc_file, 'alnwf',  i, j, old_chgres, 0)
    facsf_in  = read_NetCDF_surface_var(nc_file, 'facsf',  i, j, old_chgres, 0)
    facwf_in  = read_NetCDF_surface_var(nc_file, 'facwf',  i, j, old_chgres, 0)
    styp_in   = read_NetCDF_surface_var(nc_file, 'stype',  i, j, old_chgres, 0)
    slope_in  = read_NetCDF_surface_var(nc_file, 'slope',  i, j, old_chgres, 0)
    vtyp_in   = read_NetCDF_surface_var(nc_file, 'vtype',  i, j, old_chgres, 0)
    vfrac_in  = read_NetCDF_surface_var(nc_file, 'vfrac',  i, j, old_chgres, 0)
    shdmin_in = read_NetCDF_surface_var(nc_file, 'shdmin', i, j, old_chgres, 0)
    shdmax_in = read_NetCDF_surface_var(nc_file, 'shdmax', i, j, old_chgres, 0)
    zorlo_in  = read_NetCDF_surface_var(nc_file, 'zorl',   i, j, old_chgres, 0)
    slmsk_in  = read_NetCDF_surface_var(nc_file, 'slmsk',  i, j, old_chgres, 0)
    canopy_in = read_NetCDF_surface_var(nc_file, 'canopy', i, j, old_chgres, 0)
    hice_in   = read_NetCDF_surface_var(nc_file, 'hice',   i, j, old_chgres, 0)
    fice_in   = read_NetCDF_surface_var(nc_file, 'fice',   i, j, old_chgres, 0)
    tisfc_in  = read_NetCDF_surface_var(nc_file, 'tisfc',  i, j, old_chgres, 0)
    snwdph_in = read_NetCDF_surface_var(nc_file, 'snwdph', i, j, old_chgres, 0)
    snoalb_in = read_NetCDF_surface_var(nc_file, 'snoalb', i, j, old_chgres, 0)
    sheleg_in = read_NetCDF_surface_var(nc_file, 'sheleg', i, j, old_chgres, 0)
    f10m_in   = read_NetCDF_surface_var(nc_file, 'f10m',   i, j, old_chgres, 0)
    t2m_in    = read_NetCDF_surface_var(nc_file, 't2m',    i, j, old_chgres, 0)
    q2m_in    = read_NetCDF_surface_var(nc_file, 'q2m',    i, j, old_chgres, 0)
    ffmm_in   = read_NetCDF_surface_var(nc_file, 'ffmm',   i, j, old_chgres, 0)
    ffhh_in   = read_NetCDF_surface_var(nc_file, 'ffhh',   i, j, old_chgres, 0)
    tprcp_in  = read_NetCDF_surface_var(nc_file, 'tprcp',  i, j, old_chgres, 0)
    srflag_in = read_NetCDF_surface_var(nc_file, 'srflag', i, j, old_chgres, 0)
    sncovr_in = read_NetCDF_surface_var(nc_file, 'sncovr', i, j, old_chgres, 0)
    tsfcl_in  = read_NetCDF_surface_var(nc_file, 'tsea',   i, j, old_chgres, 0)
    zorll_in  = read_NetCDF_surface_var(nc_file, 'zorl',   i, j, old_chgres, 0)   
    zorli_in  = read_NetCDF_surface_var(nc_file, 'zorl',   i, j, old_chgres, 0)
    if (snwdph_in > 0):
        sncovr_in = 1.0
    else:
        sncovr_in = 0.0

    # present when cplwav = T
    zorlw_in = read_NetCDF_surface_var(nc_file, 'zorlw', i, j, old_chgres, 0)
    
    # NSST
    tref_in    = read_NetCDF_surface_var(nc_file, 'tref',    i, j, old_chgres, 0)
    z_c_in     = read_NetCDF_surface_var(nc_file, 'z_c',     i, j, old_chgres, 0)
    c_0_in     = read_NetCDF_surface_var(nc_file, 'c_0',     i, j, old_chgres, 0)
    c_d_in     = read_NetCDF_surface_var(nc_file, 'c_d',     i, j, old_chgres, 0)
    w_0_in     = read_NetCDF_surface_var(nc_file, 'w_0',     i, j, old_chgres, 0)
    w_d_in     = read_NetCDF_surface_var(nc_file, 'w_d',     i, j, old_chgres, 0)
    xt_in      = read_NetCDF_surface_var(nc_file, 'xt',      i, j, old_chgres, 0)
    xs_in      = read_NetCDF_surface_var(nc_file, 'xs',      i, j, old_chgres, 0)
    xu_in      = read_NetCDF_surface_var(nc_file, 'xu',      i, j, old_chgres, 0)
    xv_in      = read_NetCDF_surface_var(nc_file, 'xv',      i, j, old_chgres, 0)
    xz_in      = read_NetCDF_surface_var(nc_file, 'xz',      i, j, old_chgres, 0)
    zm_in      = read_NetCDF_surface_var(nc_file, 'zm',      i, j, old_chgres, 0)
    xtts_in    = read_NetCDF_surface_var(nc_file, 'xtts',    i, j, old_chgres, 0)
    xzts_in    = read_NetCDF_surface_var(nc_file, 'xzts',    i, j, old_chgres, 0)
    d_conv_in  = read_NetCDF_surface_var(nc_file, 'd_conv',  i, j, old_chgres, 0)
    ifd_in     = read_NetCDF_surface_var(nc_file, 'ifd',     i, j, old_chgres, 0)
    dt_cool_in = read_NetCDF_surface_var(nc_file, 'dt_cool', i, j, old_chgres, 0)
    qrain_in   = read_NetCDF_surface_var(nc_file, 'qrain',   i, j, old_chgres, 0)

    # NoahMP LSM
    snowxy_in     = read_NetCDF_surface_var(nc_file, 'snowxy',     i, j, old_chgres, 0)
    tvxy_in       = read_NetCDF_surface_var(nc_file, 'tvxy',       i, j, old_chgres, 0)
    tgxy_in       = read_NetCDF_surface_var(nc_file, 'tgxy',       i, j, old_chgres, 0)
    canicexy_in   = read_NetCDF_surface_var(nc_file, 'canicexy',   i, j, old_chgres, 0)
    canliqxy_in   = read_NetCDF_surface_var(nc_file, 'canliqxy',   i, j, old_chgres, 0)
    eahxy_in      = read_NetCDF_surface_var(nc_file, 'eahxy',      i, j, old_chgres, 0)
    tahxy_in      = read_NetCDF_surface_var(nc_file, 'tahxy',      i, j, old_chgres, 0)
    cmxy_in       = read_NetCDF_surface_var(nc_file, 'cmxy',       i, j, old_chgres, 0)
    chxy_in       = read_NetCDF_surface_var(nc_file, 'chxy',       i, j, old_chgres, 0)
    fwetxy_in     = read_NetCDF_surface_var(nc_file, 'fwetxy',     i, j, old_chgres, 0)
    sneqvoxy_in   = read_NetCDF_surface_var(nc_file, 'sneqvoxy',   i, j, old_chgres, 0)
    alboldxy_in   = read_NetCDF_surface_var(nc_file, 'alboldxy',   i, j, old_chgres, 0)
    qsnowxy_in    = read_NetCDF_surface_var(nc_file, 'qsnowxy',    i, j, old_chgres, 0)
    wslakexy_in   = read_NetCDF_surface_var(nc_file, 'wslakexy',   i, j, old_chgres, 0)
    zwtxy_in      = read_NetCDF_surface_var(nc_file, 'zwtxy',      i, j, old_chgres, 0)
    waxy_in       = read_NetCDF_surface_var(nc_file, 'waxy',       i, j, old_chgres, 0)
    wtxy_in       = read_NetCDF_surface_var(nc_file, 'wtxy',       i, j, old_chgres, 0)
    lfmassxy_in   = read_NetCDF_surface_var(nc_file, 'lfmassxy',   i, j, old_chgres, 0)
    rtmassxy_in   = read_NetCDF_surface_var(nc_file, 'rtmassxy',   i, j, old_chgres, 0)
    stmassxy_in   = read_NetCDF_surface_var(nc_file, 'stmassxy',   i, j, old_chgres, 0)
    woodxy_in     = read_NetCDF_surface_var(nc_file, 'woodxy',     i, j, old_chgres, 0)
    stblcpxy_in   = read_NetCDF_surface_var(nc_file, 'stblcpxy',   i, j, old_chgres, 0)
    fastcpxy_in   = read_NetCDF_surface_var(nc_file, 'fastcpxy',   i, j, old_chgres, 0)
    xsaixy_in     = read_NetCDF_surface_var(nc_file, 'xsaixy',     i, j, old_chgres, 0)
    xlaixy_in     = read_NetCDF_surface_var(nc_file, 'xlaixy',     i, j, old_chgres, 0)
    taussxy_in    = read_NetCDF_surface_var(nc_file, 'taussxy',    i, j, old_chgres, 0)
    smcwtdxy_in   = read_NetCDF_surface_var(nc_file, 'smcwtdxy',   i, j, old_chgres, 0)
    deeprechxy_in = read_NetCDF_surface_var(nc_file, 'deeprechxy', i, j, old_chgres, 0)
    rechxy_in     = read_NetCDF_surface_var(nc_file, 'rechxy',     i, j, old_chgres, 0)
    
    # RUC LSM
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
    
    #
    # read in profiles (3D)
    #

    # land_state
    stc_in = read_NetCDF_surface_var(nc_file, 'stc', i, j, old_chgres, missing_variable_soil_layers)
    smc_in = read_NetCDF_surface_var(nc_file, 'smc', i, j, old_chgres, missing_variable_soil_layers)
    slc_in = read_NetCDF_surface_var(nc_file, 'slc', i, j, old_chgres, missing_variable_soil_layers)
    
    # NoahMP LSM
    snicexy_in = read_NetCDF_surface_var(nc_file, 'snicexy', i, j, old_chgres, missing_variable_snow_layers)
    snliqxy_in = read_NetCDF_surface_var(nc_file, 'snliqxy', i, j, old_chgres, missing_variable_snow_layers)
    tsnoxy_in  = read_NetCDF_surface_var(nc_file, 'tsnoxy',  i, j, old_chgres, missing_variable_snow_layers)
    smoiseq_in = read_NetCDF_surface_var(nc_file, 'smoiseq', i, j, old_chgres, missing_variable_soil_layers)
    zsnsoxy_in = read_NetCDF_surface_var(nc_file, 'zsnsoxy', i, j, old_chgres, missing_variable_soil_layers + missing_variable_snow_layers)
     
    # RUC LSM
    tslb_in  = read_NetCDF_surface_var(nc_file, 'tslb',  i, j, old_chgres, missing_variable_soil_layers)
    smois_in = read_NetCDF_surface_var(nc_file, 'smois', i, j, old_chgres, missing_variable_soil_layers)
    sh2o_in  = read_NetCDF_surface_var(nc_file, 'sh2o',  i, j, old_chgres, missing_variable_soil_layers)
    smfr_in  = read_NetCDF_surface_var(nc_file, 'smfr',  i, j, old_chgres, missing_variable_soil_layers)
    flfr_in  = read_NetCDF_surface_var(nc_file, 'flfr',  i, j, old_chgres, missing_variable_soil_layers)
    
    # fractional grid
    tiice_in = read_NetCDF_surface_var(nc_file, 'tiice', i, j, old_chgres, missing_variable_ice_layers)

    # soil color (From  the UFS FV3: io/fv3atm_sfc_io.F90)
    if (slmsk_in == 1):
        scolor_in = 4
    else:
        scolor_in = 1

    #
    nc_file.close()
    
    #put data in a dictionary
    surface = {
        # Noah LSM
        "tsfco": tsfco_in,  
        "tg3": tg3_in,
        "uustar": uustar_in,
        "alvsf": alvsf_in,
        "alvwf": alvwf_in,
        "alnsf": alnsf_in,
        "alnwf": alnwf_in,
        "facsf": facsf_in,
        "facwf": facwf_in,
        "soiltyp": styp_in,
        "scolor": scolor_in,
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
        "snowd": snwdph_in,
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

########################################################################################
#
########################################################################################
def get_UFS_oro_data(dir, tile, i, j, lam, old_chgres):
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
    #fractional landmask variables
    land_frac_in  = read_NetCDF_var(nc_file, "land_frac",  i, j)
    #lake variables (optional)
    lake_frac_in  = read_NetCDF_var(nc_file, "lake_frac",  i, j)
    lake_depth_in = read_NetCDF_var(nc_file, "lake_depth", i, j)
    
    vegtype_frac_in  = read_NetCDF_surface_var(nc_file, "vegetation_type_pct",  i, j, old_chgres, missing_variable_vegetation_categories)
    soiltype_frac_in  = read_NetCDF_surface_var(nc_file, "soil_type_pct",  i, j, old_chgres, missing_variable_soil_categories)

    #
    nc_file.close()
    
    # Store data in a dictionary
    oro = {"stddev":    stddev_in,
           "convexity": convexity_in,
           "oa1":       oa1_in,
           "oa2":       oa2_in,
           "oa3":       oa3_in,
           "oa4":       oa4_in,
           "ol1":       ol1_in,
           "ol2":       ol2_in,
           "ol3":       ol3_in,
           "ol4":       ol4_in,
           "theta":     theta_in,
           "gamma":     gamma_in,
           "sigma":     sigma_in,
           "elvmax":    elvmax_in,
           "oro":       orog_filt_in,
           "oro_uf":    orog_raw_in,
           "landfrac":  land_frac_in,
           "lakefrac":  lake_frac_in,
           "lakedepth": lake_depth_in,
           "vegtype_frac": vegtype_frac_in,
           "soiltype_frac": soiltype_frac_in}

    return oro

########################################################################################
#
########################################################################################
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

########################################################################################
#
########################################################################################
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

########################################################################################
#
########################################################################################
def search_in_dict(listin,name):
    for count, dictionary in enumerate(listin):
        if dictionary["name"] == name:
            return count

def expand_neighbors_for_geostrophic_balance(dir, hist_i, hist_j, neighbors, dx, dy, nlevs, lam):
    
    # Determine UFS history file format (tiled/quilted)
    if lam:
        atm_ftag = 'atmf*.tile{0}.nc'.format(tile)
        sfc_ftag = 'sfcf*.tile{0}.nc'.format(tile)
    else:
        atm_ftag = '*atmf*.nc'
        sfc_ftag = '*sfcf*.nc'

    # Get list of UFS history files with 3D ATMospheric state variables.
    atm_filenames = []
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, atm_ftag):
          atm_filenames.append(f_name)
    if not atm_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.                \
                  format(atm_ftag,dir)
        logging.critical(message)
        raise Exception(message)
    atm_filenames = sorted(atm_filenames)
    n_filesA      = len(atm_filenames)
    
    npts  = 2*n_forcing_halo_points+1
    
    lon = np.zeros((n_filesA,npts,npts))
    lat = np.zeros((n_filesA,npts,npts))
    time = np.zeros((n_filesA))
    u_wind = np.zeros((n_filesA,nlevs,npts,npts))
    v_wind = np.zeros((n_filesA,nlevs,npts,npts))
    wind_speed = np.zeros((n_filesA,nlevs))
    wind_speed_time_mean = np.zeros((nlevs))
    
    nc_file = Dataset('{0}/{1}'.format(dir,atm_filenames[0]))
    nc_file.set_always_mask(False)
    longitude_all = np.asarray(nc_file['lon'])
    latitude_all  = np.asarray(nc_file['lat'])
    nc_file.close()
    
    #determine average wind speeds in the existing forcing stencil
    
    # Read in 3D UFS history files
    for count, filename in enumerate(atm_filenames, start=1):
        nc_file = Dataset('{0}/{1}'.format(dir,filename))
        nc_file.set_always_mask(False)
        
        time[count-1] = nc_file['time'][0]
        
        for ii in range(2*n_forcing_halo_points+1):
            for jj in range(2*n_forcing_halo_points+1):
                current_ii_index = neighbors[0,ii,jj]
                current_jj_index = neighbors[1,ii,jj]
                
                lon[count-1,ii,jj]      = nc_file['lon'][current_ii_index,current_jj_index]
                lat[count-1,ii,jj]      = nc_file['lat'][current_ii_index,current_jj_index]                
                u_wind[count-1,:,ii,jj] = nc_file['ugrd'][0,::-1,current_ii_index,current_jj_index]
                v_wind[count-1,:,ii,jj] = nc_file['vgrd'][0,::-1,current_ii_index,current_jj_index]
        
        for k in range(nlevs):
            wind_speed[count-1,k] = np.sqrt(np.mean(u_wind[count-1,k,:,:])**2 + np.mean(v_wind[count-1,k,:,:]**2))
        
        nc_file.close()
    
    critical_geos_scale_npts = np.zeros((nlevs))
    mean_lat = np.mean(lat[:,:,:])
    coriolis = 2.0*Omega*np.sin(mean_lat*deg_to_rad)
    grid_length_scale = np.sqrt(np.mean(dx)**2 + np.mean(dy)**2)*1.0E3
    for k in range(nlevs):
        wind_speed_time_mean[k] = np.mean(wind_speed[:,k])
        critical_geos_scale = wind_speed_time_mean[k]/(coriolis*critical_Rossby)
        critical_geos_scale_npts[k] = critical_geos_scale/grid_length_scale
    floor = np.floor(np.mean(critical_geos_scale_npts))
    ceiling = np.ceil(np.mean(critical_geos_scale_npts))
    if (floor%2 == 1):
        odd_pts = int(floor)
    else:
        odd_pts = int(ceiling)
        
    message = 'Need {} points for approximate geostrophic balance assuming a Rossby number of {}'.format(odd_pts, critical_Rossby)
    logging.debug(message)
    
    n_geo_halo_points = int((odd_pts - 1)/2)
    
    if hist_i < n_geo_halo_points or hist_i > latitude_all.shape[0]-1-n_geo_halo_points:
        message = 'There are not enough neighboring grid cells between the chosen point and the poles for calculation of geostrophic winds'
        logging.critical(message)
        raise Exception(message)
    else:
        #longitude = (n rows of latitude, 2*n columns of longitude), i index increases toward EAST; rows are circles of latitude, columns are meridians
        #latitude  = (n rows of latitude, 2*n columns of longitude), j index increases toward NORTH; rows are circles of latitude, columns are meridians
        new_neighbors = np.mgrid[-n_geo_halo_points:n_geo_halo_points+1,-n_geo_halo_points:n_geo_halo_points+1]
        new_neighbors[0,:,:] = new_neighbors[0,:,:] + hist_i
        new_neighbors[1,:,:] = new_neighbors[1,:,:] + hist_j
        
        
        # the closest point is near the western edge of the prime meridian, but longitude should be cyclic, so get neighbors on the eastern/western side of the prime meridian as necessary
        new_neighbors[1,:,:] = new_neighbors[1,:,:]%longitude_all.shape[1]
    
    new_neighbors = new_neighbors.astype(int) 
    #message = 'orig neighbors',neighbors,latitude_all[neighbors[0,:,0],neighbors[1,0,:]],longitude_all[neighbors[0,:,0],neighbors[1,0,:]]
    #logging.debug(message)
    #message = 'new neighbors',new_neighbors,latitude_all[new_neighbors[0,:,0],new_neighbors[1,0,:]],longitude_all[new_neighbors[0,:,0],new_neighbors[1,0,:]]
    #logging.debug(message)
        
    #calc dx, dy (n_forcing_halo_points*2 x n_forcing_halo_points*2)
    new_dx = np.zeros((n_geo_halo_points*2+1,n_geo_halo_points*2))
    new_dy = np.zeros((n_geo_halo_points*2,n_geo_halo_points*2+1))
    for ii in range(2*n_geo_halo_points+1):
         for jj in range(2*n_geo_halo_points):
             current_ii_index = new_neighbors[0,ii,jj]
             current_jj_index = new_neighbors[1,ii,jj]
             eastern_jj_index = new_neighbors[1,ii,jj+1]
             new_dx[ii,jj] = haversine_distance(latitude_all[current_ii_index,current_jj_index],longitude_all[current_ii_index,current_jj_index],latitude_all[current_ii_index,eastern_jj_index],longitude_all[current_ii_index,eastern_jj_index])
        
    for ii in range(2*n_geo_halo_points):
         for jj in range(2*n_geo_halo_points+1):
             current_ii_index  = new_neighbors[0,ii,jj]
             northern_ii_index = new_neighbors[0,ii+1,jj]
             current_jj_index  = new_neighbors[1,ii,jj]
             new_dy[ii,jj] = haversine_distance(latitude_all[current_ii_index,current_jj_index],longitude_all[current_ii_index,current_jj_index],latitude_all[northern_ii_index,current_jj_index],longitude_all[northern_ii_index,current_jj_index])
    
    return (n_geo_halo_points, new_neighbors, new_dx, new_dy)

def get_UFS_forcing_data_advective_tendency(dir, i, j, tile, neighbors, dx, dy, nlevs, lam, save_comp_data, vertical_method, geos_wind_forcing, wind_nudge, geo_pts, geo_neighbors, geo_dx, geo_dy, use_actual_geos_wind):
    
    # Determine UFS history file format (tiled/quilted)
    if lam:
        atm_ftag = '*dynf*.nc'
        sfc_ftag = '*phyf*.nc'
    else:
        atm_ftag = '*atmf*.nc'
        sfc_ftag = '*sfcf*.nc'

    # Get list of UFS history files with 3D ATMospheric state variables.
    atm_filenames = []
    for f_name in os.listdir(dir):
       if fnmatch.fnmatch(f_name, atm_ftag):
          atm_filenames.append(f_name)
    if not atm_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.                \
                  format(atm_ftag,dir)
        logging.critical(message)
        raise Exception(message)
    atm_filenames = sorted(atm_filenames)
    n_filesA      = len(atm_filenames)
    
    npts  = 2*n_forcing_halo_points+1
    npts_geo = 2*geo_pts+1
        
    lon = np.zeros((n_filesA,npts,npts))
    lat = np.zeros((n_filesA,npts,npts))
    time = np.zeros((n_filesA))
    u_wind = np.zeros((n_filesA,nlevs,npts,npts))
    v_wind = np.zeros((n_filesA,nlevs,npts,npts))
    w_wind = np.zeros((n_filesA,nlevs,npts,npts))
    temp = np.zeros((n_filesA,nlevs,npts,npts))
    qv = np.zeros((n_filesA,nlevs,npts,npts))
    pres = np.zeros((n_filesA,nlevs))
    presi = np.zeros((n_filesA,nlevs+1))
    pressfc = np.zeros((n_filesA,npts,npts))
    delz = np.zeros((n_filesA,nlevs,npts,npts))
    if(geos_wind_forcing and use_actual_geos_wind):
        temp_geo = np.zeros((n_filesA,nlevs,npts_geo,npts_geo))
        qv_geo = np.zeros((n_filesA,nlevs,npts_geo,npts_geo))
        hgtsfc = np.zeros((n_filesA,npts,npts))
    
    # Read in 3D UFS history files
    for count, filename in enumerate(atm_filenames, start=1):
        nc_file = Dataset('{0}/{1}'.format(dir,filename))
        nc_file.set_always_mask(False)
        
        time[count-1] = nc_file['time'][0]
        
        for ii in range(2*n_forcing_halo_points+1):
            for jj in range(2*n_forcing_halo_points+1):
                current_ii_index = neighbors[0,ii,jj]
                current_jj_index = neighbors[1,ii,jj]
                
                lon[count-1,ii,jj]      = nc_file['lon'][current_ii_index,current_jj_index]
                lat[count-1,ii,jj]      = nc_file['lat'][current_ii_index,current_jj_index]                
                u_wind[count-1,:,ii,jj] = nc_file['ugrd'][0,::-1,current_ii_index,current_jj_index]
                v_wind[count-1,:,ii,jj] = nc_file['vgrd'][0,::-1,current_ii_index,current_jj_index]
                w_wind[count-1,:,ii,jj] = nc_file['dzdt'][0,::-1,current_ii_index,current_jj_index]
                temp[count-1,:,ii,jj]   = nc_file['tmp'][0,::-1,current_ii_index,current_jj_index]
                qv[count-1,:,ii,jj]     = nc_file['spfh'][0,::-1,current_ii_index,current_jj_index]
                pressfc[count-1,ii,jj]  = nc_file['pressfc'][0,current_ii_index,current_jj_index]
                delz[count-1,:,ii,jj]   = -1*nc_file['delz'][0,::-1,current_ii_index,current_jj_index] #derive zh
            
        if (geos_wind_forcing and use_actual_geos_wind):
            for ii in range(npts_geo):
                for jj in range(npts_geo):
                    current_ii_index = geo_neighbors[0,ii,jj]
                    current_jj_index = geo_neighbors[1,ii,jj]
                    temp_geo[count-1,:,ii,jj]   = nc_file['tmp'][0,::-1,current_ii_index,current_jj_index]
                    qv_geo[count-1,:,ii,jj]     = nc_file['spfh'][0,::-1,current_ii_index,current_jj_index]
                    hgtsfc[count-1,ii,jj]       = nc_file['hgtsfc'][0,current_ii_index,current_jj_index]
        
        pres[count-1,:] = nc_file['pfull'][::-1]*100.0
        presi[count-1,:] = nc_file['phalf'][::-1]*100.0
    
    #check for poor time resolution (e.g. if mean number of grid points traversed by the wind between data intervals is much larger than the number of grid points used for calculating the advective terms)
    # t=0
    # for k in range(u_wind.shape[1]):
    #     mean_u = np.mean(np.sqrt(u_wind[t,k,:,:]*u_wind[t,k,:,:]))
    #     mean_dist = 43200*mean_u*1.0E-3
    #     print("k, mean_dist, dx, grid points traversed in timestep: ",k, mean_dist, np.mean(dx), mean_dist/np.mean(dx))
    
    center = n_forcing_halo_points
    smoothed_u  = np.zeros((n_filesA,nlevs,npts,npts))
    smoothed_v  = np.zeros((n_filesA,nlevs,npts,npts))
    smoothed_w  = np.zeros((n_filesA,nlevs))
    smoothed_T  = np.zeros((n_filesA,nlevs,npts,npts))
    smoothed_qv = np.zeros((n_filesA,nlevs,npts,npts))
    h_advec_u   = np.zeros((n_filesA,nlevs))
    h_advec_v   = np.zeros((n_filesA,nlevs))
    h_advec_T   = np.zeros((n_filesA,nlevs))
    h_advec_qv  = np.zeros((n_filesA,nlevs))
    v_advec_u   = np.zeros((n_filesA,nlevs))
    v_advec_v   = np.zeros((n_filesA,nlevs))
    v_advec_T   = np.zeros((n_filesA,nlevs))
    v_advec_qv  = np.zeros((n_filesA,nlevs))
    t_advec_u   = np.zeros((n_filesA,nlevs))
    t_advec_v   = np.zeros((n_filesA,nlevs))
    t_advec_T   = np.zeros((n_filesA,nlevs))
    t_advec_qv  = np.zeros((n_filesA,nlevs))
    u_g         = np.zeros((n_filesA,nlevs))
    v_g         = np.zeros((n_filesA,nlevs))
    if(geos_wind_forcing and use_actual_geos_wind):
        phii        = np.zeros((n_filesA,nlevs+1,npts_geo,npts_geo))
        phil        = np.zeros((n_filesA,nlevs,npts_geo,npts_geo))
    for t in range(n_filesA):
        #velocities gets very noisy in the UFS above the tropopause; define a linear return-to-zero profile above some pressure
        k_p_top = np.where(pres[t,:] <= top_pres_for_forcing)[0][0]
        
        #smooth velocity profile
        n_smooth_points = 5
        for ii in range(npts):
            for jj in range(npts):
                #the last argument defines the averaging "window" in grid layers; greater numbers results in greater smoothing
                smoothed_u[t,:,ii,jj]  = moving_average(u_wind[t,:,ii,jj], n_smooth_points)
                smoothed_v[t,:,ii,jj]  = moving_average(v_wind[t,:,ii,jj], n_smooth_points)
                smoothed_T[t,:,ii,jj]  = moving_average(temp  [t,:,ii,jj], n_smooth_points)
                smoothed_qv[t,:,ii,jj] = moving_average(qv    [t,:,ii,jj], n_smooth_points)
                
                #velocities gets very noisy in the UFS above the tropopause; define a linear return-to-zero profile above some pressure
                for k in range(k_p_top+1,nlevs):
                    lifrac = pres[t,k]/pres[t,k_p_top]
                    smoothed_u[t,k,ii,jj] = lifrac*smoothed_u[t,k_p_top,ii,jj]
                    smoothed_v[t,k,ii,jj] = lifrac*smoothed_v[t,k_p_top,ii,jj]

        smoothed_w[t,:] = moving_average(w_wind[t,:,center,center], n_smooth_points)
        
        for k in range(k_p_top+1,nlevs):
            lifrac = pres[t,k]/pres[t,k_p_top]
            smoothed_w[t,k] = lifrac*smoothed_w[t,k_p_top]
        
        for k in range(nlevs):
            dTdx = centered_diff_derivative(smoothed_T[t,k,center,:],dx[center,:]*1.0E3)
            dTdy = centered_diff_derivative(smoothed_T[t,k,:,center],dy[:,center]*1.0E3)
            dqdx = centered_diff_derivative(smoothed_qv[t,k,center,:],dx[center,:]*1.0E3)
            dqdy = centered_diff_derivative(smoothed_qv[t,k,:,center],dy[:,center]*1.0E3)
            dudx = centered_diff_derivative(smoothed_u[t,k,center,:],dx[center,:]*1.0E3)
            dudy = centered_diff_derivative(smoothed_u[t,k,:,center],dy[:,center]*1.0E3)
            dvdx = centered_diff_derivative(smoothed_v[t,k,center,:],dx[center,:]*1.0E3)
            dvdy = centered_diff_derivative(smoothed_v[t,k,:,center],dy[:,center]*1.0E3)
            mean_u = np.mean(smoothed_u[t,k,center,center-1:center+1+1])
            mean_v = np.mean(smoothed_v[t,k,center-1:center+1+1,center])
            
            h_advec_T[t,k]  = -mean_u*dTdx - mean_v*dTdy #K s-1
            h_advec_qv[t,k] = -mean_u*dqdx - mean_v*dqdy #kg kg-1 s-1
            h_advec_u[t,k]  = -mean_u*dudx - mean_v*dudy #m s-1 s-1
            h_advec_v[t,k]  = -mean_u*dvdx - mean_v*dvdy #m s-1 s-1        
            
        if (vertical_method == 1):
            v_advec_u[t,0] = 0.0
            v_advec_u[t,nlevs-1] = 0.0
            v_advec_v[t,0] = 0.0
            v_advec_v[t,nlevs-1] = 0.0
            v_advec_T[t,0] = 0.0
            v_advec_T[t,nlevs-1] = 0.0
            v_advec_qv[t,0] = 0.0
            v_advec_qv[t,nlevs-1] = 0.0
            for k in range(1, nlevs-1):
                 w_asc = max(smoothed_w[t,k], 0.0)
                 w_des = min(smoothed_w[t,k], 0.0)
                 #noisy wind profiles necessitate using smoothed profiles; this doesn't appear to be needed for T and qv
                 gradient_T_asc = (temp[t,k,center,center] - temp[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 gradient_T_des = (temp[t,k+1,center,center] - temp[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 #gradient_T_asc = (smoothed_T[t,k,center,center] - smoothed_T[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 #gradient_T_des = (smoothed_T[t,k+1,center,center] - smoothed_T[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 gradient_qv_asc = (qv[t,k,center,center] - qv[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 gradient_qv_des = (qv[t,k+1,center,center] - qv[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 #gradient_qv_asc = (smoothed_qv[t,k,center,center] - smoothed_qv[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 #gradient_qv_des = (smoothed_qv[t,k+1,center,center] - smoothed_qv[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 #gradient_u_asc = (u_wind[t,k,center,center] - u_wind[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 #gradient_u_des = (u_wind[t,k+1,center,center] - u_wind[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 gradient_u_asc = (smoothed_u[t,k,center,center] - smoothed_u[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 gradient_u_des = (smoothed_u[t,k+1,center,center] - smoothed_u[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 #gradient_v_asc = (v_wind[t,k,center,center] - v_wind[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 #gradient_v_des = (v_wind[t,k+1,center,center] - v_wind[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 gradient_v_asc = (smoothed_v[t,k,center,center] - smoothed_v[t,k-1,center,center])/(pres[t,k]-pres[t,k-1])
                 gradient_v_des = (smoothed_v[t,k+1,center,center] - smoothed_v[t,k,center,center])/(pres[t,k+1]-pres[t,k])
                 
                 rho = pres[t,k]/(rdgas*temp[t,k,center,center])
                 adiabatic_exp_comp_term = -(w_asc + w_des)*grav/cp
                 v_advec_u[t,k] = rho*grav*(w_asc*gradient_u_asc + w_des*gradient_u_des)
                 v_advec_v[t,k] = rho*grav*(w_asc*gradient_v_asc + w_des*gradient_v_des)
                 v_advec_T[t,k] = rho*grav*(w_asc*gradient_T_asc + w_des*gradient_T_des) + adiabatic_exp_comp_term
                 v_advec_qv[t,k] = rho*grav*(w_asc*gradient_qv_asc + w_des*gradient_qv_des)
        
        if (geos_wind_forcing):
            coriolis = 2.0*Omega*np.sin(lat[t,center,center]*deg_to_rad)
            
            if (use_actual_geos_wind):
                #calc geopotential at interface levels, starting from surface, convert to full pressure levels
                for ii in range(npts_geo):
                     for jj in range(npts_geo):
                         if (use_oro_for_geos_wind):
                             phii[t,0,ii,jj] = hgtsfc[t,ii,jj]*grav
                         else:
                             phii[t,0,ii,jj] = 0.0 #don't need to include orography -- geostrophic winds should assume geopotential over the geoid
                         for k in range(1,nlevs+1):
                             phii[t,k,ii,jj] = phii[t,k-1,ii,jj] - rdgas*temp_geo[t,k-1,ii,jj]*(1.0 + zvir*qv_geo[t,k-1,ii,jj])*np.log(presi[t,k]/presi[t,k-1])
                         
                         for k in range(0,nlevs):
                             phil[t,k,ii,jj] = 0.5*(phii[t,k,ii,jj] + phii[t,k+1,ii,jj])
            
                for k in range(1,nlevs):
                    dphidx = np.zeros((npts_geo,npts_geo))
                    dphidy = np.zeros((npts_geo,npts_geo))
                    #valid points:
                        #dphidx = np.zeros((npts_geo,npts_geo-2))
                        #dphidy = np.zeros((npts_geo-2,npts_geo))
                    for ii in range(0,npts_geo):
                        for jj in range(1,npts_geo-1):
                            dphidx[ii,jj] = centered_diff_derivative(phil[t,k,ii,jj-1:jj+2],geo_dx[ii,jj-1:jj+1]*1.0E3)
                    mean_dphidx = np.mean(dphidx[0:npts_geo,1:npts_geo-1])
                
                    for ii in range(1,npts_geo-1):
                        for jj in range(0,npts_geo):
                            dphidy[ii,jj] = centered_diff_derivative(phil[t,k,ii-1:ii+2,jj],geo_dy[ii-1:ii+1,jj]*1.0E3)
                    mean_dphidy = np.mean(dphidy[1:npts_geo-1,0:npts_geo])
                
                    u_g[t,k] = -mean_dphidy/coriolis
                    v_g[t,k] = mean_dphidx/coriolis
            else:
                for k in range(0,nlevs):
                    u_g[t,k] = np.mean(u_wind[t,k,:,:])
                    v_g[t,k] = np.mean(v_wind[t,k,:,:])
                k_pbl_top = np.where(pres[t,:] <= PBL_top_for_geos)[0][0]
                u_g[t,0] = 0.0
                v_g[t,0] = 0.0
                for k in range(1,k_pbl_top):
                    lifrac = (pres[t,k] - pres[t,0])/(pres[t,k_pbl_top] - pres[t,0])
                    u_g[t,k] = lifrac*u_g[t,k_pbl_top]
                    v_g[t,k] = lifrac*v_g[t,k_pbl_top]
            
    
    if (vertical_method == 1):
        t_advec_T[:,:] = h_advec_T[:,:] + v_advec_T[:,:]
        t_advec_qv[:,:] = h_advec_qv[:,:] + v_advec_qv[:,:]
        t_advec_u[:,:] = h_advec_u[:,:] + v_advec_u[:,:]
        t_advec_v[:,:] = h_advec_v[:,:] + v_advec_v[:,:]
    else:
        t_advec_T[:,:] = h_advec_T[:,:]
        t_advec_qv[:,:] = h_advec_qv[:,:]
        t_advec_u[:,:] = h_advec_u[:,:]
        t_advec_v[:,:] = h_advec_v[:,:]
    
    if (save_comp_data):
        # Read in 2D UFS history files
        phystends = [{"name":"dtend_temp_lw"},     {"name":"dtend_temp_sw"},  {"name":"dtend_temp_pbl"},   {"name":"dtend_temp_deepcnv"},\
                     {"name":"dtend_temp_shalcnv"},{"name":"dtend_temp_mp"},  {"name":"dtend_temp_orogwd"},{"name":"dtend_temp_phys"},   \
                     {"name":"dtend_u_pbl"},       {"name":"dtend_v_pbl"},    {"name":"dtend_u_orogwd"},   {"name":"dtend_v_orogwd"},    \
                     {"name":"dtend_u_deepcnv"},   {"name":"dtend_v_deepcnv"},{"name":"dtend_u_shalcnv"},  {"name":"dtend_v_shalcnv"},   \
                     {"name":"dtend_u_phys"},      {"name":"dtend_v_phys"},   {"name":"dtend_qv_pbl"},     {"name":"dtend_qv_deepcnv"},  \
                     {"name":"dtend_qv_shalcnv"},  {"name":"dtend_qv_mp"},    {"name":"dtend_qv_phys"},    {"name":"dtend_poop"}]
        for phystend in phystends: phystend["values"] = []

        # Variables to be added to "UFS comparison file"
        vars2d  =[{"name":"spfh2m"},    {"name":"tmp2m"},     {"name":"dswrf_ave"}, \
                  {"name":"ulwrf_ave"}, {"name":"lhtfl_ave"}, {"name":"shtfl_ave"}, \
                  {"name":"dswrf"},     {"name":"ulwrf"},     {"name":"lhtfl"},     \
                  {"name":"shtfl"},     {"name":"pwat"},      {"name":"vgrd10m"},   \
                  {"name":"ugrd10m"}]
        for var2d in vars2d: var2d["values"] = []
        
        # Get list of UFS history files with 2D fields.
        sfc_filenames = []
        for f_name in os.listdir(dir):
           if fnmatch.fnmatch(f_name, sfc_ftag):
              sfc_filenames.append(f_name)
        if not sfc_filenames:
            message = 'No filenames matching the pattern {0} found in {1}'.                \
                      format(sfc_ftag,dir)
            logging.critical(message)
            raise Exception(message)
        sfc_filenames = sorted(sfc_filenames)
        n_filesS      = len(sfc_filenames)

        if (n_filesS == n_filesA):
            n_files = n_filesA
        else:
            message = 'Number of UFS 2D/3D history files is inconsistent'
            logging.critical(message)
            raise Exception(message)
        
        center_ii_index = neighbors[0,center,center]
        center_jj_index = neighbors[1,center,center]
        for count, filename in enumerate(sfc_filenames, start=1):
            nc_file = Dataset('{0}/{1}'.format(dir,filename))
            nc_file.set_always_mask(False)
            
            for var2d in vars2d:
                data = nc_file[var2d["name"]][0,center_ii_index,center_jj_index]
                var2d["values"].append(data)
                var2d["units"]     = nc_file[var2d["name"]].getncattr(name="units")
                var2d["long_name"] = nc_file[var2d["name"]].getncattr(name="long_name")
            
            for phystend in phystends:
                if phystend["name"] in nc_file.variables.keys():
                    data = nc_file[phystend["name"]][0,::-1,center_ii_index,center_jj_index]
                    phystend["values"].append(data[:])
                    phystend["units"]     = nc_file[phystend["name"]].getncattr(name="units")
                    phystend["long_name"] = nc_file[phystend["name"]].getncattr(name="long_name")
                    phystend["missing"] = False
                else:
                    phystend["missing"] = True

            nc_file.close()
        
        # Convert to numpy arrays
        for var2d in vars2d:
            var2d["values"] = np.asarray(var2d["values"])
        for phystend in phystends:
            if (not phystend["missing"]):
                phystend["values"] = np.asarray(phystend["values"])
    sec_in_hr = 3600.
    time[0]   = 0.
    forcing = {
        "time":     time*sec_in_hr,
        "wa":       smoothed_w,
        "ps_forc":  pressfc[:,center,center],
        "pa_forc":  pres,
        "tnta_adv": t_advec_T[:,:],
        "tnqv_adv": t_advec_qv[:,:],
        "tnua_adv": t_advec_u[:,:],
        "tnva_adv": t_advec_v[:,:],
        "ug"      : u_g[:,:],
        "vg"      : v_g[:,:],
        "ua_nud"  : u_wind[:,:,center,center],
        "va_nud"  : v_wind[:,:,center,center]
    }
    
    if (save_comp_data):
        comp_data = {
            "time": time*sec_in_hr,
            "pa"  : pres,
            "ta"  : temp  [:,:,center,center],
            "qv"  : qv    [:,:,center,center],
            "ua"  : u_wind[:,:,center,center],
            "va"  : v_wind[:,:,center,center],
            "vars2d":vars2d,
            "phystends":phystends}
    else:
        comp_data = {}
    
    
    return (forcing, comp_data)

########################################################################################
#
########################################################################################
def get_UFS_forcing_data(nlevs, state_IC, location, use_nearest, forcing_dir, grid_dir,
                         tile, i, j, lam, save_comp_data, exact_mode):
    """Get the horizontal and vertical advective tendencies for the given tile and indices"""

    # Determine UFS history file format (tiled/quilted)
    if lam:
        atm_ftag = '*dynf*.nc'
        sfc_ftag = '*phyf*.nc'
    else:
        atm_ftag = '*atmf*.nc'
        sfc_ftag = '*sfcf*.nc'
    # end if
    
    # Get list of UFS history files with 3D ATMospheric state variables.
    atm_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, atm_ftag):
          atm_filenames.append(f_name)
        # end if
    # end for
    if not atm_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.                \
                  format(atm_ftag,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    # end if
    atm_filenames = sorted(atm_filenames)
    n_filesA      = len(atm_filenames)

    # Get list of UFS history files with 2D fields.
    sfc_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, sfc_ftag):
          sfc_filenames.append(f_name)
        # end if
    # end fo
    if not sfc_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.                \
                  format(sfc_ftag,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    # end if
    sfc_filenames = sorted(sfc_filenames)
    n_filesS      = len(sfc_filenames)

    if (n_filesS == n_filesA):
        n_files = n_filesA
    else:
        message = 'Number of UFS 2D/3D history files is inconsistent'
        logging.critical(message)
        raise Exception(message)
    # end if
    
    # Physical constants (used by FV3 remapping functions)
    kord_tm = -9
    kord_tr = 9
    t_min = 184.0
    q_min = 0.0
    sec_in_hr = 3600.
    
    ####################################################################################
    #
    # Read in 3D UFS history files
    #
    ####################################################################################

    # Find nearest point on UFS history file (quilted) grid.
    if use_nearest:
        (tile_jj, tile_ii, point_lon, point_lat, dist_min, theta) = find_loc_indices(location, forcing_dir, -999, lam)
        print('The closest point has indices [{0},{1}]'.format(tile_ii,tile_jj))
        print('This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat))
        print('This grid cell is approximately {0} km away from the desired location of {1} {2}'.format(dist_min/1.0E3,location[0],location[1]))
    # end if
    
    # Initialize
    ps      = []
    p_lev   = []
    p_lay   = []
    t_lay   = []
    qv_lay  = []
    u_lay   = []
    v_lay   = []
    time_hr = []

    # Get grid from UFS IC data
    (ic_grid_lon, ic_grid_lat) = get_initial_lon_lat_grid(grid_dir, tile, lam)

    # Read in 3D UFS history files
    for count, filename in enumerate(atm_filenames, start=1):
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)
        
        # Check if UFS history file grid is different than UFS initial-condition grid.
        if not use_nearest:
            try:
                data_grid_lon = nc_file['lon'][:,:]
                data_grid_lat = nc_file['lat'][:,:]
            except:
                data_grid_lon = nc_file['grid_xt'][:,:]
                data_grid_lat = nc_file['grid_yt'][:,:]
            # end try
            equal_grids = False
            if (ic_grid_lon.shape == data_grid_lon.shape and                               \
                ic_grid_lat.shape == ic_grid_lat.shape):
                if (np.equal(ic_grid_lon,data_grid_lon).all() and                          \
                    np.equal(ic_grid_lat,data_grid_lat).all()):
                    equal_grids = True
                # end if
            # end if

            # If necessary, remap history file (data_grid) to IC file (ic_grid).
            if (not equal_grids):
                print('Regridding {} onto native grid: regridding progress = {}%'.         \
                      format(filename, 100.0*count/(2*n_files)))

                grid_in   = {'lon': data_grid_lon, 'lat': data_grid_lat}
                grid_out  = {'lon': np.reshape(ic_grid_lon[j,i],(-1,1)), 'lat':            \
                             np.reshape(ic_grid_lat[j,i],(-1,1))}
                regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
                ps_data   = regridder(nc_file['pressfc'][0,:,:])
                t_data    = regridder(nc_file['tmp'][0,::-1,:,:])
                qv_data   = regridder(nc_file['spfh'][0,::-1,:,:])
                u_data    = regridder(nc_file['ugrd'][0,::-1,:,:])
                v_data    = regridder(nc_file['vgrd'][0,::-1,:,:])
                i_get     = 0
                j_get     = 0
            # Same grids for history file (data_grid) to IC file (ic_grid).
            else:
                ps_data   = nc_file['pressfc'][0,:,:]
                t_data    = nc_file['tmp'][0,::-1,:,:]
                qv_data   = nc_file['spfh'][0,::-1,:,:]
                u_data    = nc_file['ugrd'][0,::-1,:,:]
                v_data    = nc_file['vgrd'][0,::-1,:,:]
                i_get     = i
                j_get     = j
            # end if
        else:
            print('Using nearest UFS point {} progress = {}%'.format(filename, 100.0*count/(2*n_files)))
            ps_data = nc_file['pressfc'][0,:,:]
            t_data  = nc_file['tmp'][0,::-1,:,:]
            qv_data = nc_file['spfh'][0,::-1,:,:]
            u_data  = nc_file['ugrd'][0,::-1,:,:]
            v_data  = nc_file['vgrd'][0,::-1,:,:]
            j_get   = tile_jj
            i_get   = tile_ii
        # end if (use-nearest)

        # Store surface pressure.
        ps.append(ps_data[j_get,i_get])
        
        # Compute and store vertical grid information.
        ak    = getattr(nc_file, "ak")[::-1]
        bk    = getattr(nc_file, "bk")[::-1]
        nlevs = len(nc_file.dimensions['pfull'])

        # Compute and store pressure at layer-interfaces (nlev+1).
        p_interface = np.zeros(nlevs+1)
        for k in range(nlevs+1):
            p_interface[k]=ak[k]+ps[-1]*bk[k]
        # end for
        p_lev.append(p_interface)

        # Compute and store pressure at layer-centers (nlev).
        p_layer = np.zeros(nlevs)
        for k in range(nlevs):
            p_layer[k] = ((1.0/(rocp+1.0))*(p_interface[k]**(rocp+1.0) -               \
                            p_interface[k+1]**(rocp+1.0))/(p_interface[k] -            \
                            p_interface[k+1]))**(1.0/rocp)
        # end for
        p_lay.append(p_layer)

        # Store state variables.
        t_lay.append(t_data[:,j_get,i_get])
        qv_lay.append(qv_data[:,j_get,i_get])
        u_lay.append(u_data[:,j_get,i_get])
        v_lay.append(v_data[:,j_get,i_get])
        time_hr.append(nc_file['time'][0])

        # Close file
        nc_file.close()
    # end for

    # Convert from python list to numpy array
    ps      = np.asarray(ps)
    p_lev   = np.asarray(p_lev)
    p_lay   = np.asarray(p_lay)
    t_lay   = np.asarray(t_lay)
    qv_lay  = np.asarray(qv_lay)
    u_lay   = np.asarray(u_lay)
    v_lay   = np.asarray(v_lay)
    time_hr = np.asarray(time_hr)

    # Compute virtual temperature.
    tv_lay  = t_lay*(1.0 + zvir*qv_lay)

    ####################################################################################
    #
    # Read in 2D UFS history files
    #
    # Comments:
    # These files contain the 3D (diagnostic) dynamic tendencies needed to create "exact"
    # UFS forcings. These are stored in the python dictionary "dyntends".
    #
    # Additionally, a limited set of two-dimensionsal variables, python dictionary "vars2d",
    # are output to the "case comparison file".
    #
    ####################################################################################

    # Variables needed for "exact" UFS case generation.
    dyntends = [{"name":"dtend_temp_nophys", "diag_table":'"gfs_dyn", "dtend_temp_nophys", "dtend_temp_nophys", "fv3_history", "all", .false., "none", 2'},\
                {"name":"dtend_u_nophys",    "diag_table":'"gfs_dyn", "dtend_u_nophys",    "dtend_u_nophys",    "fv3_history", "all", .false., "none", 2'},\
                {"name":"dtend_v_nophys",    "diag_table":'"gfs_dyn", "dtend_v_nophys",    "dtend_v_nophys",    "fv3_history", "all", .false., "none", 2'},\
                {"name":"dtend_qv_nophys",   "diag_table":'"gfs_dyn", "dtend_qv_nophys",   "dtend_qv_nophys",   "fv3_history", "all", .false., "none", 2'}]
    for dyntend in dyntends: dyntend["values"] = []
    
    # Optional diagnostic to be copied into UFS comparison file
    phystends = [{"name":"dtend_temp_lw"},     {"name":"dtend_temp_sw"},  {"name":"dtend_temp_pbl"},   {"name":"dtend_temp_deepcnv"},\
                 {"name":"dtend_temp_shalcnv"},{"name":"dtend_temp_mp"},  {"name":"dtend_temp_orogwd"},{"name":"dtend_temp_phys"},   \
                 {"name":"dtend_u_pbl"},       {"name":"dtend_v_pbl"},    {"name":"dtend_u_orogwd"},   {"name":"dtend_v_orogwd"},    \
                 {"name":"dtend_u_deepcnv"},   {"name":"dtend_v_deepcnv"},{"name":"dtend_u_shalcnv"},  {"name":"dtend_v_shalcnv"},   \
                 {"name":"dtend_u_phys"},      {"name":"dtend_v_phys"},   {"name":"dtend_qv_pbl"},     {"name":"dtend_qv_deepcnv"},  \
                 {"name":"dtend_qv_shalcnv"},  {"name":"dtend_qv_mp"},    {"name":"dtend_qv_phys"},    {"name":"dtend_poop"}]
    for phystend in phystends: phystend["values"] = []

    # Variables to be added to "UFS comparison file"
    vars2d  =[{"name":"spfh2m"},    {"name":"tmp2m"},     {"name":"dswrf_ave"}, \
              {"name":"ulwrf_ave"}, {"name":"lhtfl_ave"}, {"name":"shtfl_ave"}, \
              {"name":"dswrf"},     {"name":"ulwrf"},     {"name":"lhtfl"},     \
              {"name":"shtfl"},     {"name":"pwat"},      {"name":"vgrd10m"},   \
              {"name":"ugrd10m"}]
    for var2d in vars2d: var2d["values"] = []

    # Read in 2D UFS history files
    for count, filename in enumerate(sfc_filenames, start=1):
        nc_file = Dataset('{0}/{1}'.format(forcing_dir,filename))
        nc_file.set_always_mask(False)
        
        # Check if UFS history file grid is different than UFS initial-condition grid.
        if not use_nearest:
            try:
                data_grid_lon = nc_file['lon'][:,:]
                data_grid_lat = nc_file['lat'][:,:]
            except:
                data_grid_lon = nc_file['grid_xt'][:,:]
                data_grid_lat = nc_file['grid_yt'][:,:]
            # end try
            equal_grids = False
            if (ic_grid_lon.shape == data_grid_lon.shape and                               \
                ic_grid_lat.shape == ic_grid_lat.shape):
                if (np.equal(ic_grid_lon,data_grid_lon).all() and                          \
                    np.equal(ic_grid_lat,data_grid_lat).all()):
                    equal_grids = True
                # end if
            # endif

            # If necessary, remap history file (data_grid) to IC file (ic_grid).
            if (not equal_grids):
                print('Regridding {} onto native grid: regridding progress = {}%'.         \
                      format(filename, 50+50.0*count/(n_files)))

                grid_in   = {'lon': data_grid_lon, 'lat': data_grid_lat}
                grid_out  = {'lon': np.reshape(ic_grid_lon[j,i],(-1,1)), 'lat':            \
                             np.reshape(ic_grid_lat[j,i],(-1,1))}
                regridder = xesmf.Regridder(grid_in, grid_out, 'bilinear')
                i_get     = 0
                j_get     = 0
            # Same grids for history file (data_grid) to IC file (ic_grid).
            else:
                i_get     = i
                j_get     = j
            # end if
        else:
            print('Using nearest UFS point {} progress = {}%'.format(filename, 50+50.0*count/(n_files)))
            j_get   = tile_jj
            i_get   = tile_ii
        # end if
        
        for var2d in vars2d:
            if not use_nearest:
                data = regridder(nc_file[var2d["name"]][0,:,:])
            else:
                data = nc_file[var2d["name"]][0,:,:]
            # end if
            var2d["values"].append(data[j_get,i_get])
            var2d["units"]     = nc_file[var2d["name"]].getncattr(name="units")
            var2d["long_name"] = nc_file[var2d["name"]].getncattr(name="long_name")
        # end for
        for dyntend in dyntends:
            if dyntend["name"] in nc_file.variables.keys():
                if not use_nearest:
                    data = regridder(nc_file[dyntend["name"]][0,::-1,:,:])
                else:
                    data = nc_file[dyntend["name"]][0,::-1,:,:]
                # end if (use-nearest)
                dyntend["values"].append(data[:,j_get,i_get])
                dyntend["units"]     = nc_file[dyntend["name"]].getncattr(name="units")
                dyntend["long_name"] = nc_file[dyntend["name"]].getncattr(name="long_name")
                dyntend["missing"] = False
            else:
                dyntend["missing"] = True
            # end if
        # end for
        for phystend in phystends:
            if phystend["name"] in nc_file.variables.keys():
                if not use_nearest:
                    data = regridder(nc_file[phystend["name"]][0,::-1,:,:])
                else:
                    data = nc_file[phystend["name"]][0,::-1,:,:]
                # end if (use-nearest)
                phystend["values"].append(data[:,j_get,i_get])
                phystend["units"]     = nc_file[phystend["name"]].getncattr(name="units")
                phystend["long_name"] = nc_file[phystend["name"]].getncattr(name="long_name")
                phystend["missing"] = False
            else:
                phystend["missing"] = True                
            # end if
        # end for
        nc_file.close()
    # end for

    # Handle the case that "exact-mode" was requested, but the UFS history files do not contains
    # the necessary fields.
    # In this scenario, create warning/message that informs the users how to rerun the UFS with
    # the correct configuration.
    # Revert to "approximate-mode"
    if exact_mode:
        missing_dtend_names = []
        missing_diag_tables = []
        for dyntend in dyntends:
            if dyntend["missing"]:
                missing_dtend_names.append(dyntend["name"])
                missing_diag_tables.append(dyntend["diag_table"])
                exact_mode = False
            # end if
        # end for
        if missing_dtend_names:
            print("#"*128)
            print("#"*59,"WARNING","#"*60)
            print("#"*128)
            print("Could not find the following fields required for ",'"exact-mode"',':')
            for missing_dtend_name in missing_dtend_names:
                print("  -",missing_dtend_name)
            # end for
            print("Using ",'"approximate mode"',' instead')
            print("")
            print("  *NOTE*")
            print("  You selected to create the UFS case using ",'"exact-mode"'," which requires additional (diagnostic)")
            print("  fields in the UFS history files.")
            print("")
            print("  These optional diagnostics tendencies can be output from any UFS forecast by:")
            print("   a) Modifying the physics namelist, gfs_phys_nml, so that:")
            print("      *) ldiag3d = .true. ")
            print("      *) qdiag3d = .true. ")
            print("   b) Addding the following lines to the UFS diag_table:")
            for missing_diag_table in missing_diag_tables:
                print('      *) ',missing_diag_table )
            # end for
            print("")
            print("#"*128)
            print("#"*128)
        # end if
    # end if (exact-mode)

    # Convert to numpy arrays
    for var2d in vars2d:
        var2d["values"] = np.asarray(var2d["values"])
    # end for
    if exact_mode:
        for dyntend in dyntends:
            dyntend["values"] = np.asarray(dyntend["values"])
        # end for
        dtend_temp_nophys = dyntends[0]["values"]
        dtend_u_nophys    = dyntends[1]["values"]
        dtend_v_nophys    = dyntends[2]["values"]
        dtend_qv_nophys   = dyntends[3]["values"]
    # end if (exact-mode)
    for phystend in phystends:
        if (not phystend["missing"]):
            phystend["values"] = np.asarray(phystend["values"])
        # end if
    # end for

    ####################################################################################
    #
    # The "total" advection, where "total" = "advective + remapping", can be computed 
    # directly by taking the difference across time. 
    # The advective part is the piece needed to force the SCM, which we can "back-out"
    # by computing the "remapping" component, and removing it from the total.
    # 
    # At initial-time (ICs)
    # - Remap IC state to vertical grid of first UFS history file.
    #
    # At subsequent time(s)
    # - Remap UFS state at current time-step to vertical grid of subsequent time-step 
    #   (used for differencing)
    #
    # This "REGRID"ed state is used to remove the impact vertical remapping on the 
    # change of state.
    #
    # *NOTE* The first UFS history file is AFTER the first physics timestep. Subsequent
    #        history files are controlled by "output_fh" in the UFS.
    #
    ####################################################################################

    # Initialize
    #nlevs      = len(p_lay[0,:])
    dummy      = np.zeros(1)
    from_p     = np.zeros([1,nlevs+1])
    to_p       = np.zeros([1,nlevs+1])
    log_from_p = np.zeros([1,nlevs+1])
    log_to_p   = np.zeros([1,nlevs+1])
    dp2        = np.zeros([1,nlevs])
    tv_rev     = np.zeros([1,nlevs])
    qv_rev     = np.zeros([1,nlevs])
    u_rev      = np.zeros([1,nlevs])
    v_rev      = np.zeros([1,nlevs])
    tv_layr    = np.zeros([n_files+1,nlevs])
    qv_layr    = np.zeros([n_files+1,nlevs])
    u_layr     = np.zeros([n_files+1,nlevs])
    v_layr     = np.zeros([n_files+1,nlevs])
    p_layr     = np.zeros([n_files+1,nlevs])
    p_levr     = np.zeros([n_files+1,nlevs+1])

    #
    # First timestep...
    # (Same for Exact-mode or Approximate mode)
    #

    # Interpolation range
    from_p[0,:]     = state_IC["pa_i"][::-1]
    to_p[0,:]       = p_lev[0,::-1]
    log_from_p[0,:] = np.log(state_IC["pa_i"][::-1])
    log_to_p[0,:]   = np.log(p_lev[0,::-1])

    # IC Virtual Temperature on vertical-grid of first UFS history file.
    tv_init      = state_IC["ta"]*(1.0 + zvir*state_IC["qv"])
    tv_init_rev  = tv_init[::-1]
    tv_rev_new   = fv3_remap.map_scalar(nlevs, log_from_p, tv_init_rev[np.newaxis, :], \
                                        dummy, nlevs, log_to_p, 0, 0, 1,               \
                                        np.abs(kord_tm), t_min)
    
    # IC Specific humidity on vertical-grid of first UFS history file.
    qv_init_rev  = state_IC["qv"][::-1]
    for k in range(0,nlevs): dp2[0,k] = from_p[0,k+1] - from_p[0,k]
    qv_rev_new   = fv3_remap.map1_q2(nlevs, from_p, qv_init_rev[np.newaxis, :],        \
                                     nlevs, to_p, dp2, 0, 0, 0, kord_tr, q_min)

    # IC Zonal wind on vertical-grid of first UFS history file.
    u_init_rev  = state_IC["ua"][::-1]
    u_rev_new   = fv3_remap.map1_ppm(nlevs, from_p, u_init_rev[np.newaxis, :], 0.0,    \
                                     nlevs, to_p, 0, 0, -1, kord_tm )

    # IC Meridional wind on vertical-grid of first UFS history file.
    v_init_rev  = state_IC["va"][::-1]
    v_rev_new   = fv3_remap.map1_ppm(nlevs, from_p, v_init_rev[np.newaxis, :], 0.0,    \
                                     nlevs, to_p, 0, 0, -1, kord_tm )
    
    # Store
    p_layr[0,:]  = p_lay[0,:]
    p_levr[0,:]  = p_lev[0,:]
    v_layr[0,:]  = v_rev_new[0,::-1]
    u_layr[0,:]  = u_rev_new[0,::-1]
    tv_layr[0,:] = tv_rev_new[0,::-1]
    qv_layr[0,:] = qv_rev_new[0,::-1]
    
    # Subsequent timestep(s). (exact-mode only)
    if exact_mode:
        for t in range(n_files-1):
            # Interpolation range
            from_p[0,:]     = p_lev[t,::-1]
            to_p[0,:]       = p_lev[t+1,::-1]
            log_from_p[0,:] = np.log(p_lev[t,::-1])
            log_to_p[0,:]   = np.log(p_lev[t+1,::-1])
            
            # Virtual Temperature @ time > 0
            tv_rev[0,:] = tv_lay[t,::-1]
            tv_rev_new  = fv3_remap.map_scalar(nlevs, log_from_p, tv_rev, dummy, nlevs,    \
                                               log_to_p, 0, 0, 1, np.abs(kord_tm), t_min)
            # Specific humidity @ time > 0
            qv_rev[0,:] = qv_lay[t,::-1]
            for k in range(0,nlevs): dp2[0,k] = to_p[0,k+1] - to_p[0,k]
            qv_rev_new = fv3_remap.map1_q2(nlevs, from_p, qv_rev, nlevs, to_p, dp2,        \
                                           0, 0, 0, kord_tr, q_min)
            # Zonal wind  @ time > 0
            u_rev[0,:] = u_lay[t,::-1]
            u_rev_new  = fv3_remap.map1_ppm(nlevs, from_p, u_rev, 0.0, nlevs, to_p,        \
                                            0, 0, -1, kord_tm )
            # Meridional wind @ time > 0
            v_rev[0,:] = v_lay[t,::-1]
            v_rev_new  = fv3_remap.map1_ppm(nlevs, from_p, v_rev, 0.0, nlevs, to_p,        \
                                            0, 0, -1, kord_tm )
            # Store
            p_layr[t+1,:]  = p_lay[t+1,:]
            p_levr[t+1,:]  = p_lev[t+1,:]
            tv_layr[t+1,:] = tv_rev_new[0,::-1]
            qv_layr[t+1,:] = qv_rev_new[0,::-1]
            u_layr[t+1,:]  = u_rev_new[0,::-1]
            v_layr[t+1,:]  = v_rev_new[0,::-1]
        # end for
    
        #
        p_layr[t+2,:]  = p_layr[t+1,:]
        p_levr[t+2,:]  = p_levr[t+1,:]
        tv_layr[t+2,:] = tv_layr[t+1,:]
        qv_layr[t+2,:] = qv_layr[t+1,:]
        u_layr[t+2,:]  = u_layr[t+1,:]
        v_layr[t+2,:]  = v_layr[t+1,:]
    # end if (exact-mode)
    
    # Temperature
    t_layr = tv_layr/(1.0 + zvir*qv_layr)

    # Save regridded initial state.
    stateInit = {"p_lay":  p_lay[0,:],  \
                 "t_lay":  t_lay[0,:],  \
                 "qv_lay": qv_lay[0,:], \
                 "u_lay":  u_lay[0,:],  \
                 "v_lay":  v_lay[0,:]}
    
    ####################################################################################
    #
    # Compute forcing.
    #
    ####################################################################################
    # Initialize tendencies.
    dtdt_adv   = np.zeros([n_files+1,nlevs])
    dqvdt_adv  = np.zeros([n_files+1,nlevs])
    dudt_adv   = np.zeros([n_files+1,nlevs])
    dvdt_adv   = np.zeros([n_files+1,nlevs])
    pres_adv   = np.zeros([n_files+1,nlevs])
    pres_i_adv = np.zeros([n_files+1,nlevs+1])
    dtdtr      = np.zeros([1,nlevs])
    dqdtr      = np.zeros([1,nlevs])
    dudtr      = np.zeros([1,nlevs])
    dvdtr      = np.zeros([1,nlevs])

    # First timestep.
    pres_adv[0,:]   = state_IC["pa"]
    pres_i_adv[0,:] = state_IC["pa_i"]
    dt              = 3600.0*(time_hr[0])
    if exact_mode:
        dqdtr[0,:]     = (qv_layr[0,:] - state_IC["qv"][:])/dt
        dqvdt_adv[0,:] = dtend_qv_nophys[0,:] - dqdtr[0,:]
        dudtr[0,:]     = (u_layr[0,:] - state_IC["ua"][:])/dt
        dudt_adv[0,:]  = (dtend_u_nophys[0,:] - dudtr[0,:])
        dvdtr[0,:]     = (v_layr[0,:] - state_IC["va"][:])/dt
        dvdt_adv[0,:]  = (dtend_v_nophys[0,:] - dvdtr[0,:])
        dtdtr[0,:]     = (t_layr[0,:] - state_IC["ta"][:])/dt
        dtdt_adv[0,:]  = (dtend_temp_nophys[0,:] - dtdtr[0,:])
    else:
        dqdtt          = qv_lay[0,:] - state_IC["qv"][:]
        dqvdt_adv[0,:] = dqdtt/dt
        dudtt          = u_lay[0,:]  - state_IC["ua"][:]
        dudt_adv[0,:]  = dudtt/dt
        dvdtt          = v_lay[0,:]  - state_IC["va"][:]
        dvdt_adv[0,:]  = dvdtt/dt
        dtdtt          = t_lay[0,:]  - state_IC["ta"][:]
        dtdt_adv[0,:]  = dtdtt/dt        

    # Subsequent timestep(s).
    dtdtr = np.zeros([n_files-1,nlevs])
    dqdtr = np.zeros([n_files-1,nlevs])
    dudtr = np.zeros([n_files-1,nlevs])
    dvdtr = np.zeros([n_files-1,nlevs])
    dtdtt = np.zeros([n_files-1,nlevs])
    dqdtt = np.zeros([n_files-1,nlevs])
    dudtt = np.zeros([n_files-1,nlevs])
    dvdtt = np.zeros([n_files-1,nlevs])
    for t in range(n_files-1):
        dt = (time_hr[t+1] - time_hr[t])*3600.
        pres_adv[t+1,:]   = p_lay[t]
        pres_i_adv[t+1,:] = p_lev[t]
        if exact_mode:
            dqdtr[t,:]       = (qv_layr[t+1,:] - qv_lay[t,:])/dt
            dqvdt_adv[t+1,:] = dtend_qv_nophys[t+1,:] - dqdtr[t,:]
            dtdtr[t,:]       = (t_layr[t+1,:] - t_lay[t,:])/dt
            dtdt_adv[t+1,:]  = dtend_temp_nophys[t+1,:] - dtdtr[t,:]
            dudtr[t,:]       = (u_layr[t+1,:] - u_lay[t,:])/dt
            dudt_adv[t+1,:]  = dtend_u_nophys[t+1,:] - dudtr[t,:]
            dvdtr[t,:]       = (v_layr[t+1,:] - v_lay[t,:])/dt
            dvdt_adv[t+1,:]  = dtend_v_nophys[t+1,:] - dvdtr[t,:]
        else:
            dqdtt            = qv_lay[t+1,:] - qv_lay[t,:]
            dqvdt_adv[t+1,:] = (dqdtt - dqdtr[t,:])/dt
            dudtt            = u_lay[t+1,:] - u_lay[t,:]
            dudt_adv[t+1,:]  = (dudtt - dudtr[t,:])/dt
            dvdtt            = v_lay[t+1,:] - v_lay[t,:]
            dvdt_adv[t+1,:]  = (dvdtt - dvdtr[t,:])/dt
            dtdtt            = t_lay[t+1,:] - t_lay[t,:]
            dtdt_adv[t+1,:]  = (dtdtt - dtdtr[t,:])/dt
        # end if (exact-mode)
    # end for (ufs history files)

    # Last timestep.
    dtdt_adv[t+2,:]   = dtdt_adv[t+1,:]
    dqvdt_adv[t+2,:]  = dqvdt_adv[t+1,:]
    dudt_adv[t+2,:]   = dudt_adv[t+1,:]
    dvdt_adv[t+2,:]   = dvdt_adv[t+1,:]
    pres_adv[t+2,:]   = pres_adv[t+1,:]
    pres_i_adv[t+2,:] = pres_i_adv[t+1,:]

    ####################################################################################
    #
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
    #
    ####################################################################################
    #time_method = 'constant_simple' #this is not implemented in the SCM code yet
    time_method = 'constant_interp'
    #time_method = 'gradient' #this produced wonky results in the SCM; avoid until investigated more
    
    if (time_method == 'constant_simple'):
        print('Forcing should not be interpolated in time. Rather, forcing should held constant at their current values until the next forcing interval is reached.')
        ntimes = n_files
        time = np.zeros(ntimes)
        
        p_s = np.zeros((ntimes),dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)        
        
        p_s[0] = ps[0]
        pressure_forc[:,0] = pres_adv[0,:]
        tot_advec_T[:,0] = dtdt_adv[0,:]
        tot_advec_qv[:,0] = dqvdt_adv[0,:]
        tot_advec_u[:,0] = dudt_adv[0,:]
        tot_advec_v[:,0] = dvdt_adv[0,:]
        
        for t in range(1,n_files):
            time[t] = sec_in_hr*time_hr[t-1]
            p_s[t] = ps[t]
            pressure_forc[:,t] = pres_adv[t,:]
            tot_advec_T[:,t] = dtdt_adv[t,:]
            tot_advec_qv[:,t] = dqvdt_adv[t,:]
            tot_advec_u[:,t] = dudt_adv[t,:]
            tot_advec_v[:,t] = dvdt_adv[t,:]
        # end for
    elif (time_method == 'constant_interp'):
        print('Forcing can be interpolated in time, but the time values are chosen such that forcing will effectively be held consant during a diagnostic time interval.')
        ntimes = 2*n_files
        
        time_setback = 1.0 #s
        
        time = np.zeros(ntimes)
        p_s = np.zeros((ntimes),dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)
            
        time[0] = 0.0
        time[1] = sec_in_hr*time_hr[0] - time_setback #forcing period should extend from beginning of diagnostic period to right BEFORE the next one
        p_s[0] = ps[0]
        p_s[1] = p_s[0]
        pressure_forc[:,0] = pres_adv[0,:]
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
            time[2*t] = sec_in_hr*time_hr[t-1]
            time[2*t+1] = sec_in_hr*time_hr[t] - time_setback
            p_s[2*t] = ps[t]
            p_s[2*t+1] = p_s[2*t]
            pressure_forc[:,2*t] = pres_adv[t,:]
            pressure_forc[:,2*t+1] = pressure_forc[:,2*t]
            tot_advec_T[:,2*t] = dtdt_adv[t,:]
            tot_advec_T[:,2*t+1] = tot_advec_T[:,2*t]
            tot_advec_qv[:,2*t] = dqvdt_adv[t,:]
            tot_advec_qv[:,2*t+1] = tot_advec_qv[:,2*t]
            tot_advec_u[:,2*t] = dudt_adv[t,:]
            tot_advec_u[:,2*t+1] = tot_advec_u[:,2*t]
            tot_advec_v[:,2*t] = dvdt_adv[t,:]
            tot_advec_v[:,2*t+1] = tot_advec_v[:,2*t]
        # end for
    elif (time_method == 'gradient'): #this produced wonky results in the SCM; avoid until investigated more
        print('Forcing can be interpolated in time since the forcing terms are assumed to follow a constant time-gradient.')
        
        ntimes = 2*n_files + 1
        time = np.zeros(ntimes)
        p_s = np.zeros((ntimes),dtype=float)
        pressure_forc = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_T = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_qv = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_u = np.zeros((nlevs,ntimes),dtype=float)
        tot_advec_v = np.zeros((nlevs,ntimes),dtype=float)
        
        p_s[0] = state_IC['ps']
        pressure_forc[:,0] = state_IC['pa']
        tot_advec_T[:,0] = 0.0
        tot_advec_qv[:,0] = 0.0
        tot_advec_u[:,0] = 0.0
        tot_advec_v[:,0] = 0.0
        
        for t in range(n_files):
            time[2*t + 1] = time[2*t] + 0.5*(sec_in_hr*time_hr[t] - time[2*t])
            time[2*t + 2] = sec_in_hr*time_hr[t]
            
            p_s[2*t+1] = ps[t]
            pressure_forc[:,2*t+1] = pres_adv[t,:]
            tot_advec_T[:,2*t+1] = dtdt_adv[t,:]
            tot_advec_qv[:,2*t+1] = dqvdt_adv[t,:]
            tot_advec_u[:,2*t+1] = dudt_adv[t,:]
            tot_advec_v[:,2*t+1] = dvdt_adv[t,:]
            
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
    # end if
    
    #
    w_ls         = np.zeros((nlevs,ntimes),dtype=float)
    omega        = np.zeros((nlevs,ntimes),dtype=float)
    rad_heating  = np.zeros((nlevs,ntimes),dtype=float)

    forcing = {
        "time":     time,
        "wa":       w_ls.swapaxes(0,1),
        "wap":      omega.swapaxes(0,1),
        "tnta_rad": rad_heating.swapaxes(0,1),
        "ps_forc":  np.ones(ntimes)*ps[0],
        "pa_forc":  pressure_forc.swapaxes(0,1),
        "tnta_adv": tot_advec_T.swapaxes(0,1),
        "tnqv_adv": tot_advec_qv.swapaxes(0,1),
        "tnua_adv": tot_advec_u.swapaxes(0,1),
        "tnva_adv": tot_advec_v.swapaxes(0,1)
    }
    
    if (save_comp_data):
        time_hr[0] = 0.0
        comp_data = {
            "time" : time_hr*sec_in_hr,
            "pa"  : p_lay,#[:,::-1],
            "ta"  : t_lay,#[:,::-1],
            "qv"  : qv_lay,#[:,::-1],
            "ua"  : u_lay,#[:,::-1],
            "va"  : v_lay,#[:,::-1],
            "vars2d":vars2d,
            "phystends":phystends}
    else:
        comp_data = {}

    return (forcing, comp_data, stateInit)

########################################################################################
#
########################################################################################
def write_SCM_case_file(state, surface, oro, forcing, init, case, date, forcing_method, vertical_method, geos_wind_forcing, wind_nudge):
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
    if (vertical_method == 2):
        nc_file.forc_wa       = forcing_on
    else:
        nc_file.forc_wa       = forcing_off
    nc_file.forc_wap          = forcing_off
    if (geos_wind_forcing):
        nc_file.forc_geo      = forcing_on
    else:
        nc_file.forc_geo      = forcing_off
    if (wind_nudge):
        nc_file.nudging_ua    = forcing_on*const_nudging_time
        nc_file.nudging_va    = forcing_on*const_nudging_time
    else:
        nc_file.nudging_ua    = forcing_off
        nc_file.nudging_va    = forcing_off
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
    vegcat_dim = nc_file.createDimension('nvegcat', len(oro["vegtype_frac"]))
    soilcat_dim =nc_file.createDimension('nsoilcat', len(oro["soiltype_frac"]))
    
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
                {"name": "theta",        "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_potential_temperature"}, \
                {"name": "thetal",       "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_liquid_potential_temperature"}, \
                {"name": "rv",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "humidity_mixing_ratio"}, \
                {"name": "rl",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "cloud_liquid_water_mixing_ratio"}, \
                {"name": "ri",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "cloud_ice_water_mixing_ratio"}, \
                {"name": "rt",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "water_mixing_ratio"}, \
                {"name": "ql",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_cloud_liquid_water_in_air"}, \
                {"name": "qi",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_cloud_ice_water_in_air", "default_value": 0.0}, \
                {"name": "qt",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "mass_fraction_of_water_in_air"}, \
                {"name": "hur",          "type":wp, "dimd": ('t0',   'lev'),    "units": "%",             "desc": "relative_humidity"}, \
                {"name": "tke",          "type":wp, "dimd": ('t0',   'lev'),    "units": "m2 s-2",        "desc": "specific_turbulen_kinetic_energy", "default_value": 0.0}, \
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
    if (forcing_method == 1 or forcing_method == 3):
        var_dict.insert(3, {"name": "ta",           "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_temperature","default_value": init["t_lay"][:], "override": True})
        var_dict.insert(10,{"name": "qv",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "specific_humidity","default_value": init["qv_lay"][:], "override": True})
        var_dict.insert(16,{"name": "ua",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "eastward_wind","default_value": init["u_lay"][:], "override": True})
        var_dict.insert(17,{"name": "va",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "northward_wind","default_value": init["v_lay"][:], "override": True})
    else:
        var_dict.insert(3, {"name": "ta",           "type":wp, "dimd": ('t0',   'lev'),    "units": "K",             "desc": "air_temperature"})
        var_dict.insert(10,{"name": "qv",           "type":wp, "dimd": ('t0',   'lev'),    "units": "kg kg-1",       "desc": "specific_humidity"})
        var_dict.insert(16,{"name": "ua",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "eastward_wind"})
        var_dict.insert(17,{"name": "va",           "type":wp, "dimd": ('t0',   'lev'),    "units": "m s-1",         "desc": "northward_wind"})
        
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
                {"name": "lakedepth",    "type":wp, "dimd": ('t0'),             "units": "none",    "desc": "lake depth", "default_value":0},
                {"name": "vegtype_frac", "type":wp, "dimd": ('t0', 'nvegcat'),  "units": "none",    "desc": "fraction of horizontal grid area occupied by given vegetation category"},
                {"name": "soiltype_frac","type":wp, "dimd": ('t0', 'nsoilcat'), "units": "none",    "desc": "fraction of horizontal grid area occupied by given soil category"}]
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
def write_comparison_file(comp_data, case_name, date, surface):
    """Write UFS history file data to netCDF file for comparison"""

    wp = np.float64
    wi = np.int32

    nlevs = comp_data["pa"].shape[1]
    ntime = comp_data["pa"].shape[0]

    start_date = datetime(date["year"],date["month"],date["day"],date["hour"],date["minute"],date["second"])
    start_date_string = start_date.strftime("%Y%m%d%H%M%S")

    loc_string = str(round(surface["lon"],2)) + "E" + str(round(surface["lat"],2)) + "N"
    case_string = 'UFS_' + start_date_string + '_' + loc_string

    com = 'mkdir -p ' + COMPARISON_DATA_DIR
    print(com)
    os.system(com)
    nc_file = Dataset(os.path.join(COMPARISON_DATA_DIR, case_name + '_comp_data.nc'), 'w', format='NETCDF3_CLASSIC')
    nc_file.case = case_string
    nc_file.title = 'UFS history file data for ' + case_string
    nc_file.reference = ''
    nc_file.author = 'Grant J. Firl'
    nc_file.version = 'Created on ' + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
    nc_file.script = os.path.basename(__file__)
    nc_file.startDate = start_date_string

    # Dimensions
    lev_dim = nc_file.createDimension('lev', size=nlevs)
    time_dim = nc_file.createDimension('time', size=ntime)
    #time_ufs_history_dim = nc_file.createDimension('time_ufs_history', size=ntime-1) #GJF why are we ignoring the first history file?
    time_ufs_history_dim = nc_file.createDimension('time_ufs_history', size=ntime)

    # Varaibles
    time_var = nc_file.createVariable('time', wp, ('time',))
    time_var.units = 'second'
    time_var.long_name = 'history file time'
    time_var[:] = comp_data['time']

    time2_var = nc_file.createVariable('time_ufs_history', wp, ('time_ufs_history',))
    time2_var.units = 'second'
    time2_var.long_name = 'UFS history file time'
    time2_var[:] = comp_data['time']

    lev_var = nc_file.createVariable('levs', wp, ('time','lev',))
    lev_var.units = 'Pa'
    lev_var.long_name = 'pressure'
    lev_var[:,:] = comp_data["pa"]

    temperature_var = nc_file.createVariable('temp', wp, ('time', 'lev',))
    temperature_var.units = 'K'
    temperature_var.long_name = 'Temperature'
    temperature_var[:,:] = comp_data["ta"]

    qv_var = nc_file.createVariable('qv', wp, ('time', 'lev',))
    qv_var.units = 'kg kg-1'
    qv_var.long_name = 'specific humidity'
    qv_var[:,:] = comp_data["qv"]

    u_var = nc_file.createVariable('u', wp, ('time', 'lev',))
    u_var.units = 'm s-1'
    u_var.long_name = 'zonal wind'
    u_var[:,:] = comp_data["ua"]

    v_var = nc_file.createVariable('v', wp, ('time', 'lev',))
    v_var.units = 'm s-1'
    v_var.long_name = 'meridional wind'
    v_var[:,:] = comp_data["va"]

    for var2d in comp_data["vars2d"]:
        tempVar           = nc_file.createVariable(var2d["name"], wp, ('time_ufs_history'))
        tempVar.units     = var2d["units"]
        tempVar.long_name = var2d["long_name"]
        tempVar[:]        = var2d["values"]

    for phystend in comp_data["phystends"]:
        if (not phystend["missing"]):
            tempVar           = nc_file.createVariable(phystend["name"], wp, ('time', 'lev',))
            tempVar.units     = phystend["units"]
            tempVar.long_name = phystend["long_name"]
            tempVar[:]        = phystend["values"]
            # end if
        # end for
    # end if
    nc_file.close()

    return

########################################################################################
#
########################################################################################
def find_date(forcing_dir, lam):

    if lam:
        atm_ftag = '*dynf*.nc'
    else:    
        atm_ftag = '*atmf*.nc'
    
    atm_filenames = []
    for f_name in os.listdir(forcing_dir):
       if fnmatch.fnmatch(f_name, atm_ftag):
          atm_filenames.append(f_name)
    if not atm_filenames:
        message = 'No filenames matching the pattern {0} found in {1}'.format(atm_ftag,forcing_dir)
        logging.critical(message)
        raise Exception(message)
    atm_filenames = sorted(atm_filenames)
    
    nc_file = Dataset('{0}/{1}'.format(forcing_dir,atm_filenames[0]))
    
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

########################################################################################
#
########################################################################################    
def main():
    setup_logging()
    
    #read in arguments
    (location, indices, date, in_dir, grid_dir, forcing_dir, tile, area, case_name, 
     old_chgres, lam, save_comp, use_nearest, forcing_method, vertical_method, geos_wind_forcing, wind_nudge) = parse_arguments()
    
    #find indices corresponding to both UFS history files and initial condition (IC) files
    (hist_i, hist_j, hist_lon, hist_lat, hist_dist_min, angle_to_hist_point, neighbors, dx, dy) = find_loc_indices_UFS_history(location, forcing_dir, lam)
    
    (IC_i, IC_j, tile, IC_lon, IC_lat, IC_dist_min, angle_to_IC_point) = find_loc_indices_UFS_IC(location, grid_dir, lam, tile, indices)
        
    #compare the locations of the found history file point and UFS IC point
    compare_hist_and_IC_points(location, hist_lon, hist_lat, IC_lon, IC_lat, hist_dist_min, angle_to_hist_point, IC_dist_min, angle_to_IC_point)
    
    #read in surface data for the intial conditions at the IC point
    surface_data = get_UFS_surface_data(in_dir, tile, IC_i, IC_j, old_chgres, lam)
    
    #when using the advective forcing methods, initial conditions are taken from the history files; check that the history file surface type is the same as the UFS IC point from which data is read
    check_IC_hist_surface_compatibility(forcing_dir, hist_i, hist_j, surface_data, lam, old_chgres, tile)
    
    #read in orographic data for the initial conditions at the IC point
    oro_data = get_UFS_oro_data(in_dir, tile, IC_i, IC_j, lam, old_chgres)
    
    #read in the initial condition profiles
    
    if (forcing_method == 1 or forcing_method == 3):
        #read initial condition profiles from UFS IC files
        (state_data, error_msg) = get_IC_data_from_UFS_ICs(in_dir, grid_dir, tile, IC_i,\
                                                                  IC_j, old_chgres, lam)
    else:
        #read initial condition profiles from UFS history files (not true "initial" conditions, but after first UFS timestep)
        state_data = get_IC_data_from_UFS_history(forcing_dir, hist_i, hist_j, lam, tile)
        error_msg = ''
    
    # get UFS IC data (TODO: flag to read in RESTART data rather than IC data and implement
    # different file reads)
    
    if (error_msg):
        print("ERROR: unknown grid orientation")
        exit()

    if not date:
        # date was not included on command line; look in atmf* file for initial date
        date = find_date(forcing_dir, lam)
    
    #get grid cell area if not given
    if not area:
        area = get_UFS_grid_area(grid_dir, tile, IC_i, IC_j, lam)
    surface_data["area"] = area
    surface_data["lon"]  = hist_lon
    surface_data["lat"]  = hist_lat
    
    # Get UFS forcing data
    stateInit = {}
    if forcing_method == 1:
        exact_mode = True
        (forcing_data, comp_data, stateInit) = get_UFS_forcing_data(state_data["nlevs"], state_data,    \
                                                                    location, use_nearest, forcing_dir, \
                                                                    grid_dir, tile, IC_i, IC_j, lam,\
                                                                    save_comp, exact_mode)
    elif forcing_method == 2:
        geo_pts = n_forcing_halo_points
        geo_neighbors = ''
        geo_dx = ''
        geo_dy = ''
        if (geos_wind_forcing and use_actual_geos_wind):
            #need to check and potentially get more neighboring data for geostrophic wind (for low enough Rossby numbers where geostrophic balance can be assumed)
            (geo_pts, geo_neighbors, geo_dx, geo_dy) = expand_neighbors_for_geostrophic_balance(forcing_dir, hist_i, hist_j, neighbors, dx, dy, state_data["nlevs"], lam)
        (forcing_data, comp_data) = get_UFS_forcing_data_advective_tendency(forcing_dir, hist_i, hist_j, tile, neighbors, dx, dy, state_data["nlevs"], lam, save_comp, vertical_method, \
                                                                            geos_wind_forcing, wind_nudge, geo_pts, geo_neighbors, geo_dx, geo_dy, use_actual_geos_wind)
    elif forcing_method == 3:
        exact_mode = False
        (forcing_data, comp_data, stateInit) = get_UFS_forcing_data(state_data["nlevs"], state_data,    \
                                                                    location, use_nearest, forcing_dir, \
                                                                    grid_dir, tile, IC_i, IC_j, lam,\
                                                                    save_comp, exact_mode)
    else:
        message = "Unsupported forcing_method chosen"
        logging.critical(message)
        raise Exception(message)
    
    
    # Write SCM case file
    fileOUT = write_SCM_case_file(state_data, surface_data, oro_data, forcing_data, stateInit, case_name, date, forcing_method, vertical_method, geos_wind_forcing, wind_nudge)

    # read in and remap the state variables to the first history file pressure profile and 
    # write them out to compare SCM output to (atmf for state variables and sfcf for physics 
    # tendencies)
    if (save_comp):
        write_comparison_file(comp_data, case_name, date, surface_data)
    
if __name__ == '__main__':
    main()
