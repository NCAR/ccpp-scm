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
import re

###############################################################################
# Global settings                                                             #
###############################################################################

#Physical constants
earth_radius = 6371000.0 #m

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
parser.add_argument('-g', '--grid_dir',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-ic','--lsm_dir',    help='directory path containing FV3 tile IC files', required=True)
parser.add_argument('-t', '--tile',       help='tile of desired point (if known - bypasses tile search if present)', type=int, choices=range(1,7))

###############################################################################
# Functions and subroutines                                                   #
###############################################################################

def parse_arguments():
    """Parse command line arguments"""
    args = parser.parse_args()
    location = args.location
    index = args.index
    grid_dir = args.grid_dir
    lsm_dir = args.lsm_dir
    tile = args.tile
    
    #validate args
    if not os.path.exists(grid_dir):
        message = 'The directory {0} does not exist'.format(grid_dir)
        logging.critical(message)
        raise Exception(message)
    if not os.path.exists(lsm_dir):
        message = 'The directory {0} does not exist'.format(lsm_dir)
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
        
    return (location, index, grid_dir, lsm_dir, tile)

def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)
    
def find_tile(loc, dir):
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
            longitude = np.array(nc_file['x']).swapaxes(0,1)
            latitude = np.array(nc_file['y']).swapaxes(0,1)
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
                    return f_name.split('tile')[1].split('.nc')[0] 
            else:
                polar_tile_filenames.append(f_name)
                
    #if the tile hasn't been found by this point, it must be contained within a polar tile
    for f_name in polar_tile_filenames:
        nc_file = Dataset('{0}/{1}'.format(dir,f_name))
        latitude = np.array(nc_file['y']).swapaxes(0,1)
        nc_file.close()
        
        #if the sign of the mean latitude of the tile is the same as that of the point, the tile has been found
        if np.sign(np.mean(latitude)) == np.sign(loc[1]):
            found_tile = True
            return f_name.split('tile')[1].split('.nc')[0]        
    return -1

def find_loc_indices(loc, dir, tile):
    """Find the nearest neighbor FV3 grid point given a lon/lat pair and the tile number"""
    #returns the indices of the nearest neighbor point in the given tile, the lon/lat of the nearest neighbor, 
    #and the distance (m) from the given point to the nearest neighbor grid cell
    
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
    lon_super = np.array(nc_file['x'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    lat_super = np.array(nc_file['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
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
    cart_loc = np.array(sph2cart(math.radians(temp_loc[0]), math.radians(temp_loc[1]), earth_radius))
    
    for i in range(len(longitude)):
        for j in range(len(longitude[i])):
            #get the Cartesian location of all grid points
            cart_cell = np.array(sph2cart(math.radians(longitude[i,j]), math.radians(latitude[i,j]), earth_radius))
            
            #calculate the euclidean distance from the given point to the current grid cell
            eucl_dist[i,j] = np.linalg.norm(cart_loc - cart_cell)
    
    #get the indices of the grid point with the minimum euclidean distance to the given point
    i,j = np.unravel_index(eucl_dist.argmin(), eucl_dist.shape)
    
    return (i,j,longitude[i,j]%360.0, latitude[i,j], eucl_dist[i,j])

def find_lon_lat_of_indices(indices, dir, tile):
    """Find the longitude and latitude of the given indices within the given tile."""
    
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
    lon_super = np.array(nc_file['x'])   #[lat,lon] or [y,x]   #.swapaxes(0,1)
    lat_super = np.array(nc_file['y'])    #[lat,lon] or [y,x]   #.swapaxes(0,1)
    #get the longitude and latitude data for the grid centers by slicing the supergrid 
    #and taking only odd-indexed values
    longitude = lon_super[1::2,1::2]
    latitude = lat_super[1::2,1::2]
    nc_file.close()
    
    return (longitude[indices[1],indices[0]], latitude[indices[1],indices[0]])

def get_UFS_surface_fix_data(fix_dir, tile, i, j):
    
    #fix_dir = dir + '/fix_sfc'
    
    filename_pattern = '*facsf.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    facsf = nc_file['facsf'][0,j,i]
    facwf = 1.0 - facsf
    nc_file.close()
    
    filename_pattern = '*maximum_snow_albedo.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    maximum_snow_albedo = nc_file['maximum_snow_albedo'][0,j,i]
    nc_file.close()
    
    filename_pattern = '*slope_type.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    slope_type = nc_file['slope_type'][0,j,i]
    nc_file.close()
    
    filename_pattern = '*snowfree_albedo.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    alvsf = nc_file['visible_black_sky_albedo'][0,j,i]
    alvwf = nc_file['visible_white_sky_albedo'][0,j,i]
    alnsf = nc_file['near_IR_black_sky_albedo'][0,j,i]
    alnwf = nc_file['near_IR_white_sky_albedo'][0,j,i]
    nc_file.close()
    
    filename_pattern = '*soil_type.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    soil_type = nc_file['soil_type'][0,j,i]
    nc_file.close()
    
    filename_pattern = '*substrate_temperature.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    substrate_temperature = nc_file['substrate_temperature'][0,j,i]
    nc_file.close()
    
    filename_pattern = '*vegetation_greenness.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    vegetation_greenness = nc_file['vegetation_greenness'][0,j,i]
    max_vegetation_greenness = np.amax(nc_file['vegetation_greenness'][0,:,:])
    min_vegetation_greenness = np.amin(nc_file['vegetation_greenness'][0,:,:])
    nc_file.close()
    
    filename_pattern = '*vegetation_type.tile{0}.nc'.format(tile)
    for f_name in os.listdir(fix_dir):
       if fnmatch.fnmatch(f_name, filename_pattern):
          filename = f_name
    if not filename:
        message = 'No filenames matching the pattern {0} found in {1}'.format(filename_pattern,fix_dir)
        logging.critical(message)
        raise Exception(message)
    
    nc_file = Dataset('{0}/{1}'.format(fix_dir,filename))
    vegetation_type = nc_file['vegetation_type'][0,j,i]
    nc_file.close()
    
    return (facsf, facwf, maximum_snow_albedo, alvsf, alvwf, alnsf, alnwf, substrate_temperature, vegetation_greenness, max_vegetation_greenness, min_vegetation_greenness, slope_type, soil_type, vegetation_type)

def sph2cart(az, el, r):
    """Calculate the Cartesian coordiates from spherical coordinates"""
    
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    
    return (x, y, z)    

def get_UFS_grid_area(dir, tile, i, j):
    """Get the horizontal grid cell area for the given tile and indices"""
    #this information is in the supergrid files
    
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
    area_in=nc_file['area'][jpt2-1:jpt2+1,ipt2-1:ipt2+1]
    
    return area_in.sum()

def main():
    setup_logging()
    
    #read in arguments
    (location, indices, grid_dir, lsm_dir, tile) = parse_arguments()
        
    #find tile containing the point using the supergrid if no tile is specified 
    if not tile:
        tile = int(find_tile(location, grid_dir))
        if tile < 0:
            message = 'No tile was found for location {0}'.format(location)
            logging.critical(message)
            raise Exception(message)
        print('Tile found: {0}'.format(tile))
    
    #find index of closest point in the tile if indices are not specified
    if not indices:
        (tile_j, tile_i, point_lon, point_lat, dist_min) = find_loc_indices(location, grid_dir, tile)
        print('The closest point in tile {0} has indices [{1},{2}]'.format(tile,tile_i,tile_j))
        print('This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat))
        print('This grid cell is approximately {0} km away from the desired location of {1} {2}'.format(dist_min/1.0E3,location[0],location[1]))
    else:
        tile_i = indices[0]
        tile_j = indices[1]
        #still need to grab the lon/lat if the tile and indices are supplied
        (point_lon, point_lat) = find_lon_lat_of_indices(indices, grid_dir, tile)
        
        print('This index has a central longitude/latitude of [{0},{1}]'.format(point_lon,point_lat))
    
    
    #get grid cell area if not given
    #area = get_UFS_grid_area(grid_dir, tile, tile_i, tile_j)
        
    (facsf, facwf, max_snow_alb, alvsf, alvwf, alnsf, alnwf, substrate_t, veg_greenness, max_veg_greenness, min_veg_greenness, slope_type, soil_type, veg_type) = get_UFS_surface_fix_data(lsm_dir, tile, tile_i, tile_j)
    
    print("facsf,facwf={0},{1}".format(facsf,facwf))
    print("maximum_snow_albedo={}".format(max_snow_alb))
    print("alvsf,alvwf,alnsf,alnwf={0},{1},{2},{3}".format(alvsf,alvwf,alnsf,alnwf))
    print("substrate_temperature={}".format(substrate_t))
    print("vegetation_greenness,max,min={0},{1},{2}".format(veg_greenness,max_veg_greenness,min_veg_greenness))
    print("slope_type={}".format(slope_type))
    print("soil_type={}".format(soil_type))
    print("vegetation_type={}".format(veg_type))
    
if __name__ == '__main__':
    main()
