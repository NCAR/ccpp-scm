#!/usr/bin/env python

from __future__ import print_function

from netCDF4 import Dataset
import numpy as np

g = 9.80665
R_dry = 287.0
c_p = 1004.0
p0 = 100000.0
L_v = 2.5E6
L_s = 2.834E6

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

def omega_to_w(omega, p, T):
    rho = p/(R_dry*T)
    w = -omega/(rho*g)
    return w

def w_to_omega(w, p, T):
    rho = p/(R_dry*T)
    omega = -w*rho*g
    return omega

def get_height_from_pres(T, p, z_sfc):
    z = np.zeros((p.size),dtype=float)

    z[0] = z_sfc
    for k in range(p.size-1):
        dphi = -R_dry*0.5*(T[k]+T[k+1])/(0.5*(p[k]+p[k+1]))*(p[k+1]-p[k])
        z[k+1] = z[k] + dphi/g
    return z

def theta_to_T(theta, p):
    T = theta*(p/p0)**(R_dry/c_p)
    return T
