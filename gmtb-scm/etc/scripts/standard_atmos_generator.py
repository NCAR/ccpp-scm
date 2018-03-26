#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np

#standard mid-latitude summer atmosphere
#data is from Anderson et al. (1986) "AFGL Atmospheric Constituent Profiles (0-120km)"

#height in m
height = 1.0E3*np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0,
    24.0, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 55.0,
    60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0])

#pressure in Pa
pressure = 100.0*np.array([1.013E3, 9.02E2, 8.02E2, 7.10E2, 6.28E2, 5.54E2, 4.87E2, 4.26E2,
    3.72E2, 3.24E2, 2.81E2, 2.43E2, 2.09E2, 1.79E2, 1.53E2, 1.30E2, 1.11E2, 9.5E1,
    8.12E1, 6.95E1, 5.95E1, 5.1E1, 4.37E1, 3.76E1, 3.22E1, 2.77E1, 1.907E1, 1.32E1,
    9.3, 6.52, 4.64, 3.33, 2.41, 1.76, 1.29, 9.51E-1, 5.15E-1, 2.72E-1, 1.39E-1,
    6.7E-2, 3.0E-2, 1.2E-2, 4.48E-3, 1.64E-3, 6.25E-4, 2.58E-4, 1.17E-4, 6.11E-5,
    3.56E-5, 2.27E-5])

#T in K
temperature = np.array([299.7, 293.7, 287.7, 283.7, 277.0, 270.3, 263.6, 257.0, 250.3,
    243.6, 237.0, 230.1, 223.6, 217.0, 210.3, 203.7, 197.0, 194.8, 198.8, 202.7, 206.7,
    210.7, 214.6, 217.0, 219.2, 221.4, 227.0, 232.3, 237.7, 243.1, 248.5, 254.0, 259.4,
    264.8, 269.6, 270.2, 263.4, 253.1, 236.0, 218.9, 201.8, 184.8, 177.1, 177.0, 184.3,
    190.7, 212.0, 241.6, 299.7, 380.0])

#r_h20 in ppmv (convert to specific humidity)
r_h2o = np.array([2.59E4, 1.95E4, 1.53E4, 8.60E3, 4.44E3, 3.35E3, 2.10E3, 1.29E3, 7.64E2,
    4.10E2, 1.91E2, 7.31E1, 2.91E1, 9.9, 6.22, 4.0, 3.0, 2.9, 2.75, 2.6, 2.6, 2.65, 2.8,
    2.9, 3.2, 3.25, 3.6, 4.0, 4.3, 4.6, 4.9, 5.2, 5.5, 5.7, 5.9, 6.0, 6.0, 6.0, 5.4, 4.5, 3.3,
    2.1, 1.3, 8.5E-1, 5.4E-1, 4.0E-1, 3.4E-1, 2.8E-1, 2.4E-1, 2.0E-1])
mol_weight_h2o = 18.01528 #g/mol
mol_weight_air = 28.97
#convert to mass mixing ratio and multiply by 1.0E-6
r_h2o = 1.0E-6*r_h2o*(mol_weight_h2o/mol_weight_air)
#convert to specific humidity (kg/kg)
q_v = r_h2o/(1.0 + r_h2o)

#r_03 in ppmv (convert to mass mixing ratio)
r_o3 = np.array([2.87E-2, 3.15E-2, 3.34E-2, 3.5E-2, 3.56E-2, 3.77E-2, 3.99E-2, 4.22E-2,
    4.47E-2, 5.0E-2, 5.6E-2, 6.61E-2, 7.82E-2, 9.29E-2, 1.05E-1, 1.26E-1, 1.44E-1,
    2.5E-1, 5.0E-1, 9.5E-1, 1.4, 1.8, 2.4, 3.4, 4.3, 5.4, 7.8, 9.3, 9.85, 9.7, 8.8,
    7.5, 5.9, 4.5, 3.45, 2.8, 1.8, 1.1, 6.5E-1, 3.0E-1, 1.8E-1, 3.3E-1, 5.0E-1, 5.2E-1,
    5.0E-1, 4.0E-1, 2.0E-1, 5.0E-2, 5.0E-3, 5.0E-4])
mol_weight_o3 = 47.997
r_o3 = 1.0E-6*r_o3*(mol_weight_o3/mol_weight_air)


writefile_fid = Dataset('../raw_case_input/mid_lat_summer_std.nc', 'w', format='NETCDF4')
writefile_fid.description = "Mid-latitude Summer Standard Atmosphere"

writefile_height_dim = writefile_fid.createDimension('height', None)
writefile_height_var = writefile_fid.createVariable('height', 'f4', ('height',))
writefile_height_var[:] = height
writefile_height_var.units = 'm'
writefile_height_var.description = 'height above sea level'

writefile_pressure_var = writefile_fid.createVariable('pressure', 'f4', ('height',))
writefile_pressure_var[:] = pressure
writefile_pressure_var.units = 'Pa'
writefile_pressure_var.description = 'atmospheric pressure'

writefile_T_var = writefile_fid.createVariable('temperature', 'f4', ('height',))
writefile_T_var[:] = temperature
writefile_T_var.units = 'K'
writefile_T_var.description = 'atmospheric absolute temperature'

writefile_qv_var = writefile_fid.createVariable('q_v', 'f4', ('height',))
writefile_qv_var[:] = q_v
writefile_qv_var.units = 'kg/kg'
writefile_qv_var.description = 'water vapor specific humidity'

writefile_o3_var = writefile_fid.createVariable('o3', 'f4', ('height',))
writefile_o3_var[:] = r_o3
writefile_o3_var.units = 'kg/kg'
writefile_o3_var.description = 'ozone mass mixing ratio'

#close file
writefile_fid.close()
