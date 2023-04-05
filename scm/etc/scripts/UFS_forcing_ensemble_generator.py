#!/usr/bin/env python

###############################################################################
# Dependencies
###############################################################################
import argparse
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import random

###############################################################################
# Argument list
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-d',    '--dir',           help='path to UFS Regression Test output', required=True)
parser.add_argument('-n',    '--case_name',     help='name of case',                       required=True)
parser.add_argument('-lonl', '--lon_limits',    help='longitude range for ensemble, separated by a space', nargs=2,   type=float, required=False)
parser.add_argument('-latl', '--lat_limits',    help='latitude range for ensemble, separated by a space',  nargs=2,   type=float, required=False)
parser.add_argument('-lons', '--lon_list',      help='longitudes for ensemble, separated by a space',      nargs='*', type=float, required=False)
parser.add_argument('-lats', '--lat_list',      help='latitudes for ensemble, separated by a space',       nargs='*', type=float, required=False)
parser.add_argument('-nens', '--nensmembers',   help='number of SCM UFS ensemble memebers to create',                 type=int,   required=False)
parser.add_argument('-dt',   '--timestep',      help='sCM timestep, in seconds',                                      type=int,   default = 3600)
parser.add_argument('-fhz',  '--fhzero',        help='UFS frequency, in hours, for emptying diagnostic buckets.',     type=int,   default = 1)
parser.add_argument('-cres', '--C_RES',         help='UFS spatial resolution',                                        type=int,   default = 96)
parser.add_argument('-sdf',  '--suite',         help='CCPP suite definition file to use for ensemble',                            default = 'SCM_GFS_v16')
parser.add_argument('-sc',   '--save_comp',     help='flag to save a file with UFS data for comparisons',                         action='store_true')
parser.add_argument('-near', '--use_nearest',   help='flag to indicate using the nearest UFS history file gridpoint, no regridding',action='store_true')

###############################################################################
# Main program
###############################################################################
def main():
    # Get command line arguments
    args  = parser.parse_args()

    if (not args.dir):
        print("ERROR: Need to provide UFS RT directory!")
        exit()
    else:
        args.dir_ic      = args.dir + "/INPUT/"
        args.dir_grid    = args.dir + "/INPUT/"
        args.dir_forcing = args.dir

    # Error checking
    if (args.lon_limits and args.lon_list):
        print("ERROR: Can't provide explicit longitude(s) AND a longitude range")
        exit()
    if (args.lat_limits and args.lat_list):
        print("ERROR: Can't provide explicit latitude(s) AND a latitude range")
        exit()
    if (args.lon_limits or args.lat_limits) and not args.nensmembers:
        print("ERROR: Longitude/Latitude range provided, but NOT ensemble count.")
        exit()

    if (args.nensmembers):
        npts = args.nensmembers
        if (args.lat_list or args.lon_list):
            print("ERROR: Can't provide explicit lon/lat range AND number of points for ensemble generation.")
            exit()
    else:
        if (args.lon_list and args.lat_list):
            if (len(args.lon_list) == len(args.lat_list)):
                npts = len(args.lon_list)
            else:
                print("ERROR: Number of longitude/latitudes are inconsistent")
                exit()

    # Make sure that SCM_WORK has been set.
    try:
        dir_scm = os.getenv('SCM_WORK')+'/'
    except:
        print("Environment variable SCM_WORK not set. Stopping.")
        exit()

    ###########################################################################
    #
    # Set longitude/latitude 
    #
    ###########################################################################
    if (args.nensmembers):
        rng1 = np.zeros((npts), dtype=float)
        rng2 = np.zeros((npts), dtype=float)
        lons = np.zeros((npts), dtype=float)
        lats = np.zeros((npts), dtype=float)
        for ipt in range(npts):
            # Here the seed is set to give the same set of points each time.
            random.seed(ipt)
            rng1[ipt] = random.randint(1,1000)*0.001
            rng2[ipt] = random.randint(2,1000)*0.001
            #
            if args.lat_limits:
                lats[ipt] = args.lat_limits[0] + (args.lat_limits[1]-args.lat_limits[0])*rng1[ipt]
            else:
                lats[ipt] = rng1[ipt]*180-90
            if args.lon_limits:
                lons[ipt] = args.lon_limits[0] + (args.lon_limits[1]-args.lon_limits[0])*rng2[ipt]
            else:
                lons[ipt] = rng2[ipt]*360
    ###########################################################################
    #
    # Use longitude and latitude provided
    #
    ###########################################################################
    else:
        lons = np.asarray(args.lon_list)
        lats = np.asarray(args.lat_list)

    ###########################################################################
    #
    # Create SCM case configuration (etc/case_config) file.
    #
    ###########################################################################

    # How many timesteps between clearing the diagnostic buckets?
    n_itt_diag = int(args.fhzero*3600/args.timestep)
    # 
    n_itt_out = int(n_itt_diag/args.fhzero)

    #
    case_config =[{"name": "input_type",  "values": str(1)},             \
                  {"name": "dt",          "values": str(args.timestep)}, \
                  {"name": "C_RES",       "values": str(args.C_RES)}]

    # What, if any, options neeed to be passsed to UFS_IC_generator.py?
    com_config = ''
    if args.save_comp:   com_config = com_config + ' -sc'
    if args.use_nearest: com_config = com_config + ' -near'

    # Create inputs to SCM
    case_list    = ""
    case_list_nf = ""
    count = 0
    run_list = []
    for pt in range(0,npts):
        # Call UFS_IC_generator.py
        case_name     = args.case_name +"_n" + str(pt).zfill(3)
        file_scminput = dir_scm+"scm/data/processed_case_input/"+case_name+"_SCM_driver.nc"
        com = "./UFS_IC_generator.py -l " +str(lons[pt]) + " " + str(lats[pt]) + \
              " -i " + args.dir_ic + " -g " + args.dir_grid + " -f " + args.dir_forcing + " -n " + case_name + com_config
        print(com)
        os.system(com)

        # Add case to ensemble list.
        case_list     = case_list + '"'+case_name+'"'
        if (count != npts-1): case_list = case_list + ', '

        # What is the surface type? (get from SCM input file)
        dataset  = xr.open_dataset(file_scminput)
        sfc_type = int(np.round_(dataset.slmsk.values[0]))

        # Create case_config file(s)
        fileOUT = dir_scm+"scm/etc/case_config/"+case_name+".nml"
        fileID  = open(fileOUT, 'w')
        fileID.write('$case_config')
        fileID.write('\n')
        fileID.write('case_name = ' + "'" + case_name + "',")
        fileID.write('\n')
        fileID.write('sfc_type = ' + str(sfc_type) + ",")
        fileID.write('\n')
        for opts in case_config:
            fileID.write(opts["name"] + ' = ' + opts["values"] + ",")
            fileID.write('\n')
        fileID.write('$end')
        fileID.write('\n')
        fileID.close()

        # Add case to dictionary to be used by run_scm.py
        run_list.append({"case": case_name, "suite": args.suite})
            
        #
        count = count + 1

    ###########################################################################
    #
    # Create "multirun file list" needed by run_scm.py
    #
    ###########################################################################
    com = "mkdir -p "+dir_scm+"scm/bin/"
    print(com)
    os.system(com)
    fileOUT = "scm_ufsens_"+case_name+".py"
    fileID  = open(dir_scm+"scm/bin/"+fileOUT, 'w')
    fileID.write('run_list = [')
    fileID.write('\n')
    for run in run_list:
        #print('            {"case": "' , run["case"] , '", "suite": "' , run["suite"] , '"},')
        fileID.write('            {"case": "' + run["case"] + '", "suite": "' + run["suite"] + '"},')
        fileID.write('\n')
    fileID.write('            ]')
    fileID.close()

    ###########################################################################
    #
    # Display run commands
    #
    ###########################################################################
    print("-------------------------------------------------------------------------------------------")
    print("Command(s) to execute in ccpp-scm/scm/bin/: ")
    print(" ")
    print("./run_scm.py --npz_type gfs --file " + fileOUT + " --n_itt_diag " + \
          str(n_itt_diag) + " --n_itt_out " + str(n_itt_out) + " --timestep "   + \
          str(args.timestep))
    print("")
    print("")
    print("-------------------------------------------------------------------------------------------")

if __name__ == '__main__':
    main()
