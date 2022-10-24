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

#
parser = argparse.ArgumentParser()
parser.add_argument('-i',   '--dir_ic',           help='Directory path containing FV3 input files',                               required=True)
parser.add_argument('-g',   '--dir_grid',         help='directory path containing FV3 tile supergrid files',                      required=True)
parser.add_argument('-f',   '--dir_forcing',      help='Directory path containing physics diag files',                            required=True)
parser.add_argument('-n',   '--case_name',        help='Name of case',                                                            required=True)
parser.add_argument('-lon', '--lon_limits',       help='Longitude range for ensemble, separated by a space', nargs=2, type=float, required=False)
parser.add_argument('-lat', '--lat_limits',       help='Latitude range for ensemble, separated by a space',  nargs=2, type=float, required=False)
parser.add_argument('-nens','--nemsmembers',      help='Number of SCM UFS ensemble memebers to create',               type=int,   default = 1)
parser.add_argument('-dt',  '--timestep',         help='SCM timestep, in seconds',                                    type=int,   default = 3600)
parser.add_argument('-fhz', '--fhzero',           help='UFS frequency, in hours, for emptying diagnostic buckets.',   type=int,   default = 1)
parser.add_argument('-cres','--UFS_resolution',   help='UFS spatial resolution',                                      type=int,   default = 96)
parser.add_argument('-sdf', '--suite',            help='CCPP suite definition file to use for ensemble',                          default = 'SCM_GFS_v16')
parser.add_argument('-sc',  '--save_comp',        help='Flag to save a file with UFS data for comparisons',                       action='store_true')
parser.add_argument('-lsm', '--add_UFS_NOAH_lsm', help='Flag to include UFS NOAH LSM surface forcing',                            action='store_true')
parser.add_argument('-ufsf','--add_UFS_dyn_tend', help='Flag to include UFS dynamic tendencies for SCM forcing',                  action='store_true')

###############################################################################
###############################################################################
def parse_arguments():
    args             = parser.parse_args()
    lon_limits       = args.lon_limits
    lat_limits       = args.lat_limits
    case_name_root   = args.case_name
    dirIC            = args.dir_ic
    dirGRID          = args.dir_grid
    dirFORCING       = args.dir_forcing
    npts             = args.nemsmembers
    suite            = args.suite
    dt               = args.timestep
    fhzero           = args.fhzero
    C_RES            = args.UFS_resolution
    save_comp        = args.save_comp
    add_UFS_NOAH_lsm = args.add_UFS_NOAH_lsm
    add_UFS_dyn_tend = args.add_UFS_dyn_tend

    return (lon_limits, lat_limits, case_name_root, npts, suite, dt, fhzero, C_RES, dirIC, dirGRID, dirFORCING, save_comp, add_UFS_NOAH_lsm, add_UFS_dyn_tend)

###############################################################################
###############################################################################
def main():

    #
    # Get command line arguments
    #
    (lon_limits, lat_limits, case_name_root, npts, suite, dt, fhzero, C_RES, dirIC, dirGRID, dirFORCING, save_comp, add_UFS_NOAH_lsm, add_UFS_dyn_tend) = parse_arguments()

    #
    # Make sure that SCM_WORK has been set.
    #
    try:
        dir_scm = os.getenv('SCM_WORK') + "/ccpp-scm/"
    except:
        print("Environment variable SCM_WORK not set. Stopping.")
        exit()


    ###########################################################################
    #
    # Generate random points for SCM ensemble
    #
    ###########################################################################
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
        if lat_limits:
            lats[ipt] = lat_limits[0] + (lat_limits[1]-lat_limits[0])*rng1[ipt]
        else:
            lats[ipt] = rng1[ipt]*180-90
        if lon_limits:
            lons[ipt] = lon_limits[0] + (lon_limits[1]-lon_limits[0])*rng2[ipt]
        else:
            lons[ipt] = rng2[ipt]*360

    ###########################################################################
    #
    # Use longitude and latitude provided (TODO)
    #
    ###########################################################################


    ###########################################################################
    #
    # Create SCM case configuration (etc/case_config) file.
    #
    ###########################################################################

    # Are we providing LSM ICs to the SCM?
    model_ics     = ".false."
    if add_UFS_NOAH_lsm:
        model_ics = ".true."

    # How many timesteps between clearing the diagnostic buckets?
    n_itt_diag = int(fhzero*3600/dt)

    # 
    n_itt_out = int(n_itt_diag/fhzero)

    # For DEPHY format, these do not change.
    input_type  = 1

    #
    case_config =[{"name": "model_ics",   "values": model_ics},       \
                  {"name": "input_type",  "values": str(input_type)}, \
                  #{"name": "dt",          "values": str(dt)},         \
                  {"name": "C_RES",       "values": str(C_RES)}]
                  #{"name": "n_itt_diag",  "values": str(n_itt_diag)}, \
                  #{"name": "n_itt_out",   "values": str(n_itt_out)}]

    #
    # What, if any, options neeed to be passsed to UFS_IC_generator.py?
    #
    com_config = ''
    if save_comp:        com_config = com_config + ' -sc'
    if add_UFS_NOAH_lsm: com_config = com_config + ' -lsm'
    if add_UFS_dyn_tend: com_config = com_config + ' -ufsf'

    #
    # Create inputs to SCM
    # (Call UFS_IC_generator and create case_config files for each case)
    #
    case_list    = ""
    case_list_nf = ""
    count = 0
    run_list = []
    run_list_nf = []
    for pt in range(0,npts):
        #
        # Call UFS_IC_generator.py
        #
        case_name     = case_name_root +"_n" + str(pt).zfill(3)
        file_scminput = dir_scm+"scm/data/processed_case_input/"+case_name+".nc"
        if add_UFS_dyn_tend:
            case_name_nf      = case_name_root +"_n" + str(pt).zfill(3)+"_noforce"
            file_scm_nf_input = dir_scm+"scm/data/processed_case_input/"+case_name+"_noforce.nc"
        #
        com = "./UFS_IC_generator.py -l " +str(lons[pt]) + " " + str(lats[pt]) + \
              " -i " + dirIC + " -g " + dirGRID + " -f " + dirFORCING + " -n " + case_name + com_config
        print(com)
        os.system(com)

        #
        # In case there are errors in UFS_IC_generator and the input is not generated.
        #
        if os.path.exists(file_scminput):
            #
            # Add case to ensemble list.
            #
            case_list     = case_list + '"'+case_name+'"'
            if (count != npts-1): case_list = case_list + ', '
            #
            if add_UFS_dyn_tend:
                case_list_nf = case_list_nf + '"'+case_name_nf+'"'
                if (count != npts-1): case_list_nf = case_list_nf + ', '

            #
            # What is the surface type? (get from SCM input file)
            #
            dataset  = xr.open_dataset(file_scminput)
            sfc_type = int(np.round_(dataset.slmsk.values[0][0][0]))

            #
            # Create case_config file(s)
            #
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
            #
            if add_UFS_dyn_tend:
                fileOUT = dir_scm+"scm/etc/case_config/"+case_name_nf+".nml"
                fileID  = open(fileOUT, 'w')
                fileID.write('$case_config')
                fileID.write('\n')
                fileID.write('case_name = ' + "'" + case_name_nf + "',")
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
            run_list.append({"case": case_name, "suite": suite})
            if add_UFS_dyn_tend:
                run_list_nf.append({"case": case_name_nf, "suite": suite})
            
            #
            count = count + 1

    #
    # Create "multirun file list" needed by run_scm.py
    #
    os.system("mkdir -p "+dir_scm+"scm/bin/")
    fileOUT = "scm_ufsens.py"
    fileID  = open(dir_scm+"scm/bin/"+fileOUT, 'w')
    fileID.write('run_list = [')
    fileID.write('\n')
    for run in run_list:
        fileID.write('            {"case": "' + run["case"] + '", "suite": "' + run["suite"] + '"},')
        fileID.write('\n')
    fileID.write('            ]')
    fileID.close()
    #
    if add_UFS_dyn_tend:
        os.system("mkdir -p "+dir_scm+"scm/bin/")
        fileOUT_nf = "scm_ufsens_nf.py"
        fileID  = open(dir_scm+"scm/bin/"+fileOUT_nf, 'w')
        fileID.write('run_list = [')
        fileID.write('\n')
        for run in run_list_nf:
            fileID.write('            {"case": "' + run["case"] + '", "suite": "' + run["suite"] + '"},')
            fileID.write('\n')
        fileID.write('            ]')
        fileID.close()

    print("-------------------------------------------------------------------------------------------")
    print("Command(s) to execute in ccpp-scm/scm/bin/: ")
    print(" ")
    print("./run_scm.py --file " + fileOUT + " --n_itt_diag " + \
          str(n_itt_diag) + " --n_itt_out " + str(n_itt_out) + " --timestep "   + \
          str(dt))
    print("")
    if add_UFS_dyn_tend:
        print("")
        print("./run_scm.py --file " + fileOUT_nf + " --n_itt_diag " + \
              str(n_itt_diag) + " --n_itt_out " + str(n_itt_out)    + " --timestep "   + \
              str(dt))
    print("")
    print("-------------------------------------------------------------------------------------------")

if __name__ == '__main__':
    main()