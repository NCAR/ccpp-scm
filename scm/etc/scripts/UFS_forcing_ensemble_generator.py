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
parser.add_argument('-lon', '--lon_limits', help='longitude range for ensemble, separated by a space', nargs=2, type=float, required=False)
parser.add_argument('-lat', '--lat_limits', help='latitude range for ensemble, separated by a space', nargs=2, type=float, required=False)
parser.add_argument('-n',   '--case_name',  help='name of case', required=True)
parser.add_argument('-i',   '--dir_ic',     help='input directory path containing FV3 input files', required=True)
parser.add_argument('-g',   '--dir_grid',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-f',   '--dir_forcing',help='directory path containing physics diag files', required=True)
parser.add_argument('-lami','--lam_ic',     help='flag to signal that the ICs are from a limited-area model run', action='store_true')
parser.add_argument('-nens','--nemsmembers',help='number of SCM UFS ensemble memeber to create', type=int,default=1)
parser.add_argument('-cfg', '--configID',   help='Tag name used to create case configuration name',required=False)
parser.add_argument('-sdf', '--suite',      help='CCPP suite definition file to use for ensemble',default='SCM_GFS_v16')
parser.add_argument('-dt',  '--deltatime',  help='', type=float, default=3600.0)
parser.add_argument('-nio', '--n_itt_out',  help='',type=int, default=1)
parser.add_argument('-nid', '--n_itt_diag', help='',type=int, default=-999)
parser.add_argument('-rr',  '--rerun',      help='Use when building additional cases for previously created case',action='store_true')
parser.add_argument('-rrid','--rerun_ID',   help='suffix to be attached when using rerun',required=False)

###############################################################################
###############################################################################
def parse_arguments():
    args           = parser.parse_args()
    lon_range      = args.lon_limits
    lat_range      = args.lat_limits
    case_name_root = args.case_name
    dirIC          = args.dir_ic
    dirGRID        = args.dir_grid
    dirFORCING     = args.dir_forcing
    lami           = args.lam_ic
    npts           = args.nemsmembers
    suite          = args.suite
    dt             = args.deltatime
    n_itt_out      = args.n_itt_out
    n_itt_diag     = args.n_itt_diag
    rerun          = args.rerun
    rerun_ID       = ''
    if (args.rerun_ID != None): rerun_ID = '_'+args.rerun_ID
    configID       = ''
    if (args.configID != None): configID = '_'+args.configID

    return (lon_range, lat_range, case_name_root, lami, npts, configID, suite, dt, n_itt_out, n_itt_diag, rerun, rerun_ID, dirIC, dirGRID, dirFORCING)

###############################################################################
###############################################################################
def main():
    #
    (lon_range, lat_range, case_name_root, lami, npts, configID, suite, dt, n_itt_out, n_itt_diag, rerun, rerun_ID, dirIC, dirGRID, dirFORCING) = parse_arguments()

    #
    fileLOG  = open('UFS_forcing_ensemble_generator.log', 'w')

    #
    dir_scm = "/glade/u/home/dswales/Projects/SCM/grantfirl/ccpp-scm/"

    #
    # Generate random points for SCM ensemble
    #
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
        if lat_range:
            lats[ipt] = lat_range[0] + (lat_range[1]-lat_range[0])*rng1[ipt]
        else:
            lats[ipt] = rng1[ipt]*180-90
        if lon_range:
            lons[ipt] = lon_range[0] + (lon_range[1]-lon_range[0])*rng2[ipt]
        else:
            lons[ipt] = rng2[ipt]*360

    model_ics           = ".true."
    input_type          = 1
    sfc_flux_spec       = ".false."
    ref_profile_choice  = 2
    C_RES               = 96
    case_config =[ \
                   {"name": "model_ics",                "values": model_ics}, \
                   {"name": "input_type",               "values": str(input_type)}, \
                   {"name": "sfc_flux_spec",            "values": sfc_flux_spec}, \
                   {"name": "reference_profile_choice", "values": str(ref_profile_choice)}, \
                   {"name": "dt",                       "values": str(dt)}, \
                   {"name": "C_RES",                    "values": str(C_RES)}, \
                   {"name": "n_itt_out",                "values": str(n_itt_out)},  \
                   {"name": "n_itt_diag",               "values": str(n_itt_diag)}]
    
    #
    # Create inputs to SCM
    # (Call UFS_IC_generator and create case_config files for each case)
    #
    case_list = ""
    count = 0
    for pt in range(0,npts):
        #
        # Call UFS_IC_generator.py
        #
        case_name     = case_name_root +"_n" + str(pt).zfill(3)
        file_scminput = dir_scm+"scm/data/processed_case_input/"+case_name+".nc"

        com = "./UFS_IC_generator.py -l " +str(lons[pt]) + " " + str(lats[pt]) + \
              " -i " + dirIC + " -g " + dirGRID + " -f " + dirFORCING + " -n " + case_name + " -lami -sc -lsm"
        if (not rerun): os.system(com)

        print(file_scminput)
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
            # What is the surface type? (get from SCM input file)
            #
            dataset  = xr.open_dataset(file_scminput)
            sfc_type = int(np.round_(dataset.slmsk.values[0][0][0]))

            fileLOG.write("Writing SCM UFS inputs at "+str(lons[pt]).zfill(6)+" "+str(lats[pt]).zfill(6)+" surface_type="+str(sfc_type))
            fileLOG.write('\n')

            #
            # Create case_config file
            #
            fileOUT = dir_scm+"scm/etc/case_config/"+case_name+".nml"
            print(fileOUT)
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
            count = count + 1

    #
    # Create "multirun file list" needed by run_scm.py
    #
    os.system("mkdir -p "+dir_scm+"scm/run/")
    fileOUT = dir_scm+"scm/run/scm_ufsens"+configID+".py"
    if rerun: fileOUT = dir_scm+"scm/run/scm_ufsens"+configID+rerun_ID+".py"
    fileID  = open(fileOUT, 'w')
    fileID.write('cases      = ['+case_list+']')
    fileID.write('\n')
    fileID.write('suites     = ["'+suite+'"]')
    fileID.write('\n')
    fileID.write('namelists  = []')
    fileID.write('\n')
    fileID.close()

    #
    fileLOG.close()

if __name__ == '__main__':
    main()
