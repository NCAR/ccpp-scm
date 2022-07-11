#!/usr/bin/env python

###############################################################################
# Dependencies
###############################################################################
import argparse
import os
from netCDF4 import Dataset
import xarray as xr
import numpy as np

#
parser = argparse.ArgumentParser()
parser.add_argument('-lon', '--lon_limits', help='longitude range for ensemble, separated by a space', nargs=2, type=float, required=False)
parser.add_argument('-lat', '--lat_limits', help='latitude range for ensemble, separated by a space', nargs=2, type=float, required=False)
parser.add_argument('-d',   '--date',       help='date corresponding to initial conditions in YYYYMMDDHHMMSS format', required=False)
parser.add_argument('-n',   '--case_name',  help='name of case', required=True)
parser.add_argument('-i',   '--dir_ic',     help='input directory path containing FV3 input files', required=True)
parser.add_argument('-g',   '--dir_grid',   help='directory path containing FV3 tile supergrid files', required=True)
parser.add_argument('-f',   '--dir_forcing',help='directory path containing physics diag files', required=True)
parser.add_argument('-lami','--lam_ic',     help='flag to signal that the ICs are from a limited-area model run', action='store_true')
parser.add_argument('-nens','--nemsmembers',help='number of SCM UFS ensemble memeber to create', type=int,default=1)
parser.add_argument('-cfg', '--configID',   help='Tag name used to create case configuration name',required=False)
parser.add_argument('-sdf', '--suite',      help='CCPP suite definition file to use for ensemble',default='SCM_GFS_v16')

###############################################################################
###############################################################################
def parse_arguments():
    args           = parser.parse_args()
    lon_range      = args.lon_limits
    lat_range      = args.lat_limits
    date           = args.date
    case_name_root = args.case_name
    dirIC          = args.dir_ic
    dirGRID        = args.dir_grid
    dirFORCING     = args.dir_forcing
    lami           = args.lam_ic
    npts           = args.nemsmembers
    suite          = args.suite
    configID       = ''
    if (args.configID != None):
        configID = '_'+args.configID

    #
    date_dict = {}
    if date:
        if len(date) != 14:
            message = 'The entered date {0} does not have the 14 characters expected in the format YYYYMMDDHHMMSS'.format(date)
            logging.critical(message)
            raise Exception(message)
        else:
            date_dict["year"]   = int(date[0:4])
            date_dict["month"]  = int(date[4:6])
            date_dict["day"]    = int(date[6:8])
            date_dict["hour"]   = int(date[8:10])
            date_dict["minute"] = int(date[10:12])
            date_dict["second"] = int(date[12:])
    #
    return (date, date_dict, lon_range, lat_range, case_name_root, lami, npts, configID, suite, dirIC, dirGRID, dirFORCING)

###############################################################################
###############################################################################
def main():
    #
    (date, date_dict, lon_range, lat_range, case_name_root, lami, npts, configID, suite, dirIC, dirGRID, dirFORCING) = parse_arguments()

    #
    fileLOG  = open('UFS_forcing_ensemble_generator.log', 'w')

    #
    dir_scm = "/glade/u/home/dswales/Projects/SCM/grantfirl/ccpp-scm/"
    
    # Generate random points for SCM ensemble
    if lat_range:
        lats = lat_range[0] + (lat_range[1]-lat_range[0])*np.random.rand(npts)
    else:
        lats = np.random.rand(npts)*180-90
    if lon_range:
        lons = lon_range[0] + (lon_range[1]-lon_range[0])*np.random.rand(npts)
    else:
        lons = np.random.rand(npts)*360
    
    runtime             = 7200
    model_ics           = ".true."
    input_type          = 1
    sfc_flux_spec       = ".false."
    ref_profile_choice  = 2
    dt                  = 600.
    C_RES               = 96
    n_itt_out           = 1#90
    n_itt_diag          = 6
    case_config =[ \
                   #{"name": "runtime",                  "values": str(runtime)}, \
                   {"name": "model_ics",                "values": model_ics}, \
                   {"name": "input_type",               "values": str(input_type)}, \
                   {"name": "sfc_flux_spec",            "values": sfc_flux_spec}, \
                   {"name": "reference_profile_choice", "values": str(ref_profile_choice)}, \
                   {"name": "dt",                       "values": str(dt)}, \
                   {"name": "C_RES",                    "values": str(C_RES)}, \
                   #{"name": "year",                     "values": str(date_dict["year"]).zfill(4)},\
                   #{"name": "month",                    "values": str(date_dict["month"]).zfill(2)},\
                   #{"name": "day",                      "values": str(date_dict["day"]).zfill(2)},\
                   #{"name": "hour",                     "values": str(date_dict["hour"]).zfill(2)}]
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
        case_list     = case_list + '"'+case_name+'"'
        file_scminput = dir_scm+"scm/data/processed_case_input/"+case_name+".nc"
        if (count != npts-1): case_list = case_list + ', '
        com = "./UFS_IC_generator.py -l " +str(lons[pt]).zfill(6) + " " + str(lats[pt]).zfill(6) + " -d " + date + \
              " -i " + dirIC + " -g " + dirGRID + " -f " + dirFORCING + " -n " + case_name + " -lami -sc"
        os.system(com)
        
        print(file_scminput)
        #
        # In case there are errors in UFS_IC_generator and the input is not generated.
        #
        if os.path.exists(file_scminput):
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
