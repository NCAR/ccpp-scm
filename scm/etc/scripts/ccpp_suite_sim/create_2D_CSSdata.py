#!/usr/bin/env python3
#################################################################################################
# Dependencies
#################################################################################################
import os
import numpy as np
import xarray as xr
import argparse as ap

##########################################################################################
# Routine to aggregate and compute hourly mean for SCM output.
##########################################################################################
def agg_scm_out(data_var, data_time, data_dt):
    # Constants
    sec_in_hr = 60.*60.*24.
    
    # Dimensions
    ntime = data_var[:,0].size
    nlev  = data_var[0,:].size

    # Local time (in seconds)
    time_local = []
    for itime in range(0,ntime-1):
        time_local.append(data_time[itime]/sec_in_hr-np.trunc(data_time[itime]/sec_in_hr))
    time_local = np.asarray(time_local)*sec_in_hr

    # Aggregate SCM output data into bins for requested frequency (data_dt)
    ntime_data = int(sec_in_hr/data_dt)
    data_out   = np.zeros((ntime_data,nlev),dtype='float')
    data_count = np.zeros((ntime_data),     dtype='int')
    t_limit    = data_dt
    t_count    = 0

    for itime in range(0,ntime-1):
        if (data_time[itime]/sec_in_hr >= 1 and t_count == 23):
            t_count = 0
            t_limit = 0
        if (time_local[itime] < t_limit):
            data_out[t_count,:] = data_out[t_count,:] + data_var[itime,:]
            data_count[t_count] = data_count[t_count] + 1
        else:
            t_limit = t_limit + data_dt
            t_count = t_count + 1
            data_out[t_count,:] = data_out[t_count,:] + data_var[itime,:]
            data_count[t_count] = data_count[t_count] + 1

    # Compute statistics
    data_mean = np.empty((ntime_data,nlev),dtype='float')
    for itime in range(0,ntime_data):
        data_mean[itime,:] = data_out[itime,:] / data_count[itime]

    return(data_mean)

#################################################################################################
# Argument list ./create_2D_CSSdata.py --cases twpice --suites SCM_GFS_v17_p8
#################################################################################################
parser = ap.ArgumentParser()
parser.add_argument('-c',   '--cases',  help='name of case(s), separated by a space', nargs='*', required=True)
parser.add_argument('-sdf', '--suites', help='CCPP SDF(s), separated by a space',     nargs='*', required=True)
parser.add_argument('-n',   '--nmls',   help='namelists, separated by a space',       nargs='*', required=False)
parser.add_argument('-dt',  '--time',   help='data interval, in seconds', type=int,    default = 3600)


def main():

    #############################################################################################
    # Get arguments
    #############################################################################################
    args      = parser.parse_args()
    case_name = args.cases
    suite     = args.suites
    data_dt   = args.time
    ncases    = len(case_name)

    if args.nmls:
        namelist  = args.nmls
    else:
        namelist = []
        for icase in range(0,ncases):
            namelist.append("")

    if (len(case_name) != len(suite) or len(case_name) != len(namelist)):
        print('ERROR: cases, suite, or namelist provided are not the same length.')
        print('ERROR: number of cases:     ',len(case_name))
        print('ERROR: number of suites:    ',len(suite))
        print('ERROR: number of namelists: ',len(case_name))
        exit()

    # Create python dictionary for cases to create.
    run_dict=[]
    for icase in range(0,ncases):
        run_dict.append({"case":case_name[icase], "suite":suite[icase], "namelist": namelist[icase]})

    #
    run_dir  = '../../../run/'
    var_dict = [{"name": "dT_dt_swrad",     "time": "time_diag"},\
                {"name": "dT_dt_lwrad",     "time": "time_diag"},\
                {"name": "dT_dt_pbl",       "time": "time_diag"},\
                {"name": "dT_dt_deepconv",  "time": "time_diag"},\
                {"name": "dT_dt_shalconv",  "time": "time_diag"},\
                {"name": "dT_dt_micro",     "time": "time_diag"},\
                {"name": "dT_dt_ogwd",      "time": "time_diag"},\
                {"name": "dT_dt_cgwd",      "time": "time_diag"},\
                {"name": "dq_dt_pbl",       "time": "time_diag"},\
                {"name": "dq_dt_deepconv",  "time": "time_diag"},\
                {"name": "dq_dt_shalconv",  "time": "time_diag"},\
                {"name": "dq_dt_micro",     "time": "time_diag"},\
                {"name": "du_dt_pbl",       "time": "time_diag"},\
                {"name": "du_dt_deepconv",  "time": "time_diag"},\
                {"name": "du_dt_shalconv",  "time": "time_diag"},\
                {"name": "du_dt_ogwd",      "time": "time_diag"},\
                {"name": "du_dt_cgwd",      "time": "time_diag"},\
                {"name": "dv_dt_pbl",       "time": "time_diag"},\
                {"name": "dv_dt_deepconv",  "time": "time_diag"},\
                {"name": "dv_dt_shalconv",  "time": "time_diag"},\
                {"name": "dv_dt_ogwd",      "time": "time_diag"},\
                {"name": "dv_dt_cgwd",      "time": "time_diag"},\
                {"name": "doz_dt_pbl",      "time": "time_diag"},\
                {"name": "doz_dt_prodloss", "time": "time_diag"},\
                {"name": "doz_dt_oz",       "time": "time_diag"},\
                {"name": "doz_dt_T",        "time": "time_diag"},\
                {"name": "doz_dt_ovhd",     "time": "time_diag"}]

    # Loop over all run(s)
    for run in run_dict:

        file_out = "data_CSS_2D."+run["case"] + "." + run["suite"]+".nc"

        # SCM data location for current run
        local_dir = run_dir + "output_" + run["case"] + "_" + run["suite"]

        # Read in SCM output and aggregate. Store in "var_dict" under "mean"
        ds = xr.open_dataset(local_dir+"/output.nc")
        for var in var_dict:
            try:
                data_var           = ds[var["name"]].values[:,:,0]
                data_time          = ds[var["time"]].values
                var["units"]       = ds[var["name"]].attrs["units"]
                var["description"] = ds[var["name"]].attrs["description"]
                var["mean"]        = agg_scm_out(data_var, data_time, data_dt)
            except:
                print(var["name"]," Not available. Skipping...")

        # Get pressure
        data_var  = ds.pres.values[:,:,0]
        data_time = ds.time_inst.values
        plev = agg_scm_out(data_var, data_time, data_dt)

        # Close file for current run
        ds.close()

        ntime = plev[:,0].size
        nlev  = plev[0,:].size

        # Create output file for CCPP suite simulator.
        time_varout = xr.Dataset({"times": (("time"),1800.+np.linspace(0,(ntime-1)*3600.,ntime))},
                                 coords = {"time": np.linspace(0,ntime,ntime)})
        time_varout["times"].attrs = {"units":"seconds"}
        lev_varout = xr.Dataset({"pressure": (("time","lev"),plev)},
                                coords = {"time": np.linspace(0,ntime,ntime),
                                          "lev": np.linspace(0,nlev,nlev)})
        lev_varout["pressure"].attrs = {"units":"Pa"}

        varsout = [time_varout, lev_varout]
        for var in var_dict:
            temp = xr.Dataset({var["name"]: (("time","lev"),var["mean"])},
                              coords = {"time": np.linspace(0,ntime,ntime),
                                        "lev": np.linspace(0,nlev,nlev)})
            temp[var["name"]].attrs = {"units": var["units"], "description": var["description"]}
            varsout.append(temp)

        xr.merge(varsout).to_netcdf(file_out)

if __name__ == '__main__':
    main()
