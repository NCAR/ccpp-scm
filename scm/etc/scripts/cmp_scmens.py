#!/usr/bin/env python
#################################################################################################
# Dependencies
#################################################################################################
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import argparse as ap
import datetime
from mpl_toolkits.basemap import Basemap

def read_SCMin(fileIN):
    #
    data = xr.open_dataset(fileIN)
    #

    surfaceType_str = data.attrs['surfaceType']
    if (surfaceType_str == 'ocean'): surfaceType = 0
    if (surfaceType_str == 'land'):  surfaceType = 1
    if (surfaceType_str == 'ice'):   surfaceType = 2

    state = {}
    state["surfaceType"] = surfaceType
    #
    return (state)
#################################################################################################
def read_SCMout(fileIN):
    #
    data = xr.open_dataset(fileIN)
    #
    state={}
    state["time_inst"]        = data["time_inst"].values
    state["time_diag"]        = data["time_diag"].values
    state["pres"]             = data["pres"].values[:,:,0]
    state["T"]                = data["T"].values[:,:,0]
    state["qv"]               = data["qv"].values[:,:,0]
    state["u"]                = data["u"].values[:,:,0]
    state["v"]                = data["v"].values[:,:,0]
    state["dtend_temp_lw"]    = data["dT_dt_lwrad"].values[:,:,0]
    state["dtend_temp_sw"]    = data["dT_dt_swrad"].values[:,:,0]
    state["dtend_temp_phys"]  = data["dT_dt_phys"].values[:,:,0]
    state["dtend_u_phys"]     = data["du_dt_phys"].values[:,:,0]
    state["dtend_v_phys"]     = data["dv_dt_phys"].values[:,:,0]
    state["dtend_qv_phys"]    = data["dq_dt_phys"].values[:,:,0]
    state["dtend_qv_deepcnv"] = data["dq_dt_deepconv"].values[:,:,0]
    state["dtend_qv_pbl"]     = data["dq_dt_pbl"].values[:,:,0]
    state["dtend_temp_pbl"]   = data["dT_dt_pbl"].values[:,:,0]
    #
    return (state)

#################################################################################################
def read_UFScomp(fileIN):
    #
    data = xr.open_dataset(fileIN)
    print(fileIN)
    #
    state = {}
    state["time"]             = data["time"].values#/1.e9
    state["pres"]             = data["lev"].values
    state["T"]                = data["temp"].values
    state["qv"]               = data["qv"].values
    state["u"]                = data["u"].values
    state["v"]                = data["v"].values
    state["dtend_qv_deepcnv"] = data["dtend_qv_deepcnv"].values
    state["dtend_qv_pbl"]     = data["dtend_qv_pbl"].values
    state["dtend_temp_pbl"]   = data["dtend_temp_pbl"].values
    #
    return (state)

#################################################################################################
cases = [1]
case_names = []
for case in cases:
    case_names.append("fv3_SCM_ensemble_"+str(case).zfill(2)+"hr_UFSforcing")

suite = "SCM_GFS_v16"
nens  = 1
ncase = len(case_names)

#
# Read in all cases
# 
init = True
for ic, case_name in enumerate(case_names):
    #
    # Read in all ensemble members for case
    #
    ensCount = 0
    scm_data_case = {}
    ufs_data_case = {}
    init = True
    for ensmember in range(nens):
        #
        # Construct filenames from current case and ensemble member
        #
        fileSCMout = "../../run/output_"+case_name+"_n"+str(ensmember).zfill(3)+"_"+suite+"/output.nc"
        fileUFSref = "../../data/comparison_data/"+case_name+"_n"+str(ensmember).zfill(3)+"_comp_data.nc"
        fileSCMin  = "../../data/processed_case_input/"+case_name+"_n"+str(ensmember).zfill(3)+".nc"
        print(fileSCMout)
        #
        # Read in data
        #
        if os.path.exists(fileSCMout):
            stateSCM = read_SCMout(fileSCMout)
            stateUFS = read_UFScomp(fileUFSref)
            stateSCMin = read_SCMin(fileSCMin)
            #
            if (init):
                ntimeUFS  = stateUFS["T"].shape[0]
                ntimeiSCM = stateSCM["T"].shape[0]
                ntimedSCM = stateSCM["dtend_temp_phys"].shape[0]
                nlev      = stateUFS["T"].shape[1]
                #
                stateUFSf = {"pres":             np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "temp":             np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "qv":               np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "u":                np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "v":                np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "dtend_qv_deepcnv": np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "dtend_qv_pbl":     np.zeros((ntimeUFS, nlev, nens),  dtype=float), \
                             "dtend_temp_pbl":   np.zeros((ntimeUFS, nlev, nens),  dtype=float)}
                stateUFSf["time"]  = stateUFS["time"]

                #
                stateSCMf = {"pres":             np.zeros((ntimeiSCM, nlev, nens), dtype=float), \
                             "temp":             np.zeros((ntimeiSCM, nlev, nens), dtype=float), \
                             "qv":               np.zeros((ntimeiSCM, nlev, nens), dtype=float), \
                             "u":                np.zeros((ntimeiSCM, nlev, nens), dtype=float), \
                             "v":                np.zeros((ntimeiSCM, nlev, nens), dtype=float), \
                             "dtend_qv_deepcnv": np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_qv_pbl":     np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_temp_pbl":   np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_temp_lw":    np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_temp_sw":    np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_temp_phys":  np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_u_phys":     np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_v_phys":     np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "dtend_qv_phys":    np.zeros((ntimedSCM, nlev, nens), dtype=float), \
                             "surfaceType":      np.zeros((                 nens), dtype=int)}
                stateSCMf["timei"] = stateSCM["time_inst"]
                stateSCMf["timed"] = stateSCM["time_diag"]
                init = False

            #
            for itime in range(0,ntimeUFS):
                stateUFSf["pres"][itime,:,ensCount]     = stateUFS["pres"][:]
            stateUFSf["temp"][:,:,ensCount]             = stateUFS["T"][:,:]
            stateUFSf["qv"][:,:,ensCount]               = stateUFS["qv"][:,:]*1000.#From kg/kg to g/kg
            stateUFSf["u"][:,:,ensCount]                = stateUFS["u"][:,:]
            stateUFSf["v"][:,:,ensCount]                = stateUFS["v"][:,:]
            stateUFSf["dtend_qv_deepcnv"][:,:,ensCount] = stateUFS["dtend_qv_deepcnv"][:,:]
            stateUFSf["dtend_qv_pbl"][:,:,ensCount]     = stateUFS["dtend_qv_pbl"][:,:]
            stateUFSf["dtend_temp_pbl"][:,:,ensCount]   = stateUFS["dtend_temp_pbl"][:,:]
            #
            stateSCMf["pres"][:,:,ensCount]             = stateSCM["pres"][:,:]
            stateSCMf["temp"][:,:,ensCount]             = stateSCM["T"][:,:]
            stateSCMf["qv"][:,:,ensCount]               = stateSCM["qv"][:,:]*1000.#From kg/kg to g/kg
            stateSCMf["u"][:,:,ensCount]                = stateSCM["u"][:,:]
            stateSCMf["v"][:,:,ensCount]                = stateSCM["v"][:,:]
            stateSCMf["dtend_qv_deepcnv"][:,:,ensCount] = stateSCM["dtend_qv_deepcnv"][:,:]
            stateSCMf["dtend_qv_pbl"][:,:,ensCount]     = stateSCM["dtend_qv_pbl"][:,:]
            stateSCMf["dtend_temp_pbl"][:,:,ensCount]   = stateSCM["dtend_temp_pbl"][:,:]

            stateSCMf["surfaceType"][ensCount]         = stateSCMin["surfaceType"]
            #
            ensCount = ensCount + 1
#    for ix in range(ntimeiSCM):
#        for iy in range(nlev):
#            print("pres: ",ix,stateSCM["pres"][ix,iy])
    #
    # Compute stats for all gridpoints, also partition results by surface type.
    #
    scm_stats = [{"name": "pres",             "dims": [nlev,ntimeiSCM], "time": "timei", "units": "Pa"},        \
                 {"name": "temp",             "dims": [nlev,ntimeiSCM], "time": "timei", "units": "K"},         \
                 {"name": "qv",               "dims": [nlev,ntimeiSCM], "time": "timei", "units": "g kg-1"},    \
                 {"name": "u",                "dims": [nlev,ntimeiSCM], "time": "timei", "units": "m s-1"},     \
                 {"name": "v",                "dims": [nlev,ntimeiSCM], "time": "timei", "units": "m s-1"},     \
                 {"name": "dtend_qv_deepcnv", "dims": [nlev,ntimedSCM], "time": "timed", "units": "g kg-1 s-1"},\
                 {"name": "dtend_qv_pbl",     "dims": [nlev,ntimedSCM], "time": "timed", "units": "g kg-1 s-1"},\
                 {"name": "dtend_temp_pbl",   "dims": [nlev,ntimedSCM], "time": "timed", "units": "K s-1"}]
    #
    ufs_stats = [{"name": "pres",             "dims": [nlev,ntimeUFS],  "time": "time",  "units": "Pa"},        \
                 {"name": "temp",             "dims": [nlev,ntimeUFS],  "time": "time",  "units": "K"},         \
                 {"name": "qv",               "dims": [nlev,ntimeUFS],  "time": "time",  "units": "g kg-1"},    \
                 {"name": "u",                "dims": [nlev,ntimeUFS],  "time": "time",  "units": "m s-1"},     \
                 {"name": "v",                "dims": [nlev,ntimeUFS],  "time": "time",  "units": "m s-1"},     \
                 {"name": "dtend_qv_deepcnv", "dims": [nlev,ntimeUFS],  "time": "time",  "units": "g kg-1 s-1"},\
                 {"name": "dtend_qv_pbl",     "dims": [nlev,ntimeUFS],  "time": "time",  "units": "g kg-1 s-1"},\
                 {"name": "dtend_temp_pbl",   "dims": [nlev,ntimeUFS],  "time": "time",  "units": "K s-1"}]


    #
    # Identify the land/ocean/ice points
    #
    isea = np.where(stateSCMf["surfaceType"] == 0)[0]
    ilnd = np.where(stateSCMf["surfaceType"] == 1)[0]
    iice = np.where(stateSCMf["surfaceType"] == 2)[0]

    #
    # Compute SCM ensemble statistics...
    #
    for data in scm_stats:
        data["mean"]      = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev"]     = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_sea"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_sea"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_lnd"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_lnd"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_ice"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_ice"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        for ilev in range(0,data["dims"][0]):
            for itime in range(0,data["dims"][1]):
                # All columns
                data["mean"][ilev,itime]  = np.mean(stateSCMf[data["name"]][itime,ilev,:])
                data["stdev"][ilev,itime] = np.std(stateSCMf[data["name"]][itime,ilev,:])
                # Ocean columns only
                if isea.size > 0:
                    data["mean_sea"][ilev,itime]  = np.mean(stateSCMf[data["name"]][itime,ilev,isea])
                    data["stdev_sea"][ilev,itime] = np.std(stateSCMf[data["name"]][itime,ilev,isea])
                # Land columns only
                if ilnd.size > 0:
                    data["mean_lnd"][ilev,itime]  = np.mean(stateSCMf[data["name"]][itime,ilev,ilnd])
                    data["stdev_lnd"][ilev,itime] = np.std(stateSCMf[data["name"]][itime,ilev,ilnd])
                # Ice columns only
                if iice.size > 0:
                    data["mean_ice"][ilev,itime]  = np.mean(stateSCMf[data["name"]][itime,ilev,iice])
                    data["stdev_ice"][ilev,itime] = np.std(stateSCMf[data["name"]][itime,ilev,iice])
#            print(data["name"]," - ",stateSCMf[data["name"]][1,ilev,:],stateUFSf[data["name"]][1,ilev,:],data["mean"][ilev,1])

    #
    # Compute UFS ensemble statistics...
    #
    for data in ufs_stats:
        data["mean"]      = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev"]     = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_sea"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_sea"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_lnd"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_lnd"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["mean_ice"]  = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        data["stdev_ice"] = np.zeros((data["dims"][0],data["dims"][1]), dtype=float)
        for ilev in range(0,data["dims"][0]):
            for itime in range(0,data["dims"][1]):
                # All columns
                data["mean"][ilev,itime]  = np.mean(stateUFSf[data["name"]][itime,ilev,:])
                data["stdev"][ilev,itime] = np.std(stateUFSf[data["name"]][itime,ilev,:])
                # Ocean columns only
                if isea.size > 0:
                    data["mean_sea"][ilev,itime]  = np.mean(stateUFSf[data["name"]][itime,ilev,isea])
                    data["stdev_sea"][ilev,itime] = np.std(stateUFSf[data["name"]][itime,ilev,isea])
                # Land columns only
                if ilnd.size > 0:
                    data["mean_lnd"][ilev,itime]  = np.mean(stateUFSf[data["name"]][itime,ilev,ilnd])
                    data["stdev_lnd"][ilev,itime] = np.std(stateUFSf[data["name"]][itime,ilev,ilnd])
                # Ice columns only
                if iice.size > 0:
                    data["mean_ice"][ilev,itime]  = np.mean(stateUFSf[data["name"]][itime,ilev,iice])
                    data["stdev_ice"][ilev,itime] = np.std(stateUFSf[data["name"]][itime,ilev,iice])

    # All three datasets are on different temporal grids. Find the correct indices for all.
    # Match the UFS comparison data times, this ensures that you are looking at the correct
    # timesteps that the tendecies are valid for.
    itUFS = [i for i in range(ntimeUFS)]
    #
    itSCMi = []
    countUFS = 0
    for count,time in enumerate(stateSCMf["timei"]):
        if (countUFS < ntimeUFS):
            if (time == stateUFSf["time"][countUFS]):
                itSCMi.append(count)
                countUFS=countUFS+1
    #
    itSCMd = []
    countUFS = 0
    for count,time in enumerate(stateSCMf["timed"]):
        if (countUFS < ntimeUFS):
            if (time == stateUFSf["time"][countUFS]):
                itSCMd.append(count)
                countUFS=countUFS+1

    #
    # This will create plots at all available UFS diagnostic times
    #
    hrs2plt = [i for i in range(cases[ic],72+1,cases[ic])]
    hrs2plt = [1,2,3,6,12,24,48,60,72,96,120,144]
    stats2plt = ["mean","mean_sea","mean_lnd","mean_ice"]
    #stats2plt = ["mean"]
    for ihour in hrs2plt:
        # Get indices for fields on different temporal grids.
        print("TIME: ",stateUFSf["time"] )
        it_ufs  = np.where(stateUFSf["time"]  == ihour*3600.)[0][0]
        it_scmi = np.where(stateSCMf["timei"] == ihour*3600.)[0][0]
        it_scmd = np.where(stateSCMf["timed"] == ihour*3600.)[0][0]
        #print(stateUFSf["time"][it_ufs],stateSCMf["timei"][it_scmi],stateSCMf["timed"][it_scmd])

        for stats in stats2plt:
            # Make figure(s)
            fig = plt.figure(figsize=(8,10))
            for count,data in enumerate(scm_stats[1::],start=1):
                #print(stats,'-',data["name"])
                if (data["time"] == "timei"): it_scm = it_scmi
                if (data["time"] == "timed"): it_scm = it_scmd
                plt.subplot(2,4,count)
                plt.plot(ufs_stats[count][stats][:,it_ufs], ufs_stats[0][stats][:,it_ufs]*0.01, color='blue')
                plt.plot(data[stats][:,it_scm],             scm_stats[0][stats][:,it_scm]*0.01, color='green')
                plt.title(data["name"])
                plt.xlabel("("+data["units"]+")")
                if count == 1 or count == 5 or count == 9: 
                    plt.ylabel("P (hPa)")
                else:
                    ax = plt.gca()
                    ax.axes.yaxis.set_ticklabels([])
                plt.ylim(1000,100)
            plt.show()
            # Save figure
            fileOUT = 'scm_ufs.'+suite+'.'+str(cases[ic]).zfill(2)+'hrdiag.'+str(ihour).zfill(3)+'fcst.'+stats+'.png'
            plt.savefig(fileOUT)
            print(fileOUT)
            #plt.title(r'$\frac{dT}{dt}shortwave$')
