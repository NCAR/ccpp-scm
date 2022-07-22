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

#################################################################################################
def read_SCMout(fileIN):
    #
    data = xr.open_dataset(fileIN)
    #
    state={}
    state["time_inst"]       = data["time_inst"].values
    state["time_diag"]       = data["time_diag"].values
    state["pres"]            = data["pres"].values[:,:,0]
    state["T"]               = data["T"].values[:,:,0]
    state["qv"]              = data["qv"].values[:,:,0]
    state["u"]               = data["u"].values[:,:,0]
    state["v"]               = data["v"].values[:,:,0]
    state["dtend_temp_lw"]   = data["dT_dt_lwrad"].values[:,:,0]
    state["dtend_temp_sw"]   = data["dT_dt_swrad"].values[:,:,0]
    state["dtend_temp_phys"] = data["dT_dt_phys"].values[:,:,0]
    state["dtend_u_phys"]    = data["du_dt_phys"].values[:,:,0]
    state["dtend_v_phys"]    = data["dv_dt_phys"].values[:,:,0]
    state["dtend_qv_phys"]   = data["dq_dt_phys"].values[:,:,0]
    #
    return (state)

#################################################################################################
def read_UFScomp(fileIN):
    #
    data = xr.open_dataset(fileIN)
    #
    state = {}
    state["time"]            = data["time"].values
    state["pres"]            = data["pres"].values[:,:,0,0]
    state["T"]               = data["T"].values[:,:,0,0]
    state["qv"]              = data["qv"].values[:,:,0,0]
    state["u"]               = data["u"].values[:,:,0,0]
    state["v"]               = data["v"].values[:,:,0,0]
    state["dtend_temp_lw"]   = data["dtend_temp_lw"].values[:,:,0,0]
    state["dtend_temp_sw"]   = data["dtend_temp_sw"].values[:,:,0,0]
    state["dtend_temp_phys"] = data["dtend_temp_phys"].values[:,:,0,0]
    state["dtend_u_phys"]    = data["dtend_u_phys"].values[:,:,0,0]
    state["dtend_v_phys"]    = data["dtend_v_phys"].values[:,:,0,0]
    state["dtend_qv_phys"]   = data["dtend_qv_phys"].values[:,:,0,0]
    state["lon"]             = data["lon"].values
    state["lat"]             = data["lat"].values
    #
    return (state)

#################################################################################################
cases = [6]#,6,3,1]
case_names = []
for case in cases:
    #case_names.append("fv3_SCM_ensemble_trop_"+str(case).zfill(2)+"hr_UFSforcing")
    case_names.append("fv3_SCM_ensemble_Atlantic_"+str(case).zfill(2)+"hr_UFSforcing")

namelist   = "SCM_GFS_v16"
nens       = 1
ncase      = len(case_names)

#
# Read in all cases
# 
init = True
case_data= []
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
        fileSCMout = "../../run/output_"+case_name+"_n"+str(ensmember).zfill(3)+"_"+namelist+"/output.nc"
        fileUFSref = "../../data/comparison_data/"+case_name+"_n"+str(ensmember).zfill(3)+"_comp_data.nc"

        print(fileSCMout)
        #
        # Read in data
        #
        if os.path.exists(fileSCMout):
            stateSCM = read_SCMout(fileSCMout)
            stateUFS = read_UFScomp(fileUFSref)
            #
            if (init):
                ntimeUFS  = stateUFS["T"].shape[0]
                ntimeiSCM = stateSCM["T"].shape[0]
                ntimedSCM = stateSCM["dtend_temp_phys"].shape[0]
                nlev      = stateUFS["T"].shape[1]
                #
                stateUFS_pres            = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_temp            = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_qv              = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_u               = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_v               = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_temp_lw   = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_temp_sw   = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_temp_phys = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_u_phys    = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_v_phys    = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_dtend_qv_phys   = np.zeros((ntimeUFS, nlev, nens),  dtype=float)
                stateUFS_lon             = np.zeros((                nens),  dtype=float)
                stateUFS_lat             = np.zeros((                nens),  dtype=float)
                #
                stateSCM_pres            = np.zeros((ntimeiSCM, nlev, nens), dtype=float)
                stateSCM_temp            = np.zeros((ntimeiSCM, nlev, nens), dtype=float)
                stateSCM_qv              = np.zeros((ntimeiSCM, nlev, nens), dtype=float)
                stateSCM_u               = np.zeros((ntimeiSCM, nlev, nens), dtype=float)
                stateSCM_v               = np.zeros((ntimeiSCM, nlev, nens), dtype=float)
                stateSCM_dtend_temp_lw   = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                stateSCM_dtend_temp_sw   = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                stateSCM_dtend_temp_phys = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                stateSCM_dtend_u_phys    = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                stateSCM_dtend_v_phys    = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                stateSCM_dtend_qv_phys   = np.zeros((ntimedSCM, nlev, nens), dtype=float)
                init = False
            #
            stateUFS_pres[:,:,ensCount]            = stateUFS["pres"][:,:]
            stateUFS_temp[:,:,ensCount]            = stateUFS["T"][:,:]
            stateUFS_qv[:,:,ensCount]              = stateUFS["qv"][:,:]*1000.#From kg/kg to g/kg
            stateUFS_u[:,:,ensCount]               = stateUFS["u"][:,:]
            stateUFS_v[:,:,ensCount]               = stateUFS["v"][:,:]
            stateUFS_dtend_temp_lw[:,:,ensCount]   = stateUFS["dtend_temp_lw"][:,:]*100000. # Scale by 1e5
            stateUFS_dtend_temp_sw[:,:,ensCount]   = stateUFS["dtend_temp_sw"][:,:]*100000. # Scale by 1e5
            stateUFS_dtend_temp_phys[:,:,ensCount] = stateUFS["dtend_temp_phys"][:,:]*100000. # Scale by 1e5
            stateUFS_dtend_u_phys[:,:,ensCount]    = stateUFS["dtend_u_phys"][:,:]*100000. # Scale by 1e5
            stateUFS_dtend_v_phys[:,:,ensCount]    = stateUFS["dtend_v_phys"][:,:]*100000. # Scale by 1e5
            stateUFS_dtend_qv_phys[:,:,ensCount]   = stateUFS["dtend_qv_phys"][:,:]*100000. # Scale by 1e5
            stateUFS_time                          = stateUFS["time"]
            stateUFS_lon[ensCount]                 = stateUFS["lon"]
            stateUFS_lat[ensCount]                 = stateUFS["lat"]
            #
            stateSCM_pres[:,:,ensCount]            = stateSCM["pres"][:,:]
            stateSCM_temp[:,:,ensCount]            = stateSCM["T"][:,:]
            stateSCM_qv[:,:,ensCount]              = stateSCM["qv"][:,:]*1000.#From kg/kg to g/kg
            stateSCM_u[:,:,ensCount]               = stateSCM["u"][:,:]
            stateSCM_v[:,:,ensCount]               = stateSCM["v"][:,:]
            stateSCM_dtend_temp_lw[:,:,ensCount]   = stateSCM["dtend_temp_lw"][:,:]*100000. # Scale by 1e5
            stateSCM_dtend_temp_sw[:,:,ensCount]   = stateSCM["dtend_temp_sw"][:,:]*100000. # Scale by 1e5
            stateSCM_dtend_temp_phys[:,:,ensCount] = stateSCM["dtend_temp_phys"][:,:]*100000. # Scale by 1e5
            stateSCM_dtend_u_phys[:,:,ensCount]    = stateSCM["dtend_u_phys"][:,:]*100000. # Scale by 1e5
            stateSCM_dtend_v_phys[:,:,ensCount]    = stateSCM["dtend_v_phys"][:,:]*100000. # Scale by 1e5
            stateSCM_dtend_qv_phys[:,:,ensCount]   = stateSCM["dtend_qv_phys"][:,:]*100000. # Scale by 1e5
            stateSCM_timei                         = stateSCM["time_inst"]
            stateSCM_timed                         = stateSCM["time_diag"]
            #
            ensCount = ensCount + 1

    #
    # Compute ensemble statistics for case
    #
    scm_mean_t               = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_std_t                = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_mean_qv              = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_std_qv               = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_mean_u               = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_std_u                = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_mean_v               = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_std_v                = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_mean_p               = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_std_p                = np.zeros((nlev,ntimeiSCM), dtype=float)
    scm_mean_dtend_temp_lw   = np.zeros((nlev,ntimedSCM), dtype=float)
    scm_std_dtend_temp_lw    = np.zeros((nlev,ntimedSCM), dtype=float)
    scm_mean_dtend_temp_sw   = np.zeros((nlev,ntimedSCM), dtype=float)
    scm_std_dtend_temp_sw    = np.zeros((nlev,ntimedSCM), dtype=float)
    scm_mean_dtend_temp_phys = np.zeros((nlev,ntimedSCM), dtype=float)
    scm_std_dtend_temp_phys  = np.zeros((nlev,ntimedSCM), dtype=float)
    #
    ufs_mean_t               = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_t                = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_qv              = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_qv               = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_u               = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_u                = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_v               = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_v                = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_p               = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_p                = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_dtend_temp_lw   = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_dtend_temp_lw    = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_dtend_temp_sw   = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_dtend_temp_sw    = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_mean_dtend_temp_phys = np.zeros((nlev,ntimeUFS),  dtype=float)
    ufs_std_dtend_temp_phys  = np.zeros((nlev,ntimeUFS),  dtype=float)

    # Make plots of ensemble pts.
    fig = plt.figure(figsize=(10,8))
    plt.scatter(stateUFS_lon[1::], stateUFS_lat[1::],color='blue')
    plt.ylim(-20,20)
    plt.xlim(0,360)
    fileOUT = "domain.png"
    plt.savefig(fileOUT)

    for ilev in range(0,nlev):
        for itime in range(0,ntimeiSCM):
            scm_mean_t[ilev,itime]               = np.mean(stateSCM_temp[itime,ilev,:])
            scm_std_t[ilev,itime]                = np.std( stateSCM_temp[itime,ilev,:])
            scm_mean_p[ilev,itime]               = np.mean(stateSCM_pres[itime,ilev,:])
            scm_std_p[ilev,itime]                = np.std( stateSCM_pres[itime,ilev,:])
            scm_mean_qv[ilev,itime]              = np.mean(stateSCM_qv[itime,ilev,:])
            scm_std_qv[ilev,itime]               = np.std( stateSCM_qv[itime,ilev,:])
            scm_mean_v[ilev,itime]               = np.mean(stateSCM_v[itime,ilev,:])
            scm_std_v[ilev,itime]                = np.std( stateSCM_v[itime,ilev,:])
            scm_mean_u[ilev,itime]               = np.mean(stateSCM_u[itime,ilev,:])
            scm_std_u[ilev,itime]                = np.std( stateSCM_u[itime,ilev,:])
        for itime in range(0,ntimedSCM):
            scm_mean_dtend_temp_lw[ilev,itime]   = np.mean(stateSCM_dtend_temp_lw[itime,ilev,:])
            scm_std_dtend_temp_lw[ilev,itime]    = np.std( stateSCM_dtend_temp_lw[itime,ilev,:])
            scm_mean_dtend_temp_sw[ilev,itime]   = np.mean(stateSCM_dtend_temp_sw[itime,ilev,:])
            scm_std_dtend_temp_sw[ilev,itime]    = np.std( stateSCM_dtend_temp_sw[itime,ilev,:])
            scm_mean_dtend_temp_phys[ilev,itime] = np.mean(stateSCM_dtend_temp_phys[itime,ilev,:])
            scm_std_dtend_temp_phys[ilev,itime]  = np.std( stateSCM_dtend_temp_phys[itime,ilev,:])
        for itime in range(0,ntimeUFS):
            ufs_mean_t[ilev,itime]               = np.mean(stateUFS_temp[itime,ilev,:])
            ufs_std_t[ilev,itime]                = np.std( stateUFS_temp[itime,ilev,:])
            ufs_mean_p[ilev,itime]               = np.mean(stateUFS_pres[itime,ilev,:])
            ufs_std_p[ilev,itime]                = np.std( stateUFS_pres[itime,ilev,:])
            ufs_mean_qv[ilev,itime]              = np.mean(stateUFS_qv[itime,ilev,:])
            ufs_std_qv[ilev,itime]               = np.std( stateUFS_qv[itime,ilev,:])
            ufs_mean_v[ilev,itime]               = np.mean(stateUFS_v[itime,ilev,:])
            ufs_std_v[ilev,itime]                = np.std( stateUFS_v[itime,ilev,:])
            ufs_mean_u[ilev,itime]               = np.mean(stateUFS_u[itime,ilev,:])
            ufs_std_u[ilev,itime]                = np.std( stateUFS_u[itime,ilev,:])
            ufs_mean_dtend_temp_lw[ilev,itime]   = np.mean(stateUFS_dtend_temp_lw[itime,ilev,:])
            ufs_std_dtend_temp_lw[ilev,itime]    = np.std( stateUFS_dtend_temp_lw[itime,ilev,:])
            ufs_mean_dtend_temp_sw[ilev,itime]   = np.mean(stateUFS_dtend_temp_sw[itime,ilev,:])
            ufs_std_dtend_temp_sw[ilev,itime]    = np.std( stateUFS_dtend_temp_sw[itime,ilev,:])
            ufs_mean_dtend_temp_phys[ilev,itime] = np.mean(stateUFS_dtend_temp_phys[itime,ilev,:])
            ufs_std_dtend_temp_phys[ilev,itime]  = np.std( stateUFS_dtend_temp_phys[itime,ilev,:])
    case_data.append({"name":        case_name,\
                      "scm_mean_t":  scm_mean_t,  "scm_std_t":  scm_std_t, \
                      "scm_mean_p":  scm_mean_p,  "scm_std_p":  scm_std_p, \
                      "scm_mean_qv": scm_mean_qv, "scm_std_qv": scm_std_qv,\
                      "scm_mean_u":  scm_mean_u,  "scm_std_u":  scm_std_u, \
                      "scm_mean_v":  scm_mean_v,  "scm_std_v":  scm_std_v, \
                      "scm_mean_dtend_temp_lw":   scm_mean_dtend_temp_lw,  \
                      "scm_std_dtend_temp_lw":    scm_std_dtend_temp_lw,   \
                      "scm_mean_dtend_temp_sw":   scm_mean_dtend_temp_sw,  \
                      "scm_std_dtend_temp_sw":    scm_std_dtend_temp_sw,   \
                      "scm_mean_dtend_temp_phys": scm_mean_dtend_temp_phys,\
                      "scm_std_dtend_temp_phys":  scm_std_dtend_temp_phys, \
                      "ufs_mean_t":  ufs_mean_t,  "ufs_std_t":  ufs_std_t, \
                      "ufs_mean_p":  ufs_mean_p,  "ufs_std_p":  ufs_std_p, \
                      "ufs_mean_qv": ufs_mean_qv, "ufs_std_qv": ufs_std_qv,\
                      "ufs_mean_u":  ufs_mean_u,  "ufs_std_u":  ufs_std_u, \
                      "ufs_mean_v":  ufs_mean_v,  "ufs_std_v":  ufs_std_v, \
                      "ufs_mean_dtend_temp_lw":   ufs_mean_dtend_temp_lw,  \
                      "ufs_std_dtend_temp_lw":    ufs_std_dtend_temp_lw,   \
                      "ufs_mean_dtend_temp_sw":   ufs_mean_dtend_temp_sw,  \
                      "ufs_std_dtend_temp_sw":    ufs_std_dtend_temp_sw,   \
                      "ufs_mean_dtend_temp_phys": ufs_mean_dtend_temp_phys,\
                      "ufs_std_dtend_temp_phys":  ufs_std_dtend_temp_phys })

    # All three datasets are on different temporal grids. Find the correct indices for all.
    # Match the UFS comparison data times, this ensures that you are looking at the correct
    # timesteps that the tendecies are valid for.

    itUFS = [i for i in range(ntimeUFS)]
    #
    itSCMi = []
    countUFS = 0
    for count,time in enumerate(stateSCM_timei):
        if (countUFS < ntimeUFS):
            if (time == stateUFS_time[countUFS]):
                itSCMi.append(count)
                countUFS=countUFS+1
    #
    itSCMd = []
    countUFS = 0
    for count,time in enumerate(stateSCM_timed):
        if (countUFS < ntimeUFS):
            if (time == stateUFS_time[countUFS]):
                itSCMd.append(count)
                countUFS=countUFS+1

        
    if (cases[ic] == 3):
        hrs2plt = [3,6,9,12,15,18,21,24]
    elif (cases[ic] == 6):
        hrs2plt = [6,12,18,24,36,48,60,72]
    elif (cases[ic] == 12):
        hrs2plt = [12,24]
    elif (cases[ic] == 1):
        hrs2plt = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
        hrs2plt = [1,3,6,9,12,15,18,21,24,36,48,60,72]

    for ihour in hrs2plt:
        it_ufs  = np.where(stateSCM_timei == ihour*3600.)[0][0]
        it_scmi = np.where(stateUFS_time  == ihour*3600.)[0][0]
        it_scmd = np.where(stateSCM_timed == ihour*3600.)[0][0]

        #
        # Make some plots
        #

        # Plot State variables (T, qv, u, v)
        fig = plt.figure(figsize=(10,8))
        #
        plt.subplot(2,4,1)
        plt.plot(ufs_mean_t[:,it_scmi], ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.plot(scm_mean_t[:,it_ufs],  ufs_mean_p[:,it_scmi]*0.01,color='green')
        plt.ylim(1000,100)
        plt.ylabel("P (hPa)")
        plt.title("Temperature")
        #
        plt.legend(["UFS","SCM_ensemble"],loc='lower left',fontsize='x-small')
        #
        plt.subplot(2,4,2)
        plt.plot(ufs_mean_dtend_temp_phys[:,it_scmi], ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.plot(scm_mean_dtend_temp_phys[:,it_scmd], ufs_mean_p[:,it_scmi]*0.01,color='green')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.title(r'$\frac{dT}{dt}physics$')
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        plt.subplot(2,4,3)
        plt.plot(ufs_mean_dtend_temp_lw[:,it_scmi], ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.plot(scm_mean_dtend_temp_lw[:,it_scmd], ufs_mean_p[:,it_scmi]*0.01,color='green')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.title(r'$\frac{dT}{dt}longwave$')
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        plt.subplot(2,4,4)
        plt.plot(ufs_mean_dtend_temp_sw[:,it_scmi], ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.plot(scm_mean_dtend_temp_sw[:,it_scmd], ufs_mean_p[:,it_scmi]*0.01,color='green')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.title(r'$\frac{dT}{dt}shortwave$')
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        plt.subplot(2,4,5)
        plt.plot(np.zeros((nlev),dtype=float),ufs_mean_p[:,it_scmi]*0.01,color='black')
        plt.plot(ufs_mean_t[:,it_scmi] - scm_mean_t[:,it_ufs],ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.ylabel("P (hPa)")
        plt.xlabel("(K)")
        #
        plt.subplot(2,4,6)
        plt.plot(np.zeros((nlev),dtype=float),ufs_mean_p[:,it_scmi]*0.01,color='black')
        plt.plot(ufs_mean_dtend_temp_phys[:,it_scmi] - scm_mean_dtend_temp_phys[:,it_scmd],ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.xlabel("(K/s)*(1e5)")
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        plt.subplot(2,4,7)
        plt.plot(np.zeros((nlev),dtype=float),ufs_mean_p[:,it_scmi]*0.01,color='black')
        plt.plot(ufs_mean_dtend_temp_lw[:,it_scmi] - scm_mean_dtend_temp_lw[:,it_scmd],ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.xlabel("(K/s)*(1e5)")
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        plt.subplot(2,4,8)
        plt.plot(np.zeros((nlev),dtype=float),ufs_mean_p[:,it_scmi]*0.01,color='black')
        plt.plot(ufs_mean_dtend_temp_sw[:,it_scmi] - scm_mean_dtend_temp_sw[:,it_scmd],ufs_mean_p[:,it_scmi]*0.01,color='blue')
        plt.ylim(1000,100)
        plt.xlim(-5,5)
        plt.xlabel("(K/s)*(1e5)")
        ax = plt.gca()
        ax.axes.yaxis.set_ticklabels([])
        #
        fileOUT = 'scm_ufs.'+namelist+'.'+str(cases[ic]).zfill(2)+'hrdiag.'+str(ihour).zfill(2)+'fcst.png'
        print(fileOUT)
        plt.savefig(fileOUT)
        #plt.show()

