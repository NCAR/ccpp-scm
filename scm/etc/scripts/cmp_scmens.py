#!/usr/bin/env python
################################################################
# Dependencies
################################################################
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import argparse as ap

################################################################
def read_SCMout(fileIN):
    var_list = ["pres","T","qv","u","v"]
    data = xr.open_dataset(fileIN)

    state={}
    for var in var_list:
        state[var] = data[var].values[:,:,0]

    state["time"]=data["time_inst"].values
    return (state)

################################################################
def read_UFScomp(fileIN,dephy):
    var_list = ["pres","T","qv","u","v"]
    data = xr.open_dataset(fileIN)

    state = {}
    for var in var_list:
        if dephy:
            state[var] = data[var].values[:,:,0,0]
        else:
            if (var == "pres"):
                state[var] = data[var].values[:]
            else:
                state[var] = data[var].values[:,:]

    state["time"]=data["time"].values
    return (state)

################################################################
DEPHY = True
case_name = "fv3_SCM_ensemble_12hr_UFSforcing"
namelist  = "SCM_GFS_v16"
ensmember = 0

#
# Construct filenames from configurations information provided.
#
fileSCMout = "../../run/output_"+case_name+"_n"+str(ensmember).zfill(3)+"_"+namelist+"/output.nc"
fileUFSref = "../../data/comparison_data/"+case_name+"_n"+str(ensmember).zfill(3)+"_comp_data.nc"

#
# Read in data
#
stateSCM   = read_SCMout(fileSCMout)
stateUFS   = read_UFScomp(fileUFSref,DEPHY)

time = 0
#
# Make some plots
#

for ij in range(0,stateUFS["T"].shape[1]):
    print(ij,stateUFS["pres"][time,ij],stateSCM["pres"][time,ij])

# Plot State variables (T, qv, u, v)
for time in range(1,7):
    fig = plt.figure(figsize=(10,8))
    #
    plt.subplot(1,4,1)
    plt.plot(stateUFS["T"][time,:],stateUFS["pres"][time,:]*0.01,color='blue')
    plt.plot(stateSCM["T"][time,:],stateUFS["pres"][time,:]*0.01,color='green')
    plt.ylim(1000,0)
    plt.ylabel("P (hPa)")
    plt.xlabel("(K)")
    plt.title("Temperature")
    #
    plt.subplot(1,4,2)
    plt.plot(stateUFS["qv"][time,:],stateUFS["pres"][time,:]*0.01,color='blue')
    plt.plot(stateSCM["qv"][time,:],stateUFS["pres"][time,:]*0.01,color='green')
    plt.ylim(1000,0)
    plt.xlabel("kg/kg")
    plt.title("Specific Humidity")
    #
    plt.subplot(1,4,3)
    plt.plot(stateUFS["u"][time,:],stateUFS["pres"][time,:]*0.01,color='blue')
    plt.plot(stateSCM["u"][time,:],stateUFS["pres"][time,:]*0.01,color='green')
    plt.ylim(1000,0)
    plt.xlabel("m/s")
    plt.title("Zonal-wind")
    #
    plt.subplot(1,4,4)
    plt.plot(stateUFS["v"][time,:],stateUFS["pres"][time,:]*0.01,color='blue')
    plt.plot(stateSCM["v"][time,:],stateUFS["pres"][time,:]*0.01,color='green')
    plt.ylim(1000,0)
    plt.xlabel("m/s")
    plt.title("Meridional-wind")
    #
    plt.show()

