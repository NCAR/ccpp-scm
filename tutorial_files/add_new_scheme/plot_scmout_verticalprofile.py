#!/usr/bin/env python
######################################################################
# Dependencies
######################################################################
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

if os.getenv("SCM_WORK"):
    dirSCM = os.getenv("SCM_WORK")
else:
    print("ERROR: System variable $SCM_WORK not set")
    exit()

######################################################################
# Plot configuration
######################################################################

# Where is data located?
datadir = dirSCM+"/ccpp-scm/scm/run/"

# What variable(s) to compare?
vars  = ["dT_dt_phys"]
nvars = len(vars)

# What hours to include in temporal mean?
t0 = 0
tf = 150

# Which case(s) to compare?
caseORG = "twpice"
caseNEW = "twpice"

# Which suite(s) to compare?
suiteORG = "SCM_RRFS_v1beta"
suiteNEW = "SCM_RRFS_v1beta_sas_sfcmod"

######################################################################
# No changes below
######################################################################
# Data location
dirORG = datadir + "output_"+caseORG+"_"+suiteORG+"/"
dirNEW = datadir + "output_"+caseNEW+"_"+suiteNEW+"/"

# Create plot for each variable requested.
for var in vars:

    # Read in data
    dataset = xr.open_dataset(dirORG+"output.nc")
    layer   = dataset.pres.values[1::,:]/100.
    var_org = dataset[var]
    dataset = xr.open_dataset(dirNEW+"output.nc")
    var_new = dataset[var]

    # Get variable attributes
    description = dataset[var].description
    units       = dataset[var].units

    # Convert tendencies from (K s-1) to (K day-1)
    if (units == "K s-1"):
        var_org = var_org*3600.*24.
        var_new = var_new*3600.*24.
        units   = "K day-1"

    # Compute mean profile from t0->tf
    var_org = np.mean(var_org[t0:tf,:],axis=0)
    var_new = np.mean(var_new[t0:tf,:],axis=0)
    layer   = np.mean(layer[t0:tf,:],axis=0)

    # Plot data
    fig  = plt.figure(figsize=(10,8))
    minX = min(np.concatenate((var_org,var_new)))
    maxX = max(np.concatenate((var_org,var_new)))
    #
    plt.plot(var_org, layer, color='red')
    plt.plot(var_new, layer, color='green')
    plt.title(var+" ("+str(t0)+"-"+str(tf)+"z)")
    plt.ylabel("(hPa)")
    plt.xlabel(description+" ("+units+")")
    plt.ylim([1000,0])
    plt.xlim([minX,maxX])
    plt.legend([suiteORG,suiteNEW],loc='upper right',fontsize='x-small')

    #
    #plt.show()
    plotOUT = 'cmp.'+suiteORG+'.'+suiteNEW+'.'+var+'.png'
    plt.savefig(plotOUT)
