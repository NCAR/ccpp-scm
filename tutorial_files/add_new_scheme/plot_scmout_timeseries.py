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
vars  =["lhf","shf"]
nvars = len(vars)

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
    time    = dataset.time_inst.values/3600.
    var_org = dataset[var]
    dataset = xr.open_dataset(dirNEW+"output.nc")
    var_new = dataset[var]

    # Get variable attributes
    description = dataset[var].description
    units       = dataset[var].units
    
    # Plot data
    fig = plt.figure(figsize=(10,8))
    #
    plt.plot(time,var_org,color='red')
    plt.plot(time,var_new,color='green')
    plt.ylabel(description+" ("+units+")")
    plt.xlabel("Forecast Time (hours)")
    plt.ylim([0,max(var_org)])
    plt.xlim([20,150])
    plt.legend([suiteORG,suiteNEW],loc='upper left',fontsize='x-small')

    #
    #plt.show()
    plotOUT = 'cmp.'+suiteORG+'.'+suiteNEW+'.'+var+'.png'
    plt.savefig(plotOUT)
