#!/usr/bin/env python
#################################################################################################
# Dependencies
#################################################################################################
import os
import xarray as xr
import matplotlib.pyplot as plt
import argparse as ap

#################################################################################################
def read_SCMout2d(fileIN,vars2plt):
    #
    data = xr.open_dataset(fileIN)

    # Get list of available variables.
    varsavail = list(data.keys())

    # Check that requested variables (vars2plt) exist.
    for var2plt in vars2plt:
        if (not var2plt in varsavail):
            print('ERROR: ',var2plt, ' is not in ',fileIN,'.')
            exit()

    # Read in data, store in python dictonary "state", return.
    state={}
    state["time"] = data["time_inst"].values/3600.  # sec -> hour
    state["pres"] = data["pres"].values[:,:,0]*0.01 # Pa -> hPa
    for var2plt in vars2plt:
        state[var2plt]          = data[var2plt]
        state[var2plt+"_units"] = data[var2plt].units
        state[var2plt+"_time"]  = data[data[var2plt].dims[0][0:len(data[var2plt].dims[0])-4]].values/3600.

    #
    return (state)

#################################################################################################
# Argument list
# ./plt_scmout_2d.py -n twpice -sdf SCM_GFS_v17_p8 -nmls input_GFS_v17_p8_simA input_GFS_v17_p8_simB -vars sfc_dwn_sw sfc_up_sw sfc_dwn_lw
#################################################################################################
parser = ap.ArgumentParser()
parser.add_argument('-n',    '--case_name', help='name of case', required=True)
parser.add_argument('-sdf',  '--suite',     help='CCPP suite definition file',required=True)
parser.add_argument('-nmls', '--nml_list',  help='namelists, separated by a space', nargs='*')
parser.add_argument('-vars', '--var_list',  help='varaibles to plot, separated by a space', nargs='*', required=True)

def main():

    #############################################################################################
    # Get arguments
    #############################################################################################
    args      = parser.parse_args()
    case_name = args.case_name
    suite     = args.suite
    have_nml  = True
    if (args.nml_list):
        namelist = args.nml_list
    else:
        have_nml = False
        namelist = ['']
    vars2plt  = args.var_list

    #############################################################################################
    # Get data
    #############################################################################################

    # Get files to plot
    nfiles  = 0
    fileSCM = []
    for nml in namelist:
        if have_nml:
            fileSCM.append("../../../../scm/run/output_" + case_name + "_" + suite+"_"+nml+"/output.nc")
        else:
            fileSCM.append("../../../../scm/run/output_" + case_name + "_" + suite+"/output.nc")

        # Check that file exists (Exit if not)
        if (not os.path.exists(fileSCM[nfiles])):
            print('ERROR: ',fileSCM[nfiles],' does not exist.')
            exit()
        nfiles = nfiles + 1

    # Read in data
    stateSCM = []
    for ifile in range(0,nfiles):
        stateSCM.append(read_SCMout2d(fileSCM[ifile],vars2plt))

    #############################################################################################
    # Make plots
    #############################################################################################
    fontsize = 7
    fig      = plt.figure(figsize=(10,7))

    colors = ["black","red","turquoise","magenta","green"]
    colors.append(colors)

    # Plot full temporal range, unless tlimit is provided (ToDo)
    tlimit   = max(stateSCM[0]["time"])

    #
    varcount = 1
    for var2plt in vars2plt:

        # Top row) Absolute plot
        if (nfiles  > 1): plt.subplot(2,len(vars2plt),varcount)
        if (nfiles == 1): plt.subplot(1,len(vars2plt),varcount)
        plt.title(var2plt, fontsize=fontsize*1.5)
        ccount = 0
        for ifile in range(0,nfiles):
            plt.plot(stateSCM[ifile][var2plt+"_time"], stateSCM[ifile][var2plt][:], color=colors[ccount])
            ccount = ccount + 1
        plt.xlim(0,tlimit)
        plt.yticks(fontsize=fontsize)
        plt.ylabel("("+stateSCM[ifile][var2plt+"_units"]+")")

        # Bottom row) Difference plot (assumes first element is reference)
        if (nfiles > 1):
            plt.subplot(2,len(vars2plt),len(vars2plt)+varcount)
            ccount = 1
            for ifile in range(1,nfiles):
                plt.plot(stateSCM[ifile][var2plt+"_time"], stateSCM[ifile][var2plt][:] - stateSCM[0][var2plt][:], colors[ccount])
                plt.plot([0,tlimit],[0,0],color='grey',linestyle='--')
                ccount = ccount+ 1
            plt.xlim(0,tlimit)
            ylimit = max(stateSCM[ifile][var2plt][:] - stateSCM[0][var2plt][:], key=abs)*1.05
            plt.ylim(-ylimit,ylimit)
            plt.yticks(fontsize=fontsize)
            plt.ylabel("("+stateSCM[ifile][var2plt+"_units"]+")")
            plt.xlabel("(hours)")
        #
        varcount = varcount + 1
    plt.show()

if __name__ == '__main__':
    main()
