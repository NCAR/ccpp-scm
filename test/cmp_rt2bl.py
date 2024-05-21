#!/usr/bin/env python

##############################################################################
#
# This script compares SCM RT output to baselines.
#
##############################################################################
import os
import sys
from rt_test_cases import run_list
from os.path import exists
import argparse
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#
parser = argparse.ArgumentParser()
parser.add_argument('-b',   '--build_type',  help='SCM build type')
parser.add_argument('-drt', '--dir_rt',      help='Directory containing SCM RT output')
parser.add_argument('-dbl', '--dir_bl',      help='Directory containing SCM RT baselines')

#
def parse_args():
    args       = parser.parse_args()
    build_type = args.build_type
    dir_rt     = args.dir_rt 
    dir_bl     = args.dir_bl

    return (build_type,dir_rt,dir_bl)

#
def plot_results(file_BL,file_RT):
    # List of SCM output fields to plot (This everything)
    vars2plot = ["time_inst","time_diag","time_swrad","time_lwrad","time_rad","pres",    \
                 "pres_i","sigma","sigma_i","pres_s","qv","T","u","v","ql","qi","qc",    \
                 "qv_force_tend","T_force_tend","u_force_tend","v_force_tend","w_ls",    \
                 "u_g","v_g","dT_dt_rad_forc","h_advec_thil","h_advec_qt",               \
                 "v_advec_thil","v_advec_qt","T_s","lhf","shf","tprcp_inst",             \
                 "tprcp_rate_inst","t2m","q2m","ustar","tsfc","tau_u","tau_v","upd_mf",  \
                 "dwn_mf","det_mf","sfc_up_lw_land","sfc_up_lw_ice","sfc_up_lw_water",   \
                 "sfc_up_sw_dir_nir","sfc_up_sw_dif_nir","sfc_up_sw_dir_vis",            \
                 "sfc_up_sw_dif_vis","sfc_dwn_sw_dir_nir","sfc_dwn_sw_dif_nir",          \
                 "sfc_dwn_sw_dir_vis","sfc_dwn_sw_dif_vis","mp_prcp_inst",               \
                 "dcnv_prcp_inst","scnv_prcp_inst","pwat","dT_dt_lwrad","dT_dt_swrad",   \
                 "dT_dt_pbl","dT_dt_deepconv","dT_dt_shalconv","dT_dt_micro",            \
                 "dT_dt_ogwd","dT_dt_cgwd","dT_dt_phys","dT_dt_nonphys","dq_dt_pbl",     \
                 "dq_dt_deepconv","dq_dt_shalconv","dq_dt_micro","dq_dt_phys",           \
                 "dq_dt_nonphys","doz_dt_pbl","doz_dt_prodloss","doz_dt_oz","doz_dt_T",  \
                 "doz_dt_ovhd","doz_dt_phys","doz_dt_nonphys","du_dt_pbl","du_dt_ogwd",  \
                 "du_dt_deepconv","du_dt_cgwd","du_dt_shalconv","du_dt_phys",            \
                 "du_dt_nonphys","dv_dt_pbl","dv_dt_ogwd","dv_dt_deepconv","dv_dt_cgwd", \
                 "dv_dt_shalconv","dv_dt_phys","dv_dt_nonphys","sfc_dwn_sw","sfc_up_sw", \
                 "sfc_net_sw","sfc_dwn_lw","gflux","u10m","v10m","hpbl","tprcp_accum",   \
                 "ice_accum","snow_accum","graupel_accum","conv_prcp_accum",             \
                 "tprcp_rate_accum","ice_rate_accum","snow_rate_accum",                  \
                 "graupel_rate_accum","conv_prcp_rate_accum","max_cloud_fraction",       \
                 "toa_total_albedo","vert_int_lwp_mp","vert_int_iwp_mp",                 \
                 "vert_int_lwp_cf","vert_int_iwp_cf"]

    # Open SCM datasets
    SCM_BL = Dataset(file_BL)
    SCM_RT = Dataset(file_RT)

    plot_files = []
    for var in SCM_BL.variables.keys():
        if (var in vars2plot):
            # Handle temporal axis.
            # There are 4 different dimensions in the SCM output, identified by the suffix "_dim".
            # Here the suffix is stripped and used to identify the temporal dimenesion (index 0 in netcdf file)
            timeD = SCM_BL[var].dimensions[0]
            timeD = timeD[0:len(timeD)-4]
            x1    = SCM_BL[timeD][:].squeeze()/3600. #seconds -> hours
            x2    = SCM_RT[timeD][:].squeeze()/3600. #seconds - >hours
            # Is this a 2D (time, x) variable? (Really 1D since x=1 in SCM)
            is2D  = False
            if (len(SCM_BL[var].dimensions)==2):
                is2D  = True
            # endif
            # one/two-dimensional variables
            if (len(SCM_BL[var].shape) != 3):
                if (is2D):
                    y1 = SCM_BL[var][:,0].squeeze()
                    y2 = SCM_RT[var][:,0].squeeze()
                else:
                    y1 = SCM_BL[var][:]
                    y2 = SCM_RT[var][:]
                # endif
                
                # Make figure
                fig = plt.figure(figsize=(13,10))
                # Baselines and SCM RTs on same plot
                plt.subplot(2,1,1)
                plt.title(SCM_BL[var].description)
                plt.plot(x1, y1,  color='blue')
                plt.plot(x2, y2,  color='black')
                plt.ylabel('('+SCM_BL[var].units+')')
                plt.xlabel('(hours)')
                # Difference (Baseline-SCMRT)
                plt.subplot(2,1,2)
                plt.title("Difference (blue - black)")
                plt.plot(x1, y1 - y2,  color='red')
                plt.plot(x1, np.zeros(len(y1)), color='grey',linestyle='dashed')
                plt.ylabel('('+SCM_RT[var].units+')')
                plt.xlabel('(hours)')
                #
                fileOUT = 'scm.' + var +'.png'
                plt.savefig(fileOUT)
                #
                plot_files.append(fileOUT)
            # three-dimensional variables
            elif len(SCM_BL[var].shape) == 3:
                z1 = np.transpose(SCM_BL[var][:,:,0]).squeeze()
                z2 = np.transpose(SCM_RT[var][:,:,0]).squeeze()

                # vertical axis
                y1 = SCM_BL["pres"][0,:].squeeze()*0.01
                y2 = SCM_RT["pres"][0,:].squeeze()*0.01
                nlev = SCM_BL[var][:,:,0].shape[1]
                # Layer (nlev) quantities are the default, so check for case where we have an
                # interface (nlev+1) variable to plot.
                if (SCM_BL[var][:,:,0].shape[1] > len(y1)):
                    y1 = SCM_BL["pres_i"][0,:].squeeze()*0.01
                    y2 = SCM_RT["pres_i"][0,:].squeeze()*0.01
                # endif

                # Comppute differences and determine valid plot range(s).
                dz = z1-z2
                if np.min(z1) != np.max(z1):
                    clev  = np.arange(np.min(z1),np.max(z1),(np.max(z1)-np.min(z1))*0.05)
                    if np.count_nonzero(dz) > 0:
                        clevd = np.arange(np.min(dz),np.max(dz),(np.max(dz)-np.min(dz))*0.05)
                    else:
                        clevd = 0
                    # end if
                else:
                    clev  = 0
                    clevd = 0
                # end if

                # Finally, make figure.
                fig = plt.figure(figsize=(13,10))
                # Baselines
                plt.subplot(3,1,1)
                plt.title(SCM_BL[var].description, fontsize=12)
                print("var: ",var)
                print("x1:  ",x1)
                print("y1:  ",y1)
                print("z1:  ",z1)

                plt.contourf(x1, y1, z1, clev, cmap='YlGnBu')
                plt.ylim(1000,200)
                plt.ylabel('(Pa)')
                cbr = plt.colorbar()
                cbr.set_label('('+SCM_RT[var].units+')')
                # SCM RTs
                plt.subplot(3,1,2)
                plt.contourf(x2, y2, z2, clev, cmap='YlGnBu')
                plt.ylim(1000,200)
                plt.ylabel('(Pa)')
                plt.xlabel('(hours)')
                cbr = plt.colorbar()
                cbr.set_label('('+SCM_RT[var].units+')')
                # Only plot differences if they exist.
                if np.count_nonzero(dz) > 0:
                    plt.subplot(3,1,3)
                    plt.title("Difference (top - middle)", fontsize=8)
                    plt.contourf(x2, y2, dz, clevd, cmap='bwr')
                    plt.ylim(1000,200)
                    plt.ylabel('(Pa)')
                    plt.xlabel('(hours)')
                    cbr = plt.colorbar()
                    cbr.set_label('('+SCM_RT[var].units+')')
                # end if (no differences exist)
                fileOUT = 'scm.' + var +'.png'
                plt.savefig(fileOUT)
                #
                plot_files.append(fileOUT)
            # end if (fields exist?)
        # end if     (field requested?)
    # end for        (fields in file)

    return(plot_files)

#
def main():
    #
    (build_type, dir_rt, dir_bl) = parse_args()

    #
    error_count = 0
    for run in run_list:
        file_rt = dir_rt + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        file_bl = dir_bl + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        if exists(file_rt) and exists(file_bl):
            com = "cmp "+file_rt+" "+file_bl+" > logfile.txt"
            result = os.system(com)
            if (result != 0):
                print("Output for "+run["case"]+"_"+run["suite"]+ " DIFFERS from baseline")
                error_count = error_count + 1
            else:
                print("Output for "+run["case"]+"_"+run["suite"]+ " is IDENTICAL to baseline")
            # end if

            # Create plots between RTs and baselines
            plot_files = plot_results(file_rt,file_bl)

            # Setup output directories for plots.
            result = os.system("mkdir -p scm_rt_out/"+run["case"]+"/"+run["suite"])

            # Archive plots.
            com = "mv"
            for plot_file in plot_files:
                com = com + " " + plot_file
            # end if
            com = com + " scm_rt_out/" + run["case"] + "/" + run["suite"]
            result = os.system(com)
        else:
            if not exists(file_rt):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from output")
            # end if
            if not exists(file_bl):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from baseline")
            # end if
            error_count = error_count + 1
        # end if
    # end for

    #
    if error_count == 0:
        print("ALL TESTS PASSED, OUTPUT IS IDENTICAL.")
    else:
        print("ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE.")
        #1/0
    # end if

#
if __name__ == '__main__':
    main()
