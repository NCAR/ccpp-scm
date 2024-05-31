#!/usr/bin/env python

import argparse
##############################################################################
#
# This script compares SCM RT output to baselines.
#
##############################################################################
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
#
def plot_results(file_BL, file_RT, plot_RT_only):
    # List of SCM output fields to plot (This is everything)
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
    # Use subset of available SCM output for plots.
    vars2plot = ["qv","T","u","v","ql","qi","qc","sfc_dwn_sw","sfc_up_sw",               \
                 "sfc_net_sw","sfc_dwn_lw", "u10m","v10m","hpbl","gflux",                \
                 "qv_force_tend","T_force_tend","u_force_tend","v_force_tend","w_ls",    \
                 "u_g","v_g","dT_dt_rad_forc","h_advec_thil","h_advec_qt",               \
                 "v_advec_thil","v_advec_qt","T_s","lhf","shf","tprcp_inst",             \
                 "tprcp_rate_inst","t2m","q2m","ustar","tsfc","tau_u","tau_v"]

    # Open SCM datasets
    SCM_BL = Dataset(file_BL)
    SCM_RT = Dataset(file_RT)

    plot_diff = True
    if plot_RT_only:
        plot_diff = False
    # end if
    
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
                if (np.size(x1) > 1):
                    fig = plt.figure(figsize=(13,10))
                    # Baselines and SCM RTs on same plot
                    if plot_diff:
                        plt.subplot(2,1,1)
                    # end if
                    plt.title(SCM_BL[var].description)
                    if plot_RT_only:
                        plt.plot(x1, y1,  color='blue')
                    else:
                        plt.plot(x1, y1,  color='blue')
                        plt.plot(x2, y2,  color='black')
                    # end if
                    plt.ylabel('('+SCM_BL[var].units+')')
                    plt.xlabel('(hours)')
                    
                    # Difference (Baseline-SCMRT)
                    if plot_diff:
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

                # Finally, make figure.
                if (np.size(x1) > 1):
                    fig = plt.figure(figsize=(13,10))
                    if plot_RT_only:
                        plt.contourf(x2, y2, z2, 20, cmap='YlGnBu')
                        plt.ylim(1000,200)
                        plt.ylabel('(Pa)')
                        plt.xlabel('(hours)')
                        cbr = plt.colorbar()
                        cbr.set_label('('+SCM_RT[var].units+')')
                    else:
                        # Baselines
                        plt.subplot(3,1,1)
                        plt.title(SCM_BL[var].description, fontsize=12)
                        plt.contourf(x1, y1, z1, 20, cmap='YlGnBu')
                        plt.ylim(1000,200)
                        plt.ylabel('(Pa)')
                        cbr = plt.colorbar()
                        cbr.set_label('('+SCM_RT[var].units+')')
                        # SCM RTs
                        plt.subplot(3,1,2)
                        plt.contourf(x2, y2, z2, 20, cmap='YlGnBu')
                        plt.ylim(1000,200)
                        plt.ylabel('(Pa)')
                        plt.xlabel('(hours)')
                        cbr = plt.colorbar()
                        cbr.set_label('('+SCM_RT[var].units+')')
                        # Only plot differences if they exist.
                        if np.count_nonzero(dz) > 0:
                            plt.subplot(3,1,3)
                            plt.title("Difference (top - middle)", fontsize=8)
                            plt.contourf(x2, y2, dz, 20, cmap='bwr')
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
                # end if (Have enought pts to plot?)
            # end if     (fields exist?)
        # end if         (field requested?)
    # end for            (fields in file)

    return(plot_files)

if __name__ == "__main__":

    #Parse arguments
    parser = argparse.ArgumentParser(
                     description="Script for setting up a forecast and creating a workflow"\
                     "according to the parameters specified in the config file\n")

    parser.add_argument('-bl', '--baseline', type=str, required=True,
                        help='Baseline file')
    parser.add_argument('-rt', '--rt', type=str, required=True,
                        help='Regression test file')
    pargs = parser.parse_args()

    plot_results(pargs.baseline, pargs.rt, False)
