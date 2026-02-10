#!/usr/bin/env python

##############################################################################
#
# This script compares SCM RT output to baselines.
#
##############################################################################
import os
import sys
from rt_test_cases_supported import run_list as run_list_supported
from rt_test_cases_legcy     import run_list as run_list_legacy
from rt_test_cases_dev       import run_list as run_list_dev
from os.path import exists
import argparse
from plot_scm_out import plot_results

#
parser = argparse.ArgumentParser()
parser.add_argument('-drt',  '--dir_rt',   help='Directory containing SCM RT output',              required=True)
parser.add_argument('-dbl',  '--dir_bl',   help='Directory containing SCM RT baselines',           required=True)
parser.add_argument('-np',   '--no_plots', help='flag to turn off generation of difference plots', required=False, action='store_true')
parser.add_argument('-s',    '--sdfs',     help='SCM SDFs and cases')
#
def parse_args():
    args      = parser.parse_args()
    dir_rt    = args.dir_rt 
    dir_bl    = args.dir_bl
    no_plots  = args.no_plots
    sdfs      = args.sdfs
    return (dir_rt, dir_bl, no_plots, sdfs)

#
def main():
    #
    (dir_rt, dir_bl, no_plots, sdfs) = parse_args()

    if (sdfs == 'supported'): run_list = run_list_supported
    if (sdfs == 'legacy'):    run_list = run_list_legacy
    if (sdfs == 'dev'):       run_list = run_list_dev
    #
    error_count = 0
    for run in run_list:
        file_rt = dir_rt + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        file_bl = dir_bl + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        if exists(file_rt) and exists(file_bl):
            com = "cmp "+file_rt+" "+file_bl+" > logfile.txt"
            result = os.system(com)
            if (result != 0):
                message = "Output for "+run["case"]+"_"+run["suite"]+ " DIFFERS from baseline."
                if (not no_plots):
                    message += " Difference plots will be created."
                print(message)
                error_count = error_count + 1
            else:
                print("Output for "+run["case"]+"_"+run["suite"]+ " is IDENTICAL to baseline")
            # end if

            # Create plots between RTs and baselines (only if differences exist)
            if (result != 0 and not no_plots):
                plot_files = plot_results(file_bl, file_rt)

                # Setup output directories for plots.
                result = os.system("mkdir -p scm_rt_out/"+run["case"]+"/"+run["suite"])

                # Archive plots.
                com = "mv"
                for plot_file in plot_files:
                    com = com + " " + plot_file
                # end if
                com = com + " scm_rt_out/" + run["case"] + "/" + run["suite"]
                result = os.system(com)
            # end if
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

    # Create tarball with plots.
    if (not no_plots):
        result = os.system('tar -cvf scm_rt_out.tar scm_rt_out/*')
    
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
