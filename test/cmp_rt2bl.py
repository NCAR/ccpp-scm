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
from plot_scm_out import plot_results

#
parser = argparse.ArgumentParser()
parser.add_argument('-drt',  '--dir_rt',  help='Directory containing SCM RT output',      required=True)
parser.add_argument('-dbl',  '--dir_bl',  help='Directory containing SCM RT baselines',   required=True)
parser.add_argument('-prt',  '--plt_rt',  help='If true, create plots of dir_rt',         action='store_true')
parser.add_argument('-pbl',  '--plt_bl',  help='If true, create plots of dir_bl',         action='store_true')
parser.add_argument('-pall', '--plt_all', help='If true, create plots of all SCM fields', action='store_true')
parser.add_argument('-dbg',  '--debug',   help='Debug mode ',                             action='store_true')

#
def parse_args():
    args    = parser.parse_args()
    dir_rt  = args.dir_rt 
    dir_bl  = args.dir_bl
    plt_rt  = args.plt_rt
    plt_bl  = args.plt_bl
    plt_all = args.plt_all
    debug   = args.debug
    return (dir_rt, dir_bl, plt_rt, plt_bl, plt_all, debug)

#
def main():
    #
    (dir_rt, dir_bl, plt_rt, plt_bl, plt_all, debug) = parse_args()

    #
    error_count = 0
    for run in run_list:
        file_rt = dir_rt + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        file_bl = dir_bl + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        if exists(file_rt) and exists(file_bl):
            com = "cmp "+file_rt+" "+file_bl+" > logfile.txt"
            result = os.system(com)
            if (result != 0):
                print("Output for "+run["case"]+"_"+run["suite"]+ " DIFFERS from baseline. Difference plots will be created")
                error_count = error_count + 1
            else:
                print("Output for "+run["case"]+"_"+run["suite"]+ " is IDENTICAL to baseline")
            # end if

            # Create plots between RTs and baselines (only if differences exist)
            if (result != 0):
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
