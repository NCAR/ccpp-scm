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

#
parser = argparse.ArgumentParser()
parser.add_argument('-b',   '--build_type',  help='SCM build type')
parser.add_argument('-drt', '--dir_rt',      help='Directory containing SCM RT output')
parser.add_argument('-dbl', '--dir_bl',      help='Directory containing SCM RT baselines')

def parse_args():
    args       = parser.parse_args()
    build_type = args.build_type
    dir_rt     = args.dir_rt 
    dir_bl     = args.dir_bl

    return (build_type,dir_rt,dir_bl)

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
        else:
            if not exists(file_rt):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from output")
            if not exists(file_bl):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from baseline")
            error_count = error_count + 1

    #
    if error_count == 0:
        print("ALL TESTS PASSED, OUTPUT IS IDENTICAL.")
    else:
        print("ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE.")
        #1/0

#
if __name__ == '__main__':
    main()
