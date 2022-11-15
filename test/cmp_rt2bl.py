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
    for run in run_list:
        file_rt = dir_rt + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        file_bl = dir_bl + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        print(file_rt)
        print(file_bl)
        if exists(file_rt) and exists(file_bl):
            com = "cmp "+file_rt+" "+file_bl
            result = os.system(com)
            print(com)
            if (result != 0):
                print("ERROR: Results do not agree with baselines for"+run["case"]+"_"+run["suite"])
        else:
            print("FAIL: "+file_rt+" "+file_bl)
            exit()

#
if __name__ == '__main__':
    main()
