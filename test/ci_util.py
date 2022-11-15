#!/usr/bin/env python

##############################################################################
#
# This script gathers the output from the $SCM_ROOT/ccpp-scm/scm/run directory
# for the SCM CI Regression Tests. 
#
##############################################################################
import os
import sys
from rt_test_cases import run_list
from os.path import exists
import argparse

#
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--build_type',  help='SCM build type')

def parse_args():
    args       = parser.parse_args()
    build_type = args.build_type
    return (build_type)

def main():
    
    (build_type) = parse_args()

    #
    for run in run_list:
        file_out = "../scm/run/" + "output_"+run["case"]+"_"+run["suite"]+"/output.nc"
        if exists(file_out):
            os.system("mkdir -p artifact-"+build_type+"/"+run["case"]+"_"+run["suite"]+"/")
            os.system("cp " + file_out + " artifact-"+build_type+"/"+run["case"]+"_"+run["suite"]+"/output.nc")
        else:
            print("FAIL: Could not copy output to for baseline generation")
            exit()

if __name__ == '__main__':
    main()
