#!/usr/bin/env python

##############################################################################
#
# This script gathers the output from the $SCM_ROOT/ccpp-scm/scm/run directory
# for the SCM CI Regression Tests. 
#
##############################################################################
import os
import sys
from rt_test_cases
from os.path import exists
import argparse

#
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--build_type', help='SCM build type')
parser.add_argument('-s', '--sdfs',       help='SCM SDFs and cases')

def parse_args():
    args       = parser.parse_args()
    build_type = args.build_type
    sdfs       = args.sdfs
    return (build_type, sdfs)

def main():
    
    (build_type, sdfs) = parse_args()

    if (sdfs == 'supported'): run_list = run_list_supported
    if (sdfs == 'legacy'):    run_list = run_list_legacy
    if (sdfs == 'dev'):       run_list = run_list_dev
    if (sdfs == 'sp'):        run_list = run_list_sp
    if (sdfs == 'nvhpc'):     run_list = run_list_nvhpc
    #
    errmsgs=[]
    for run in run_list:
        case_tag = run["case"]+"_"+run["suite"]
        file_out = "../scm/run/" + "output_"+case_tag+"/output.nc"
        if exists(file_out):
            os.system("mkdir -p artifact-"+build_type+"/"+case_tag+"/")
            os.system("cp " + file_out + " artifact-"+build_type+"/"+case_tag+"/output.nc")
        else:
            errmsgs.append("Could not copy output for baseline generation "+ case_tag)
        # end if
    # end for

    for errmsg in errmsgs:
        print(errmsg)
    # end for

if __name__ == '__main__':
    main()
