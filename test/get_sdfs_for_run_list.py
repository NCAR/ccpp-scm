#!/usr/bin/env python

##############################################################################
#
# This script gets SDFs needed for a run list.
#
##############################################################################
import os
import sys
import rt_test_cases
from os.path import exists
import argparse

#
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--sdf_list',   help='List of SCM SDFs')

def parse_args():
    args       = parser.parse_args()
    sdf_list   = args.sdf_list
    return (sdf_list)

def main():
    
    (sdf_list_name) = parse_args()
    sdf_list = getattr(rt_test_cases, sdf_list_name)
    list_out = ''
    for count,sdf in enumerate(sdf_list):
        if (count < len(sdf_list)-1):
            list_out = list_out+sdf+','
        else:
            list_out = list_out+sdf
        # endif
    # end for
    print(list_out)
# end def

if __name__ == '__main__':
    main()
