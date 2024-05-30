#!/usr/bin/env python

##############################################################################
#
# This script compares SCM RT output to baselines.
#
##############################################################################
import os
import sys
from os.path import exists
import argparse
from plot_scm_out import plot_results

#
parser = argparse.ArgumentParser()
parser.add_argument('-fbl',  '--file_bl', help='File containing SCM RT baselines', required=True)
parser.add_argument('-frt',  '--file_rt', help='File containing SCM RT output')

#
def parse_args():
    args    = parser.parse_args()
    file_rt = args.file_rt
    file_bl = args.file_bl
    return (file_bl, file_rt)

#
def main():
    #
    (file_bl, file_rt) = parse_args()

    plot_files = plot_results(file_bl, file_rt)

    # Put plots in local directory (scm_plots/)
    result = os.system("mkdir -p scm_plots/")
    com = "mv"
    for plot_file in plot_files:
        com = com + " " + plot_file
    # end for
    result = os.system(com+" scm_plots/")

    # Housekeeping
    if file_rt is not None:
        result = os.system('tar -cvf scm_out_abs.tar scm_plots/*.png')
        result = os.system('tar -cvf /scratch1/data_untrusted/Dustin.Swales/scm_out_abs.tar scm_plots/*.png')
    else:
        result = os.system('tar -cvf scm_out_diff.tar scm_plots/*.png')
        result = os.system('tar -cvf /scratch1/data_untrusted/Dustin.Swales/scm_out_diff.tar scm_plots/*.png')
    # end if
    result = os.system('rm -rf scm_plots/')
#
if __name__ == '__main__':
    main()
