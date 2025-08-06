#!/usr/bin/env python
"""
-----------------------------------------------------------------------
 Description:  Example script to submit a list of jobs through the HPC batch system.

 Usage: ./scm_slurm_example.py -d ${SCM_ROOT} -a <ACCOUNT> -c [case or list of cases]  # from the bin directory
 Usage: ./scm_slurm_example.py -d ${PWD}/../../ -a gmtb -c twpice bomex # from the bin directory
-----------------------------------------------------------------------
"""

from subprocess import Popen, PIPE
import os
import time
import argparse

###############################################################################
# Argument list
###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-c',  '--case_name',  help='name of CCPP SCM case',          default='twpice', nargs='+')
parser.add_argument('-a',  '--account',    help='account name for job scheduler', default='gmtb')
parser.add_argument('-d',  '--scm_dir',    help='CCPP SCM directory',             required=True)

###############################################################################
# Main program
###############################################################################
def main():
    # Get command line arguments
    args     = parser.parse_args()
    cases    = args.case_name
    account  = args.account
    scm_dir  = args.scm_dir

    if (scm_dir is None):
        print("ERROR: No scm diretory <scm_dir> provided")
        exit()
    # endif

    # Process arguments
    ncase = len(cases)
    print(str(ncase)," cases provided: ", cases)

    for scm_case in cases:

        run_dir = scm_dir+"/ccpp-scm/scm/run_"+scm_case
        
        # Open process
        PROC = Popen('sbatch', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

        ### User-editable section ### 
        JOB_NAME = "test_"+scm_case
        ACCOUNT = account
        WALLTIME = "20"
        PROCESSORS = "1"
        QUEUE = "batch"
        COMMAND = "./run_scm.py -c " + scm_case + " --run_dir "+run_dir
        SERIAL_MEM = "512M"
        WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
        ### End User-editable section ###

        # Create job script for current case
        JOB_STRING = """#!/bin/bash
#SBATCH -J %s
#SBATCH -A %s
#SBATCH --time=%s
#SBATCH -n %s
#SBATCH -o ./%s.out
#SBATCH -q %s
#SBATCH --mem=%s
#SBATCH -D%s
%s""" % (JOB_NAME, ACCOUNT, WALLTIME, PROCESSORS, JOB_NAME, QUEUE,
         SERIAL_MEM, WORKING_DIR, COMMAND)

        # Send job script (JOB_STRING) to qsub
        PROC.stdin.write(JOB_STRING.encode())
        OUT, ERR = PROC.communicate()

        # Print your job and the system response to the screen as it's submitted
        print(JOB_STRING)
        print(OUT.decode())
        print(ERR.decode())

        time.sleep(0.1)
    # end for (scm case loop)
# end def

###############################################################################
if __name__ == '__main__':
    main()
