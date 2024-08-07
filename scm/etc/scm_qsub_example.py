#!/usr/bin/env python
"""
-----------------------------------------------------------------------
 Description:  Example script to submit a job through the HPC batch system

 Assumptions: For use on Derecho. This script must be copied to the bin directory.
              The user should edit the JOB_NAME, ACCOUNT, etc.
              Before running this script the environment must be setup.
              The following are instructions for that:
              $ module purge
              $ module use scm/etc/modules
              $ module load derecho_gnu
                or
              $ module load derecho_intel

 Command line arguments: none

 Usage: ./scm_qsub_example.py   # from the bin directory
-----------------------------------------------------------------------
"""

from subprocess import Popen, PIPE
import os
import time

USER = os.getenv('USER')
STRINGS = [USER, 'ucar.edu']
MY_EMAIL = '@'.join(STRINGS)
PROC = Popen('qsub -V', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

### Begin User-editable section ###
JOB_NAME = "test_job"
ACCOUNT = "ENTER_ACCOUNT_NUMBER"
WALLTIME = "walltime=00:20:00"
PROCESSORS = "select=1:ncpus=1"
QUEUE = "develop"
COMMAND = "./run_scm.py -c twpice"
EMAIL_ADDR = MY_EMAIL
SERIAL_MEM = "512M"
WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
### End User-editable section ###


JOB_STRING = """#!/bin/bash
#PBS -N %s
#PBS -A %s
#PBS -l %s
#PBS -q %s
#PBS -j oe
#PBS -k eod
#PBS -l %s
#PBS -m abe
#PBS -M %s
%s""" % (JOB_NAME, ACCOUNT, WALLTIME, QUEUE, PROCESSORS, EMAIL_ADDR, COMMAND)

# Send job_string to qsub
PROC.stdin.write(JOB_STRING.encode())
OUT, ERR = PROC.communicate()

# Print your job and the system response to the screen as it's submitted
print(JOB_STRING)
print(OUT.decode())
print(ERR.decode())

time.sleep(0.1)
