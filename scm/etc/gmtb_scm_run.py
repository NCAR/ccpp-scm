#!/usr/bin/env python

from subprocess import Popen, PIPE
import time

proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

### User-editable section -- See https://theiadocs.rdhpcs.noaa.gov/wikis/theiadocs/doku.php?id=running_and_monitoring_jobs for appropriate values ###
job_name = "test_job"
account = "gmtb"
walltime = "0:20:00"
processors = "procs=1"
queue = "batch"
command = "./gmtb_scm_ens.py"
email_addr = "Grant.Firl@noaa.gov"
email_occasion = "abe"
serial_mem = "512M"
### End User-editable section ###


job_string = """#!/bin/bash
#PBS -N %s
#PBS -A %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o ./%s.out
#PBS -e ./%s.err
#PBS -q %s
#PBS -M %s
#PBS -m %s
#PBS -l vmem=%s
#PBS -V
cd $PBS_O_WORKDIR
module load intel
module load netcdf
%s""" % (job_name, account, walltime, processors, job_name, job_name, queue, email_addr, email_occasion, serial_mem, command)

# Send job_string to qsub
proc.stdin.write(job_string)
out, err = proc.communicate()

# Print your job and the system response to the screen as it's submitted
print job_string
print out
print err

time.sleep(0.1)
