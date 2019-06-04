#!/usr/bin/env python

from subprocess import Popen, PIPE
import time

proc = Popen('sbatch', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)

### User-editable section -- See https://theiadocs.rdhpcs.noaa.gov/wikis/theiadocs/doku.php?id=running_and_monitoring_jobs for appropriate values ###
job_name = "test_job"
account = "gmtb"
walltime = "20"
processors = "1"
queue = "batch"
command = "./run_gmtb_scm.py -c twpice"
email_addr = "Grant.Firl@noaa.gov"
email_occasion = "BEGIN,END,FAIL"
serial_mem = "512M"
working_dir = "../bin"
### End User-editable section ###


job_string = """#!/bin/bash
#SBATCH -J %s
#SBATCH -A %s
#SBATCH --time=%s
#SBATCH -n %s
#SBATCH -o ./%s.out
#SBATCH -e ./%s.err
#SBATCH -q %s
#SBATCH --mail-user=%s
#SBATCH --mail-type=%s
#SBATCH --mem=%s
#SBATCH -D%s
%s""" % (job_name, account, walltime, processors, job_name, job_name, queue, email_addr, email_occasion, serial_mem, working_dir, command)

# Send job_string to qsub
proc.stdin.write(job_string)
out, err = proc.communicate()

# Print your job and the system response to the screen as it's submitted
print job_string
print out
print err

time.sleep(0.1)
