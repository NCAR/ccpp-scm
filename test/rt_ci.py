#!/usr/bin/env python

import os
import sys
from rt_test_cases import run_list
from os.path import exists

#
for run in run_list:
    file_out = "../scm/run/" + "output_"+run["case"]+"_"+run["suite"]+"/output.nc"
    if exists(file_out):
        os.system("mkdir -p artifact/"+run["case"]+"_"+run["suite"]+"/")
        os.system("cp " + file_out + " artifact/"+run["case"]+"_"+run["suite"]+"/output.nc")
    else:
        print("FAIL:  Could not copy output to for baseline generation")
        exit()
