#!/usr/bin/env python

from __future__ import print_function
import shlex
import subprocess
import os
import shutil
import f90nml
from Queue import Queue
from threading import Thread
import sys
import multiprocessing as mp

def run_SCM(case_config_file):
    cmd = './gmtb_scm ' + case_config_file
    args = shlex.split(cmd)
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try: run_SCM(*args, **kwargs)
    except Exception as e:
        print("error: %s run(*%r, **%r)" % (e, args, kwargs))

def run(f):
    real = f[f.find('ens')+3:]
    print('Realization %s' % real)
    cmd = './gmtb_scm ' + f
    args = shlex.split(cmd)
    with open(os.devnull, "w") as f:
        subprocess.check_call(args, stdout=f, stderr=f)

def worker(queue):
    """Process files from the queue."""
    for args in iter(queue.get, None):
        try:
            run(*args)
        except Exception as e: # catch exceptions to avoid exiting the
                               # thread prematurely
            print('%r failed: %s' % (args, e,), file=sys.stderr)

# populate files
# ws = r'D:\Data\Users\jbellino\Project\stJohnsDeepening\model\xsec_a'
# workdir = os.path.join(ws, r'fieldgen\reals')
# files = ((os.path.join(workdir, f), ws)
#          for f in os.listdir(workdir) if f.endswith('.npy'))

### ensemble of forcing files ###

#prepare the case config and model_config files (assuming ensemble forcing files already exist and are named case_name_ensXX)

base_case_config = 'twpice_gf_test_t1534_ctl'
case_config_dir = 'case_config'
model_config_dir = 'model_config'
num_files = 100
#q = Queue()
ens_case_configs = []

#read in the base case config namelist
base_nml_filename = case_config_dir + '/' + base_case_config
base_nml = f90nml.read(base_nml_filename + '.nml')
if base_nml['case_config']:
    case_name = base_nml['case_config']['case_name']
    output_dir = base_nml['case_config']['output_dir']
    output_file = base_nml['case_config']['output_file']
    for i in range(num_files):

        patch_nml = {'case_config': {'case_name': case_name + '_ens' + str(i).zfill(2), 'output_dir': output_dir + '/ens', 'output_file': output_file + '_ens' + str(i).zfill(2)}}
        f90nml.patch(base_nml_filename + '.nml', patch_nml, base_nml_filename + '_ens' + str(i).zfill(2) + '.nml')
        ens_case_configs.append(base_case_config + '_ens' + str(i).zfill(2))
        #q.put_nowait((base_case_config + '_ens' + str(i).zfill(2),))

        #copy over the same model_config
        shutil.copy(model_config_dir + '/' + base_case_config + '.nml', model_config_dir + '/' + base_case_config + '_ens' + str(i).zfill(2) + '.nml')
else:
    print('Cannot make an ensemble out of a single timestep case config file')

os.chdir('bin')

# # start threads
# threads = [Thread(target=worker, args=(q,)) for _ in range(mp.cpu_count())]
# for t in threads:
#     t.daemon = True # threads die if the program dies
#     t.start()
#
# for _ in threads: q.put_nowait(None) # signal no more files
# for t in threads: t.join() # wait for completion

# def main():
#     # populate files
#     # ws = r'D:\Data\Users\jbellino\Project\stJohnsDeepening\model\xsec_a'
#     # workdir = os.path.join(ws, r'fieldgen\reals')
#     # files = ((os.path.join(workdir, f), ws)
#     #          for f in os.listdir(workdir) if f.endswith('.npy'))
#
#     # start processes
#     pool = mp.Pool() # use all available CPUs
#     pool.map(safe_run, ens_case_configs)
#
# if __name__=="__main__":
#     mp.freeze_support() # optional if the program is not frozen
#     main()

for i in range(len(ens_case_configs)):
    cmd = './gmtb_scm ' + ens_case_configs[i]
    args = shlex.split(cmd)
    print(args)
    with open(os.devnull, "w") as f:
        #process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process = subprocess.Popen(args, stdout=f, stderr=f)
    out, err = process.communicate()
    exit_code = process.wait()
    while exit_code != 0:
        with open(os.devnull, "w") as f:
            #process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process = subprocess.Popen(args, stdout=f, stderr=f)
        out, err = process.communicate()
        exit_code = process.wait()
    print(out, err, exit_code)

# os.chdir(output_dir + '/ens')
#
# cmd = 'ncecat -O ' + output_file + '_ens??.nc ' + output_file + '_ens.nc'
# args = shlex.split(cmd)
# print(os.getcwd())
# print(cmd, args)
# process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
