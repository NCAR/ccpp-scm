#!/usr/bin/env python

import sys
import argparse
import importlib
import os
import logging
import subprocess
from supported_suites import suites
from supported_cases import cases
import timeit, functools

# Name of the python runscript executable to run, including path (relative to run dir)
RUN_SCRIPT = './run_scm.py'

# number of realizations to time if timer is used
timer_iterations = 1

# error strings to check in the output
ERROR_STRINGS = ['systemexit','sigsegv','backtrace']

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()#required=True)
group.add_argument('-c', '--case',       help='name of case to run',)
group.add_argument('-s', '--suite',      help='name of suite to use',)
group.add_argument('-f', '--file',       help='name of file where SCM runs are defined',)
parser.add_argument('-v', '--verbose',   help='once: set logging level to debug; twice: set logging level to debug '\
                                              'and write log to file', action='count', default=0)
parser.add_argument('-t', '--timer',     help='set to time each subprocess', action='store_true', default=False)
parser.add_argument('-d', '--docker',    help='include if scm is being run in a docker container to mount volumes', action='store_true', default=False)
parser.add_argument('--runtime',         help='set the runtime in the namelists', action='store', required=False)
parser.add_argument('-l', '--levels',    help='number of vertical levels', required=False)
parser.add_argument('--npz_type',        help='type of FV3 vertical grid to produce (see scm_vgrid.F90 for valid values)', required=False)
parser.add_argument('--vert_coord_file', help='filename with coefficients to produce a vertical grid', required=False)
parser.add_argument('--case_data_dir',   help='directory containing the case input data netCDF file', required=False)

# Results are recorded in this global list (to avoid complications with getting return values from the partial functions used below)
RESULTS = []

def setup_logging(verbose):
    """Sets up the logging module."""
    # print out debug messages (logs and output from subprocesses) if verbose argument is set
    if verbose:
        LOG_LEVEL = logging.DEBUG
    else:
        LOG_LEVEL = logging.INFO
    LOG_FILE = 'multi_run_scm.log'
    LOG_FORMAT = '%(levelname)s: %(message)s'

    logging.basicConfig(format=LOG_FORMAT, level=LOG_LEVEL)

    # write out a log file if verbosity is set twice (-vv)
    if verbose > 1:
        fh = logging.FileHandler(LOG_FILE, mode='w')
        logger = logging.getLogger()
        logger.addHandler(fh)

# function to spawn subprocesses directly (default) or through a timer
def spawn_subprocess(command, timer):
    if timer:
        t = timeit.Timer(functools.partial(subprocess_work, command))
        #timer_iterations determines how many realizations are executed in order to calculate the elapsed time
        elapsed_time = t.timeit(timer_iterations)
        return elapsed_time
    else:
        subprocess_work(command)

# function to actually do the work of executing a subprocess and adding the stderr and stdout to the log
# at the DEBUG level
def subprocess_work(command):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    (output, err) = p.communicate()
    p_status = p.wait()
    logging.debug(output)
    exit_code = p.returncode
    if not exit_code == 0:
        message = '####### The subprocess started using the command ({0}) exited with code {1}. #######\n'\
                  'Run the command ({0}) by itself again or use the -v or -vv options for more details.'.format(command, exit_code)
        logging.critical(message)
    else:
        found_error_strings = {x for x in ERROR_STRINGS if x in str(output).lower()}
        if found_error_strings:
             message = '####### The subprocess started using the command ({0}) exited with a normal exit code, but\n'\
             'the terminal output indicated the following strings occured in the output: ({1}). #######\n'\
             'Run the command ({0}) by itself again or use the -v or -vv options for more details.'.format(command, found_error_strings)
             logging.critical(message)
             exit_code = -999 #unknown
    RESULTS.append([command, exit_code])

def main():
    args = parser.parse_args()

    setup_logging(args.verbose)
    
    passthrough_args = ''
    if args.levels:
        passthrough_args += ' -l ' + args.levels
    if args.npz_type:
        if args.npz_type == 'input':
            if args.vert_coord_file:
                passthrough_args += ' --npz_type ' + args.npz_type + ' --vert_coord_file ' + args.vert_coord_file
            else:
                message = 'The npz_type was set to \'input\' but no file was specified via --vert_coord_file. Please specify the name of the file via --vert_coord_file name_of_file'
                logging.critical(message)
                raise Exception(message)
        else:
            passthrough_args += ' --npz_type ' + args.npz_type
    if args.case_data_dir:
        passthrough_args += ' --case_data_dir ' + args.case_data_dir
    if args.docker:
        passthrough_args += ' -d'
    if args.runtime:
        passthrough_args += ' --runtime ' + args.runtime
    
    # if the case argument is specified, run through all supported suites with the specified case
    if args.case:
        logging.info('Running all supported suites with case {0} and passthrough args: {1}'.format(args.case, passthrough_args))
        for i, suite in enumerate(suites,1):
            command = RUN_SCRIPT + ' -c ' + args.case + ' -s ' + suite + passthrough_args
            logging.info('Executing process {0} of {1} ({2})'.format(i, len(suites), command))
            elapsed_time = spawn_subprocess(command, args.timer)
            if args.timer:
               logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                  .format(elapsed_time/timer_iterations, args.case, suite))

    # if the suite argument is specified, run through all supported cases with the specified suite
    if args.suite:
        logging.info('Running all supported cases with suite {0} and passthrough args: {1}'.format(args.suite, passthrough_args))
        for i, case in enumerate(cases,1):
            command = RUN_SCRIPT + ' -c ' + case + ' -s ' + args.suite + passthrough_args
            logging.info('Executing process {0} of {1} ({2})'.format(i, len(cases), command))
            elapsed_time = spawn_subprocess(command, args.timer)
            if args.timer:
              logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                 .format(elapsed_time/timer_iterations, case, args.suite))

    # For maximum flexibility, run the SCM as specified from an external file where cases, suites, and physics namelists
    # are all specified. This file must contain python lists called 'cases','suites', and 'namelists'. The suites and
    # namelists lists can be empty ([]) if necessary.
    #The following rules apply:
    # 1. The case list in the file must not be empty.
    # 2. If only a case list is specified, the cases are run with the default suite specified in run_scm.py with
    #       the default namelists specified in default_namelists.py.
    # 3. If a case list and suite list is provided without a namelist list, all permutations of cases and suites will
    #       be run using default namelists specified in default_namelists.py.
    # 4. If a case list and suite list is provided with a namelist:
    # 4a. If only one suite is specified, it can be run with any number of namelists.
    # 4b. If more than one suite is specified, the number of namelists must match, and each case is run with each
    #       (suite,namelist) pair, by order specified in the lists.
    # 5. If a case list and namelist list are specified without a suite list, each case is run with the default suite
    #       specified in run_scm.py using the supplied namelists.
    if args.file:
        logging.info('Importing {0} to run requested combinations'.format(args.file))
        try:
            dirname, basename = os.path.split(args.file)
            sys.path.append(dirname)
            module_name = os.path.splitext(basename)[0]
            scm_runs = importlib.import_module(module_name)
            sys.path.pop()
        except ImportError:
            message = 'There was a problem loading {0}. Please check that the path exists.'.format(args.file)
            logging.critical(message)
            raise Exception(message)

        if not scm_runs.cases:
            message = 'The cases list in {0} must not be empty'.format(args.file)
            logging.critical(message)
            raise Exception(message)

        if scm_runs.cases and not scm_runs.suites and not scm_runs.namelists:
            logging.info(
                'Only cases were specified in {0}, so running all cases with the default suite and passthrough args: {1}'.format(args.file, passthrough_args))
            for i, case in enumerate(scm_runs.cases,1):
                command = RUN_SCRIPT + ' -c ' + case + passthrough_args
                logging.info('Executing process {0} of {1} ({2})'.format(i, len(scm_runs.cases), command))
                elapsed_time = spawn_subprocess(command, args.timer)
                if args.timer:
                   logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                      .format(elapsed_time/timer_iterations, case, args.file))

        if scm_runs.cases and scm_runs.suites:
            if scm_runs.namelists:
                if len(scm_runs.suites) == 1:
                    logging.info('Cases and namelists were specified with 1 suite in {0}, so running all cases with '\
                        'the suite {1} for all specified namelists and passthrough args: {2}'.format(args.file, scm_runs.suites[0], passthrough_args))
                    for i, case in enumerate(scm_runs.cases):
                        for j, namelist in enumerate(scm_runs.namelists,1):
                            command = RUN_SCRIPT + ' -c ' + case + ' -s ' + scm_runs.suites[0] + ' -n ' + namelist + passthrough_args
                            logging.info('Executing process {0} of {1} ({2})'.format(
                                len(scm_runs.namelists)*i+j, len(scm_runs.cases)*len(scm_runs.namelists), command))
                            elapsed_time = spawn_subprocess(command, args.timer)
                            if args.timer:
                               logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                                  .format(elapsed_time/timer_iterations, case, scm_runs.suites[0]))
                elif len(scm_runs.suites) == len(scm_runs.namelists):
                    logging.info('Cases, suites, and namelists were specified in {0}, so running all cases with all '\
                        'suites, matched with namelists by order and passthrough args: {1}'.format(args.file, passthrough_args))
                    for i, case in enumerate(scm_runs.cases):
                        for j, suite in enumerate(scm_runs.suites,1):
                            command = RUN_SCRIPT + ' -c ' + case + ' -s ' + suite + ' -n ' + scm_runs.namelists[j-1] + passthrough_args
                            logging.info('Executing process {0} of {1} ({2})'.format(
                                len(scm_runs.suites)*i+j, len(scm_runs.cases)*len(scm_runs.suites), command))
                            elapsed_time = spawn_subprocess(command, args.timer)
                            if args.timer:
                               logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                                  .format(elapsed_time/timer_iterations, case, suite))
                else:
                    message = 'The number of suites and namelists specified in {0} is incompatible. Either use one '\
                        'suite with many namelists or the number of suites must match the number of namelists '\
                        'provided.'.format(args.file)
                    logging.critical(message)
                    raise Exception(message)
            else:
                logging.info('Cases and suites specified in {0}, so running all cases with all suites using default '\
                    'namelists for each suite and passthrough args: {1}'.format(args.file, passthrough_args))
                for i, case in enumerate(scm_runs.cases):
                    for j, suite in enumerate(scm_runs.suites,1):
                        command = RUN_SCRIPT + ' -c ' + case + ' -s ' + suite + passthrough_args
                        logging.info('Executing process {0} of {1} ({2})'.format(
                            len(scm_runs.suites)*i+j, len(scm_runs.cases)*len(scm_runs.suites), command))
                        elapsed_time = spawn_subprocess(command, args.timer)
                        if args.timer:
                           logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                              .format(elapsed_time/timer_iterations, case, suite))

        if scm_runs.cases and not scm_runs.suites and scm_runs.namelists:
            logging.info('Cases and namelists were specified in {0}, so running all cases with the default suite '\
                'using the list of namelists and passthrough args: {1}'.format(args.file, passthrough_args))
            for i, case in enumerate(scm_runs.cases):
                for j, namelist in enumerate(scm_runs.namelists,1):
                    command = RUN_SCRIPT + ' -c ' + case + ' -n ' + namelist + passthrough_args
                    logging.info('Executing process {0} of {1} ({2})'.format(
                        len(scm_runs.namelists)*i+j, len(scm_runs.cases)*len(scm_runs.namelists), command))
                    elapsed_time = spawn_subprocess(command, args.timer)
                    if args.timer:
                      logging.info('elapsed time: {0} s for case: {1} namelist: {2}'\
                         .format(elapsed_time/timer_iterations, case, namelist))

    # If running the script with no arguments, run all supported (case,suite) permutations.
    if not args.case and not args.suite and not args.file:
        logging.info('Since no arguments were specified, running through all permuatations of supported cases and '\
            'suites with passthrough args: {0}').format(passthrough_args)
        for i, case in enumerate(cases):
            for j, suite in enumerate(suites,1):
                command = RUN_SCRIPT + ' -c ' + case + ' -s ' + suite + passthrough_args
                logging.info('Executing process {0} of {1} ({2})'.format(
                    len(suites)*i+j, len(cases)*len(suites), command))
                elapsed_time = spawn_subprocess(command, args.timer)
                if args.timer:
                   logging.info('elapsed time: {0} s for case: {1} suite: {2}'\
                      .format(elapsed_time/timer_iterations, case, suite))

    # Generate report at the end of the log file when verbose flag is set
    if args.verbose > 0:
        for (command, exit_code) in RESULTS:
            if exit_code == 0:
                logging.info('Process "{}" completed successfully'.format(command))
            else:
                logging.error('Process "{}" exited with code {}'.format(command, exit_code))

        logging.info('Done'.format(command))

if __name__ == '__main__':
    main()
