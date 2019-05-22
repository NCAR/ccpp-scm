#!/usr/bin/env python

import argparse
import atexit
import f90nml
import logging
import os
import shutil
import subprocess
import sys
import time

###############################################################################
# Global settings                                                             #
###############################################################################

# Name of the Fortran executable to run, including path (relative to run dir)
EXECUTABLE = './gmtb_scm'

# Path to the directory containing experiment namelists (relative to run dir)
EXPERIMENT_NAMELIST_DIR = '../etc/experiment_config'

# Standard name of experiment namelist in run directory
STANDARD_EXPERIMENT_NAMELIST = 'input_experiment.nml'

# Path to the directory containing physics namelists (relative to run dir)
PHYSICS_NAMELIST_DIR = '../../ccpp/physics_namelists'

# Path to the directory containing physics namelists (relative to run dir)
PHYSICS_SUITE_DIR = '../../ccpp/suites'

# Default output directory (relative to run dir), must match default value in gmtb_scm_input.f90
DEFAULT_OUTPUT_DIR = 'output'

# Default settings and filenames of input data for ozone physics;
# these must match the default settings in GFS_typedefs.F90.
DEFAULT_OZ_PHYS      = True
DEFAULT_OZ_PHYS_2015 = False
OZ_PHYS_TARGET       = 'global_o3prdlos_orig.f77'
OZ_PHYS_2015_TARGET  = 'ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77'
OZ_PHYS_LINK         = 'global_o3prdlos.f77'

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-e', '--experiment', help='name of experiment to run', required=True)
parser.add_argument('-g', '--gdb',        help='invoke gmtb_scm through gdb', action='store_true', default=False)

###############################################################################
# Functions and subroutines                                                   #
###############################################################################

def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)

def execute(cmd):
    """Runs a local command in a shell. Waits for completion and
    returns status, stdout and stderr."""
    logging.debug('Executing "{0}"'.format(cmd))
    p = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE, shell = True)
    (stdout, stderr) = p.communicate()
    status = p.returncode
    if status == 0:
        message = 'Execution of "{0}" returned with exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        logging.debug(message)
    else:
        message = 'Execution of command "{0}" failed, exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.rstrip('\n'))
        logging.debug(message)
    return (status, stdout.rstrip('\n'), stderr.rstrip('\n'))

def parse_arguments():
    """Parse command line arguments"""
    args = parser.parse_args()
    experiment = args.experiment
    gdb = args.gdb
    return (experiment, gdb)

def find_gdb():
    """Detect gdb, abort if not found"""
    logging.info('Searching for gdb ...')
    cmd = 'which gdb'
    (status, stdout, stderr) = execute(cmd)
    if status==1:
        message = 'gdb not found'
        logging.critical(message)
        raise Exception(message)
    gdb = stdout.strip()
    logging.info('Found {0}'.format(gdb))
    return gdb

class Experiment(object):

    def __init__(self, name):
        """Initialize experiment. This routine does most of the work,
        including setting and checking the experiment configuration
        (namelist)."""
        self._name = name
        self._namelist = os.path.join(EXPERIMENT_NAMELIST_DIR, self._name + '.nml')
        if not os.path.isfile(self._namelist):
            message = 'Experiment {0} with namelist {1} not found'.format(self._name, self._namelist)
            logging.critical(message)
            raise Exception(message)

    @property
    def name(self):
        """Get the name of the experiment."""
        return self._name

    @name.setter
    def name(self, value):
        """Set the name of the experiment."""
        self._name = value

    @property
    def namelist(self):
        """Get the namelist of the experiment."""
        return self._namelist

    @namelist.setter
    def namelist(self, value):
        """Set the namelist of the experiment."""
        self._namelist = value

    def setup_rundir(self):
        """Set up run directory for this experiment."""
        # Link namelist to run directory with name STANDARD_EXPERIMENT_NAMELIST
        logging.info('Linking experiment namelist {0} to {1} in run directory'.format(self._namelist, STANDARD_EXPERIMENT_NAMELIST))
        if os.path.isfile(STANDARD_EXPERIMENT_NAMELIST):
            os.remove(STANDARD_EXPERIMENT_NAMELIST)
        cmd = "ln -sf {0} {1}".format(self._namelist, STANDARD_EXPERIMENT_NAMELIST)
        execute(cmd)

        # Parse STANDARD_EXPERIMENT_NAMELIST and extract
        # - physics namelist
        # - physics suites (list)
        # - output directory
        logging.info('Parsing experiment namelist {0}'.format(STANDARD_EXPERIMENT_NAMELIST))
        nml = f90nml.read(STANDARD_EXPERIMENT_NAMELIST)
        physics_namelist = nml['physics_config']['physics_nml']
        physics_suites = [ 'suite_' + suite + '.xml' for suite in nml['physics_config']['physics_suite'].split(',')]
        try:
            output_dir = nml['experiment_config']['output_dir']
        except KeyError:
            output_dir = DEFAULT_OUTPUT_DIR
        del nml

        # Link physics namelist to run directory with its original name
        logging.info('Linking physics namelist {0} to run directory'.format(physics_namelist))
        if os.path.isfile(physics_namelist):
            os.remove(physics_namelist)
        if not os.path.isfile(os.path.join(PHYSICS_NAMELIST_DIR, physics_namelist)):
            message = 'Physics namelist {0} not found in directory {1}'.format(physics_namelist, PHYSICS_NAMELIST_DIR)
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(PHYSICS_NAMELIST_DIR, physics_namelist), physics_namelist)
        execute(cmd)

        # Link physics suites to run directory with its original name
        for physics_suite in physics_suites:
            logging.info('Linking physics suite {0} to run directory'.format(physics_suite))
            if os.path.isfile(physics_suite):
                os.remove(physics_suite)
            if not os.path.isfile(os.path.join(PHYSICS_SUITE_DIR, physics_suite)):
                message = 'Physics suite {0} not found in directory {1}'.format(physics_suite, PHYSICS_SUITE_DIR)
                logging.critical(message)
                raise Exception(message)
            cmd = "ln -sf {0} {1}".format(os.path.join(PHYSICS_SUITE_DIR, physics_suite), physics_suite)
            execute(cmd)

        # Parse physics namelist and extract
        # - oz_phys
        # - oz_phys_2015
        logging.info('Parsing physics namelist {0}'.format(physics_namelist))
        nml = f90nml.read(physics_namelist)
        # oz_phys
        try:
            oz_phys = nml['gfs_physics_nml']['oz_phys']
        except KeyError:
            oz_phys = DEFAULT_OZ_PHYS
        # oz_phys_2015
        try:
            oz_phys_2015 = nml['gfs_physics_nml']['oz_phys_2015']
        except KeyError:
            oz_phys_2015 = DEFAULT_OZ_PHYS_2015
        # Make sure that only one of the two ozone physics options is activated
        if oz_phys_2015 and oz_phys:
            message = 'Logic error, both oz_phys and oz_phys_2015 are set to true in the physics namelist'
            logging.critical(message)
            raise Exception(message)

        # Link input data for oz_phys or oz_phys_2015
        if os.path.exists(OZ_PHYS_LINK):
            os.remove(OZ_PHYS_LINK)
        if oz_phys:
            logging.info('Linking input data for oz_phys')
            cmd = 'ln -sf {0} {1}'.format(OZ_PHYS_TARGET, OZ_PHYS_LINK)
            execute(cmd)
        elif oz_phys_2015:
            logging.info('Linking input data for oz_phys_2015')
            cmd = 'ln -sf {0} {1}'.format(OZ_PHYS_2015_TARGET, OZ_PHYS_LINK)
            execute(cmd)

        # Create output directory (delete existing directory)
        logging.info('Creating output directory {0} in run directory'.format(output_dir))
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)

def launch_executable(use_gdb, gdb):
    """Configure model run command and pass control to shell/gdb"""
    if use_gdb:
        cmd = '{gdb} {executable}'.format(gdb=gdb, executable=EXECUTABLE)
    else:
        cmd = '{executable}'.format(executable=EXECUTABLE)
    logging.info('Passing control to "{0}"'.format(cmd))
    time.sleep(2)
    os.system(cmd)

def main():
    # Basics
    setup_logging()
    (experiment_name, use_gdb) = parse_arguments()
    # Experiment
    logging.info('Setting up experiment {0}'.format(experiment_name))
    exp = Experiment(experiment_name)
    exp.setup_rundir()
    # Debugger
    if use_gdb:
        gdb = find_gdb()
    else:
        gdb = None
    # Launch model on exit
    atexit.register(launch_executable, use_gdb, gdb)

if __name__ == '__main__':
    main()