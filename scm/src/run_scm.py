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
from default_namelists import default_physics_namelists
from default_tracers import default_tracers
from netCDF4 import Dataset

###############################################################################
# Global settings                                                             #
###############################################################################

# Name of the Fortran executable to run, including path (relative to run dir)
EXECUTABLE = './scm'

# Path to the directory containing experiment namelists (relative to run dir)
CASE_NAMELIST_DIR = '../etc/case_config'

# Path to the directory containing tracer configurations (relative to run dir)
TRACERS_DIR = '../etc/tracer_config'
TRACERS_LINK = 'tracers.txt'

# Standard name of experiment namelist in run directory, must match value in scm_input.f90
STANDARD_EXPERIMENT_NAMELIST = 'input_experiment.nml'

# Path to the directory containing physics namelists (relative to run dir)
PHYSICS_NAMELIST_DIR = '../../ccpp/physics_namelists'

# Path to the directory containing physics namelists (relative to run dir)
PHYSICS_SUITE_DIR = '../../ccpp/suites'

# Default suite to use if none is specified
DEFAULT_SUITE = 'SCM_GFS_v15p2'

# Path to physics data files
PHYSICS_DATA_DIR = '../data/physics_input_data'

# Path to analysis script
SCM_ANALYSIS_SCRIPT_DIR = '../etc/scripts'

# Path to analysis script configuration files
SCM_ANALYSIS_CONFIG_DIR = '../etc/scripts/plot_configs'

# Default settings and filenames of input data for ozone physics;
# these must match the default settings in GFS_typedefs.F90.
DEFAULT_OZ_PHYS      = True
DEFAULT_OZ_PHYS_2015 = False
OZ_PHYS_TARGET       = 'global_o3prdlos_orig.f77'
OZ_PHYS_2015_TARGET  = 'ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77'
OZ_PHYS_LINK         = 'global_o3prdlos.f77'

#Default filename/targets for UGWPv1
DEFAULT_DO_UGWP_V1   = False
TAU_TARGET           = 'ugwp_c384_tau.nc'
TAU_LINK             = 'ugwp_limb_tau.nc'

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--case',       help='name of case to run', required=True)
parser.add_argument('-g', '--gdb',        help='invoke scm through gdb', action='store_true', default=False)
parser.add_argument('-s', '--suite',      help='name of suite to use', default=DEFAULT_SUITE)
parser.add_argument('-n', '--namelist',   help='physics namelist to use')
parser.add_argument('-t', '--tracers',    help='tracer configuration to use')
parser.add_argument('-d', '--docker',     help='include if scm is being run in a docker container to mount volumes', action='store_true', default=False)

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
        message += '    stdout: "{0}"\n'.format(stdout.decode(encoding='ascii', errors='ignore').rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.decode(encoding='ascii', errors='ignore').rstrip('\n'))
        logging.debug(message)
    else:
        message = 'Execution of command "{0}" failed, exit code {1}\n'.format(cmd, status)
        message += '    stdout: "{0}"\n'.format(stdout.decode(encoding='ascii', errors='ignore').rstrip('\n'))
        message += '    stderr: "{0}"'.format(stderr.decode(encoding='ascii', errors='ignore').rstrip('\n'))
        logging.debug(message)
    return (status, stdout.decode(encoding='ascii', errors='ignore').rstrip('\n'), stderr.decode(encoding='ascii', errors='ignore').rstrip('\n'))

def parse_arguments():
    """Parse command line arguments"""
    args = parser.parse_args()
    case = args.case
    gdb = args.gdb
    suite = args.suite
    namelist = args.namelist
    docker = args.docker
    tracers = args.tracers
    return (case, gdb, suite, namelist, docker, tracers)

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

    def __init__(self, case, suite, physics_namelist, tracers):
        """Initialize experiment. This routine does most of the work,
        including setting and checking the experiment configuration
        (namelist)."""
        
        self._case = case
        self._suite = suite
        self._name = case + '_' + suite
        
        #if a physics namelist is specified (entire filename), it will be used; 
        #otherwise, a default physics namelist for the given suite is used from default_namelists.py
        if physics_namelist:
            self._physics_namelist = physics_namelist
        else:
            if self._suite in default_physics_namelists:
                self._physics_namelist = default_physics_namelists.get(self._suite)
            else:
                message = 'A default physics namelist for suite {0} is not found in default_namelists.py'.format(self._suite)
                logging.critical(message)
                raise Exception(message)
        
        #check to see that the physics namelists exists in the right dir
        if not os.path.isfile(os.path.join(PHYSICS_NAMELIST_DIR, self._physics_namelist)):
            message = 'The physics namelist {0} was not found'.format(os.path.join(PHYSICS_NAMELIST_DIR, self._physics_namelist))
            logging.critical(message)
            raise Exception(message)
        
        #if a tracer configuration is specified (entire filename), it will be used; 
        #otherwise, a default tracer configuration for the given suite is used from default_tracers.py
        if tracers:
            self._tracers = tracers
        else:
            if self._suite in default_tracers:
                self._tracers = default_tracers.get(self._suite)
            else:
                message = 'A default tracer configuration for suite {0} is not found in default_tracers.py'.format(self._suite)
                logging.critical(message)
                raise Exception(message)
        
        #check to see that the physics namelists exists in the right dir
        if not os.path.isfile(os.path.join(TRACERS_DIR, self._tracers)):
            message = 'The tracer configuration {0} was not found'.format(os.path.join(TRACERS_DIR, self._tracers))
            logging.critical(message)
            raise Exception(message)
                        
        #check to see if the case namelists exists in the right dir
        self._namelist = os.path.join(CASE_NAMELIST_DIR, self._case + '.nml')
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
        """Get the case namelist of the experiment."""
        return self._namelist

    @namelist.setter
    def namelist(self, value):
        """Set the case namelist of the experiment."""
        self._namelist = value
    
    @property
    def case(self):
        """Get the case of the experiment."""
        return self._case
    
    @name.setter
    def case(self, value):
        """Set the case of the experiment."""
        self._case = value
    
    @property
    def suite(self):
        """Get the suite of the experiment."""
        return self._suite
    
    @suite.setter
    def suite(self, value):
        """Set the suite of the experiment."""
        self._suite = value
    
    @property
    def physics_namelist(self):
        """Get the physics namelist of the experiment."""
        return self._physics_namelist

    @physics_namelist.setter
    def physics_namelist(self, value):
        """Set the physics namelist of the experiment."""
        self._physics_namelist = value

    def setup_rundir(self):
        """Set up run directory for this experiment."""
        
        # Parse case configuration namelist and extract
        # - output directory
        # - surface_flux_spec
        logging.info('Parsing case configuration namelist {0}'.format(self._namelist))
        case_nml = f90nml.read(self._namelist)
        # look for the output_dir variable in the case configuration namelist and use it if it does; 
        # if it doesn't exist, create a default output directory name (from the case and suite names) and create a namelist patch
        try:
            output_dir = case_nml['case_config']['output_dir']
            custom_output_dir = True
        except KeyError:
            # If using the default namelist, no need to include it in the output directory name; if not, need to use write custom namelist in output dir name in case running multiple experiments with the same case and suite but different namelists
            if self._physics_namelist == default_physics_namelists.get(self._suite):
                output_dir = 'output_' + self._case + '_' + self._suite
            else:
                output_dir = 'output_' + self._case + '_' + self._suite + '_' + os.path.splitext(self._physics_namelist)[0]
            output_dir_patch_nml = {'case_config':{'output_dir':output_dir}}
            custom_output_dir = False
        
        #if using the DEPHY format, need to also check the case data file for the surfaceForcing global attribute for 'Flux' or 'surfaceFlux', which denotes prescribed surface fluxes
        try:
            input_type = case_nml['case_config']['input_type']
            if input_type == 1:
                #open the case data file and read the surfaceForcing global attribute
                case_data_dir = case_nml['case_config']['case_data_dir']
                nc_fid = Dataset(case_data_dir + '/' + self._case + '_SCM_driver.nc' , 'r')
                surfaceForcing = nc_fid.getncattr('surfaceForcing')
                nc_fid.close()
                if (surfaceForcing.lower() == 'flux' or surfaceForcing.lower() == 'surfaceflux'):
                    surface_flux_spec = True
        except KeyError:
            # if not using DEPHY format, check to see if surface fluxes are specified in the case configuration file (default is False)
            try:
                surface_flux_spec = case_nml['case_config']['sfc_flux_spec']
            except KeyError:
                surface_flux_spec = False
        
        # If surface fluxes are specified for this case, use the SDF modified to use them
        if surface_flux_spec:
            logging.info('Specified surface fluxes are used for case {0}. Switching to SDF {1} from {2}'.format(self._case,'suite_' + self._suite + '_ps' + '.xml','suite_' + self._suite + '.xml'))
            self._suite = self._suite + '_ps'
                
        # Create physics_config namelist for experiment configuration file
        physics_config = {"physics_suite":self._suite,
                          "physics_nml":self._physics_namelist,}
        physics_config_dict = {"physics_config":physics_config}
        physics_config_nml = f90nml.namelist.Namelist(physics_config_dict)
        
        # Create STANDARD_EXPERIMENT_NAMELIST in the run directory with the case configuration and physics configuration namelists
        logging.info('Creating experiment configuration namelist {0} in the run directory from {1} using {2} and {3} '.format(STANDARD_EXPERIMENT_NAMELIST,self._namelist,self._suite,self._physics_namelist))
        
        with open(STANDARD_EXPERIMENT_NAMELIST, "w+") as nml_file:
            case_nml.write(nml_file)
        
        with open(STANDARD_EXPERIMENT_NAMELIST, "a") as nml_file:
            physics_config_nml.write(nml_file)
        
        # if using the default output dir name created in this script, patch the experiment namelist with the new output_dir variable
        if(not custom_output_dir):
            # GJF TODO: this implementation is clunky; newer versions of f90nml can handle this better, but this works with v0.19 so no need to require newer version
            f90nml.patch(STANDARD_EXPERIMENT_NAMELIST, output_dir_patch_nml, 'temp.nml')
            cmd = "mv {0} {1}".format('temp.nml', STANDARD_EXPERIMENT_NAMELIST)
            execute(cmd)
        
        # Link physics namelist to run directory with its original name
        logging.info('Linking physics namelist {0} to run directory'.format(self._physics_namelist))
        if os.path.isfile(self._physics_namelist):
            os.remove(self._physics_namelist)
        if not os.path.isfile(os.path.join(PHYSICS_NAMELIST_DIR, self._physics_namelist)):
            message = 'Physics namelist {0} not found in directory {1}'.format(self._physics_namelist, PHYSICS_NAMELIST_DIR)
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(PHYSICS_NAMELIST_DIR, self._physics_namelist), self._physics_namelist)
        if surface_flux_spec:
            #check for optional prescribed-surface-specific physics namelist and link it instead (if present)
            opt_ps_nml_filename = os.path.splitext(os.path.join(PHYSICS_NAMELIST_DIR,self._physics_namelist))[0] + '_ps.nml'
            if os.path.isfile(opt_ps_nml_filename):
                logging.info('Found optional prescribed surface physics namelist {0}; linking it to run directory'.format(opt_ps_nml_filename))
                cmd = "ln -sf {0} {1}".format(opt_ps_nml_filename, self._physics_namelist)
        execute(cmd)
        
        # Link tracer configuration to run directory with its original name
        logging.info('Linking tracer configuration {0} to run directory'.format(self._tracers))
        if os.path.isfile(self._tracers):
            os.remove(self._tracers)
        if not os.path.isfile(os.path.join(TRACERS_DIR, self._tracers)):
            message = 'Tracer configuration {0} not found in directory {1}'.format(self._tracers, TRACERS_DIR)
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(TRACERS_DIR, self._tracers), TRACERS_LINK)
        execute(cmd)
        
        # Link physics SDF to run directory
        physics_suite = 'suite_' + self._suite + '.xml'
        logging.info('Linking physics suite {0} to run directory'.format(physics_suite))
        if os.path.isfile(physics_suite):
            os.remove(physics_suite)
        if not os.path.isfile(os.path.join(PHYSICS_SUITE_DIR, physics_suite)):
            message = 'Physics suite {0} not found in directory {1}'.format(physics_suite, PHYSICS_SUITE_DIR)
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(PHYSICS_SUITE_DIR, physics_suite), physics_suite)
        execute(cmd)
        
        # Link physics data needed for schemes to run directory
        logging.info('Linking physics input data from {0} into run directory'.format(PHYSICS_DATA_DIR))
        for entry in os.listdir(PHYSICS_DATA_DIR):
            if os.path.isfile(os.path.join(PHYSICS_DATA_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(PHYSICS_DATA_DIR, entry), entry)
                    execute(cmd)
        
        # Parse physics namelist and extract
        # - oz_phys
        # - oz_phys_2015
        logging.info('Parsing physics namelist {0}'.format(self._physics_namelist))
        nml = f90nml.read(self._physics_namelist)
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
        
        # Look for do_ugwp_v1
        try:
            do_ugwp_v1 = nml['gfs_physics_nml']['do_ugwp_v1']
        except KeyError:
            do_ugwp_v1 = DEFAULT_DO_UGWP_V1
        
        # Link the tau file if do_ugwp_v1
        if do_ugwp_v1:
            if os.path.exists(TAU_LINK):
                os.remove(TAU_LINK)
            logging.info('Linking input data for UGWP_v1')
            cmd = 'ln -sf {0} {1}'.format(TAU_TARGET, TAU_LINK)
            execute(cmd)
        
        # Link scripts needed to run SCM analysis
        logging.info('Linking analysis scripts from {0} into run directory'.format(SCM_ANALYSIS_SCRIPT_DIR))
        analysis_script_files = ['scm_analysis.py','configspec.ini']
        for entry in analysis_script_files:
            if os.path.isfile(os.path.join(SCM_ANALYSIS_SCRIPT_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ANALYSIS_SCRIPT_DIR, entry), entry)
                    execute(cmd)
        
        # Link plot configuration files needed to run SCM analysis
        logging.info('Linking plot configuration files from {0} into run directory'.format(SCM_ANALYSIS_CONFIG_DIR))
        for entry in os.listdir(SCM_ANALYSIS_CONFIG_DIR):
            if os.path.isfile(os.path.join(SCM_ANALYSIS_CONFIG_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ANALYSIS_CONFIG_DIR, entry), entry)
                    execute(cmd)
        
        # Create output directory (delete existing directory)
        logging.info('Creating output directory {0} in run directory'.format(output_dir))
        if os.path.isdir(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        
        # Write experiment configuration file to output directory
        logging.info('Writing experiment configuration {0}.nml to output directory'.format(self._name))
        cmd = 'cp {0} {1}'.format(STANDARD_EXPERIMENT_NAMELIST, os.path.join(output_dir,self._name + '.nml'))
        execute(cmd)
        
        return output_dir

def launch_executable(use_gdb, gdb):
    """Configure model run command and pass control to shell/gdb"""
    if use_gdb:
        cmd = '{gdb} {executable}'.format(gdb=gdb, executable=EXECUTABLE)
    else:
        cmd = '{executable}'.format(executable=EXECUTABLE)
    logging.info('Passing control to "{0}"'.format(cmd))
    time.sleep(2)
    return_code = os.system(cmd)
    if return_code != 0:
        sys.exit(return_code)
    
def copy_outdir(exp_dir):
    """Copy output directory to /home for this experiment."""
    home_output_dir = '/home/'+exp_dir
    if os.path.isdir(home_output_dir):
        shutil.rmtree(home_output_dir)
    shutil.copytree(exp_dir, home_output_dir)

def main():
    (case, use_gdb, suite, namelist, docker, tracers) = parse_arguments()
    
    setup_logging()
    
    #Experiment
    if namelist:
        if tracers:
            logging.info('Setting up experiment {0} with suite {1} using namelist {2} and tracers {3}'.format(case,suite,namelist,tracers))
        else:
            logging.info('Setting up experiment {0} with suite {1} using namelist {2} using default tracers for the suite'.format(case,suite,namelist))
    else:
        if tracers:
            logging.info('Setting up experiment {0} with suite {1} using the default namelist for the suite and tracers {2}'.format(case,suite,tracers))
        else:
            logging.info('Setting up experiment {0} with suite {1} using the default namelist and tracers for the suite'.format(case,suite))
    exp = Experiment(case, suite, namelist, tracers)
    exp_dir = exp.setup_rundir()
    # Debugger
    if use_gdb:
        gdb = find_gdb()
    else:
        gdb = None
    # Launch model on exit
    if docker:
        #registering this function first should mean that it executes last, which is what we want
        atexit.register(copy_outdir, exp_dir)
    atexit.register(launch_executable, use_gdb, gdb)
    
    
if __name__ == '__main__':
    main()
