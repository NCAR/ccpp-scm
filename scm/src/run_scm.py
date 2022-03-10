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
from suite_info import suite, suite_list
from netCDF4 import Dataset

###############################################################################
# Global settings                                                             #
###############################################################################

SCM_ROOT = ''
SCM_RUN = ''
SCM_BIN = ''
EXECUTABLE = ''

# Name of the Fortran executable to run, including path
EXECUTABLE_NAME = 'scm'

# Path to default run directory (relative to scm_root)
DEFAULT_RUN_DIR = 'scm/run'

# Path to default bin directory (relative to scm_root)
DEFAULT_BIN_DIR = 'scm/bin'

# Copy executable to run directory if true (otherwise it will be linked)
COPY_EXECUTABLE = False

# Default output periods
DEFAULT_OUTPUT_PERIOD = 1
DEFAULT_DIAG_PERIOD = 6

# Path to the directory containing experiment namelists (relative to scm_root)
CASE_NAMELIST_DIR = 'scm/etc/case_config'

# Path to the directory containing case input data files (relative to scm_root)
DEFAULT_CASE_DATA_DIR = 'scm/data/processed_case_input'

# Path to the directory containing vertical coordinate files (relative to scm_root)
VERT_COORD_DATA_DIR = 'scm/data/vert_coord_data'

# Path to reference profiles
REFERENCE_PROFILE_DIR = DEFAULT_CASE_DATA_DIR

# Reference profile file list
REFERENCE_PROFILE_FILE_LIST = ['McCProfiles.dat','mid_lat_summer_std.nc']

# Path to the directory containing tracer configurations (relative to scm_root)
TRACERS_DIR = 'scm/etc/tracer_config'
TRACERS_LINK = 'tracers.txt'

# Standard name of experiment namelist in run directory, must match value in scm_input.f90
STANDARD_EXPERIMENT_NAMELIST = 'input_experiment.nml'

# Path to the directory containing physics namelists (relative to scm_root)
PHYSICS_NAMELIST_DIR = 'ccpp/physics_namelists'

# Path to the directory containing physics namelists (relative to scm_root)
PHYSICS_SUITE_DIR = 'ccpp/suites'

# Default suite to use if none is specified
DEFAULT_SUITE = 'SCM_GFS_v15p2'

# Path to physics data files (relative to scm_root)
PHYSICS_DATA_DIR = 'scm/data/physics_input_data'

# Path to analysis script (relative to scm_root)
SCM_ANALYSIS_SCRIPT_DIR = 'scm/etc/scripts'

# Path to analysis script configuration files (relative to scm_root)
SCM_ANALYSIS_CONFIG_DIR = 'scm/etc/scripts/plot_configs'

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
parser.add_argument('--runtime',          help='set the runtime in the namelists', action='store', type=int, required=False)
parser.add_argument('--runtime_mult',     help='multiply the existing runtime in the namelist by some factor', action='store', type=float, required=False)
parser.add_argument('-d', '--docker',     help='include if scm is being run in a docker container to mount volumes', action='store_true', default=False)
parser.add_argument('-l', '--levels',     help='number of vertical levels', required=False, type=int)
parser.add_argument('--npz_type',         help='type of FV3 vertical grid to produce (see scm_vgrid.F90 for valid values)', required=False)
parser.add_argument('--vert_coord_file',  help='filename with coefficients to produce a vertical grid', required=False)
parser.add_argument('--run_dir',          help='path for the run directory', required=False)
parser.add_argument('--case_data_dir',    help='directory containing the case input data netCDF file', required=False)
parser.add_argument('--n_itt_out',        help='period of instantaneous output (number of timesteps)', required=False, type=int)
parser.add_argument('--n_itt_diag',       help='period of diagnostic output (number of timesteps)', required=False, type=int)
parser.add_argument('-dt', '--timestep',  help='timestep (s)', required=False, type=float)

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
    runtime = args.runtime
    runtime_mult = args.runtime_mult
    levels = args.levels
    npz_type = args.npz_type
    vert_coord_file = args.vert_coord_file
    case_data_dir = args.case_data_dir
    n_itt_out = args.n_itt_out
    n_itt_diag = args.n_itt_diag
    run_dir = args.run_dir
    timestep = args.timestep
    
    return (case, gdb, suite, namelist, docker, tracers, runtime, runtime_mult, levels, npz_type, vert_coord_file, case_data_dir, n_itt_out, n_itt_diag, run_dir, timestep)

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

    def __init__(self, case, suite, runtime, runtime_mult, levels, npz_type, vert_coord_file, case_data_dir, n_itt_out, n_itt_diag):
        """Initialize experiment. This routine does most of the work,
        including setting and checking the experiment configuration
        (namelist)."""
        
        self._case = case
        self._suite_obj = suite
        self._suite = suite._name
        self._name = case + '_' + suite._name
        
        self._physics_namelist = suite.namelist
                
        #check to see that the physics namelists exists in the right dir
        if not os.path.isfile(os.path.join(SCM_ROOT, PHYSICS_NAMELIST_DIR, self._physics_namelist)):
            message = 'The physics namelist {0} was not found'.format(os.path.join(SCM_ROOT, PHYSICS_NAMELIST_DIR, self._physics_namelist))
            logging.critical(message)
            raise Exception(message)
        
        self._tracers = suite.tracers
        
        #check to see that the tracers exists in the right dir
        if not os.path.isfile(os.path.join(SCM_ROOT, TRACERS_DIR, self._tracers)):
            message = 'The tracer configuration {0} was not found'.format(os.path.join(TRACERS_DIR, self._tracers))
            logging.critical(message)
            raise Exception(message)
                        
        #check to see if the case namelists exists in the right dir
        self._namelist = os.path.join(SCM_ROOT, CASE_NAMELIST_DIR, self._case + '.nml')
        if not os.path.isfile(self._namelist):
            message = 'Experiment {0} with namelist {1} not found'.format(self._name, self._namelist)
            logging.critical(message)
            raise Exception(message)

        if runtime:
            self._runtime = runtime
            message = 'Namelist runtime adjustment {0} IS applied'.format(self._runtime)
            logging.info(message)
        else:
            self._runtime = None
            message = 'Namelist runtime adjustment {0}  IS NOT applied'.format(self._runtime)
            logging.info(message)
        
        if runtime_mult:
            self._runtime_mult = runtime_mult
            message = 'Existing case namelist runtime multipled by {0}'.format(self._runtime_mult)
            logging.info(message)
        
        if levels:
            self._levels = levels
            message = 'The number of vertical levels is set to {0}'.format(self._levels)
            logging.info(message)
        else:
            self._levels = None
            message = 'The number of vertical levels contained in the case configuration file is used if present, otherwise the default value in scm_input.F90 is used.'
            logging.info(message)
        
        if npz_type:
            self._npz_type = npz_type
            message = 'The npz_type of vertical levels is set to {0}'.format(self._npz_type)
            logging.info(message)
            if npz_type == 'input':
                if vert_coord_file:
                    self._vert_coord_file = vert_coord_file
                    #check to see that the vert coords file exists in the right dir
                    if not os.path.isfile(os.path.join(SCM_ROOT, VERT_COORD_DATA_DIR, self._vert_coord_file)):
                        message = 'The vertical coordinate file {0} was not found'.format(os.path.join(SCM_ROOT, VERT_COORD_DATA_DIR, self._vert_coord_file))
                        logging.critical(message)
                        raise Exception(message)
                else:
                    message = 'The npz_type was set to \'input\' but no file was specified via --vert_coord_file. Please specify the name of the file via --vert_coord_file name_of_file'
                    logging.critical(message)
                    raise Exception(message)
            else:
                self._vert_coord_file = None
        else:
            self._npz_type = None
            self._vert_coord_file = None
            message = 'The npz_type contained in the case configuration file is used if present, otherwise the default value in scm_input.F90 is used.'
            logging.info(message)
        
        if case_data_dir:
            self._case_data_dir = case_data_dir
        else:
            self._case_data_dir = DEFAULT_CASE_DATA_DIR
        #check to see that the provided case data directory exists
        if not os.path.isdir(os.path.join(SCM_ROOT, self._case_data_dir)):
            message = 'The case data directory {0} was not found'.format(self._case_data_dir)
            logging.critical(message)
            raise Exception(message)
        
        if n_itt_out:
            self._n_itt_out = n_itt_out
        else:
            self._n_itt_out = DEFAULT_OUTPUT_PERIOD
        
        if n_itt_diag:
            self._n_itt_diag = n_itt_diag
        else:
            self._n_itt_diag = DEFAULT_DIAG_PERIOD
        
        if suite.timestep is not None:
            self._timestep = suite.timestep
        
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
    
    @property
    def tracers(self):
        """Get the tracer file for the experiment."""
        return self._tracers

    @tracers.setter
    def tracers(self, value):
        """Set the tracer file for the experiment."""
        self._tracers = value
        
    @property
    def levels(self):
        """Get the number of vertical levels for the experiment."""
        return self._levels

    @levels.setter
    def levels(self, value):
        """Set the number of vertical levels for the experiment."""
        self._levels = value
    
    @property
    def npz_type(self):
        """Get the vertical level type for the experiment."""
        return self._npz_type

    @npz_type.setter
    def npz_type(self, value):
        """Set the vertical level type for the experiment."""
        self._npz_type = value
    
    @property
    def vert_coord_file(self):
        """Get the file containing vertical levels for the experiment."""
        return self._vert_coord_file

    @vert_coord_file.setter
    def vert_coord_file(self, value):
        """Set the file containing vertical levels for the experiment."""
        self._vert_coord_file = value
    
    @property
    def case_data_dir(self):
        """Get the case data directory for the experiment."""
        return self._case_data_dir

    @case_data_dir.setter
    def case_data_dir(self, value):
        """Set the case data directory for the experiment."""
        self._case_data_dir = value
    
    @property
    def n_itt_out(self):
        """Get the output period (in timesteps) for the experiment."""
        return self._n_itt_out

    @n_itt_out.setter
    def n_itt_out(self, value):
        """Set the output period (in timesteps) for the experiment."""
        self._n_itt_out = value
    
    @property
    def n_itt_diag(self):
        """Get the diagnostic period (in timesteps) for the experiment."""
        return self._n_itt_diag

    @n_itt_diag.setter
    def n_itt_diag(self, value):
        """Set the diagnostic period (in timesteps) for the experiment."""
        self._n_itt_diag = value
    
    def setup_rundir(self):
        """Set up run directory for this experiment."""
        
        # Parse case configuration namelist and extract
        # - output directory
        # - surface_flux_spec
        logging.info('Parsing case configuration namelist {0}'.format(self._namelist))
        case_nml = f90nml.read(self._namelist)
        # If running the regression test, reduce the runtime
        if self._runtime:
           case_nml['case_config']['runtime'] = self._runtime
        if self._runtime_mult:
            try:
                old_runtime = case_nml['case_config']['runtime']
                case_nml['case_config']['runtime'] = old_runtime*self._runtime_mult
            except KeyError:
                logging.info('The runtime multiplier argument was set, but the runtime is not set in {0} '.format(self._namelist))
        # If the number of levels is specified, set the namelist value
        if self._levels:
           case_nml['case_config']['n_levels'] = self._levels
        # If the npz_type is specified, set the namelist value
        if self._npz_type:
           case_nml['case_config']['npz_type'] = self._npz_type
        if self._vert_coord_file:
           case_nml['case_config']['vert_coord_file'] = self._vert_coord_file
        if self._n_itt_out:
           case_nml['case_config']['n_itt_out'] = self._n_itt_out
        if self._n_itt_diag:
           case_nml['case_config']['n_itt_diag'] = self._n_itt_diag
        if self._timestep:
           case_nml['case_config']['dt'] = self._timestep
        # look for the output_dir variable in the case configuration namelist and use it if it does; 
        # if it doesn't exist, create a default output directory name (from the case and suite names) and create a namelist patch
        try:
            output_dir = case_nml['case_config']['output_dir']
            custom_output_dir = True
        except KeyError:
            # If using the default namelist, no need to include it in the output directory name; if not, need to use write custom namelist in output dir name in case running multiple experiments with the same case and suite but different namelists
            if self._physics_namelist == self._suite_obj._default_namelist:
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
                nc_fid = Dataset(os.path.join(SCM_ROOT, case_data_dir) + '/' + self._case + '_SCM_driver.nc' , 'r')
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
        
        with open(os.path.join(SCM_RUN, STANDARD_EXPERIMENT_NAMELIST), "w+") as nml_file:
            case_nml.write(nml_file)
        
        with open(os.path.join(SCM_RUN, STANDARD_EXPERIMENT_NAMELIST), "a") as nml_file:
            physics_config_nml.write(nml_file)
        
        # if using the default output dir name created in this script, patch the experiment namelist with the new output_dir variable
        if(not custom_output_dir):
            # GJF TODO: this implementation is clunky; newer versions of f90nml can handle this better, but this works with v0.19 so no need to require newer version
            f90nml.patch(os.path.join(SCM_RUN, STANDARD_EXPERIMENT_NAMELIST), output_dir_patch_nml, 'temp.nml')
            cmd = "mv {0} {1}".format('temp.nml', os.path.join(SCM_RUN, STANDARD_EXPERIMENT_NAMELIST))
            execute(cmd)
        
        # Link physics namelist to run directory with its original name
        logging.info('Linking physics namelist {0} to run directory'.format(self._physics_namelist))
        if os.path.isfile(os.path.join(SCM_RUN, self._physics_namelist)):
            os.remove(os.path.join(SCM_RUN,self._physics_namelist))
        if not os.path.isfile(os.path.join(SCM_ROOT, PHYSICS_NAMELIST_DIR, self._physics_namelist)):
            message = 'Physics namelist {0} not found in directory {1}'.format(SCM_ROOT, self._physics_namelist, PHYSICS_NAMELIST_DIR)
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(SCM_ROOT, PHYSICS_NAMELIST_DIR, self._physics_namelist), os.path.join(SCM_RUN, self._physics_namelist))
        if surface_flux_spec:
            #check for optional prescribed-surface-specific physics namelist and link it instead (if present)
            opt_ps_nml_filename = os.path.splitext(os.path.join(SCM_ROOT, PHYSICS_NAMELIST_DIR,self._physics_namelist))[0] + '_ps.nml'
            if os.path.isfile(opt_ps_nml_filename):
                logging.info('Found optional prescribed surface physics namelist {0}; linking it to run directory'.format(opt_ps_nml_filename))
                cmd = "ln -sf {0} {1}".format(opt_ps_nml_filename, os.path.join(SCM_RUN, self._physics_namelist))
        execute(cmd)
        
        # Link tracer configuration to run directory with standard name
        logging.info('Linking tracer configuration {0} to run directory'.format(self._tracers))
        if os.path.isfile(os.path.join(SCM_RUN, self._tracers)):
            os.remove(os.path.join(SCM_RUN, self._tracers))
        if not os.path.isfile(os.path.join(SCM_ROOT, TRACERS_DIR, self._tracers)):
            message = 'Tracer configuration {0} not found in directory {1}'.format(self._tracers, os.path.join(SCM_ROOT, TRACERS_DIR))
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(SCM_ROOT, TRACERS_DIR, self._tracers), os.path.join(SCM_RUN, TRACERS_LINK))
        execute(cmd)
        
        # Link case data file to run directory with original name
        try:
            input_type = case_nml['case_config']['input_type']
            if input_type == 1:
                case_data_netcdf_file = self._case + '_SCM_driver.nc'
            else:
                case_data_netcdf_file = self._case + '.nc'
        except KeyError:
            case_data_netcdf_file = self._case + '.nc'
        logging.info('Linking case input data file {0} to run directory'.format(case_data_netcdf_file))
        if os.path.isfile(os.path.join(SCM_RUN, case_data_netcdf_file)):
            os.remove(os.path.join(SCM_RUN, case_data_netcdf_file))
        if not os.path.isfile(os.path.join(SCM_ROOT, self._case_data_dir, case_data_netcdf_file)):
            message = 'Case input data file {0} not found in directory {1}'.format(case_data_netcdf_file, os.path.join(SCM_ROOT, self._case_data_dir))
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(SCM_ROOT, self._case_data_dir, case_data_netcdf_file), os.path.join(SCM_RUN, case_data_netcdf_file))
        execute(cmd)
        
        # Link vertical coordinate file to run directory with its original name
        if (self._npz_type == 'input'):
            logging.info('Linking vertical coordinate file {0} to run directory'.format(self._vert_coord_file))
            if os.path.isfile(os.path.join(SCM_RUN, self._vert_coord_file)):
                os.remove(os.path.join(SCM_RUN, self._vert_coord_file))
            if not os.path.isfile(os.path.join(SCM_ROOT, VERT_COORD_DATA_DIR, self._vert_coord_file)):
                message = 'Vertical coordinate file {0} not found in directory {1}'.format(self._vert_coord_file, os.path.join(SCM_ROOT, VERT_COORD_DATA_DIR))
                logging.critical(message)
                raise Exception(message)
            cmd = "ln -sf {0} {1}".format(os.path.join(SCM_ROOT, VERT_COORD_DATA_DIR, self._vert_coord_file), os.path.join(SCM_RUN, self._vert_coord_file))
            execute(cmd)
        
        # Link physics SDF to run directory
        physics_suite = 'suite_' + self._suite + '.xml'
        logging.info('Linking physics suite {0} to run directory'.format(physics_suite))
        if os.path.isfile(os.path.join(SCM_RUN, physics_suite)):
            os.remove(os.path.join(SCM_RUN, physics_suite))
        if not os.path.isfile(os.path.join(SCM_ROOT, PHYSICS_SUITE_DIR, physics_suite)):
            message = 'Physics suite {0} not found in directory {1}'.format(physics_suite, os.path.join(SCM_ROOT, PHYSICS_SUITE_DIR))
            logging.critical(message)
            raise Exception(message)
        cmd = "ln -sf {0} {1}".format(os.path.join(SCM_ROOT, PHYSICS_SUITE_DIR, physics_suite), os.path.join(SCM_RUN, physics_suite))
        execute(cmd)
        
        # Link physics data needed for schemes to run directory
        logging.info('Linking physics input data from {0} into run directory'.format(os.path.join(SCM_ROOT, PHYSICS_DATA_DIR)))
        for entry in os.listdir(os.path.join(SCM_ROOT, PHYSICS_DATA_DIR)):
            if os.path.isfile(os.path.join(SCM_ROOT, PHYSICS_DATA_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ROOT, PHYSICS_DATA_DIR, entry), os.path.join(SCM_RUN, entry))
                    execute(cmd)
        
        # Link reference profile data to run directory
        logging.info('Linking reference profile data from {0} into run directory'.format(os.path.join(SCM_ROOT, REFERENCE_PROFILE_DIR)))
        for entry in REFERENCE_PROFILE_FILE_LIST:
            if os.path.isfile(os.path.join(SCM_ROOT, REFERENCE_PROFILE_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ROOT, REFERENCE_PROFILE_DIR, entry), os.path.join(SCM_RUN, entry))
                    execute(cmd)
        
        # Parse physics namelist and extract
        # - oz_phys
        # - oz_phys_2015
        logging.info('Parsing physics namelist {0}'.format(os.path.join(SCM_RUN, self._physics_namelist)))
        nml = f90nml.read(os.path.join(SCM_RUN, self._physics_namelist))
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
        if os.path.exists(os.path.join(SCM_RUN, OZ_PHYS_LINK)):
            os.remove(os.path.join(SCM_RUN, OZ_PHYS_LINK))
        if oz_phys:
            logging.info('Linking input data for oz_phys')
            cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_RUN, OZ_PHYS_TARGET), os.path.join(SCM_RUN, OZ_PHYS_LINK))
            execute(cmd)
        elif oz_phys_2015:
            logging.info('Linking input data for oz_phys_2015')
            cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_RUN, OZ_PHYS_2015_TARGET), os.path.join(SCM_RUN, OZ_PHYS_LINK))
            execute(cmd)
        
        # Look for do_ugwp_v1
        try:
            do_ugwp_v1 = nml['gfs_physics_nml']['do_ugwp_v1']
        except KeyError:
            do_ugwp_v1 = DEFAULT_DO_UGWP_V1
        
        # Link the tau file if do_ugwp_v1
        if do_ugwp_v1:
            if os.path.exists(os.path.join(SCM_RUN, TAU_LINK)):
                os.remove(os.path.join(SCM_RUN, TAU_LINK))
            logging.info('Linking input data for UGWP_v1')
            cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_RUN, TAU_TARGET), os.path.join(SCM_RUN, TAU_LINK))
            execute(cmd)
        
        # Link scripts needed to run SCM analysis
        logging.info('Linking analysis scripts from {0} into run directory'.format(os.path.join(SCM_ROOT, SCM_ANALYSIS_SCRIPT_DIR)))
        analysis_script_files = ['scm_analysis.py','configspec.ini']
        for entry in analysis_script_files:
            if os.path.isfile(os.path.join(SCM_ROOT, SCM_ANALYSIS_SCRIPT_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ROOT, SCM_ANALYSIS_SCRIPT_DIR, entry), os.path.join(SCM_RUN, entry))
                    execute(cmd)
        
        # Link plot configuration files needed to run SCM analysis
        logging.info('Linking plot configuration files from {0} into run directory'.format(os.path.join(SCM_ROOT, SCM_ANALYSIS_CONFIG_DIR)))
        for entry in os.listdir(os.path.join(SCM_ROOT, SCM_ANALYSIS_CONFIG_DIR)):
            if os.path.isfile(os.path.join(SCM_ROOT, SCM_ANALYSIS_CONFIG_DIR, entry)):
                if not os.path.exists(entry):
                    logging.debug('Linking file {0}'.format(entry))
                    cmd = 'ln -sf {0} {1}'.format(os.path.join(SCM_ROOT, SCM_ANALYSIS_CONFIG_DIR, entry), os.path.join(SCM_RUN, entry))
                    execute(cmd)
        
        # Create output directory (delete existing directory)
        logging.info('Creating output directory {0} in run directory'.format(output_dir))
        if os.path.isdir(os.path.join(SCM_RUN, output_dir)):
            shutil.rmtree(os.path.join(SCM_RUN, output_dir))
        os.makedirs(os.path.join(SCM_RUN, output_dir))
        
        # Write experiment configuration file to output directory
        logging.info('Writing experiment configuration {0}.nml to output directory'.format(self._name))
        cmd = 'cp {0} {1}'.format(os.path.join(SCM_RUN, STANDARD_EXPERIMENT_NAMELIST), os.path.join(SCM_RUN, output_dir,self._name + '.nml'))
        execute(cmd)
        
        # Move executable to run dir
        if COPY_EXECUTABLE:
            logging.info('Copying executable to run directory')
            cmd = 'cp {0} {1}'.format(os.path.join(SCM_ROOT, SCM_BIN, EXECUTABLE_NAME), os.path.join(SCM_RUN, EXECUTABLE_NAME))
            execute(cmd)
        else:
            logging.info('Linking executable to run directory')
            cmd = 'ln -s {0} {1}'.format(os.path.join(SCM_ROOT, SCM_BIN, EXECUTABLE_NAME), os.path.join(SCM_RUN, EXECUTABLE_NAME))
            execute(cmd)
        
        return output_dir

def launch_executable(use_gdb, gdb):
    """Configure model run command and pass control to shell/gdb"""
    if use_gdb:
        cmd = '(cd {scm_run} && exec {gdb} {executable})'.format(scm_run=SCM_RUN, gdb=gdb, executable=EXECUTABLE)
    else:
        cmd = '(cd {scm_run} && exec {executable})'.format(scm_run=SCM_RUN, executable=EXECUTABLE)
    logging.info('Passing control to "{0}"'.format(cmd))
    time.sleep(2)
    failure = os.system(cmd)
    if failure:
        logging.info('Execution of {0} failed!\n'.format(cmd))
        sys.exit(1)
    
def copy_outdir(exp_dir):
    """Copy output directory to /home for this experiment."""
    home_output_dir = '/home/'+exp_dir
    if os.path.isdir(home_output_dir):
        shutil.rmtree(home_output_dir)
    shutil.copytree(exp_dir, home_output_dir)

def main():
    (case, use_gdb, suite_name, namelist, docker, tracers, runtime, runtime_mult, levels, npz_type, vert_coord_file, case_data_dir, n_itt_out, n_itt_diag, run_dir, timestep) = parse_arguments()
    
    global SCM_ROOT
    SCM_ROOT = os.getenv('SCM_ROOT')
    if SCM_ROOT is None:
        message = 'The SCM_ROOT environment variable is not set. Please set the SCM_ROOT environment variable to the top-level path where the model was cloned.'
        logging.critical(message)
        raise Exception(message)
    
    global SCM_BIN
    SCM_BIN = os.path.join(SCM_ROOT, DEFAULT_BIN_DIR)
    
    global SCM_RUN
    if run_dir:
        SCM_RUN = run_dir
    else:
        SCM_RUN = os.path.join(SCM_ROOT, DEFAULT_RUN_DIR)
    if not os.path.isdir(SCM_RUN):
        os.makedirs(SCM_RUN)
    
    global EXECUTABLE
    EXECUTABLE = os.path.join(SCM_RUN, EXECUTABLE_NAME)
    
    setup_logging()
    
    #Experiment
    if namelist:
        if tracers:
            logging.info('Setting up experiment {0} with suite {1} using namelist {2} and tracers {3}'.format(case,suite_name,namelist,tracers))
        else:
            logging.info('Setting up experiment {0} with suite {1} using namelist {2} using default tracers for the suite'.format(case,suite_name,namelist))
    else:
        if tracers:
            logging.info('Setting up experiment {0} with suite {1} using the default namelist for the suite and tracers {2}'.format(case,suite_name,tracers))
        else:
            logging.info('Setting up experiment {0} with suite {1} using the default namelist and tracers for the suite'.format(case,suite_name))
    
    active_suite = None
    for s in suite_list:
        if suite_name == s._name:
            active_suite = s
            break
        
    if (active_suite is None):
        if (namelist and tracers):
            if timestep:
                active_suite = suite(suite_name, tracers, namelist, timestep, -1, False)
            else:
                active_sutie = suite(suite_name, tracers, namelist, -1, -1, False)
        else:
            message = 'The given suite ({0}), does not have defaults set in suite_info.py and either the tracers file or physics namelist file (or both) were not provided.'.format(suite_name)
            logging.critical(message)
            raise Exception(message)
    else:
        if namelist:
            active_suite.namelist = namelist
        if tracers:
            active_suite.tracers = tracers
        if timestep:
            active_suite.timestep = timestep
    
    exp = Experiment(case, active_suite, runtime, runtime_mult, levels, npz_type, vert_coord_file, case_data_dir, n_itt_out, n_itt_diag)
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
