#!/usr/bin/env python

import argparse
import logging
from netCDF4 import Dataset
import os
import f90nml
import numpy as np
from scipy import interpolate
import forcing_file_common as ffc
import shutil

###############################################################################
# Global settings                                                             #
###############################################################################

SCM_ROOT = ''

# Path to default run directory (relative to scm_root)
DEFAULT_RUN_DIR = 'scm/run'

DEFAULT_CASE_INPUT_DIR = 'scm/data/processed_case_input'

DEFAULT_CASE_CONFIG_DIR = 'scm/etc/case_config'

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--case',       help='name of case to modify forcing', required=True)
parser.add_argument('--lwrad',            help='add longwave radiation tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--swrad',            help='add shortwave radiation tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--pbl',              help='add plentary boundary layer tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--deepconv',         help='add deep convective tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--shalconv',         help='add shallow convective tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--micro',            help='add microphysics tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--ogwd',             help='add orographic gravity wave drag tendency to forcing terms', action='store_true', default=False)
parser.add_argument('--cgwd',             help='add convective gravity wave drag tendency to forcing terms', action='store_true', default=False)
parser.add_argument('-v', '--verbose',    help='once: set logging level to debug; twice: set logging level to debug '\
                                               'and write log to file', action='count', default=0)

###############################################################################
# Functions and subroutines                                                   #
###############################################################################
def setup_logging(verbose):
    """Sets up the logging module."""
    # print out debug messages (logs and output from subprocesses) if verbose argument is set
    if verbose==2:
        LOG_LEVEL = logging.DEBUG
    elif verbose==1:
        LOG_LEVEL = logging.INFO
    else:
        LOG_LEVEL = logging.WARNING
    LOG_FILE = 'mod_forcing.log'
    LOG_FORMAT = '%(levelname)s: %(message)s'

    logging.basicConfig(format=LOG_FORMAT, level=LOG_LEVEL)

    # write out a log file if verbosity is set twice (-vv)
    if verbose > 1:
        fh = logging.FileHandler(LOG_FILE, mode='w')
        logger = logging.getLogger()
        logger.addHandler(fh)

def parse_arguments():
    args = parser.parse_args()
    case = args.case
    verbose = args.verbose
    
    terms = []
    if args.lwrad:
        terms.append('lwrad')
    if args.swrad:
        terms.append('swrad')
    if args.pbl:
        terms.append('pbl')
    if args.deepconv:
        terms.append('deepconv')
    if args.shalconv:
        terms.append('shalconv')
    if args.micro:
        terms.append('micro')
    if args.ogwd:
        terms.append('ogwd')
    if args.cgwd:
        terms.append('cgwd')
    
    return (case, terms, verbose)

def read_tendency_terms(files, terms):
    
    time_inst_all_files = []
    time_diag_all_files = []
    pres_l_all_files = []
    
    dT_dt_all_files = []
    dq_dt_all_files = []
    du_dt_all_files = []
    dv_dt_all_files = []
    
    for file in files:
        
        dT_dt_terms = []
        dq_dt_terms = []
        du_dt_terms = []
        dv_dt_terms = []
        
        nc_fid = Dataset(os.path.join(SCM_ROOT, DEFAULT_RUN_DIR, file), 'r')
        nc_fid.set_auto_mask(False)
    
        time_inst_all_files.append(nc_fid.variables['time_inst'][:])
        time_diag_all_files.append(nc_fid.variables['time_diag'][:])
        pres_l_all_files.append(nc_fid.variables['pres'][:])  #on instantaneous time interval
        
        for term in terms:
            var_label = 'dT_dt_' + term
            try:
                data = (nc_fid.variables[var_label][:])
                dT_dt_terms.append(np.swapaxes(data, 0, 1))
                logging.debug('{0} variable found in {1}'.format(var_label,file))
            except KeyError:
                logging.debug('{0} variable not found in {1}'.format(var_label,file))
            
            var_label = 'dq_dt_' + term
            try:
                data = nc_fid.variables[var_label][:]
                dq_dt_terms.append(np.swapaxes(data, 0, 1))
                logging.debug('{0} variable found in {1}'.format(var_label,file))
            except KeyError:
                logging.debug('{0} variable not found in {1}'.format(var_label,file))
            
            var_label = 'du_dt_' + term
            try:
                data = nc_fid.variables[var_label][:]
                du_dt_terms.append(np.swapaxes(data, 0, 1))
                logging.debug('{0} variable found in {1}'.format(var_label,file))
            except KeyError:
                logging.debug('{0} variable not found in {1}'.format(var_label,file))
            
            var_label = 'dv_dt_' + term
            try:
                data = nc_fid.variables[var_label][:]
                dv_dt_terms.append(np.swapaxes(data, 0, 1))
                logging.debug('{0} variable found in {1}'.format(var_label,file))
            except KeyError:
                logging.debug('{0} variable not found in {1}'.format(var_label,file))
        nc_fid.close()
        
        dT_dt_all_files.append(dT_dt_terms)
        dq_dt_all_files.append(dq_dt_terms)
        du_dt_all_files.append(du_dt_terms)
        dv_dt_all_files.append(dv_dt_terms)
    
    #check whether all files have the same number of terms
    it = iter(dT_dt_all_files)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('Not all output files have the same number of temperature tendency terms in the set {}'.format(terms))
    
    it = iter(dq_dt_all_files)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('Not all output files have the same number of water vapor tendency terms in the set {}'.format(terms))
    
    it = iter(du_dt_all_files)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('Not all output files have the same number of u-wind tendency terms in the set {}'.format(terms))
    
    it = iter(dv_dt_all_files)
    the_len = len(next(it))
    if not all(len(l) == the_len for l in it):
        raise ValueError('Not all output files have the same number of v-wind tendency terms in the set {}'.format(terms))
    
    return (time_inst_all_files, time_diag_all_files, pres_l_all_files, dT_dt_all_files, dq_dt_all_files, du_dt_all_files, dv_dt_all_files)

def read_forcing_terms(case):
    case_nml = f90nml.read(os.path.join(SCM_ROOT, DEFAULT_CASE_CONFIG_DIR, case + '.nml'))
    
    try:
        input_type = case_nml['case_config']['input_type']
        if input_type == 1:
            logging.debug('{} uses the DEPHY case data format.')
            message = 'Editing DEPHY-formatted cases is not supported yet. Exiting.'
            logging.critical(message)
            raise Exception(message)
    except KeyError:
        logging.debug('{} uses the proprietary case data format.'.format(case))
        logging.warning('Since {} uses the proprietary case data format, straightforward application of momentum tendencies is not possible, so momentum tendency terms will be ignored.'.format(case))
        input_type = 0
    
    if input_type == 0:
        try:
            thermo_forcing_type = case_nml['case_config']['thermo_forcing_type']
        except KeyError:
            thermo_forcing_type = 2 #default in scm_input.F90
    
        nc_fid = Dataset(os.path.join(SCM_ROOT, DEFAULT_CASE_INPUT_DIR, case + '.nc'), 'r')
        nc_fid.set_auto_mask(False)
    
        time = nc_fid.variables['time'][:]
        levels = nc_fid.variables['levels'][:]
        h_advec_thil = nc_fid.groups['forcing'].variables['h_advec_thetail'][:] #levels, time
        h_advec_qt = nc_fid.groups['forcing'].variables['h_advec_qt'][:]
        
        nc_fid.close()
    
    return (time, levels, h_advec_thil, h_advec_qt)

def calc_new_forcing(time_inst, time_diag, pres, t_tend, qv_tend, u_tend, v_tend, forcing_time, forcing_levels):
    num_tend_datasets = len(time_inst)
    
    t_tend_final = np.zeros((forcing_levels.shape[0],forcing_time.shape[0],num_tend_datasets))
    thil_tend_final = np.zeros((forcing_levels.shape[0],forcing_time.shape[0],num_tend_datasets))
    qv_tend_final = np.zeros((forcing_levels.shape[0],forcing_time.shape[0],num_tend_datasets))
    u_tend_final = np.zeros((forcing_levels.shape[0],forcing_time.shape[0],num_tend_datasets))
    v_tend_final = np.zeros((forcing_levels.shape[0],forcing_time.shape[0],num_tend_datasets))
    
    for i in range(num_tend_datasets):
        #collect all tendency terms into one array for interpolation
        if(len(t_tend[i]) > 0):
            t_tend_all = np.zeros((t_tend[i][0].shape[0],t_tend[i][0].shape[1])) #levels, time for this file
            for tend in t_tend[i]:
                t_tend_all = t_tend_all + tend[:,:,0]
        
        if(len(qv_tend[i]) > 0):
            qv_tend_all = np.zeros((qv_tend[i][0].shape[0],qv_tend[i][0].shape[1])) #levels, time for this file
            for tend in qv_tend[i]:
                qv_tend_all = qv_tend_all + tend[:,:,0]
        
        if(len(t_tend[i]) > 0):
            u_tend_all = np.zeros((u_tend[i][0].shape[0],u_tend[i][0].shape[1])) #levels, time for this file
            for tend in u_tend[i]:
                u_tend_all = u_tend_all + tend[:,:,0]
        
        if(len(v_tend[i]) > 0):
            v_tend_all = np.zeros((v_tend[i][0].shape[0],v_tend[i][0].shape[1])) #levels, time for this file
            for tend in v_tend[i]:
                v_tend_all = v_tend_all + tend[:,:,0]
        
        if ('t_tend_all' in locals()):
            t_tend_all_on_forcing_levels = np.zeros((forcing_levels.shape[0],t_tend[i][0].shape[1]))
            thil_tend_all_on_forcing_levels = np.zeros((forcing_levels.shape[0],t_tend[i][0].shape[1]))
        
        if ('qv_tend_all' in locals()):
            qv_tend_all_on_forcing_levels = np.zeros((forcing_levels.shape[0],qv_tend[i][0].shape[1]))
        
        if ('u_tend_all' in locals()):
            u_tend_all_on_forcing_levels = np.zeros((forcing_levels.shape[0],u_tend[i][0].shape[1]))
        
        if ('v_tend_all' in locals()):
            v_tend_all_on_forcing_levels = np.zeros((forcing_levels.shape[0],v_tend[i][0].shape[1]))
                
        #interpolate the output tendencies in the vertical to the forcing levels
        for t_index, t in enumerate(time_diag[i]):
            #find the instantaneous index corresponding to t from time_diag[i]
            inst_time_index = np.where(time_inst[i] == t)[0][0]
            
            #find the pressure levels for tendency terms for interpolation
            diag_levels = pres[i][inst_time_index,:,0]
            
            if ('t_tend_all' in locals()): 
                t_tend_all_f = interpolate.interp1d(diag_levels, t_tend_all[:,t_index], bounds_error=False, fill_value=0.0)
                t_tend_all_on_forcing_levels[:,t_index] = t_tend_all_f(forcing_levels.tolist())
                thil_tend_all_on_forcing_levels[:,t_index] = (ffc.p0/forcing_levels)**(ffc.R_dry/ffc.c_p)*t_tend_all_on_forcing_levels[:,t_index]
                
            if ('qv_tend_all' in locals()): 
                qv_tend_all_f = interpolate.interp1d(diag_levels, qv_tend_all[:,t_index], bounds_error=False, fill_value=0.0)
                qv_tend_all_on_forcing_levels[:,t_index] = qv_tend_all_f(forcing_levels.tolist())
            
            if ('u_tend_all' in locals()): 
                u_tend_all_f = interpolate.interp1d(diag_levels, u_tend_all[:,t_index], bounds_error=False, fill_value=0.0)
                u_tend_all_on_forcing_levels[:,t_index] = u_tend_all_f(forcing_levels.tolist())
            
            if ('v_tend_all' in locals()): 
                v_tend_all_f = interpolate.interp1d(diag_levels, v_tend_all[:,t_index], bounds_error=False, fill_value=0.0)
                v_tend_all_on_forcing_levels[:,t_index] = v_tend_all_f(forcing_levels.tolist())
                
        for t_index, force_t in enumerate(forcing_time):
            if t_index == 0:
                #find how much state variables have changed between forcing_time[t] and forcing_time[t+1]
                begin_forcing_time = force_t
                end_forcing_time = forcing_time[t_index+1]
                            
                begin_forcing_interval_diag_index = np.where(time_diag[i] >= begin_forcing_time)[0][0]
                end_forcing_interval_diag_index = np.where(time_diag[i] >= end_forcing_time)[0][0]
                
                T_change = np.zeros(forcing_levels.shape[0])
                thil_change = np.zeros(forcing_levels.shape[0])
                qv_change = np.zeros(forcing_levels.shape[0])
                u_change = np.zeros(forcing_levels.shape[0])
                v_change = np.zeros(forcing_levels.shape[0])
                accum_diag_time = 0.0
                for t in range(begin_forcing_interval_diag_index,end_forcing_interval_diag_index+1):
                    if t == 0:
                        if ('t_tend_all' in locals()):
                            T_change = t_tend_all_on_forcing_levels[:,t]*(time_diag[i][t])
                            thil_change = thil_tend_all_on_forcing_levels[:,t]*(time_diag[i][t])
                        if ('qv_tend_all' in locals()):
                            qv_change = qv_tend_all_on_forcing_levels[:,t]*(time_diag[i][t])
                        if ('u_tend_all' in locals()):
                            u_change = u_tend_all_on_forcing_levels[:,t]*(time_diag[i][t])
                        if ('v_tend_all' in locals()):
                            v_change = v_tend_all_on_forcing_levels[:,t]*(time_diag[i][t])
                        accum_diag_time = time_diag[i][t]
                    else:
                        if ('t_tend_all' in locals()):
                            T_change = T_change + t_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                            thil_change = thil_change + thil_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                        if ('qv_tend_all' in locals()):
                            qv_change = qv_change + qv_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                        if ('u_tend_all' in locals()):
                            u_change = u_change + u_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                        if ('v_tend_all' in locals()):
                            v_change = v_change + v_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                        accum_diag_time = accum_diag_time + (time_diag[i][t]-time_diag[i][t-1])
                
                if ('t_tend_all' in locals()):
                    t_tend_final[:,t_index,i] = T_change/accum_diag_time
                    thil_tend_final[:,t_index,i] = thil_change/accum_diag_time
                if ('qv_tend_all' in locals()):
                    qv_tend_final[:,t_index,i] = qv_change/accum_diag_time
                if ('u_tend_all' in locals()):
                    u_tend_final[:,t_index,i] = u_change/accum_diag_time
                if ('v_tend_all' in locals()):
                    v_tend_final[:,t_index,i] = v_change/accum_diag_time
            elif t_index < forcing_time.shape[0]-1:
                begin_forcing_time = force_t
                end_forcing_time = forcing_time[t_index+1]
                
                begin_forcing_interval_diag_index = np.where(time_diag[i] > begin_forcing_time)[0][0]
                end_forcing_interval_diag_index = np.where(time_diag[i] >= end_forcing_time)[0][0]
                
                T_change = np.zeros(forcing_levels.shape[0])
                thil_change = np.zeros(forcing_levels.shape[0])
                qv_change = np.zeros(forcing_levels.shape[0])
                u_change = np.zeros(forcing_levels.shape[0])
                v_change = np.zeros(forcing_levels.shape[0])
                accum_diag_time = 0.0
                for t in range(begin_forcing_interval_diag_index,end_forcing_interval_diag_index+1):
                    if ('t_tend_all' in locals()):
                        T_change = T_change + t_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                        thil_change = thil_change + thil_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                    if ('qv_tend_all' in locals()):
                        qv_change = qv_change + qv_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                    if ('u_tend_all' in locals()):
                        u_change = u_change + u_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                    if ('v_tend_all' in locals()):
                        v_change = v_change + v_tend_all_on_forcing_levels[:,t]*(time_diag[i][t]-time_diag[i][t-1])
                    accum_diag_time = accum_diag_time + (time_diag[i][t]-time_diag[i][t-1])
                                    
                if ('t_tend_all' in locals()):
                    t_tend_final[:,t_index,i] = 2.0*T_change/accum_diag_time - t_tend_final[:,t_index-1,i]
                    thil_tend_final[:,t_index,i] = 2.0*thil_change/accum_diag_time - thil_tend_final[:,t_index-1,i]
                if ('qv_tend_all' in locals()):
                    qv_tend_final[:,t_index,i] = 2.0*qv_change/accum_diag_time - qv_tend_final[:,t_index-1,i]
                if ('u_tend_all' in locals()):
                    u_tend_final[:,t_index,i] = 2.0*u_change/accum_diag_time - u_tend_final[:,t_index-1,i]
                if ('v_tend_all' in locals()):
                    v_tend_final[:,t_index,i] = 2.0*v_change/accum_diag_time - v_tend_final[:,t_index-1,i]
            else:
                if ('t_tend_all' in locals()):
                    t_tend_final[:,t_index,i] = t_tend_final[:,t_index-1,i]
                    thil_tend_final[:,t_index,i] = thil_tend_final[:,t_index-1,i]
                if ('qv_tend_all' in locals()):
                    qv_tend_final[:,t_index,i] = qv_tend_final[:,t_index-1,i]
                if ('u_tend_all' in locals()):
                    u_tend_final[:,t_index,i] = u_tend_final[:,t_index-1,i]
                if ('v_tend_all' in locals()):
                    v_tend_final[:,t_index,i] = v_tend_final[:,t_index-1,i]
                
    new_t_forcing = np.mean(t_tend_final, axis=2)
    new_thil_forcing = np.mean(thil_tend_final, axis=2)
    new_qv_forcing = np.mean(qv_tend_final, axis=2)
    new_u_forcing = np.mean(u_tend_final, axis=2)
    new_v_forcing = np.mean(v_tend_final, axis=2)
    
    return (new_thil_forcing, new_qv_forcing, new_u_forcing, new_v_forcing)
    
#def calc_new_forcing_dephy():

def write_new_forcing(case, new_thil_forcing, new_qv_forcing, new_u_forcing, new_v_forcing):
        shutil.copy2(os.path.join(SCM_ROOT, DEFAULT_CASE_INPUT_DIR, case + '.nc'),os.path.join(SCM_ROOT, DEFAULT_CASE_INPUT_DIR, case + '_mod.nc'))
        
        nc_fid = Dataset(os.path.join(SCM_ROOT, DEFAULT_CASE_INPUT_DIR, case + '_mod.nc'), 'r+')
        
        nc_fid.groups['forcing'].variables['h_advec_thetail'][:] = nc_fid.groups['forcing'].variables['h_advec_thetail'][:] + new_thil_forcing
        nc_fid.groups['forcing'].variables['h_advec_qt'][:] = nc_fid.groups['forcing'].variables['h_advec_qt'][:] + new_qv_forcing
                
        nc_fid.close()
        
        logging.warning('The new forcing was written to {}'.format(os.path.join(SCM_ROOT, DEFAULT_CASE_INPUT_DIR, case + '_mod.nc')))

def main():
    global SCM_ROOT
    SCM_ROOT = os.getenv('SCM_ROOT')
    if SCM_ROOT is None:
        message = 'The SCM_ROOT environment variable is not set. Please set the SCM_ROOT environment variable to the top-level path where the model was cloned.'
        logging.critical(message)
        raise Exception(message)
    
    (case, terms, verbose) = parse_arguments()
    
    setup_logging(verbose)
    
    output_files = ['output_' + case + '_SCM_GFS_v16/output.nc', 'output_' + case + '_SCM_GFS_v17_p8/output.nc', 'output_' + case + '_SCM_RAP/output.nc']
    
    if len(terms) > 0:
        (time_inst, time_diag, pres, t_tend, qv_tend, u_tend, v_tend) = read_tendency_terms(output_files, terms)
    else:
        logging.warning('No tendencies were specified in the command line arguments to add to the existing forcing of {}'.format(case))
        logging.warning('No new file forcing file was generated.')
        exit()
        
    (forcing_time, forcing_levels, forcing_h_advec_thil, forcing_h_advec_qt) = read_forcing_terms(case)
    
    (new_thil_forcing, new_qv_forcing, new_u_forcing, new_v_forcing) = calc_new_forcing(time_inst, time_diag, pres, t_tend, qv_tend, u_tend, v_tend, forcing_time, forcing_levels)
    
    write_new_forcing(case, new_thil_forcing, new_qv_forcing, new_u_forcing, new_v_forcing)
    
if __name__ == '__main__':
    main()
