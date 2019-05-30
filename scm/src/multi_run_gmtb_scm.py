#!/usr/bin/env python

import argparse
import imp
import os
import logging
import run_gmtb_scm
from supported_suites import suites
from supported_cases import cases

# For developers: set logging level to DEBUG for additional output
#LOGLEVEL = logging.DEBUG
LOGLEVEL = logging.INFO

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()#required=True)
group.add_argument('-c', '--case',       help='name of case to run',)
group.add_argument('-s', '--suite',      help='name of suite to use',)
group.add_argument('-f', '--file',      help='name of file where SCM runs are defined',)

def setup_logging():
    """Sets up the logging module."""
    logging.basicConfig(format='%(levelname)s: %(message)s', level=LOGLEVEL)

def main():
    args = parser.parse_args()
    
    setup_logging()
    
    # if the case argument is specified, run through all supported suites with the specified case
    if args.case:
        logging.info('Running all supported suites with case {0}'.format(args.case))
        for suite in suites:
            run_gmtb_scm.main(case=args.case, suite=suite)
    
    # if the suite argument is specified, run through all supported cases with the specified suite
    if args.suite:
        logging.info('Running all supported cases with suite {0}'.format(args.suite))
        for case in cases:
            run_gmtb_scm.main(case=case, suite=args.suite)
    
    # For maximum flexibility, run the SCM as specified from an external file where cases, suites, and physics namelists are all specified. 
    # This file must contain python lists called 'cases','suites', and 'namelists'. The suites and namelists lists can be empty ([]) if necessary.
    #The following rules apply:
    # 1. The case list in the file must not be empty.
    # 2. If only a case list is specified, the cases are run with the default suite specified in run_gmtb_scm.py with the default namelists specified in default_namelists.py.
    # 3. If a case list and suite list is provided without a namelist list, all permutations of cases and suites will be run using default namelists specified in default_namelists.py.
    # 4. If a case list and suite list is provided with a namelist:
    # 4a. If only one suite is specified, it can be run with any number of namelists.
    # 4b. If more than one suite is specified, the number of namelists must match, and each case is run with each (suite,namelist) pair, by order specified in the lists.
    # 5. If a case list and namelist list are specified without a suite list, each case is run with the default suite specified in run_gmtb_scm.py using the supplied namelists.
    if args.file:
        logging.info('Importing {0} to run requested combinations'.format(args.file))
        try:
            scm_runs = imp.load_source(os.path.splitext(args.file)[0], args.file)
        except ImportError:
            message = 'There was a problem loading {0}. Please check that the path exists.'.format(args.file)
            logging.critical(message)
            raise
                
        if not scm_runs.cases:
            message = 'The cases list in {0} must not be empty'.format(args.file)
            logging.critical(message)
            raise Exception(message)
        
        if scm_runs.cases and not scm_runs.suites and not scm_runs.namelists:
            logging.info('Only cases were specified in {0}, so running all cases with the default suite'.format(args.file))
            for case in cases:
                run_gmtb_scm.main(case=case)
        
        if scm_runs.cases and scm_runs.suites:
            if scm_runs.namelists:
                if len(scm_runs.suites) == 1:
                    logging.info('Cases and namelists were specified with 1 suite in {0}, so running all cases with the suite {1} for all specified namelists'.format(args.file, scm_runs.suites[0]))
                    for case in scm_runs.cases:
                        for namelist in scm_runs.namelists:
                            run_gmtb_scm.main(case=case, suite=scm_runs.suites[0], namelist=namelist)
                elif len(scm_runs.suites) == len(scm_runs.namelists):
                    logging.info('Cases, suites, and namelists were specified in {0}, so running all cases with all suites, matched with namelists by order'.format(args.file))
                    for case in scm_runs.cases:
                        for num, suite in enumerate(scm_runs.suites):
                            run_gmtb_scm.main(case=case, suite=suite, namelist=scm_runs.namelists[num])
                else:
                    message = 'The number of suites and namelists specified in {0} is incompatible. Either use one suite with many namelists or the number of suites must match the number of namelists provided.'.format(args.file)
                    logging.critical(message)
                    raise Exception(message)
            else:
                logging.info('Cases and suites specified in {0}, so running all cases with all suites using default namelists for each suite'.format(args.file))
                for case in scm_runs.cases:
                    for suite in scm_runs.suites:
                        run_gmtb_scm.main(case=case, suite=suite)
                        
        if scm_runs.cases and not scm_runs.suites and scm_runs.namelists:
            logging.info('Cases and namelists were specified in {0}, so running all cases with the default suite using the list of namelists'.format(args.file))
            for case in scm_runs.cases:
                for namelist in scm_runs.namelists:
                    run_gmtb_scm.main(case=case, namelist=namelist)
        
    # If running the script with no arguments, run all supported (case,suite) permutations.
    if not args.case and not args.suite and not args.file:
        logging.info('Since no arguments were specified, running through all permuatations of supported cases and suites')
        for case in cases:
            for suite in suites:
                run_gmtb_scm.main(case=case, suite=suite)
        
if __name__ == '__main__':
    main()