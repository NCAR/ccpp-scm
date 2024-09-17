#!/usr/bin/env python

import sys, os

#DEFAULT_SUITE_BEHAVIOR = 'supported'
DEFAULT_SUITE_BEHAVIOR = 'regression_test'

class suite(object):
  
    DEFAULT_MAX_TIMESTEP = 1800.0
  
    def __init__(self, name, tracers, namelist, timestep, max_timestep, supported):
        DEFAULT_MAX_TIMESTEP = 1800.0
        self._name = name                     #should remain unchanged after init
        self._default_tracers = tracers       #should remain unchanged after init
        self._default_namelist = namelist     #should remain unchanged after init
        self.tracers = self._default_tracers   #can be modified after init
        self.namelist = self._default_namelist #can be modified after init
        self._supported = supported
      
        if max_timestep > 0:
            self._max_timestep = max_timestep     #should remain unchanged after init
        else:
            self._max_timestep = DEFAULT_MAX_TIMESTEP
        
        if timestep <= self._max_timestep and timestep > 0:
            self._default_timestep = timestep #should remain unchanged after init
            self.timestep = self._default_timestep #can be modified after init
        elif timestep <= 0:
            self.timestep = None #use the default dt in the SCM code
        else:
            message = 'The timestep for suite {0} cannot be set greater than the max_timestep of {1}'.format(self._name, self._max_timestep)
            raise Exception(message)
      
        @property
        def timestep(self):
            """Get the timestep for the given suite."""
            return self.timestep

        @timestep.setter
        def timestep(self, value):
            """Set the timestep for the given suite."""
            if value <= self._max_timestep:
                self.timestep = value
            else:
                message = 'The timestep for suite {0} cannot be set greater than the max_timestep of {1}'.format(self._name, self._max_timestep)
                raise Exception(message)
      
suite_list = []
suite_list.append(suite('SCM_GFS_v16',           'tracers_GFS_v16.txt',                  'input_GFS_v16.nml',                 600.0, 1800.0, True ))
suite_list.append(suite('SCM_GFS_v17_p8_ugwpv1', 'tracers_GFS_v17_p8_ugwpv1.txt',        'input_GFS_v17_p8_ugwpv1.nml',       600.0, 600.0,  True ))
suite_list.append(suite('SCM_RAP',               'tracers_RAP.txt',                      'input_RAP.nml',                     600.0, 600.0 , True ))
suite_list.append(suite('SCM_HRRR_gf',           'tracers_HRRR_gf.txt',                  'input_HRRR_gf.nml',                 600.0, 600.0 , True ))
suite_list.append(suite('SCM_WoFS_v0',           'tracers_WoFS_v0.txt',                  'input_WoFS_v0.nml',                 600.0, 600.0 , True ))
suite_list.append(suite('SCM_GFS_v16_RRTMGP',    'tracers_GFS_v16.txt',                  'input_GFS_v16_RRTMGP.nml',          600.0, 1800.0, True ))

suite_list.append(suite('SCM_GFS_v17_p8',        'tracers_GFS_v17_p8.txt',               'input_GFS_v17_p8.nml',              600.0, 600.0,  False))
suite_list.append(suite('SCM_HRRR',              'tracers_HRRR.txt',                     'input_HRRR.nml',                    600.0, 600.0 , False))
suite_list.append(suite('SCM_RRFS_v1beta',       'tracers_RRFS_v1beta.txt',              'input_RRFS_v1beta.nml',             600.0, 600.0 , False))
suite_list.append(suite('SCM_GFS_v15p2',         'tracers_GFS_v15p2.txt',                'input_GFS_v15p2.nml',               600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_RRTMGP',  'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_RRTMGP.nml',        600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_no_nsst', 'tracers_GFS_v15p2.txt',                'input_GFS_v15p2.nml',               600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_noahmp',  'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_noahmp.nml',        600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_MYJ',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_MYJ.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_YSU',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_YSU.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_saYSU',   'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_saYSU.nml',         600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_ACM',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_ACM.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v16_no_nsst',   'tracers_GFS_v16.txt',                  'input_GFS_v16.nml',                 600.0, 1800.0, False))
suite_list.append(suite('HAFS_v0_hwrf',          'tracers_HAFS_v0_hwrf.txt',             'input_HAFS_v0_hwrf.nml',            600.0, 1800.0, False))
suite_list.append(suite('HAFS_v0_hwrf_thompson', 'tracers_HAFS_v0_hwrf_thompson.txt',    'input_HAFS_v0_hwrf_thompson.nml',   600.0, 600.0 , False))
suite_list.append(suite('SCM_GSD_v1nssl',        'tracers_gsd_nssl.txt',                 'input_GSD_v1nssl.nml',              600.0, 600.0 , False))
suite_list.append(suite('SCM_GSD_v1',            'tracers_gsd.txt',                      'input_GSD_v1.nml',                  600.0, 600.0 , False))
suite_list.append(suite('SCM_RRFS_v1nssl',       'tracers_RRFS_v1nssl_nohail_noccn.txt', 'input_RRFS_v1nssl_nohailnoccn.nml', 600.0, 600.0 , False))
suite_list.append(suite('SCM_csawmg',            'tracers_csawmg.txt',                   'input_csawmg.nml',                  600.0, 1800.0, False))

def main():
    
    #print supported suites separated by commas
    suite_string = ''
    
    if DEFAULT_SUITE_BEHAVIOR == 'regression_test':
        dir_path = os.path.dirname(os.path.realpath(__file__))
        sys.path.insert(1, dir_path + '/../../test/')
        
        rt_suite_list = []
        
        import rt_test_cases
        import rt_test_cases_sp
        import rt_test_cases_nvidia
        
        for item in rt_test_cases.run_list:
            rt_suite_list.append(item.get("suite"))
        
        for item in rt_test_cases_sp.run_list:
            rt_suite_list.append(item.get("suite"))
        
        for item in rt_test_cases_nvidia.run_list:
            rt_suite_list.append(item.get("suite"))
        
        unique_suite_list = list(set(rt_suite_list))
        
        for s in unique_suite_list:    
            suite_string += s + ',' + s + '_ps' + ','
    else:
        for s in suite_list:
            if s._supported:
                suite_string += s._name + ',' + s._name + '_ps' + ','  
    print(suite_string[:-1])

if __name__ == '__main__':
    main()  
 
