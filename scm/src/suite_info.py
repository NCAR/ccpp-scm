class suite(object):
  
    DEFAULT_MAX_TIMESTEP = 1800.0
  
    def __init__(self, name, tracers, namelist, timestep, max_timestep, supported):
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
suite_list.append(suite('SCM_GFS_v15p2',         'tracers_GFS_v15p2.txt',                'input_GFS_v15p2.nml',               600.0, 1800.0, True ))
suite_list.append(suite('SCM_GFS_v16',           'tracers_GFS_v16.txt',                  'input_GFS_v16.nml',                 600.0, 1800.0, True ))
suite_list.append(suite('SCM_GSD_v1',            'tracers_gsd.txt',                      'input_GSD_v1.nml',                  600.0, 600.0 , True ))
suite_list.append(suite('SCM_RRFS_v1alpha',      'tracers_RRFS_v1alpha.txt',             'input_RRFS_v1alpha.nml',            600.0, 600.0 , True ))
suite_list.append(suite('SCM_csawmg',            'tracers_csawmg.txt',                   'input_csawmg.nml',                  600.0, 1800.0, True ))

suite_list.append(suite('SCM_GFS_v15p2_RRTMGP',  'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_RRTMGP.nml',        600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_no_nsst', 'tracers_GFS_v15p2.txt',                'input_GFS_v15p2.nml',               600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_noahmp',  'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_noahmp.nml',        600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_MYJ',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_MYJ.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_YSU',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_YSU.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_saYSU',   'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_saYSU.nml',         600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v15p2_ACM',     'tracers_GFS_v15p2.txt',                'input_GFS_v15p2_ACM.nml',           600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v16_RRTMGP',    'tracers_GFS_v16.txt',                  'input_GFS_v16_RRTMGP.nml',          600.0, 1800.0, False))
suite_list.append(suite('SCM_GFS_v16_no_nsst',   'tracers_GFS_v16.txt',                  'input_GFS_v16.nml',                 600.0, 1800.0, False))
suite_list.append(suite('HAFS_v0_hwrf',          'tracers_HAFS_v0_hwrf.txt',             'input_HAFS_v0_hwrf.nml',            600.0, 1800.0, False))
suite_list.append(suite('HAFS_v0_hwrf_thompson', 'tracers_HAFS_v0_hwrf_thompson.txt',    'input_HAFS_v0_hwrf_thompson.nml',   600.0, 600.0 , False))
suite_list.append(suite('SCM_GSD_v1nssl',        'tracers_gsd_nssl.txt',                 'input_GSD_v1nssl.nml',              600.0, 600.0 , False))
suite_list.append(suite('SCM_RRFS_v1nssl',       'tracers_RRFS_v1nssl_nohail_noccn.txt', 'input_RRFS_v1nssl_nohailnoccn.nml', 600.0, 600.0 , False))
  
  
 
