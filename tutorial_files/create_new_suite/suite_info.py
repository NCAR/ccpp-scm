#!/usr/bin/env python

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
suite_list.append(suite('SCM_GFS_v16',           'tracers_GFS_v16.txt',                  'input_GFS_v16.nml',                 600.0, 1800.0, True ))
suite_list.append(suite('SCM_GFS_v17_p8',        'tracers_GFS_v17_p8.txt',               'input_GFS_v17_p8.nml',              600.0, 600.0,  True ))
suite_list.append(suite('SCM_RAP',               'tracers_RAP.txt',                      'input_RAP.nml',                     600.0, 600.0 , True ))
suite_list.append(suite('SCM_RRFS_v1beta',       'tracers_RRFS_v1beta.txt',              'input_RRFS_v1beta.nml',             600.0, 600.0 , True ))
suite_list.append(suite('SCM_WoFS_v0',           'tracers_WoFS_v0.txt',                  'input_WoFS_v0.nml',                 600.0, 600.0 , True ))
suite_list.append(suite('SCM_HRRR',              'tracers_HRRR.txt',                     'input_HRRR.nml',                    600.0, 600.0 , True ))
suite_list.append(suite('SCM_RRFS_v1beta_sas',        'tracers_RRFS_v1beta.txt',         'input_RRFS_v1beta_sas.nml',         600.0, 600.0 , False))
suite_list.append(suite('SCM_RRFS_v1beta_sas_sfcmod', 'tracers_RRFS_v1beta.txt',         'input_RRFS_v1beta_sas.nml',         600.0, 600.0 , False))

def main():
    
    #print supported suites separated by commas
    suite_string = ''
    for s in suite_list:
        if s._supported:
            suite_string += s._name + ',' + s._name + '_ps' + ','  
    print(suite_string[:-1])

if __name__ == '__main__':
    main()  
 
