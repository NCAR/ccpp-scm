#!/usr/bin/env python

###############################################################################
# Dependencies
###############################################################################
import argparse
import os
import sys
sys.path.append('ccpp/config')
from ccpp_prebuild_config import SCHEME_FILES, TYPEDEFS_NEW_METADATA, VARIABLE_DEFINITION_FILES

###############################################################################
# Main program
###############################################################################
def main():
    # This list of files is not compared since they are known to fail.
    files_known_to_not_pass_test = [
        'ccpp/physics/physics/MP/GFDL/fv_sat_adj.F90',
        'ccpp/physics/physics/Radiation/RRTMG/rrtmg_lw_post.F90',
        'ccpp/physics/physics/SFC_Layer/UFS/sfc_diff.f',
        'ccpp/physics/physics/SFC_Models/Lake/CLM/clm_lake.f90',
        'ccpp/physics/physics/smoke_dust/rrfs_smoke_wrapper.F90']

    # Comparison script from framework
    script_dir = 'ccpp/framework/scripts/'
    f2m_script = script_dir+'/ccpp_fortran_to_metadata.py'

    # Create list of known DDTs (using prebuild config script)
    ddts = []
    for ct, typedefs in enumerate(TYPEDEFS_NEW_METADATA):
        for typedef in TYPEDEFS_NEW_METADATA[typedefs]:
            ddts.append(typedef)
        # end for
    # end for

    # Create string from DDT list.
    init = True
    known_ddts = ''
    for ddt in ddts:
        if (init):
            known_ddts = ddt
            init = False
        else:
            known_ddts = known_ddts+','+ddt
        # end if
    # end for

    # Run Fortran-to-metadata comparision script
    error_count = 0
    error_file  = []
    com = f2m_script+' --ddt-names '+known_ddts+' '
    for scheme_file in SCHEME_FILES:
        if (scheme_file not in files_known_to_not_pass_test):
            result = os.system(com+scheme_file)
            if (result != 0):
                error_count = error_count + 1
                error_file.append(scheme_file)
            # end if
        else:
            print("skipping ",scheme_file)
        # end if
    # end for

    # Display results
    print('------------------------------------------------------------------------------------')
    if (error_count == 0):
        print("SUCCESS! All metadata and Fortran files are consistent")
        sys.exit(0)
    else:
        print("ERROR! There are differences between the metadata and Fortran source files")
        print("      ",str(error_count)+' files have issues:')
        for ierror in error_file:
            print("        ",ierror)
        # end for
        print("   See messages above for more details")
        sys.exit(1)
    # end if
    print('------------------------------------------------------------------------------------')

###############################################################################
###############################################################################
if __name__ == '__main__':
    main()
