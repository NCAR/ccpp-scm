#!/usr/bin/env python

###############################################################################
# Dependencies
###############################################################################
import argparse
import os
import sys
# Using Prebuild
sys.path.append('ccpp/config')
from ccpp_prebuild_config import SCHEME_FILES, TYPEDEFS_NEW_METADATA, VARIABLE_DEFINITION_FILES
# Using Capgen
sys.path.append('ccpp/framework/scripts')
from metadata_table import register_ddts

###############################################################################
# Argument list
###############################################################################
parser = argparse.ArgumentParser()

###############################################################################
# Main program
###############################################################################
def main():
    # Get command line arguments (N/A)
    args  = parser.parse_args()

    #
    script_dir = 'ccpp/framework/scripts/'
    f2m_script = script_dir+'/ccpp_fortran_to_metadata.py'

    # Create list of known DDTs (using prebuild config script)
    ddts_prebuild = []
    for ct, typedefs in enumerate(TYPEDEFS_NEW_METADATA):
        for typedef in TYPEDEFS_NEW_METADATA[typedefs]:
            ddts_prebuild.append(typedef)
        # end for
    # end for
    #print(ddts_prebuild)

    # Create list of known DDTs (using capgen register_ddts)
    # ######################################################################## 
    # ########################################################################
    # DJS: It's hard to know how Capgen is implemented, but to find the DDTs,
    #      we need names for all the source files (scheme and host). Here I
    #      just used the files lists used by the CCPP prebuild configuration.
    # ########################################################################
    # ######################################################################## 
    # First, convert source file suffixes from fortran to .meta
    # Scheme files
    scheme_metadata = []
    for scheme_file in SCHEME_FILES:
        scheme_metadata.append(os.path.splitext(scheme_file)[0]+'.meta')
    # end for
    # Host files
    host_metadata = []
    # DJS: The last file in VARIABLE_DEFINITION_FILES is not a valid CCPP metadata file. Omit.
    nhost = len(VARIABLE_DEFINITION_FILES)
    for host_file in VARIABLE_DEFINITION_FILES[0:nhost-1]:
        host_metadata.append(os.path.splitext(host_file)[0]+'.meta')
    # end for

    # Finally, register the DDTs
    # NOTE: This returns a list of (lowercase) ddt names.
    ddts_capgen = register_ddts(scheme_metadata+host_metadata)
    #print(ddts_capgen)

    # Create string from DDT list.
    init = True
    known_ddts = ''
    for ddt in ddts_prebuild:
        if (init):
            known_ddts = ddt
            init = False
        else:
            known_ddts = known_ddts+','+ddt
        # end if
    # end for

    #
    # Run Fortran-to-metadata comparision script
    #
    # DJS: This fails if you use the known_ddts from Capgen, which are all lowercase.
    error_count = 0
    error_file  = []
    com = f2m_script+' --ddt-names '+known_ddts+' '
    for scheme_file in SCHEME_FILES:
        result = os.system(com+scheme_file)
        if (result != 0):
            error_count = error_count + 1
            error_file.append(scheme_file)
        # end if
    # end for

    print('------------------------------------------------------------------------------------')
    if (error_count == 0):
        print("SUCCESS! All metadata and Fortran files are consistent")
    else:
        print("ERROR! There are differences between the metadata and Fortran source files")
        print("      ",str(error_count)+' files have issues:')
        for ierror in error_file:
            print("        ",ierror)
        # end for
        print("   See messages above for more details")
    # end if
    print('------------------------------------------------------------------------------------')

###############################################################################
###############################################################################
if __name__ == '__main__':
    main()
