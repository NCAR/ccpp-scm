#!/usr/bin/env python

# CCPP Capgen config for CCPP Single Column Model (SCM)

# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'scm/bin'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE     = '{build_dir}/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE    = '{build_dir}/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE   = '{build_dir}/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE      = '{build_dir}/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE     = '{build_dir}/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE    = '{build_dir}/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE         = '{build_dir}/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE        = '{build_dir}/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE       = '{build_dir}/physics/CCPP_CAPS.sh'

# Auto-generated makefile/cmakefile snippets that contain all working precisions
KINDS_MAKEFILE        = '{build_dir}/physics/CCPP_KINDS.mk'
KINDS_CMAKEFILE       = '{build_dir}/physics/CCPP_KINDS.cmake'
KINDS_SOURCEFILE      = '{build_dir}/physics/CCPP_KINDS.sh'

# CCPP Static API.
STATIC_API_MAKEFILE   = '{build_dir}/../../src/CCPP_STATIC_API.mk'
STATIC_API_CMAKEFILE  = '{build_dir}/../../src/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = '{build_dir}/../../src/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = 'ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = 'ccpp/physics/CCPP_VARIABLES_SCM.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_SCM.tex'
