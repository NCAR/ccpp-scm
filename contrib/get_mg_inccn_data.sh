#!/bin/bash

# Function to display help message
print_help() {
    echo "get_mg_inccn_data.sh: contrib/get_mg_inccn_data.sh [-v,--verbose]"
    echo "    Script for downloading/extracting the Morrison-Gettelman data."
    echo ""
    echo "Options:"
    echo "    -v, --verbose    Turn on wget verbose output."
    echo "    --help           Show this help message and exit."
}

verbose="-nv"
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --help)
            print_help
            exit 0
            ;;
        -v|--verbose)
            verbose="-v"
            ;;
        *)
            echo "Unknown option: $1"
            print_help
            exit 1
            ;;
    esac
    shift
done

set -ex

# Directory where this script is located
if [[ $(uname -s) == Darwin ]]; then
  if [[ $(sw_vers -productVersion) < 12.3 ]]; then
    MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
  else
    MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
  fi
else
  MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
BASEDIR=$MYDIR/..

# Change to directory containing the physics input data, download and extract archive
cd $BASEDIR/scm/data/physics_input_data/
wget ${verbose} https://github.com/NCAR/ccpp-scm/releases/download/v7.0.0/MG_INCCN_data.tar.gz
tar -xvf MG_INCCN_data.tar.gz
rm -f MG_INCCN_data.tar.gz
cd $BASEDIR/
