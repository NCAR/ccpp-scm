#!/bin/bash

# Function to display help message
print_help() {
    echo "get_aerosol_climo.sh: contrib/get_aerosol_climo.sh [-v,--verbose]"
    echo "    Script for downloading/extracting the GOCART climatological aerosol data."
    echo ""
    echo "Options:"
    echo "    -v, --verbose    Turn on wget verbose output."
    echo "    --help           Show this help message and exit."
}

verbose="-v"
wget_options=(
    --tries=10
    --waitretry=10
    --retry-connrefused
    --retry-on-http-error=429,500,502,503,504
)
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
data_files=("FV3_aeroclim1" "FV3_aeroclim2" "FV3_aeroclim3" "FV3_aeroclim_optics")

mkdir -p $BASEDIR/scm/data/physics_input_data/
cd $BASEDIR/scm/data/physics_input_data/
for file in "${data_files[@]}"; do
    echo "Retrieving $file.tar.gz"
    wget ${verbose} "${wget_options[@]}" https://github.com/NCAR/ccpp-scm/releases/download/v7.0.0/${file}.tar.gz
    tar -zxf ${file}.tar.gz
    rm -f ${file}.tar.gz
done

cd $BASEDIR/
