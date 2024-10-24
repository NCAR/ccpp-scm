#!/bin/bash

# Function to display help message
print_help() {
    echo "get_all_static_data.sh: contrib/get_all_static_data.sh [-v,--verbose]"
    echo "    Script for downloading/extracting the processed SCM case data."
    echo ""
    echo "Options:"
    echo "    -v, --verbose    Turn on wget verbose output."
    echo "    --help           Show this help message and exit."
}

verbose="-q"
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

#set -ex

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
data_files=("comparison_data" "physics_input_data" "processed_case_input" "raw_case_input")

for file in "${data_files[@]}"; do
    mkdir -p $BASEDIR/scm/data/$file
    cd $BASEDIR/scm/data/$file
    echo "Retrieving $file"
    wget ${verbose} https://github.com/NCAR/ccpp-scm/releases/download/v7.0.0/${file}.tar.gz
    tar -xf ${file}.tar.gz
    rm -f ${file}.tar.gz
done

cd $BASEDIR/
