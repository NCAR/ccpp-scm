#!/bin/bash

#set -ex

# Directory where this script is located
if [[ $(uname -s) == Darwin ]]; then
  MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
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
#    wget https://github.com/NCAR/ccpp-scm/releases/download/v5.1.0/${file}.tar.gz
#    tar -xf ${file}.tar.gz
#    rm -f ${file}.tar.gz
done

cd $MYDIR

#Legacy static data have their own scripts
echo "retrieving mg_inccn_data"
source get_mg_inccn_data.sh
cd $MYDIR
echo "retrieving thompson_tables"
source get_thompson_tables.sh

cd $BASEDIR/

