#!/bin/bash

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
data_files=("FV3_aeroclim1" "FV3_aeroclim2" "FV3_aeroclim3" "FV3_aeroclim_optics")

cd $BASEDIR/scm/data/physics_input_data/
for file in "${data_files[@]}"; do
    echo "Retrieving $file.tar.gz"
    wget https://github.com/NCAR/ccpp-scm/releases/download/v7.0.0-beta/${file}.tar.gz
    tar -xvf ${file}.tar.gz
    rm -f ${file}.tar.gz
done

cd $BASEDIR/

