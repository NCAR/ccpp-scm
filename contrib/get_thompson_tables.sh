#!/bin/bash

set -ex

# Directory where this script is located
if [[ $(uname -s) == Darwin ]]; then
  MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
BASEDIR=$MYDIR/..

# Change to directory containing the physics input data, download and extract archive
cd $BASEDIR/scm/data/physics_input_data/
wget https://github.com/NCAR/gmtb-scm/releases/download/v4.0.0/thompson_tables.tar
tar -xvf thompson_tables.tar
rm -f thompson_tables.tar
