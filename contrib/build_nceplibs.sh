#!/bin/bash

set -e

if [ "$#" -ne 1 ]; then
  echo "Illegal number of parameters, need to specify install directory for NCEPLIBS."
  exit 1
fi

NCEPLIBS_DIR=$1
if [ -d $NCEPLIBS_DIR ]; then
  while true; do
    read -p "Warning, destination $NCEPLIBS_DIR already exists. Proceed [y/n]? " yn
    case $yn in
      [Yy]* ) break;;
      [Nn]* ) exit;;
      * ) echo "Please answer yes or no.";;
    esac
  done
fi
NCEPLIBS_SRC=$NCEPLIBS_DIR/src
mkdir -p $NCEPLIBS_SRC

cd $NCEPLIBS_SRC
git clone -b v2.4.1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-bacio
cd NCEPLIBS-bacio
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR ..
make VERBOSE=1
make install

cd ../..
git clone -b v2.3.3 --recursive https://github.com/NOAA-EMC/NCEPLIBS-sp
cd NCEPLIBS-sp
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR ..
make VERBOSE=1
make install

cd ../..
git clone -b v2.4.1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-w3nco
cd NCEPLIBS-w3nco
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR ..
make VERBOSE=1
make install

echo " "
echo " "
echo "Set the following environment variables for building the single column model:"
echo "============================================================================="
echo "# for bash"
echo "export bacio_ROOT=$NCEPLIBS_DIR"
echo "export sp_ROOT=$NCEPLIBS_DIR"
echo "export w3nco_ROOT=$NCEPLIBS_DIR"
echo "# for csh"
echo "setenv bacio_ROOT $NCEPLIBS_DIR"
echo "setenv sp_ROOT $NCEPLIBS_DIR"
echo "setenv w3nco_ROOT $NCEPLIBS_DIR"
echo "============================================================================="
echo " "
echo " "
