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
git clone -b release/public-v1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-bacio
cd NCEPLIBS-bacio
BACIO_VERSION=`cat VERSION`
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR ..
make VERBOSE=1
make install

cd $NCEPLIBS_SRC
git clone -b release/public-v1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-sp
cd NCEPLIBS-sp
SP_VERSION=`cat VERSION`
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR ..
make VERBOSE=1
make install

cd $NCEPLIBS_SRC
git clone -b release/public-v1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-w3nco
cd NCEPLIBS-w3nco
W3NCO_VERSION=`cat VERSION`
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
echo "export BACIO_LIB4=$NCEPLIBS_DIR/lib/libbacio_v${BACIO_VERSION}_4.a"
echo "export SP_LIBd=$NCEPLIBS_DIR/lib/libsp_v${SP_VERSION}_d.a"
echo "export W3NCO_LIBd=$NCEPLIBS_DIR/lib/libw3nco_v${W3NCO_VERSION}_d.a"
echo "# for csh"
echo "setenv BACIO_LIB4 $NCEPLIBS_DIR/lib/libbacio_v${BACIO_VERSION}_4.a"
echo "setenv SP_LIBd $NCEPLIBS_DIR/lib/libsp_v${SP_VERSION}_d.a"
echo "setenv W3NCO_LIBd $NCEPLIBS_DIR/lib/libw3nco_v${W3NCO_VERSION}_d.a"
echo "============================================================================="
echo " "
echo " "
