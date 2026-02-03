#!/bin/bash
# Configure Intel compilers to use appropriate GCC version
# This script populates the .cfg files for Intel compilers to ensure they use
# a compatible GCC version based on the Intel compiler version requested.

set -e

ONEAPI_VERSION=${1:?}

# Determine appropriate GCC version based on Intel compiler version
# Intel 2024.x series (oneAPI only) supports GCC 13.x+
determine_gcc_version() {
  local intel_version=$1
  local major=$(echo $intel_version | cut -d. -f1)
  local minor=$(echo $intel_version | cut -d. -f2)
  
  if [ "$major" = "2023" ]; then
    echo "12"
  elif [ "$major" -ge "2024" ]; then
    echo "13"
  else
    # Default to GCC 13 for unknown versions
    echo "13"
  fi
}

GCC_VERSION_ONEAPI=$(determine_gcc_version $ONEAPI_VERSION)

# Check if the GCC version is available, install if needed
echo "Installing GCC ..."
sudo apt-get update
sudo apt-get install -y gcc-${GCC_VERSION_ONEAPI} g++-${GCC_VERSION_ONEAPI} gfortran-${GCC_VERSION_ONEAPI}


# Create/update icx.cfg
icx_cfg_path=$(which icx).cfg
echo "  Modifying icx.cfg to use GCC ${GCC_VERSION_ONEAPI}:"
echo "--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/${GCC_VERSION_ONEAPI}" | sudo tee $icx_cfg_path

# Create/update icpx.cfg
icpx_cfg_path=$(which icpx).cfg
echo "  Modifying icpx.cfg to use GCC ${GCC_VERSION_ONEAPI}:"
echo "--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/${GCC_VERSION_ONEAPI}" | sudo tee $icpx_cfg_path

# Create/update ifx.cfg
ifx_config_path=$(which ifx).cfg
echo "  Modifying ifx.cfg to use GCC ${GCC_VERSION_ONEAPI}:"
echo "-gcc-name=gcc-${GCC_VERSION_ONEAPI}" | sudo tee $ifx_config_path

