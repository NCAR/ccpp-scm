FROM dtcenter/common-community-container:gnu9

MAINTAINER Michelle Harrold <harrold@ucar.edu> or Grant Firl <grantf@ucar.edu> or Michael Kavulich <kavulich@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

# Note: The common community container image contains the following packages, which are used to build the SCM:
# gfortran, gcc, cmake, netCDF, HDF5, ZLIB, SZIP, Python, and libxml2
# To access the common community container repository: https://github.com/NCAR/Common-Community-Container

# Obtain CCPP SCM source code
RUN cd /comsoftware \
  && git clone --recursive -b release/public-v6 https://github.com/NCAR/ccpp-scm

# Obtain static data that was previously stored in repository
RUN cd /comsoftware/ccpp-scm/ \
  && . contrib/get_all_static_data.sh

# Obtain the pre-computed look-up tables for running with Thompson microphysics
RUN cd /comsoftware/ccpp-scm/ \
  && . contrib/get_thompson_tables.sh
  
# Run the machine setup script to set environment variables
ENV CC=/opt/rh/devtoolset-9/root/usr/bin/gcc
ENV CXX=/opt/rh/devtoolset-9/root/usr/bin/g++
ENV F77=/opt/rh/devtoolset-9/root/usr/bin/gfortran
ENV F90=/opt/rh/devtoolset-9/root/usr/bin/gfortran
ENV FC=/opt/rh/devtoolset-9/root/usr/bin/gfortran

ENV NETCDF=/comsoftware/libs/netcdf

RUN cd /comsoftware/ccpp-scm/contrib \
  && wget https://raw.githubusercontent.com/NCAR/ccpp-scm/3f501aa8af0fb00ff124d8301c932292d1d0abf3/contrib/build_nceplibs.sh \
  && chmod +x build_nceplibs.sh \
  && cd .. \
  && ./contrib/build_nceplibs.sh $PWD/nceplibs

ENV bacio_ROOT /comsoftware/ccpp-scm/nceplibs
ENV sp_ROOT /comsoftware/ccpp-scm/nceplibs
ENV w3nco_ROOT /comsoftware/ccpp-scm/nceplibs

# Create your own link from python -> python3 
# This works without setting the system PATH env var
# since /usr/local/bin is before /usr/bin in the search path.
USER root
RUN ln -s /usr/bin/python3 /usr/local/bin/python
USER comuser

# Invoke cmake on the source code to build
RUN cd /comsoftware/ccpp-scm/scm \
  && mkdir bin \
  && cd bin \
  && cmake ../src \
  && make -j4

# The analysis scripts have options for using LaTeX when making figure labels.
# If you would like to install LaTeK, uncomment the section below. 
# Note: This will increase the image size by 1 GB.
#USER root
#RUN yum -y update \
#  && yum -y install texlive-*
#USER comuser

# Set working directory
WORKDIR /comsoftware/ccpp-scm/scm/bin
ENV SCM_WORK=/comsoftware
ENV SCM_ROOT=/comsoftware/ccpp-scm/
