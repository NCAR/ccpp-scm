FROM dtcenter/common-community-container:gnu9

MAINTAINER Michelle Harrold <harrold@ucar.edu> or Grant Firl <grantf@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

# Note: The common community container image contains the following packages, which are used to build the SCM:
# gfortran, gcc, cmake, netCDF, HDF5, ZLIB, SZIP, Python, and libxml2
# To access the common community container repository: https://github.com/NCAR/Common-Community-Container

# Obtain CCPP SCM source code
RUN cd /comsoftware \
  && git clone --recursive -b v5.0.0 https://github.com/NCAR/ccpp-scm

# Obtain the pre-computed look-up tables for running with Thompson microphysics
RUN cd /comsoftware/ccpp-scm/ \
  && . contrib/get_thompson_tables.sh
  
# Obtain the IN-CCN data for use with Morrison-Gettelman microphysics
RUN cd /comsoftware/ccpp-scm/ \
  && . contrib/get_mg_inccn_data.sh

# Run the machine setup script to set environment variables
RUN cd /comsoftware/ccpp-scm/scm/ \
  && . etc/CENTOS_docker_setup.sh

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
