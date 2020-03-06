FROM dtcenter/common-community-container:latest
MAINTAINER Michelle Harrold <harrold@ucar.edu> or Grant Firl <grantf@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

# Note: The common community container image contains the following packages, which are used to build the SCM:
# gfortran, gcc, cmake, netCDF, HDF5, ZLIB, SZIP, Python, and libxml2
# To access the common community container repository: https://github.com/NCAR/Common-Community-Container

# Obtain CCPP SCM source code
RUN cd /comsoftware \
  && git clone --recursive -b release/public-v4 https://github.com/grantfirl/gmtb-scm

# Obtain the pre-computed look-up tables for running with Thompson microphysics
RUN cd /comsoftware/gmtb-scm/scm/data/physics_input_data \
  && curl -L -O https://dtcenter.org/GMTB/freezeH2O.dat \
  && curl -L -O https://dtcenter.org/GMTB/qr_acr_qg.dat \
  && curl -L -O https://dtcenter.org/GMTB/qr_acr_qs.dat

# Upgrade cmake for NCEPLIBS
RUN cd /comsoftware/gmtb-scm/ \
  && mkdir -p cmake \
  && source /opt/rh/devtoolset-8/enable \
  && curl -L -O https://github.com/Kitware/CMake/releases/download/v3.16.5/cmake-3.16.5.tar.gz \
  && tar -xf cmake-3.16.5.tar.gz \
  && cd cmake-3.16.5 \
  && cmake -DCMAKE_USE_OPENSSL=OFF -DCMAKE_INSTALL_PREFIX=/comsoftware/gmtb-scm/cmake  . \
  && make install

ENV PATH /comsoftware/gmtb-scm/cmake/bin:$PATH

# Run the machine setup script to set environment variables
RUN cd /comsoftware/gmtb-scm/scm/ \
  && . etc/CENTOS_docker_setup.sh

ENV BACIO_LIB4 /comsoftware/gmtb-scm/nceplibs/lib/libbacio_v2.2.0_4.a
ENV SP_LIBd /comsoftware/gmtb-scm/nceplibs/lib/libsp_v2.1.0_d.a
ENV W3NCO_LIBd /comsoftware/gmtb-scm/nceplibs/lib/libw3nco_v2.1.0_d.a

# Invoke cmake on the source code to build
RUN cd /comsoftware/gmtb-scm/scm \
  && mkdir bin \
  && cd bin \
  && cmake ../src \
  && make

# The analysis scripts have options for using LaTeX when making figure labels.
# If you would like to install LaTeK, uncomment the section below. 
# Note: This will increase the image size by 1 GB.
#USER root
#RUN yum -y update \
#  && yum -y install texlive-*
#USER comuser

# Set working directory
WORKDIR /comsoftware/gmtb-scm/scm/bin
