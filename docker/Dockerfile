FROM debian:12
MAINTAINER Michael Kavulich <kavulich@ucar.edu>

# Set up base OS environment

RUN apt-get -y update

# Get "essential" tools and libraries
RUN apt-get -y install build-essential \
  && apt-get -y install cmake curl git file gfortran-12 ksh m4 python3 tcsh time wget vim \
  && apt-get -y install libnetcdf-pnetcdf-19 libnetcdff7 libnetcdf-dev libnetcdff-dev libxml2 \
  && apt-get -y install python3-pip python3.11-venv

MAINTAINER Grant Firl <grantf@ucar.edu> or Michael Kavulich <kavulich@ucar.edu>

#Compiler environment variables
ENV CC  /usr/bin/gcc
ENV FC  /usr/bin/gfortran
ENV CXX /usr/bin/g++
ENV F77 /usr/bin/gfortran
ENV F90 /usr/bin/gfortran

# Other necessary environment variables
ENV LD_LIBRARY_PATH /usr/lib/

# Set up unpriviledged user account, set up user home space and make sure user has permissions on all stuff in /comsoftware
RUN groupadd comusers -g 9999 \
 && useradd -u 9999 -g comusers -M -s /bin/bash -c "Unpriviledged user account" -d /home comuser \
 && mkdir /comsoftware \
 && chown -R comuser:comusers /home \
 && chmod 6755 /home \
 && chown -R comuser:comusers /comsoftware \
 && chmod -R 6755 /comsoftware

# Link version-specific aliases (python3 will be created later with virtual environment)
RUN ln -s ~comuser/.venv/bin/python3 /usr/local/bin/python
RUN ln -s /usr/bin/gfortran-12 /usr/bin/gfortran

# all root steps completed above, now continue below as regular userID comuser
USER comuser
WORKDIR /home

# Build NCEP libraries we need for SCM

ENV NCEPLIBS_DIR /comsoftware/nceplibs

RUN mkdir -p $NCEPLIBS_DIR/src && cd $NCEPLIBS_DIR/src \
  && git clone -b v2.4.1 --recursive https://github.com/NOAA-EMC/NCEPLIBS-bacio \
  && mkdir NCEPLIBS-bacio/build && cd NCEPLIBS-bacio/build \
  && cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR .. \
  && make VERBOSE=1 \
  && make install 

RUN cd $NCEPLIBS_DIR/src \
  && git clone -b v2.3.3 --recursive https://github.com/NOAA-EMC/NCEPLIBS-sp \
  && mkdir NCEPLIBS-sp/build && cd NCEPLIBS-sp/build \
  && cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR .. \
  && make VERBOSE=1 \
  && make install

RUN cd $NCEPLIBS_DIR/src \
  && git clone -b v2.11.0 --recursive https://github.com/NOAA-EMC/NCEPLIBS-w3emc \
  && mkdir NCEPLIBS-w3emc/build && cd NCEPLIBS-w3emc/build \
  && cmake -DCMAKE_INSTALL_PREFIX=$NCEPLIBS_DIR .. \
  && make VERBOSE=1 \
  && make install

ENV bacio_ROOT /comsoftware/nceplibs
ENV sp_ROOT /comsoftware/nceplibs
ENV w3emc_ROOT /comsoftware/nceplibs

# Obtain CCPP SCM source code and static data, build code
RUN cd /comsoftware \
  && git clone --recursive -b main https://github.com/NCAR/ccpp-scm \
  && cd /comsoftware/ccpp-scm/ \
  && ./contrib/get_all_static_data.sh \
  && ./contrib/get_thompson_tables.sh \
  && cd /comsoftware/ccpp-scm/scm \
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

# For interactive use, vim mouse settings are infuriating
RUN echo "set mouse=" > ~/.vimrc

# Set up python virtual environment and install needed packages
ENV VIRTUAL_ENV=~/.venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip3 install f90nml==1.4.4 netcdf4==1.6.5
