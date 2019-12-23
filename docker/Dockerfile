FROM base:latest
MAINTAINER Michelle Harrold <harrold@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

# Obtain NCEPlibs source code
# Change location to cd /comsoftware/libs
RUN source /opt/rh/devtoolset-8/enable \
 && cd /comsoftware \
 && git clone https://github.com/NCAR/NCEPlibs \
 && cd NCEPlibs \
 && ./make_ncep_libs.sh -s linux -c gnu -d /comsoftware/NCEPlibs -o 1 -m 1 <<< $'yes\r1\r'

# Set environment for NCEPLibs
ENV NCEPLIBS_DIR /comsoftware/NCEPlibs

# Obtain CCPP SCM source code
RUN cd /comsoftware \
  && git clone --recursive -b v3.0 https://github.com/NCAR/gmtb-scm

# Run the CCPP prebuild script
RUN cd /comsoftware/gmtb-scm \
  && ./ccpp/framework/scripts/ccpp_prebuild.py --config=./ccpp/config/ccpp_prebuild_config.py

# Run the machine setup script to set environment variables
RUN cd /comsoftware/gmtb-scm/scm/etc \
  && cat CENTOS_setup.sh | sed 's/CC=\/usr\/local\/bin\/gcc/CC=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gcc/g; s/CXX=\/usr\/local\/bin\/g++/CXX=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/g++/g; s/F77=\/usr\/local\/bin\/gfortran/F77=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/F90=\/usr\/local\/bin\/gfortran/F90=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/FC=\/usr\/local\/bin\/gfortran/FC=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/NETCDF=\/usr\/local/NETCDF=\/comsoftware\/libs\/netcdf/g' > new_CENTOS.sh \
  && mv new_CENTOS.sh CENTOS_setup.sh \
  && . CENTOS_setup.sh

# Invoke cmake on the source code to build
RUN cd /comsoftware/gmtb-scm/scm \
  && mkdir bin \
  && cd bin \
  && cmake ../src \
  && make > make.log 2>&1

# Set working directory
WORKDIR /comsoftware/gmtb-scm/scm/bin
