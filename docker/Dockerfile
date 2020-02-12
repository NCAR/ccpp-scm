#FROM dtcenter/common-community-container:latest
FROM base:latest
MAINTAINER Michelle Harrold <harrold@ucar.edu> or Grant Firl <grantf@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

# Note: The common community container image contains the following packages, which are used to build the SCM:
# gfortran, gcc, cmake, netCDF, HDF5, ZLIB, SZIP, Python, and libxml2
# To access the common community container repository: https://github.com/NCAR/Common-Community-Container

# Obtain NCEPlibs source code
RUN source /opt/rh/devtoolset-8/enable \
 && cd /comsoftware/libs \
 && git clone https://github.com/NCAR/NCEPlibs \
 && cd NCEPlibs \
 && ./make_ncep_libs.sh -s linux -c gnu -d /comsoftware/libs/NCEPlibs -o 1 -m 1 <<< $'yes\r1\r'

# Set environment for NCEPLibs
ENV NCEPLIBS_DIR /comsoftware/libs/NCEPlibs

# Obtain CCPP SCM source code
RUN cd /comsoftware \
  && git clone --recursive -b v3.0 https://github.com/NCAR/gmtb-scm \
  && cd gmtb-scm/ccpp/physics \
  && git clone https://github.com/NCAR/ccpp-physics.git \
  && cd ../framework \
  && git clone https://github.com/NCAR/ccpp-framework.git

# Obtain the pre-computed look-up tables for running with Thompson microphysics
RUN cd /comsoftware/gmtb-scm/scm/data/physics_input_data \
  && curl -L -O https://dtcenter.org/GMTB/freezeH2O.dat \
  && curl -L -O https://dtcenter.org/GMTB/qr_acr_qg.dat \
  && curl -L -O https://dtcenter.org/GMTB/qr_acr_qs.dat

# Run the CCPP prebuild script
RUN cd /comsoftware/gmtb-scm \
  && ./ccpp/framework/scripts/ccpp_prebuild.py --config=./ccpp/config/ccpp_prebuild_config.py

# Run the machine setup script to set environment variables
RUN cd /comsoftware/gmtb-scm/scm/etc \
  && cat CENTOS_setup.sh | sed 's/CC=\/usr\/local\/bin\/gcc/CC=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gcc/g; s/CXX=\/usr\/local\/bin\/g++/CXX=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/g++/g; s/F77=\/usr\/local\/bin\/gfortran/F77=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/F90=\/usr\/local\/bin\/gfortran/F90=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/FC=\/usr\/local\/bin\/gfortran/FC=\/opt\/rh\/devtoolset-8\/root\/usr\/bin\/gfortran/g; s/NETCDF=\/usr\/local/NETCDF=\/comsoftware\/libs\/netcdf/g' > new_CENTOS.sh \
  && mv new_CENTOS.sh CENTOS_setup.sh \
  && . CENTOS_setup.sh

# To allow for functional volume mounts from "local" space to Docker space, modifications must be made to run_gmtb_scm.py.
# For the running the CCPP SCM in Docker, output will be copied to the /home directory and mounted from there.
RUN cd /comsoftware/gmtb-scm/scm/src \
  && cat run_gmtb_scm.py | sed "/def launch_executable(use_gdb, gdb):/i \       \ return output_dir\n\ndef copy_outdir(exp_dir):\n    \"\"\"Copy output directory to \/home for this experiment.\"\"\"\n    home_output_dir = \'\/home\/\'+exp_dir\n    if os.path.isdir(home_output_dir):\n        shutil.rmtree(home_output_dir)\n    shutil.copytree(exp_dir, home_output_dir)\n" | sed "/atexit.register/i \   \ atexit.register(copy_outdir, exp_dir)" | sed "s/exp.setup_rundir()/exp_dir = exp.setup_rundir()/g" > run_gmtb_scm_sed.py \
  && mv run_gmtb_scm_sed.py run_gmtb_scm.py \
  && chmod +x run_gmtb_scm.py

# Invoke cmake on the source code to build
RUN cd /comsoftware/gmtb-scm/scm \
  && mkdir bin \
  && cd bin \
  && cmake ../src \
  && make > make.log 2>&1

# The analysis scripts have options for using LaTeX when making figure labels.
# If you would like to install LaTeK, uncomment the section below. 
# Note: This will increase the image size 1 GB.
#USER root
#RUN yum -y update \
#  && yum -y install texlive-*
#USER comuser

# Set working directory
WORKDIR /comsoftware/gmtb-scm/scm/bin
