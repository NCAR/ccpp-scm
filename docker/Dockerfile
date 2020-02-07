FROM base:latest
#FROM dtcenter/base_image:latest
MAINTAINER Michelle Harrold <harrold@ucar.edu>

#
# Dockerfile for building CCPP SCM container
#

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

RUN cd /comsoftware/gmtb-scm/scm/src \
#  && cat run_gmtb_scm.py | sed "s/output_'/\/home\/output_'/g" > new_run_gmtb_scm.py \
#  && cp new_run_gmtb_scm.py run_gmtb_scm.py
#  && cat run_gmtb_scm.py | sed "s/if os.path.isdir/#if os.path.isdir/g; s/shutil.rmtree/#shutil.rmtree/g; s/os.makedirs/os.mkdir/g" > new_run_gmtb_scm.py \
#  && cp new_run_gmtb_scm.py run_gmtb_scm.py
#  && cat run_gmtb_scm.py | sed "s/shutil.rmtree(output_dir)/files = glob.glob(output_dir+\"\\/*\")\n            for f in files:\n                 os.remove(f)\n        if not os.path.exists(output_dir):/g; s/os.makedirs/    os.makedirs/g" > new_run_gmtb_scm.py \
#  && cp new_run_gmtb_scm.py run_gmtb_scm.py
  && cat run_gmtb_scm.py | sed "s/os.makedirs(output_dir)/os.makedirs(output_dir)\n        os.symlink(os.path.join('\/comsoftware\/gmtb-scm\/scm\/bin\/' \+ output_dir),os.path.join('\/home\/' \+ output_dir))/g" > new_run_gmtb_scm.py \
  && cp new_run_gmtb_scm.py run_gmtb_scm.py

# Invoke cmake on the source code to build
RUN cd /comsoftware/gmtb-scm/scm \
  && mkdir bin \
  && cd bin \
  && cmake ../src \
  && make > make.log 2>&1

# Set working directory
WORKDIR /comsoftware/gmtb-scm/scm/bin
