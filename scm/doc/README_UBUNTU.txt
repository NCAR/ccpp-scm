# Dom Heinzeller, 08/28/2018

In order to build and run SCM-CCPP v2 on Ubuntu Linux, the following installation steps are suggested:

(tested on Ubuntu 18.04.1 LTS with gcc/gfortran 7.3.0)

1. Install default UBUNTU system from ubuntu-18.04.1-desktop-amd64

2. As root (sudo su), install additional packages:
    apt-get update
    apt install synaptic
    apt install gfortran
    apt install libnetcdf-dev
    apt install libnetcdff-dev
    apt install git
    apt install cmake
    apt install libxml2-dev
    # the following two packages are for convenience
    apt install medit
    apt install ncview

3. As standard user, proceed with checking out the code and building as described in the quick start guide
