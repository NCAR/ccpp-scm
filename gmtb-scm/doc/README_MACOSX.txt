# Dom Heinzeller, 03/23/2018

In order to build and run SCM-CCPP v1 on Mac OS X, the following installation steps are suggested:

1. install homebrew (enter sudo password when requested)
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

2. Install gcc-7.2.0, gfortran-7.2.0
    brew install -v gcc --verbose --without-multilib

3. Install clang-5.0.0 with openmp support
    brew install -v llvm

4. Install mpich-3.2.1 (note: MPI itself not required to run SCM-CCPP, but installing mpich fixes issues with detecting OpenMP for the clang compilers)
    brew install -v mpich

5. Install netCDF library
    brew install -v netcdf

6. Install realpath for MacOSX
    mkdir -p /usr/local/src
    cd /usr/local/src
    git clone https://github.com/harto/realpath-osx.git
    cd realpath-osx
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd /usr/local/src
    rm -fr realpath-osx

7. Install ncview for viewing netCDF files
    brew install -v ncview
