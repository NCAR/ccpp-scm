# Dom Heinzeller, 04/05/2018

In order to build and run SCM-CCPP v1 on CentOS Linux, the following installation steps can be used:

1. Install "Gnome Desktop / Development Tools" CentOS system from CentOS-7-x86_64-Everything-1708.iso with network enabled

2. As root (su), install additional packages:
    yum --enablerepo=extras install epel-release
    yum update
    yum install gmp-devel
    yum install mpfr-devel
    yum install libmpc-devel
    yum install autogen
    yum install libxml2-devel
    yum install libaec-devel
    yum install cmake
    yum install dejagnu
    yum install texinfo
    yum install ncview

3. As root (su), install modern GNU compilers
    curl https://ftp.gnu.org/gnu/gcc/gcc-7.3.0/gcc-7.3.0.tar.gz -O
    tar xvf gcc-7.3.0.tar.gz
    mkdir gcc-7.3.0-build
    cd gcc-7.3.0-build
    ../gcc-7.3.0/configure --disable-multilib 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install

4. As root (su), install netCDF library/headers for same compiler

    export LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH"

    cd /usr/local/src

    # Download the following src files from the web to /usr/local/src
    hdf5-1.8.20.tar.gz
    netcdf-4.5.0.tar.gz
    netcdf-cxx4-4.3.0.tar.gz
    netcdf-fortran-4.4.4.tar.gz
    parallel-netcdf-1.8.1.tar.gz
    szip-2.1.1.tar.gz
    zlib-1.2.11.tar.gz

    # zlib-1.2.11
    tar -xvf zlib-1.2.11.tar.gz 
    cd zlib-1.2.11/
    ./configure 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr zlib-1.2.11

    # szip-2.1.1
    gunzip szip-2.1.1.tar.gz 
    tar -xvf szip-2.1.1.tar 
    gzip szip-2.1.1.tar 
    cd szip-2.1.1/
    ./configure 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr szip-2.1.1

    # hdf5-1.8.20
    tar -xvzf hdf5-1.8.20.tar.gz 
    cd hdf5-1.8.20/
    ./configure --with-szlib=/usr/local --with-zlib=/usr/local --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr hdf5-1.8.20

    # netcdf-4.5.0
    tar -xvf netcdf-4.5.0.tar.gz 
    cd netcdf-4.5.0/
    CFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib -lhdf5 -lhdf5_hl -lsz -lz" ./configure 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr netcdf-4.5.0

    # netcdf-fortran-4.4.4
    tar -xvf netcdf-fortran-4.4.4.tar.gz 
    cd netcdf-fortran-4.4.4/
    FFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib -lnetcdf -lhdf5 -lhdf5_hl -lsz -lz" ./configure 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    sync
    cd ..
    rm -fr netcdf-fortran-4.4.4

5. As standard user, add to .bashrc

    export PATH="/usr/local/bin:$PATH"
    export LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH"
