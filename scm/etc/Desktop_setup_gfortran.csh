#!/bin/tcsh

echo "Setting environment variables for CCPP-SCM on Desktop (MacOS) with gcc/gfortran"

echo "Setting CC/CXX/FC environment variables"
setenv CC /opt/local/bin/gcc-mp-10
setenv CXX /opt/local/bin/g++-mp-10
setenv FC gfortran-mp-10

echo "Setting location of NETCDF"
setenv NETCDF /opt/local
setenv LDFLAGS "-I${NETCDF}/include -L${NETCDF}/lib -Wl,-rpath,${NETCDF}/lib"

echo "Setting location of NCEPLIBS libraries"
setenv BACIO_LIB4 /Users/$USER/NCEPLIBS/lib/libbacio_v2.2.0_4.a
setenv SP_LIBd /Users/$USER/NCEPLIBS/lib/libsp_v2.1.0_d.a
setenv W3NCO_LIBd /Users/$USER/NCEPLIBS/lib/libw3nco_v2.1.0_d.a

#check to see if CMake is installed locally
echo "Checking if CMake is installed"
cmake --version

if ( $? != 0 )  then
        echo "CMake not found; installing CMake"
        pip install cmake
else
        echo "CMake is installed"
endif

#check to see if f90nml is installed locally
echo "Checking if f90nml python module is installed"
python -c "import f90nml"

if ( $? != 0 ) then
        echo "Not found; installing f90nml"
        pip install --no-cache-dir f90nml==0.19
else
        echo "f90nml is installed"
endif
