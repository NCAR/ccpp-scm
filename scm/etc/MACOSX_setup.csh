#!/bin/tcsh

echo "Setting environment variables for SCM-CCPP on MACOSX with clang/gfortran"

setenv CC /usr/local/bin/mpicc
setenv CXX /usr/local/bin/mpicxx
setenv F77 /usr/local/bin/mpif77
setenv F90 /usr/local/bin/mpif90
setenv FC /usr/local/bin/mpif90
setenv CPP "/usr/local/bin/mpif90 -E -x f95-cpp-input"

setenv LDFLAGS "-L/usr/local/opt/zlib/lib -L/usr/local/opt/llvm/lib"
setenv CPPFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv CFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv CXXFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv FFLAGS "-I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"
setenv FCFLAGS " -I/usr/local/opt/zlib/include -I/usr/local/opt/llvm/include"

if (! $?PATH) then
  setenv PATH "/usr/local/opt/llvm/bin"
else
  setenv PATH "/usr/local/opt/llvm/bin:$PATH"
endif
if (! $?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib"
else
  setenv LD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib:$LD_LIBRARY_PATH"
endif
if (! $?DYLD_LIBRARY_PATH) then
  setenv DYLD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib"
else
  setenv DYLD_LIBRARY_PATH "/usr/local/opt/zlib/lib:/usr/local/opt/llvm/lib:$DYLD_LIBRARY_PATH"
endif

setenv NETCDF /usr/local
