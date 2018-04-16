SET(sp_source_code ${SPLIB_SRC}/splat.F
                    ${SPLIB_SRC}/lapack_gen.F

)

if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(f_flags -O3 -fdefault-real-8 -fconvert=big-endian -cpp -DLINUX -fPIC)
  # -fp-model strict -real-size 32
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(f_flags -O3 -auto -i4 -r8 -convert big_endian -assume byterecl -fp-model strict -fpp -DLINUX -fPIC)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set(f_flags -O3 -Mrecursive -i4 -r8 -byteswapio -Kieee -Mpreprocess -DLINUX -fPIC)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message ("This program has only been compiled with gfortran, pgf90 and ifort. If another compiler is needed, the appropriate flags must be added in ${SPLIB_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#add OpenMP
set(c_flags ${c_flags} ${OpenMP_C_FLAGS})
set(f_flags ${f_flags} ${OpenMP_Fortran_FLAGS})

if(${CMAKE_VERSION} LESS 3.3)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_SOURCE_FILES_PROPERTIES(${sp_source_code} PROPERTIES COMPILE_FLAGS ${f_flags_str})
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

ADD_LIBRARY(sp STATIC ${sp_source_code})
