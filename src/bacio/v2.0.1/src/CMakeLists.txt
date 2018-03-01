SET(bacio_source_code ${BACIOLIB_SRC}/baciof.f
                      ${BACIOLIB_SRC}/bacio.c
)

#set Fortran compiler flags
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(f_flags -O3 -fbacktrace -fPIC)
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(f_flags -O3 -xHOST -traceback -free)
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${CMAKE_Fortran_COMPILER_ID})
  message ("This program has only been compiled with gfortran and ifort. If another compiler is needed, the appropriate flags must be added in ${BACIOLIB_SRC}/CMakeLists.txt")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#set C compiler flags
if (${CMAKE_C_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  set(c_flags -O3 -DUNDERSCORE -DMACOSX -fPIC)
elseif (${CMAKE_C_COMPILER_ID} MATCHES "Intel")
  set(c_flags -O3 -DUNDERSCORE -DLINUX)
else (${CMAKE_C_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_C_COMPILER})
  message ("Fortran compiler: " ${CMAKE_C_COMPILER_ID})
  message ("This program has only been compiled with GNU and Intel. If another compiler is needed, the appropriate flags must be added in ${BACIOLIB_SRC}/CMakeLists.txt")
endif (${CMAKE_C_COMPILER_ID} MATCHES "GNU" OR ${CMAKE_C_COMPILER_ID} MATCHES "Clang")

if(${CMAKE_VERSION} LESS 3.3)
  file(GLOB f_files *.f)
  string (REPLACE ";" " " f_flags_str "${f_flags}")
  SET_SOURCE_FILES_PROPERTIES(${f_files} PROPERTIES COMPILE_FLAGS ${f_flags_str})
  file(GLOB c_files *.c)
  string (REPLACE ";" " " c_flags_str "${c_flags}")
  SET_SOURCE_FILES_PROPERTIES(${c_files} PROPERTIES COMPILE_FLAGS ${c_flags_str})
else(${CMAKE_VERSION} LESS 3.3)
  add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${f_flags}>")
  add_compile_options("$<$<COMPILE_LANGUAGE:C>:${c_flags}>")
endif (${CMAKE_VERSION} LESS 3.3)

ADD_LIBRARY(bacio STATIC ${bacio_source_code})
