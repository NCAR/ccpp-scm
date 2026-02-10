#!/bin/bash
# Verify MPI is installed and configured correctly
# Usage: verify-mpi.sh [oneapi|classic]

set +e

COMPILER_TYPE="${1:-oneapi}"

echo "=== Checking MPI Installation (${COMPILER_TYPE}) ==="

if [ "$COMPILER_TYPE" = "ifx" ]; then
    for compiler in mpiicx mpiicpx mpiifx; do
        which $compiler || { echo "ERROR: $compiler not found"; }
    done
    mpiicx --version
    echo "CC: $CC , should be mpiicx"
    echo "CXX: $CXX , should be mpiicpx"
    echo "FC: $FC , should be mpiifx"
    test "$CC" = "mpiicx"   || { echo "ERROR: CC should be mpiicx but is $CC"; exit 1; }
    test "$CXX" = "mpiicpx" || { echo "ERROR: CXX should be mpiicpx but is $CXX"; exit 1; }
    test "$FC" = "mpiifx"   || { echo "ERROR: FC should be mpiifx but is $FC"; exit 1; }
fi
if [ "$COMPILER_TYPE" = "gfortran" ]; then
    for compiler in mpicc mpicxx mpif90; do
        which $compiler || { echo "ERROR: $compiler not found"; }
    done
    mpicc --version
    echo "CC: $CC , should be mpicc"
    echo "CXX: $CXX , should be mpicxx"
    echo "FC: $FC , should be mpif90"
    test "$CC" = "mpicc"   || { echo "ERROR: CC should be mpicc but is $CC"; exit 1; }
    test "$CXX" = "mpicxx" || { echo "ERROR: CXX should be mpicxx but is $CXX"; exit 1; }
    test "$FC" = "mpif90"  || { echo "ERROR: FC should be mpif90 but is $FC"; exit 1; }    
fi
if [ "$COMPILER_TYPE" = "nvfortran" ]; then
    for compiler in mpicc mpic++ mpifort; do
        which $compiler || { echo "ERROR: $compiler not found"; }
    done
    mpicc --version
    echo "CC: $CC , should be mpicc"
    echo "CXX: $CXX , should be mpic++"
    echo "FC: $FC , should be mpifort"
    test "$CC" = "mpicc"   || { echo "ERROR: CC should be mpicc but is $CC"; exit 1; }
    test "$CXX" = "mpic++" || { echo "ERROR: CXX should be mpic++ but is $CXX"; exit 1; }
    test "$FC" = "mpifort" || { echo "ERROR: FC should be mpifort but is $FC"; exit 1; }
fi

which mpirun
mpirun --version

echo "✓ MPI installation verified successfully"
