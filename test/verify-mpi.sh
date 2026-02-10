#!/bin/bash
# Verify Intel MPI is installed and configured correctly
# Usage: verify-mpi.sh [oneapi|classic]

set +e

COMPILER_TYPE="${1:-oneapi}"

echo "=== Checking MPI Installation (${COMPILER_TYPE}) ==="

if [ "$COMPILER_TYPE" = "oneapi" ]; then
    for compiler in mpiicx mpiicpx mpiifx; do
        which $compiler || { echo "ERROR: $compiler not found"; }
    done
    mpiicx --version
    echo "CC: $CC , should be mpiicx"
    echo "CXX: $CXX , should be mpiicpx"
    echo "FC: $FC , should be mpiifx"
    test "$CC" = "mpiicx" || { echo "ERROR: CC should be mpiicx but is $CC"; exit 1; }
    test "$CXX" = "mpiicpx" || { echo "ERROR: CXX should be mpiicpx but is $CXX"; exit 1; }
    test "$FC" = "mpiifx" || { echo "ERROR: FC should be mpiifx but is $FC"; exit 1; }
elif ["$COMPILER_TYPE" = "gfortran"]; then
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
else
    for compiler in mpiicc mpiicx mpiifort; do
        which $compiler || { echo "ERROR: $compiler not found"; }
    done
    mpiicc --version
    echo "CC: $CC , should be mpiicc"
    echo "CXX: $CXX , should be mpiicpc"
    echo "FC: $FC , should be mpiifort"
    test "$CC" = "mpiicc" || { echo "ERROR: CC should be mpiicc but is $CC"; exit 1; }
    test "$CXX" = "mpiicpc" || { echo "ERROR: CXX should be mpiicpc but is $CXX"; exit 1; }
    test "$FC" = "mpiifort" || { echo "ERROR: FC should be mpiifort but is $FC"; exit 1; }
fi

which mpirun
mpirun --version

echo "✓ MPI installation verified successfully"
