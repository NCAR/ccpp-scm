#!/bin/bash

echo "Setting environment variables for SCM-CCPP on MACOSX with clang/gfortran"

export CFLAGS="-fmax-stack-var-size=30000 -I /Volumes/Cluster/opt/local/include/"
export FFLAGS="-fmax-stack-var-size=30000 -fbounds-check -fcheck=all -gdwarf-4 -fvar-tracking-assignments -fbacktrace -finit-real=nan -fcheck=bounds -ffree-line-length-none -cpp -fdefault-real-8 -fdec"