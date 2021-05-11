#!/bin/bash
#=======================================================================
# Description:  This script is called from the CCPP SCM regression test
#               script rt.sh as each run completes and parses the contents.
#
# Assumptions:
#
# Input arguments: input_file  File to be parsed
#                  output_file File to store output
#
# Usage: ./summarize.sh $input_file $output_file
# set -eux    # Uncomment for debugging
#=======================================================================
function summary_usage() {
  echo
  echo "Usage: $0 input_file output_file | -h"
  echo
  echo "       input_file     [required]"
  echo "       output_file    [required]"
  echo "       -h            display this help"
  echo
  exit 1
}

[[ $# -eq 0 ]] && summary_usage
if [ "$1" = "-h" ] ; then summary_usage ; fi

file=$1
TEST_OUTPUT=$2

run_dir=$(dirname "${file}")
echo "Run directory is ${run_dir}" >> ${TEST_OUTPUT}
#-----------------------------------------------------------------------
# Get total number of tests run
#-----------------------------------------------------------------------
num_tests=$(grep "Executing process" ${file} | head -n 1 | cut -d" " -f6)
echo "Total number of tests is ${num_tests}" >> ${TEST_OUTPUT}

#-----------------------------------------------------------------------
# Count processes executed, completed successfully, and exited with an error code
#-----------------------------------------------------------------------
echo "Process summary:" >> ${TEST_OUTPUT}
echo "================" >> ${TEST_OUTPUT}

num_exec_processes=$(grep -c "Executing process" ${file}) || true
num_completed_processes=$(grep -c "completed successfully" ${file}) || true
num_exited_processes=$(grep -c "exited with code" ${file}) || true

#-----------------------------------------------------------------------
# Determine PASS/FAIL criteria
#-----------------------------------------------------------------------
if [ "${num_exited_processes}" == "0" ] ; then
  echo "PASS:  Number of processes exited = ${num_exited_processes}" >> ${TEST_OUTPUT}
else
  echo "FAIL:  Number of processes exited = ${num_exited_processes}" >> ${TEST_OUTPUT}
fi

if [ "${num_exec_processes}" == "${num_completed_processes}" ] ; then
  echo "PASS: All processes completed successfully ${num_completed_processes}" >> ${TEST_OUTPUT}
else
  echo "FAIL: Number processes executed ${num_exec_processes} /= \
    processes completed successfully ${num_completed_processes}" >> ${TEST_OUTPUT}
fi

#-----------------------------------------------------------------------
# Count output_*/output.nc files, should equal num_tests
#-----------------------------------------------------------------------
num_output_nc_files=$(ls -l ${run_dir}/output_*/output.nc | wc -l)
echo "Number of processes with output.nc files = ${num_output_nc_files}" >> ${TEST_OUTPUT}

if [ "${num_output_nc_files}" == "${num_tests}" ] ; then
  echo "PASS: Number output.nc files ${num_output_nc_files} = ${num_tests}" >> ${TEST_OUTPUT}
else
  echo "FAIL: Number output.nc files ${num_output_nc_files} /= ${num_tests}" >> ${TEST_OUTPUT}
fi

echo " " >> ${TEST_OUTPUT}
echo "Process details:" >> ${TEST_OUTPUT}
echo "================" >> ${TEST_OUTPUT}
if grep -q "Executing process" ${file}; then grep "Executing process" ${file} >> ${TEST_OUTPUT}; fi
if grep -q "completed successfully" ${file}; then grep "completed successfully" ${file} >> ${TEST_OUTPUT}; fi
if grep -q "exited with code" ${file}; then grep "exited with code" ${file} >> ${TEST_OUTPUT}; fi
echo "Summary complete" >> ${TEST_OUTPUT}
