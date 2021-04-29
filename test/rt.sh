#!/bin/bash
#=======================================================================
# Description:  This script builds and runs the CCPP SCM
#               The executable gmtb_scm is built for each compiler and
#               build type on the platform (machine name).
#
# Assumptions:
#
# Command line arguments: machine name (hera or cheyenne)
#
# Usage: ./rt.sh $machine >& test.out &      # Run the regression tests
# set -eux    # Uncomment for debugging
#=======================================================================

fail() { echo -e "\n$1\n" >> ${TEST_OUTPUT} && exit 1; }

function usage() {
  echo
  echo "Usage: $0 machine | -h"
  echo
  echo "       machine       [required] is one of: ${machines[@]}"
  echo "       -h            display this help"
  echo
  exit 1
}

function wait_for_criteria() {
 
# Function that iteratively checks the return code of a given Unix command
# and returns 0 when the command is successful or 1 if not

  n_sleep_sec=$1    ; shift    # Number of seconds to sleep before retry
  max_iterations=$1 ; shift    # Max number of times to try
  msg=$1            ; shift    # Message for output
  TEST_OUTPUT=$1    ; shift    # File for output
  cmd=$@                       # Unix command to test for criteria
  echo "Command to test: $cmd" >> ${TEST_OUTPUT}

  local criteria=0
  local counter=0

  while [ ${criteria} -eq 0 ] ; do
    while [ ${counter} -le ${max_iterations} ] ; do
      criteria=$(eval $cmd)
      if [ "${criteria}" -eq 1 ] ; then break ; fi    # criteria met
      echo "criteria = ${criteria}, counter = ${counter}, sleeping..." >> ${TEST_OUTPUT}
      sleep ${n_sleep_sec}
      let counter=counter+1
    done
    break   # criteria not met
  done
  if [ ${counter} -ge ${max_iterations} ] ; then
    echo "ERROR: ${msg} not found"  >> ${TEST_OUTPUT}
    exit 1
  else
    echo "SUCCESS: ${msg} found" >> ${TEST_OUTPUT}
  fi
}

machines=( hera cheyenne )
executable_name=gmtb_scm

[[ $# -eq 0 ]] && usage
if [ "$1" = "-h" ] ; then usage ; fi

export machine=${1}
machine=$(echo "${machine}" | tr '[A-Z]' '[a-z]')  # Make machine name all lower case
machine=$(echo "${machine^}")                      # Capitalize first letter for setup script name
build_it=0                                         # Set to 1 to skip build (for testing run, pass/fail criteria)
run_it=0                                           # Set to 1 to skip run (for testing pass/fail criteria)

build_types=( Release Debug )                      # Set all instances of CMAKE_BUILD_TYPE

if [ "${machine}" == "Cheyenne" ] ; then
  users=( $USER@ucar.edu )
  compilers=( intel gnu )
  job_submission_script=gmtb_scm_qsub_example.py
  memory="512M"
  walltime="walltime=00:40:00"
  max_iterations=240                                 # Max # of iterations to wait for batch output file to
                                                     # be created and to job to complete
  n_sleep_sec=30                                     # Interval in seconds between interation
else
  users=( $USER@noaa.gov )
  compilers=( intel )
  job_submission_script=gmtb_scm_slurm_example.py
  memory="1G"
  walltime="40"
  max_iterations=120
  n_sleep_sec=30
fi

#-----------------------------------------------------------------------
# Set some directories
#-----------------------------------------------------------------------
PID=$$
TEST_DIR=$( pwd )                   # Directory with this script
TEST_OUTPUT=${TEST_DIR}/rt_summary${PID}.out
TOP_DIR=$TEST_DIR/..
SCM_DIR=$TOP_DIR/scm
ETC_DIR=$TOP_DIR/scm/etc
SRC_DIR=$TOP_DIR/scm/src
PHYS_DATA_DIR=$TOP_DIR/scm/data/physics_input_data
#TODO:  Add script that gets cam5_4_143_NAAI_monclimo2.nc and cam5_4_143_NPCCN_monclimo2.nc files
phys_data_files=( CCN_ACTIVATE.BIN freezeH2O.dat qr_acr_qgV2.dat qr_acr_qsV2.dat )
batch_out=test_job.o*           # Without PID extension

#-----------------------------------------------------------------------
# Check if lookup tables and other large datasets have been downloaded
# If any of the phys_data_files are missing, run the download script
#-----------------------------------------------------------------------
for phys_data_file in "${phys_data_files[@]}"; do
  if [ ! -f "${PHYS_DATA_DIR}/${phys_data_file}" ]; then
    ${TOP_DIR}/contrib/get_thompson_tables.sh
    break
  else
    echo "Data file ${phys_data_file} exists..." >> ${TEST_OUTPUT}
  fi
done

#-----------------------------------------------------------------------
# Create the output file if it doesn't exist
#-----------------------------------------------------------------------
if [ ! -f "$TEST_OUTPUT" ]; then
   touch ${TEST_OUTPUT}
fi

#-----------------------------------------------------------------------
# Set up the build environment and run the build script.
#-----------------------------------------------------------------------
for compiler in "${compilers[@]}"; do
  for build_type in "${build_types[@]}"; do

    echo "-------------------------------------------------------" >> ${TEST_OUTPUT}
    echo "Starting ${compiler} ${build_type} build" >> ${TEST_OUTPUT}
    echo "-------------------------------------------------------" >> ${TEST_OUTPUT}

    build_type_lc=$(echo "${build_type}" | tr '[A-Z]' '[a-z]')  # Make build_type all lower case
    BIN_DIR=$TOP_DIR/scm/bin_${compiler}_${build_type_lc}       # Assign a bin and run dir
    RUN_DIR=$TOP_DIR/scm/run_${compiler}_${build_type_lc}       # for each build/run in test
    BUILD_OUTPUT=${BIN_DIR}/build.out
    if [ "${build_type}" == "Debug" ] ; then
      # Add --runtime ${runtime} to multi_run_gmtb_scm.py to reduce runtime for tests
      test_run_cmd="${RUN_DIR}/multi_run_gmtb_scm.py -f ${TEST_DIR}/rt_test_cases.py -v --runtime 86400" # 1 day
    else
      test_run_cmd="${RUN_DIR}/multi_run_gmtb_scm.py -f ${TEST_DIR}/rt_test_cases.py -v"
    fi

    . ${ETC_DIR}/${machine}_setup_${compiler}.sh

    if [ $build_it -eq 0 ] ; then
#-----------------------------------------------------------------------
# Build the SCM
#-----------------------------------------------------------------------
      if [ -d "${BIN_DIR}" ] ; then rm -rf ${BIN_DIR}; fi
      mkdir ${BIN_DIR}
      cd ${BIN_DIR}
      cmake -DCMAKE_BUILD_TYPE=${build_type} ../src >& log.cmake || fail "cmake FAILED, see ${BIN_DIR}/log.cmake"
      make -j4 >& ${BUILD_OUTPUT} || fail "Build ${machine} ${compiler} FAILED, see ${BUILD_OUTPUT}"
    fi    # End of skip build for testing

    exec_file=${BIN_DIR}/${executable_name}
    if [ -f "${exec_file}" ]; then
      echo "SUCCESS: ${machine} ${compiler} ${build_type} executable file ${exec_file} exists" >> ${TEST_OUTPUT}
    else
      echo "FAIL: ${machine} ${compiler} ${build_type} executable file ${exec_file} does NOT exist" >> ${TEST_OUTPUT}
    fi

#-----------------------------------------------------------------------
# Set up the run directory for the build
#-----------------------------------------------------------------------
    if [ $run_it -eq 0 ] ; then
      if [ -d "${RUN_DIR}" ] ; then rm -rf ${RUN_DIR}; fi
      mkdir ${RUN_DIR}
      cd ${RUN_DIR}
      ln -s ${BIN_DIR}/${executable_name} ${executable_name}
      ln -s ${SRC_DIR}/multi_run_gmtb_scm.py multi_run_gmtb_scm.py
      ln -s ${SRC_DIR}/run_gmtb_scm.py run_gmtb_scm.py
      cp ${ETC_DIR}/${job_submission_script} ${job_submission_script}.tmp
#-----------------------------------------------------------------------
# Substitute COMMAND in job_submission script to use multi_run
#-----------------------------------------------------------------------
      sed "s,^.*COMMAND \= .*,COMMAND \= \"${test_run_cmd}\"," ${job_submission_script}.tmp > ${job_submission_script}.tmp2
      sed "s,^.*WALLTIME \= .*,WALLTIME \= \"${walltime}\"," ${job_submission_script}.tmp2 > ${job_submission_script}.tmp3
      sed "s,^.*SERIAL_MEM \= .*,SERIAL_MEM \= \"${memory}\"," ${job_submission_script}.tmp3 > ${job_submission_script}
      chmod +x ${job_submission_script}
      rm -rf ${job_submission_script}.tmp
      rm -rf ${job_submission_script}.tmp2
      rm -rf ${job_submission_script}.tmp3

#-----------------------------------------------------------------------
# Submit job_submission to queue
#-----------------------------------------------------------------------
      ./${job_submission_script}

    fi # End of skip run for testing

  done # End build type loop
done   # End compiler loop

#-----------------------------------------------------------------------
# Wait for job to finish, search for "Done" at the end of the file
#-----------------------------------------------------------------------
n_tests=0
for compiler in "${compilers[@]}"; do
  for build_type in "${build_types[@]}"; do
  
    echo "-------------------------------------------------------" >> ${TEST_OUTPUT}
    echo "Monitoring ${compiler} ${build_type} run..." >> ${TEST_OUTPUT}
    echo "-------------------------------------------------------" >> ${TEST_OUTPUT}

    build_type_lc=$(echo "${build_type}" | tr '[A-Z]' '[a-z]')  # Make build_type all lower case
    RUN_DIR=$TOP_DIR/scm/run_${compiler}_${build_type_lc}      

#-----------------------------------------------------------------------
# Wait until the batch output file exists
#-----------------------------------------------------------------------
    msg="Batch output file for ${compiler} ${build_type}"
    cmd="ls -1 ${RUN_DIR}/${batch_out} | wc -l"

    wait_for_criteria ${n_sleep_sec} ${max_iterations} "${msg}" ${TEST_OUTPUT} $cmd
    ret_val=$?

    if [ ${ret_val} -eq 0 ] ; then
      batch_out_file=`ls -1 ${RUN_DIR}/${batch_out}`
#-----------------------------------------------------------------------
# Wait until the run completes by searching for "Done"
#-----------------------------------------------------------------------
      msg=" String Done for ${compiler} ${build_type}"
      cmd="grep -c "Done" ${batch_out_file}"

      wait_for_criteria ${n_sleep_sec} ${max_iterations} "${msg}" ${TEST_OUTPUT} $cmd
      ret_val=$?
      if [ ${ret_val} -eq 0 ] ; then
        echo "SUCCESS:  ${batch_out_file} completed, creating summary..." >> ${TEST_OUTPUT}
        echo "-------------------------------------------------------" >> ${TEST_OUTPUT}
        echo "SUMMARY for ${compiler} ${build_type}:" >> ${TEST_OUTPUT}
        echo "-------------------------------------------------------" >> ${TEST_OUTPUT}
        . ${TEST_DIR}/summarize.sh ${batch_out_file} ${TEST_OUTPUT}
      else
        echo "ERROR: ${batch_out_file} did not complete in allowed time, continuing to next run..." >> ${TEST_OUTPUT}
      fi
    else
      echo "ERROR: File ${batch_out_file} not found, continuing to next run..." >> ${TEST_OUTPUT}
    fi
    ((n_tests=n_tests+1))
  done #End build type loop
done   #End compiler loop

#-----------------------------------------------------------------------
# Determine PASS/FAIL of all tests
#-----------------------------------------------------------------------
msg="????"
n_fail=$(grep -ci "fail" ${TEST_OUTPUT}) || true
n_error=$(grep -ci "error" ${TEST_OUTPUT}) || true
n_completed=$(grep -c "PASS: All processes completed successfully" ${TEST_OUTPUT}) || true

if [[ $n_completed -eq $n_tests && $n_fail -eq 0 && $n_error -eq 0 ]] ; then
  echo "ALL TESTS SUCCEEDED" >> ${TEST_OUTPUT}
  msg="PASS"
else
  echo "TEST(S) FAILED" >> ${TEST_OUTPUT}
  msg="FAIL"
fi

#-----------------------------------------------------------------------
# Send email with PASS/FAIL message
#-----------------------------------------------------------------------
cat ${TEST_OUTPUT} | mail -s "[ ${msg} ] ccpp-scm test $1" ${users}
