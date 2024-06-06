#!/bin/bash
#=======================================================================
# Description:  This script builds and runs the CCPP SCM.
#               The executable scm is built for each compiler and
#               build type on the platform (machine name).
#               Command line options are available for baseline
#               generation and comparison.
#
# Assumptions:
#
# Command line arguments: machine name (hera or derecho)
#                         -g <dir> run tests, generate baseline in <dir>
#                         -c <dir> run tests, compare to baseline in <dir>
#                         -h display usage
#
# Usage: ./rt.sh $machine >& test.out &      # Run the regression tests
#        ./rt.sh $machine -g /top_level/path/for/generated/baseline >& test.out &
#        ./rt.sh $machine -c /top_level/path/for/baseline/comparison >& test.out &
#
# set -eux    # Uncomment for debugging
#=======================================================================

fail() { echo -e "\n$1\n" >> ${TEST_OUTPUT} && exit 1; }

function usage() {
  echo
  echo "Usage: $0 machine [-g <dir>] [-c <dir>] | -h"
  echo
  echo "       machine       [required] is one of: ${machines[@]}"
  echo "       -g <dir>      [optional] run tests, write generated baseline to <dir>"
  echo "       -c <dir>      [optional] run tests, compare to baseline in <dir>"
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

machines=( hera derecho desktop )
executable_name=scm

[[ $# -eq 0 ]] && usage                     # Display usage if no command-line arguments

#-----------------------------------------------------------------------
# Parse command-line arguments
#-----------------------------------------------------------------------
gen_baseline=false
cmp_baseline=false
while [ $# -gt 0 ]; do
  while getopts "g:c:h" opt; do            # Parse parameters, exit on non-option parameter
    case "$opt" in
      g) gen_baseline_dir=$OPTARG; gen_baseline=true; shift;;
      c) cmp_baseline_dir=$OPTARG; cmp_baseline=true; shift;;
      h) usage;;
      \?) echo unknown Option;;
      :) echo missing required parameter for Option $opt;;
    esac
    shift
    OPTIND=1
  done
  if [ $# -gt 0 ]; then                    # Save required parameter(s) in array
    REQUIRED_ARGS=(${REQUIRED_ARGS[@]} $1)
    shift
    OPTIND=1
  fi
done

export machine=$REQUIRED_ARGS

#Check that only one required parameter has been entered
if [ ${#REQUIRED_ARGS[@]} -gt 1 ] ; then
  echo "ERROR: Number of required parameters entered is > 1"
  usage
  exit 1
fi

machine=$(echo "${machine}" | tr '[A-Z]' '[a-z]')  # Make machine name all lower case

#-----------------------------------------------------------------------
#Check that machine is valid
#-----------------------------------------------------------------------
if [[ "${machines[@]}" =~ "$machine" ]]; then
  echo "machine ${machine} is valid"
else
  echo "ERROR: machine ${machine} is NOT valid"
  exit 1
fi

#-----------------------------------------------------------------------
#Check that -g and -c are not both set
#-----------------------------------------------------------------------
if [[ "$gen_baseline" == "true" ]] && [[ "$cmp_baseline" == "true" ]]; then
  echo "ERROR: Baseline cannot be generated and compared in same run, set either -g OR -c"
  exit 1
fi
if [[ "$gen_baseline" = "true" ]]; then echo "Generated baseline is in ${gen_baseline_dir}"; fi
if [[ "$cmp_baseline" = "true" ]]; then
  if [ ! -d "${cmp_baseline_dir}" ]; then
    fail "FAIL: Directory for baseline comparison ${cmp_baseline_dir} does not exist"
  fi
  echo "Comparison baseline is in ${cmp_baseline_dir}"
fi

machine="$(tr '[:lower:]' '[:upper:]' <<< ${machine:0:1})${machine:1}"  # Capitalize first letter for setup script name
build_it=0                                         # Set to 1 to skip build (for testing run, pass/fail criteria)
run_it=0                                           # Set to 1 to skip run (for testing pass/fail criteria)

build_types=( Release Debug )                      # Set all instances of CMAKE_BUILD_TYPE

if [ "${machine}" == "Derecho" ] ; then
  users=( $USER@ucar.edu )
  compilers=( intel gnu )
  use_batch_system=true
  job_submission_script=scm_qsub_example.py
  memory="512M"
  walltime="walltime=00:40:00"
  max_iterations=240                                 # Max # of iterations to wait for batch output file to
                                                     # be created and to job to complete
  n_sleep_sec=30                                     # Interval in seconds between interation
elif [ "${machine}" == "Hera" ] ; then
  users=( $USER@noaa.gov )
  compilers=( intel )
  use_batch_system=true
  job_submission_script=scm_slurm_example.py
  memory="1G"
  walltime="40"
  max_iterations=120
  n_sleep_sec=30
elif [ "${machine}" == "Desktop" ] ; then
  users=( $USER@ucar.edu )
  compilers=( gfortran )
  use_batch_system=false
  max_iterations=10
  n_sleep_sec=2
fi

#-----------------------------------------------------------------------
# Set some directories
#-----------------------------------------------------------------------
PID=$$
TEST_DIR=$( pwd )                   # Directory with this script
TEST_OUTPUT=${TEST_DIR}/rt_summary${PID}.out
TEST_LOGFILE=${TEST_DIR}/rt_log${PID}.out
TOP_DIR=$TEST_DIR/..
SCM_DIR=$TOP_DIR/scm
ETC_DIR=$TOP_DIR/scm/etc
SRC_DIR=$TOP_DIR/scm/src
PHYS_DATA_DIR=$TOP_DIR/scm/data/physics_input_data
phys_data_files=( CCN_ACTIVATE.BIN freezeH2O.dat qr_acr_qgV2.dat qr_acr_qsV2.dat )
mg_inccn_data_files=( cam5_4_143_NAAI_monclimo2.nc cam5_4_143_NPCCN_monclimo2.nc )
job_prefix=test_job             # Batch job and std out file prefix
suites_for_RTs="SCM_GFS_v15p2,SCM_GFS_v16,SCM_GFS_v17_p8,SCM_HRRR,SCM_RRFS_v1beta,SCM_RAP,SCM_WoFS_v0"

#-----------------------------------------------------------------------
# Create the output and log files if they doesn't exist
#-----------------------------------------------------------------------
if [ ! -f "$TEST_OUTPUT" ]; then
   touch ${TEST_OUTPUT}
fi
if [ ! -f "$TEST_LOGFILE" ]; then
   touch ${TEST_LOGFILE}
fi

#-----------------------------------------------------------------------
# Check if lookup tables and other large datasets have been downloaded
# If any of the phys_data_files are missing, run the download script
#-----------------------------------------------------------------------
for phys_data_file in "${phys_data_files[@]}"; do
  if [ ! -f "${PHYS_DATA_DIR}/${phys_data_file}" ]; then
    ${TOP_DIR}/contrib/get_thompson_tables.sh
    break
  else
    echo "Data file ${phys_data_file} exists..." >> ${TEST_LOGFILE}
  fi
done

for mg_inccn_data_file in "${mg_inccn_data_files[@]}"; do
  if [ ! -f "${PHYS_DATA_DIR}/${mg_inccn_data_file}" ]; then
    ${TOP_DIR}/contrib/get_mg_inccn_data.sh
    break
  else
    echo "Data file ${mg_inccn_data_file} exists..." >> ${TEST_LOGFILE}
  fi
done

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
      # Add --runtime ${runtime} to run_scm.py to reduce runtime for tests
      test_run_cmd="${BIN_DIR}/run_scm.py --file ${TEST_DIR}/rt_test_cases.py -v --runtime_mult 0.1 --run_dir ${RUN_DIR} --bin_dir ${BIN_DIR}"
    else
      test_run_cmd="${BIN_DIR}/run_scm.py --file ${TEST_DIR}/rt_test_cases.py --runtime_mult 0.1 --run_dir ${RUN_DIR} --bin_dir ${BIN_DIR}"
    fi

    . ${ETC_DIR}/${machine}_setup_${compiler}.sh

    if [ $build_it -eq 0 ] ; then
#-----------------------------------------------------------------------
# Build the SCM
#-----------------------------------------------------------------------
      if [ -d "${BIN_DIR}" ] ; then rm -rf ${BIN_DIR}; fi
      mkdir ${BIN_DIR}
      cd ${BIN_DIR}
      cmake -DCMAKE_BUILD_TYPE=${build_type} -DCCPP_SUITES=${suites_for_RTs} ../src >& log.cmake || fail "cmake FAILED, see ${BIN_DIR}/log.cmake"
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
      ln -s ${SRC_DIR}/run_scm.py run_scm.py
      job_name=${job_prefix}_${compiler}_${build_type_lc}
      if ${use_batch_system} ; then
        cp ${ETC_DIR}/${job_submission_script} ${job_submission_script}.tmp
#-----------------------------------------------------------------------
# Substitute COMMAND in job_submission script to use multi_run
#-----------------------------------------------------------------------
        sed "s,^.*COMMAND \= .*,COMMAND \= \"${test_run_cmd}\"," ${job_submission_script}.tmp > ${job_submission_script}.tmp2
        sed "s,^.*WALLTIME \= .*,WALLTIME \= \"${walltime}\"," ${job_submission_script}.tmp2 > ${job_submission_script}.tmp3
        sed "s,^.*SERIAL_MEM \= .*,SERIAL_MEM \= \"${memory}\"," ${job_submission_script}.tmp3 > ${job_submission_script}.tmp4
        sed "s,^.*JOB_NAME \= .*,JOB_NAME \= \"${job_name}\"," ${job_submission_script}.tmp4 > ${job_submission_script}
        chmod +x ${job_submission_script}
        rm -rf ${job_submission_script}.tmp
        rm -rf ${job_submission_script}.tmp2
        rm -rf ${job_submission_script}.tmp3
        rm -rf ${job_submission_script}.tmp4

#-----------------------------------------------------------------------
# Submit job_submission to queue
#-----------------------------------------------------------------------
        ./${job_submission_script}
      else
        ${test_run_cmd} >& ${job_name}.o
      fi

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
    batch_out=${job_prefix}_${compiler}_${build_type_lc}.o*           # Without PID extension

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
#-----------------------------------------------------------------------
# Generate baseline if option is set
#-----------------------------------------------------------------------
  if [[ "$gen_baseline" == "true" ]]; then
    echo "Generating baselines in ${gen_baseline_dir}..." >> ${TEST_OUTPUT}
    if [ -d "${gen_baseline_dir}" ]; then rm -rf ${gen_baseline_dir}; fi
    mkdir -p ${gen_baseline_dir} || \
           fail "FAIL: Could not make directory ${gen_baseline_dir} for generated baseline"

    for compiler in "${compilers[@]}"; do
      for build_type in "${build_types[@]}"; do
        build_type_lc=$(echo "${build_type}" | tr '[A-Z]' '[a-z]')  # Make build_type all lower case
        run_dir_ext=run_${compiler}_${build_type_lc}
        RUN_DIR=$TOP_DIR/scm/${run_dir_ext}
        echo "Copying output from run diretory ${RUN_DIR}..." >> ${TEST_OUTPUT}
        for run_path in ${RUN_DIR}/output_*; do
          [ -d "${run_path}" ] || continue                          # if not a directory, skip
          output_dir="$(basename "${run_path}")"
          gen_baseline_path=${gen_baseline_dir}/${run_dir_ext}/${output_dir}
          if [ ! -d "${gen_baseline_path}" ]; then mkdir -p ${gen_baseline_path}; fi
          cp -r ${run_path}/output.nc ${gen_baseline_path}/output.nc || \
            fail "FAIL:  Could not copy ${run_path} output to ${gen_baseline_path} for baseline generation"
        done
      done #End build type loop for baseline generation
    done   #End compiler loop for baseline generation
  fi
#-----------------------------------------------------------------------
# Compare baseline if option is set
#-----------------------------------------------------------------------
  if [[ "$cmp_baseline" == "true" ]]; then
    echo "Comparing output to baselines in ${cmp_baseline_dir}..." >> ${TEST_OUTPUT}

    n_identical=0
    n_differs=0
    n_comparisons=0
    for compiler in "${compilers[@]}"; do
      for build_type in "${build_types[@]}"; do
        build_type_lc=$(echo "${build_type}" | tr '[A-Z]' '[a-z]')  # Make build_type all lower case
        run_dir_ext=run_${compiler}_${build_type_lc}
        RUN_DIR=$TOP_DIR/scm/${run_dir_ext}
        echo "Comparing output from run diretory ${RUN_DIR}..." >> ${TEST_OUTPUT}
        for run_path in ${RUN_DIR}/output_*; do
          [ -d "${run_path}" ] || continue                          # if not a directory, skip
          output_dir="$(basename "${run_path}")"
          cmp_baseline_path=${cmp_baseline_dir}/${run_dir_ext}/${output_dir}
          ((n_comparisons=n_comparisons+1))
          if cmp ${run_path}/output.nc ${cmp_baseline_path}/output.nc ; then
            echo "Output for ${run_dir_ext}/${output_dir} is IDENTICAL to baseline" >> ${TEST_LOGFILE}
            ((n_identical=n_identical+1))
          else
            echo "Output for ${run_dir_ext}/${output_dir} DIFFERS from baseline" >> ${TEST_LOGFILE}
            ((n_differs=n_differs+1))
          fi
        done
      done #End build type loop for baseline comparison
    done   #End compiler loop for baseline comparison
#-----------------------------------------------------------------------
# PASS = All tests succeeded and output is identical
#-----------------------------------------------------------------------
    if [[ $n_identical -eq $n_comparisons && $n_differs -eq 0 ]] ; then
      echo "ALL TESTS PASSED, OUTPUT IS IDENTICAL." >> ${TEST_OUTPUT}
      msg="PASS"
    else
      echo "ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE." >> ${TEST_OUTPUT}
      msg="FAIL"
    fi
  fi    # End of if cmp_baseline=true
else
  echo "TEST(S) FAILED" >> ${TEST_OUTPUT}
  msg="FAIL"
fi

echo "More output from this test can be found in ${TEST_LOGFILE}" >> ${TEST_OUTPUT}

#-----------------------------------------------------------------------
# Send email with PASS/FAIL message
#-----------------------------------------------------------------------
cat ${TEST_OUTPUT} | mail -s "[ ${msg} ] ccpp-scm test $1" ${users}
