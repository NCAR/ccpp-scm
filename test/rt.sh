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
set -eux    # Uncomment for debugging
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

machines=( hera cheyenne )
users=( schramm@ucar.edu )
executable_name=gmtb_scm

[[ $# -eq 0 ]] && usage
if [ "$1" = "-h" ] ; then usage ; fi

export machine=${1}
machine=`echo "${machine}" | tr '[A-Z]' '[a-z]'`  # Make machine name all lower case
machine=`echo "${machine^}"`                      # Capitalize first letter for setup script name
build_it=0                                        # Set to 1 to skip build (for testing pass/fail criteria)

build_types=( Release Debug )                     # Set all instances of CMAKE_BUILD_TYPE

if [ "${machine}" == "Cheyenne" ] ; then
  compilers=( intel gnu )
  job_submission_script=gmtb_scm_qsub_example.py
else
  compilers=( intel )
  job_submission_script=gmtb_scm_slurm_example.py
fi

#-----------------------------------------------------------------------
# Set some directories
#-----------------------------------------------------------------------
PID=$$
TEST_DIR=$( pwd )                   # Directory with this script
TEST_OUTPUT=${TEST_DIR}/rt${PID}.out
TOP_DIR=$TEST_DIR/..
SCM_DIR=$TOP_DIR/scm
ETC_DIR=$TOP_DIR/scm/etc
SRC_DIR=$TOP_DIR/scm/src
PHYS_DATA_DIR=$TOP_DIR/scm/data/physics_input_data
phys_data_files=( CCN_ACTIVATE.BIN freezeH2O.dat qr_acr_qg.dat qr_acr_qs.dat )

#-----------------------------------------------------------------------
# Check if lookup tables and other large datasets have been downloaded
# If any of the phys_data_files are missing, run the download script
#-----------------------------------------------------------------------
for phys_data_file in "${phys_data_files[@]}"; do
  if [ ! -f "${PHYS_DATA_DIR}/${phys_data_file}" ]; then
    ${TOP_DIR}/contrib/get_thompson_tables.sh
    break
  else
    echo " Data file ${phys_data_file} exists..."
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

    echo "Starting ${compiler} ${build_type} build"

    build_type_lc=`echo "${build_type}" | tr '[A-Z]' '[a-z]'`  # Make build_type all lower case
    BIN_DIR=$TOP_DIR/scm/bin_${compiler}_${build_type_lc}      # Assign a bin and run dir
    RUN_DIR=$TOP_DIR/scm/run_${compiler}_${build_type_lc}      # for each build/run in test
    BUILD_OUTPUT=${BIN_DIR}/build.out
    test_run_cmd="${RUN_DIR}/multi_run_gmtb_scm.py -f ${TEST_DIR}/rt_test_cases.py"
    walltime="walltime=00:40:00"

    . ${ETC_DIR}/${machine}_setup_${compiler}.sh

    if [ $build_it -eq 0 ] ; then
#-----------------------------------------------------------------------
# Build the SCM
#-----------------------------------------------------------------------
      if [ -d "${BIN_DIR}" ] ; then rm -rf ${BIN_DIR}; fi
      mkdir ${BIN_DIR}
      cd ${BIN_DIR}
      cmake -DCMAKE_BUILD_TYPE=${build_type} ../src >& log.cmake
      if [ "${build_type_lc}" == "release" ] ; then
        make -j4 >& ${BUILD_OUTPUT} || fail "Build ${machine} ${compiler} FAILED"
      else   # Intel debug fails with -j4
        make >& ${BUILD_OUTPUT} || fail "Build ${machine} ${compiler} FAILED"
      fi
    fi    # End of skip build for testing

    exec_file=${BIN_DIR}/${executable_name}
    if [ -f "${exec_file}" ]; then
      echo "SUCCEED: ${machine} ${compiler} executable file ${exec_file} exists" >> ${TEST_OUTPUT}
    else
      echo "FAIL: ${machine} ${compiler} executable file ${exec_file} does NOT exist" >> ${TEST_OUTPUT}
    fi

#-----------------------------------------------------------------------
# Set up the run directory for the build
#-----------------------------------------------------------------------
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
    sed "s,^.*WALLTIME \= .*,WALLTIME \= \"${walltime}\"," ${job_submission_script}.tmp2 > ${job_submission_script}
    chmod +x ${job_submission_script}
    rm -rf ${job_submission_script}.tmp
    rm -rf ${job_submission_script}.tmp2

#-----------------------------------------------------------------------
# Submit job_submission to queue
#-----------------------------------------------------------------------
    ./${job_submission_script}

  done #End build type loop
done   #End compiler loop
