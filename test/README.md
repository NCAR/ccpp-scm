# CCPP SCM Regresstion Test

## Description

This regression test builds and runs the CCPP SCM for all supported compilers, build types, cases
and suites.  It consists of the following scripts:

* ``rt.sh`` - Driver for the regresion test that builds, runs and calls the summarize.sh script
* ``rt_test_cases.py``- list of all supported cases and suites
* ``summarize.sh`` - Called by ``rt.sh`` to parse the output of each test to determine pass/fail		

Currently, the following configurations are supported:

Machine     | Cheyenne       | Hera           |
------------| ---------------|----------------|
Compiler(s) | Intel, GNU     | Intel          |
Build Types | Release, Debug | Release, Debug |

The executable for each build is created in its own directory under the ``scm`` directory:

```
bin_gnu_debug
bin_gnu_release
bin_intel_debug
bin_intel_release
```
For each build, all supported cases and suites are run in its own directory:

```
run_gnu_debug
run_gnu_release
run_intel_debug
run_intel_release
```

The debug tests use a reduced runtime for faster turnaround and to conserve computer resources.

## Usage

On Cheyenne:

```
cd test
./rt.sh cheyenne >& test.out &
```

On Hera:

```
cd test
./rt.sh hera >& test.out &
```

Upon completion, an email summary will be sent to ``$USER@ucar.edu`` or ``$USER.noaa.gov`` depending on the platform.  A summary of the tests will be in the file ``rt_summary$PID.out`` under the ``test`` directory.  More detailed output for each run can be found in the standard output file under each run directory, for example: ``scm/run_intel_release/test_job*``.
