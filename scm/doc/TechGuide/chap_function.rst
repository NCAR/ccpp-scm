.. _`chapter: algorithm`:

Algorithm
=========

Algorithm Overview
------------------

Like most SCMs, the algorithm for the CCPP SCM is quite simple. In a
nutshell, the SCM code performs the following:

-  Read in an initial profile and the forcing data.

-  Create a vertical grid and interpolate the initial profile and
   forcing data to it.

-  Initialize the physics suite.

-  Perform the time integration, applying forcing and calling the
   physics suite each time step.

-  Output the state and physics data.

In this chapter, it will briefly be described how each of these tasks is
performed.

Reading input
-------------

The following steps are performed at the beginning of program execution:

#. Call in the module to read in the ` <#subsection: case config>`__ and
   ` <#subsection: physics config>`__ namelists. This subroutine also
   sets some variables within the derived type from the data that was
   read.

#. Call (or if using the DEPHY format) in the module to read in the
   `case input data file <#subsection: case input>`__. This subroutine
   also sets some variables within the derived type from the data that
   was read.

#. Call in the module to read in the reference profile data. This
   subroutine also sets some variables within the derived type from the
   data that was read. At this time, there is no “standard” format for
   the reference profile data file. There is a statement within the
   subroutine that reads in differently-formatted data. If adding a new
   reference profile, it will be required to add a section that reads
   its data in this subroutine.

Setting up vertical grid and interpolating input data
-----------------------------------------------------

The CCPP SCM uses pressure for the vertical coordinate (lowest index is
the surface). The pressure levels are calculated using the surface
pressure and coefficients (:math:`a_k` and :math:`b_k`), which are taken
directly from FV3 code (). For vertical grid options, inspect for valid
values of . The default vertical coordinate uses 127 levels and sets to
the empty string. Alternatively, one can specify the (:math:`a_k` and
:math:`b_k`) coefficients via an external file in the directory and pass
it in to the SCM via the argument of the run script. If changing the
number of vertical levels or the algorithm via the or run script
arguments, be sure to check and that the vertical coordinate is as
inteneded.

After the vertical grid has been set up, the state variable profiles
stored in the derived data type are interpolated from the input and
reference profiles in the subroutine of the module.

.. _`section: physics init`:

Physics suite initialization
----------------------------

With the CCPP framework, initializing a physics suite is a 3-step
process:

#. Initialize variables needed for the suite initialization routine. For
   suites originating from the GFS model, this involves setting some
   values in a derived data type used in the initialization subroutine.
   Call the suite initialization subroutine to perform suite
   initialization tasks that are not already performed in the routines
   of the CCPP-compliant schemes (or associated initialization stages
   for groups or suites listed in the suite definition file). Note: As
   of this release, this step will require another suite intialization
   subroutine to be coded for a non-GFS-based suite to handle any
   initialization that is not already performed within CCPP-compliant
   scheme initialization routines.

#. Associate the variables with the appropriate pointers in the derived
   data type.

#. Call with the derived data type as input. This call executes the
   initialization stages of all schemes, groups, and suites that are
   defined in the suite definition file.

.. _`section: time integration`:

Time integration
----------------

The CCPP SCM uses the simple forward Euler scheme for time-stepping.

During each step of the time integration, the following sequence occurs:

#. Update the elapsed model time.

#. Calculate the current date and time given the initial date and time
   and the elapsed time.

#. Call the subroutine in the module to interpolate the forcing data in
   space and time.

#. Recalculate the pressure variables (pressure, Exner function,
   geopotential) in case the surface pressure has changed.

#. Call in the module. Within this subroutine:

   -  Call the appropriate subroutine from the module.

   -  For each column, call to call all physics schemes within the suite
      (this assumes that all suite parts are called sequentially without
      intervening code execution)

#. Check to see if output should be written during the current time step
   and call in the module if necessary.

Writing output
--------------

Output is accomplished via writing to a NetCDF file. If not in the
initial spin-up period, a subroutine is called to determine whether data
needs to be added to the output file during every timestep. Variables
can be written out as instantaneous or time-averaged and there are 5
output periods:

#. one associated with how often instantaneous variables should be
   written out (controlled by the run script variable).

#. one associated with how often diagnostic (either instantaneous or
   time-averaged) should be written out (controlled by the run script
   variable)

#. one associated with the shortwave radiation period (controlled by
   variable in the physics namelist)

#. one associated with the longwave radiation period (controlled by the
   variable in the physics namelist)

#. one associated with the minimum of the shortwave and longwave
   radiation intervals (for writing output if any radiation is called)

Further, which variables are output and on each interval are controlled
via the source file. Of course, any changes to this file must result in
a recompilation to take effect. There are several subroutines for
initializing the output file () and for appending to it () that are
organized according to their membership in physics derived data types.
See the source file to understand how to change output variables.
