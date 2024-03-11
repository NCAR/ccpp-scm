.. _`chapter: ccpp_interface`:

CCPP Interface
==============

Chapter 6 of the CCPP v6 Technical Documentation
(https://ccpp-techdoc.readthedocs.io/en/v6.0.0/) provides a wealth of
information on the overall process of connecting a host model to the
CCPP framework for calling physics. This chapter describes the
particular implementation within this SCM, including how to set up,
initialize, call, and change a physics suite using the CCPP framework.

Setting up a suite
------------------

Setting up a physics suite for use in the CCPP SCM with the CCPP
framework involves three steps: preparing data to be made available to
physics through the CCPP, running the script to reconcile SCM-provided
variables with physics-required variables, and preparing a suite
definition file.

Preparing data from the SCM
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As described in sections 6.1 and 6.2 of the `CCPP Technical
Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__ a host
model must allocate memory and provide metadata for variables that are
passed into and out of the schemes within the physics suite. As of this
release, in practice this means that a host model must do this for all
variables needed by all physics schemes that are expected to be used
with the host model. For this SCM, all variables needed by the physics
schemes are allocated and documented in the file and are contained
within the derived data type. This derived data type initializes its
component variables in a type-bound procedure. As mentioned in section
6.2 of the `CCPP Technical
Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__, files
containing all required metadata was constructed for describing all
variables in the derived data type. These files are , and . Further,
exists to provide metadata for derived data type definitions and their
instances, which is needed by the ccpp-framework to properly reference
the data. The standard names of all variables in this table must match
with a corresponding variable within one or more of the physics schemes.
A list of all standard names used can be found in .

Editing and running 
~~~~~~~~~~~~~~~~~~~~

General instructions for configuring and running the script can be found
in chapter 8 of the `CCPP Technical
Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__. The
script expects to be run with a host-model-dependent configuration file,
passed as argument . Within this configuration file are variables that
hold paths to the variable definition files (where metadata tables can
be found on the host model side), the scheme files (a list of paths to
all source files containing scheme entry points), the auto-generated
physics schemes makefile snippet, the auto-generated physics scheme caps
makefile snippet, and the directory where the auto-generated physics
caps should be written out to. As mentioned in section
`[section: compiling] <#section: compiling>`__, this script must be run
to reconcile data provided by the SCM with data required by the physics
schemes before compilation – this is done automatically by .

Preparing a suite definition file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The suite definition file is a text file read by the model at compile
time. It is used to specify the physical parameterization suite, and
includes information about the number of parameterization groupings,
which parameterizations that are part of each of the groups, the order
in which the parameterizations should be run, and whether subcycling
will be used to run any of the parameterizations with shorter timesteps.

In addition to the six or so major parameterization categories (such as
radiation, boundary layer, deep convection, resolved moist physics,
etc.), the suite definition file can also have an arbitrary number of
additional interstitial schemes in between the parameterizations to
prepare or postprocess data. In many models, this interstitial code is
not known to the model user but with the suite definition file, both the
physical parameterizations and the interstitial processing are listed
explicitly.

For this release, supported suite definition files used with this SCM
are found in and have default namelist, tracer configuration, and
timesteps attached in . For all of these suites, the physics schemes
have been organized into 3 groupings following how the physics are
called in the UFS Atmosphere model, although no code is executed in the
SCM time loop between execution of the grouped schemes. Several
“interstitial” schemes are included in the suite definition file to
execute code that previously was part of a hard-coded physics driver.
Some of these schemes may eventually be rolled into the schemes
themselves, improving portability.

Initializing/running a suite
----------------------------

The process for initializing and running a suite in this SCM is
described in sections
`[section: physics init] <#section: physics init>`__ and
`[section: time integration] <#section: time integration>`__,
respectively. A more general description of the process for performing
suite initialization and running can also be found in sections 6.4 and
6.5 of the `CCPP Technical
Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__.

Changing a suite
----------------

Replacing a scheme with another
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prior to being able to swap a scheme within a suite, one must first add
a CCPP-compliant scheme to the pool of available schemes in the CCPP
physics repository. This process is described in chapter 2 of the `CCPP
Technical
Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__.

Once a CCPP-compliant scheme has been added to the CCPP physics
repository, the process for modifying an existing suite should take the
following steps into account:

-  Examine and compare the arguments of the scheme being replaced and
   the replacement scheme.

   -  Are there any new variables that the replacement scheme needs from
      the host application? If so, these new variables must be added to
      the host model cap. For the SCM, this involves adding a component
      variable to the derived data type and a corresponding entry in the
      metadata table. The new variables must also be allocated and
      initialized in the type-bound procedure.

   -  Do any of the new variables need to be calculated in an
      interstitial scheme? If so, one must be written and made
      CCPP-compliant itself. The `CCPP Technical
      Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__
      will help in this endeavor, and the process outlined in its
      chapter 2 should be followed.

   -  Do other schemes in the suite rely on output variables from the
      scheme being replaced that are no longer being supplied by the
      replacement scheme? Do these output variables need to be
      derived/calculated in an interstitial scheme? If so, see the
      previous bullet about adding one.

-  Examine existing interstitial schemes related to the scheme being
   replaced.

   -  There may be scheme-specific interstitial schemes (needed for one
      specific scheme) and/or type-generic interstitial schemes (those
      that are called for all schemes of a given type, i.e. all PBL
      schemes). Does one need to write analogous scheme-specific
      interstitial schemes for the replacement?

   -  Are the type-generic interstitial schemes relevant or do they need
      to be modified?

-  Depending on the answers to the above considerations, edit the suite
   definition file as necessary. Typically, this would involve finding
   the elements associated with the scheme to be replaced and its
   associated interstitial elements and simply replacing the scheme
   names to reflect their replacements. See chapter 4 of the `CCPP
   Technical
   Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/>`__ for
   further details.

Modifying “groups” of parameterizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The concept of grouping physics in the suite definition file (currently
reflected in the elements) enables “groups” of parameterizations to be
called with other computation (perhaps related to the dycore, I/O, etc.)
in between. In the suite definition file included in this release, three
groups are specified, but currently no computation happens between calls
for these groups. However, one can edit the groups to suit the needs of
the host application. For example, if a subset of physics schemes needs
to be more tightly connected with the dynamics and called more
frequently, one could create a group consisting of that subset and place
a call in the appropriate place in the host application. The remainder
of the parameterizations groups could be called using calls in a
different part of the host application code.

Subcycling parameterizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The suite definition file allows subcycling of schemes, or calling a
subset of schemes at a smaller time step than others. The element in the
suite definition file controls this function. All schemes within such an
element are called times during one call. An example of this is found in
the suite definition file, where the surface schemes are executed twice
for each timestep (implementing a predictor/corrector paradigm). Note
that no time step information is included in the suite definition file.
**If subcycling is used for a set of parameterizations, the smaller time
step must be an input argument for those schemes. This is not handled
automatically by the ccpp-framework yet.**

Adding variables
----------------

.. _adding_physics_only_variable:

Adding a physics-only variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose that one wants to add the variable to a scheme that spans the
depth of the column and that this variable is internal to physics, not
part of the SCM state or subject to external forcing. Here is how one
would do so:

#. First, add the new variable to the derived data type definition in .
   Within the definition, you’ll notice that there are nested derived
   data types (which contain most of the variables needed by the physics
   and are used for mainly legacy reasons) and several other
   integers/reals/logicals. One could add the new variable to one of the
   nested GFS derived data types if the variable neatly fits inside one
   of them, but it is suggested to bypass the GFS derived data types and
   add a variable directly to the type definition:

   .. code:: fortran

      real(kind=kind_phys), allocatable :: foo(:,:)

#. Second, within the subroutine, add an allocate and initialization
   statement.

   .. code:: fortran

      allocate(foo(n_columns, n_levels))
      physics%foo = 0.0

   Note that even though foo only needs to have the vertical dimension,
   it is also allocated with the dimension as the first dimension since
   this model is intended to be used with multiple independent columns.
   Also, the initialization in this creation subroutine can be
   overwritten by an initialization subroutine associated with a
   particular scheme.

#. At this point, these changes are enough to allocate the new variable
   ( is called in the main subroutine of ), although this variable
   cannot be used in a physics scheme yet. For that, you’ll need to add
   an entry in the corresponding metadata file. See section 2.2 of the
   `CCPP Technical
   Documentation <https://ccpp-techdoc.readthedocs.io/en/v6.0.0/CompliantPhysicsParams.html#metadata-table-rules>`__
   for more information regarding the format.

#. On the physics scheme side, there will also be a metadata file entry
   for . For example, say that scheme uses . If is further initialized
   in ’s subroutine, a metadata entry for must be found in the
   corresponding section in the metadata file. If it is used in ’s run
   subroutine, a metadata entry for foo must also appear in the metadata
   file section for . The metadata entry on the physics scheme side has
   the same format as the one on the host model side described above.
   The standard name, rank, type, and kind must match the entry from the
   host model table. Others attributes (local name, units (assuming that
   an automatic conversion exists in the ccpp-framework), long_name,
   intent) can differ. The local name corresponds to the name of the
   variable used within the scheme subroutine, and the intent attribute
   should reflect how the variable is actually used within the scheme.

   Note: In addition to the metadata file, the argument list for the
   scheme subroutine must include the new variable (i.e., must actually
   be in the argument list for and be declared appropriately in regular
   Fortran).

If a variable is declared following these steps, it can be used in any
CCPP-compliant physics scheme and it will retain its value from timestep
to timestep. A variable will ONLY be zeroed out (either every timestep
or periodically) if it is in the or data types. So, if one needs the new
variable to be ‘prognostic’, one would need to handle updating its value
within the scheme, something like:

.. math:: \text{foo}^{t+1} = \text{foo}^t + \Delta t*\text{foo\_tendency}

Technically, the host model can “see” foo between calls to physics
(since the host model allocated its memory at initialization), but it
will not be touching it.

Adding a prognostic SCM variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following instructions are valid for adding a passive, prognostic
tracer to the SCM. Throughout these instructions, the new tracer is
called ‘smoke’.

#. Add a new tracer to the SCM state. In do the following:

   -  Add an index for the new tracer in the definition.

   -  Do the following in the subroutine:

      -  Increment

      -  Set =

      -  Set =

      -  Note: is initialized to zero in this subroutine already, so
         there is no need to do so again.

#. Initialize the new tracer to something other than zero (from an input
   file).

   -  Edit an existing input file (in ): add a field in the ‘initial’
      group of the NetCDF file(s) (with vertical dimension in pressure
      coordinates) with an appropriate name in one (or all) of the input
      NetCDF files and populate with whatever values are necessary to
      initialize the new tracer.

   -  Create a new input variable to read in the initialized values. In
      :

      -  Add a new input variable in

         .. code:: fortran

            real(kind=dp), allocatable              :: input_smoke(:)

      -  In , allocate and initialize the new variable to 0.

   -  Read in the input values to initialize the new tracer. In :

      -  Add a variable under the initial profile section:

         .. code:: fortran

            real(kind=dp), allocatable  :: input_smoke(:) !< smoke profile (fraction)

      -  Add the new input variable to the allocate statement.

      -  Read the values in from the file:

         .. code:: fortran

            call check(NF90_INQ_VARID(grp_ncid,"smoke",varID))
                    call check(NF90_GET_VAR(grp_ncid,varID,input_smoke))

      -  set =

   -  Interpolate the input values to the model grid. Edit :

      -  Add a loop over the columns to call that puts on grid levels in

         .. code:: fortran

            do i=1, scm_state%n_cols
                            call interpolate_to_grid_centers(scm_input%input_nlev, scm_input%input_pres, scm_input%input_smoke, scm_state%pres_l(i,1,:), &
                                    scm_state%n_levels, scm_state%state_tracer(i,1,:,scm_state%smoke_index,1), last_index_init, 1)
                    end do

   -  At this point, you have a new tracer initialized to values
      specified in the input file on the model vertical grid, but it is
      not connected to any physics or changed by any forcing.

#. For these instructions, we’ll assume that the tracer is not subject
   to any external forcing (e.g., horizontal advective forcing, sources,
   sinks). If it is, further work is required to:

   -  One needs to provide data on how tracer is forced in the input
      file, similar to specifying its initial state, as above.

   -  Create, allocate, and read in the new variable for forcing
      (similar to above).

   -  Add to (similar to above, but interpolates the forcing to the
      model grid and model time).

   -  Add statements to time loop to handle the first time step and
      time-advancing.

   -  Edit in .

#. In order to connect the new tracer to the CCPP physics, perform steps
   1-4 in section `1.4.1 <#adding_physics_only_variable>`__ for adding a
   physics variable. In addition, do the following in order to associate
   the variable with variables used in the physics through a pointer:

   -  Point the new physics variable to in .

#. There may be additional steps depending on how the tracer is used in
   the physics and how the physics scheme is integrated with the current
   GFS physics suite. For example, the GFS physics has two tracer
   arrays, one for holding tracer values before the physics timestep ()
   and one for holding tracer values that are updated during/after the
   physics (). If the tracer needs to be part of these arrays, there are
   a few additional steps to take. If you need help, please post on the
   support forum at:
   https://dtcenter.org/forum/ccpp-user-support/ccpp-single-column-model.
