.. _`chapter: introduction`:

Introduction
============

A single-column model (SCM) can be a valuable tool for diagnosing the
performance of a physics suite, from validating that schemes have been
integrated into a suite correctly to deep dives into how physical
processes are being represented by the approximating code. This SCM has
the advantage of working with and being developed alongside the Common Community Physics Package
(CCPP), a library of physical parameterizations for atmospheric
numerical models (CCPP physics) and the associated framework for connecting potentially
any atmospheric model to physics suites constructed from its member
parameterizations (CCPP framework). In fact, this SCM serves as perhaps the simplest
example for using the CCPP and its framework in an atmospheric model.
This version contains all parameterizations of NOAA’s evolved
operational GFS v16 suite (implemented in 2021), plus additional
developmental schemes. The schemes are grouped in five supported suites
described in detail in the `CCPP Scientific
Documentation <https://dtcenter.ucar.edu/GMTB/v8.0.0/sci_doc/>`__
(GFS_v16, GFS_v16_RRTMGP, GFS_v17_p8_ugwpv1, HRRR_gf, and WoFS_v0).

This document serves as both the User and Technical Guides for this
model. It contains a Quick Start Guide with instructions for obtaining
the code, compiling, and running a sample test case; an explanation for
what is included in the repository; a brief description of the operation
of the model; a description of how cases are set up and run, and
finally, an explanation for how the model interfaces with physics
through the CCPP infrastructure.

| Please refer to the release web page for further documentation and
  user notes:
| https://dtcenter.org/community-code/common-community-physics-package-ccpp/download

Version Notes
-------------


The CCPP SCM v8.0.0 contains the following changes since v7.0.

CCPP Framework

-  Update SCM to use new CCPP framework, capgen-v1

   -  Flipping of arrays, e.g. for cases in which the physics is
      oriented bottom to top and the host model top to bottom.
   -  Kind transformations e.g., double to single precision, creating
      flexibility in connections between host and physics.
   -  Superior memory management: suite-owned variables do not need to
      be allocated by the host model.
   -  Introspection. Capgen provides a data file containing all
      necessary scheme and host files, required dependencies, and all
      variable lists needed by host-model. This can be used directly by
      the build-system, or for debugging purposes.
   -  Integer state machine. Allows physics working precision to be
      build-time configurable. Currently this is hardcoded and included
      as part of the CCPP Physics library.
   -  Registered scalar-index dimensions for threading and ensembles.
      This allows for host configurability for multi-threaded and
      multi-instance implementations.
   -  Write-if-changed build integration. Facilitates quicker model
      build time and debugging.
   -  Enhanced error-checking: a Fortran-vs-metadata validator tool
      identifies inconsistencies and prompts the developer to address it
      at build-time.
   -  Handling of constituents array for the flexible management of
      atmospheric tracers (like moisture, aerosols, or chemical species).
   -  Ease of code inspection, since the code does not involve a static
      driver with conditional statements.
   -  Improved computational performance.

CCPP Physics

-  Newer code base for the parameterizations and suites supported in CCPP
   v7
  
CCPP SCM

-  SCM code cleanup: Removed unsupported code and eliminated
   nonfunctional or unimplemented components.

-  New SCM Case: Atlantic Tradewind Ocean-Atmosphere Mesoscale
   Interaction Campaign (ATOMIC) for January/February 2020 marine shallow
   convection. The ATOMIC test case consists of seven 33-hour periods.

Updated technical documentation

Limitations
~~~~~~~~~~~

This release bundle has some known limitations:

-  Using the HRRR_gf and WoFS_v0 suites for cases where deep
   convection is expected to be active will likely produce
   strange/unreliable results, unless the forcing has been modified to
   account for the deep convection. This is because forcing for existing
   cases assumes a horizontal scale for which deep convection is
   subgrid-scale and is expected to be parameterized. The suites without
   convection are intended for use with regional UFS simulations with
   horizontal scale small enough not to need a deep convection
   parameterization active, and it does not contain a deep convective
   scheme. Nevertheless, these suites are included with the SCM as-is
   for research purposes.

-  The provided cases over land points cannot use an LSM at this time
   due to the lack of initialization data for the LSMs. Therefore, for
   the provided cases over land points (ARM_SGP_summer_1997\_\* and
   LASSO\_\*, where sfc_type = 1 is set in the case configuration file),
   prescribed surface fluxes must be used:

   -  surface sensible and latent heat fluxes must be provided in the
      case data file

   -  sfc_flux_spec must be set to true in the case configuration file

   -  the surface roughness length in cm must be set in the case
      configuration file

   -  the suite defintion file used (physics_suite variable in the case
      configuration file) must have been modified to use prescribed
      surface fluxes rather than an LSM.

   -  NOTE: If one can develop appropriate initial conditions for the
      LSMs for the supplied cases over land points, there should be no
      technical reason why they cannot be used with LSMs, however.

-  Using the SCM over a land point with an LSM is
   possible through the use of UFS initial conditions (see
   :numref:`Section %s <UFScasegen>`).

-  There are several capabilities of the developmental code that have
   not been tested sufficiently to be considered part of the supported
   release. Those include additional parameterizations and the CCPP
   Suite Simulator. Users that want to use experimental capabilities
   should refer to :numref:`Subsection %s <development_code>`.
