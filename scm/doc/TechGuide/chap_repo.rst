.. _`chapter: repository`:

Repository
==========

What is included in the repository?
-----------------------------------

The repository contains all code required to build the CCPP SCM and
scripts that can be used to obtain data to run it (e.g. downloading
large initialization tables for the Thompson microphysics schemes
discussed in subsection
`[subsection: singlerunscript] <#subsection: singlerunscript>`__ and
processed case data). It is functionally separated into 3 subdirectories
representing the SCM model infrastructure ( directory), the CCPP
infrastructure ( directory), and the CCPP physics schemes ( directory).
The entire repository resides on Github’s NCAR space, and the and
directories are git submodules that point to repositories and on the
same space. The structure of the entire repository is represented below.
Note that the repository also contains files needed for using the CCPP
with the UFS Atmosphere host model that uses the Finite-Volume
Cubed-Sphere (FV3) dynamical core.
