Preface
=======

Meaning of typographic changes and symbols
------------------------------------------

The table below describes the type changes and symbols used in this document.

.. _scheme_suite_table:

.. list-table:: *Type changes and symbols used in this guide.*
   :header-rows: 1

   * - Typeface or Symbol
     - Meaning
     - Examples
   * - ``AaBbCc123``
     - 
         * The names of commands, files, and directories
         * On-screen terminal output
     - 
         * Edit your ``.bashrc`` file 
         * Use ``ls -a`` to list all files. 
         * ``host$ You have mail!``
   * - *AaBbCc123*
     - 
         * The names of SCM- and CCPP-specific terms, subroutines, etc.
         * Captions for figures, tables, etc.
     - 
         * Each scheme must include at least one of the following subroutines: ``{schemename}_timestep_init``, ``{schemename}_init``, ``{schemename}_run``, ``{schemename}_finalize``, and ``{schemename}_timestep_finalize``.
         * *Listing 2.1: Fortran template for a CCPP-compliant scheme showing the* _run *subroutine.*
   * - **AaBbCc123**
     - Words or phrases requiring particular emphasis
     - Fortran77 code should **not** be used
