!> \file gmtb_scm_type_defs.f90
!!  Contains type definitions for SCM-related variables and physics-related variables

module gmtb_scm_type_defs

  use gmtb_scm_kinds, only : sp, dp, qp

  implicit none

  integer, parameter :: character_length = 80

  character(len = character_length) :: clear_char = ''

  type scm_state_type

    character(len=character_length)                 :: experiment_name !> name of model configuration file
    character(len=character_length)                 :: model_name !< name of "host" model (must be "GFS" for prototype)
    character(len=character_length)                 :: output_dir !< name of output directory to place netCDF file
    character(len=character_length)                 :: physics_suite_dir !< location of the physics suite XML files for the IPD (relative to the executable path)
    character(len=character_length)                 :: output_file !< name of output file (without the file extension)
    character(len=character_length)                 :: case_name !< name of case initialization and forcing to use (different than experiment name, which names the model run (as a control, experiment_1, etc.))

    contains
      procedure :: create  => scm_state_create

  end type scm_state_type

  contains

  subroutine scm_state_create(scm_state)
    class(scm_state_type)             :: scm_state

    scm_state%experiment_name = clear_char
    scm_state%model_name = clear_char
    scm_state%output_dir = clear_char
    scm_state%physics_suite_dir = clear_char
    scm_state%output_file = clear_char
    scm_state%case_name = clear_char

  end subroutine scm_state_create

end module gmtb_scm_type_defs
