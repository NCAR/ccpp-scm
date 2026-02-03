!> \file ccpp_config.F90
!!  Contains CCPP configuration

module ccpp_config
  use mpi_f08, only: MPI_Comm
  implicit none
!> \section arg_table_ccpp_config
!! \htmlinclude ccpp_config.html
!!
  !> @var The default loop counter indicating outside of a subcycle loop
  integer, parameter :: CCPP_DEFAULT_LOOP_CNT = -999
  integer, parameter :: CCPP_DEFAULT_LOOP_MAX = -999

  !> @var The default values for block and chunk numbers indicating invalid data
  integer, parameter :: CCPP_DEFAULT_CHUNK_NUMBER  = -999

  !> @var The default maximum number of threads for CCPP
  integer, parameter :: CCPP_DEFAULT_THREAD_NUMBER = -999
  integer, parameter :: CCPP_DEFAULT_THREAD_COUNT  = -999

!! \section arg_table_ty_ccpp_config
!! \htmlinclude ty_ccpp_config.html
!!
  type :: ty_ccpp_config
     ! CCPP-internal variables for physics schemes
     integer            :: loop_cnt      = CCPP_DEFAULT_LOOP_CNT
     integer            :: loop_max      = CCPP_DEFAULT_LOOP_MAX
     integer            :: thrd_no       = CCPP_DEFAULT_THREAD_NUMBER
     integer            :: thrd_cnt      = CCPP_DEFAULT_THREAD_COUNT
     integer            :: chunk_no      = CCPP_DEFAULT_CHUNK_NUMBER
     integer            :: chunk_cnt     = 1
     integer            :: ccpp_instance = 1
     integer            :: ccpp_errflg   = 0
     character(len=512) :: ccpp_errmsg   = ''
   contains
     procedure :: initialized  => ccpp_cfg_initialized
  end type ty_ccpp_config

contains

  function ccpp_cfg_initialized(ccpp_d) result(initialized)
    implicit none
    class(ty_ccpp_config) :: ccpp_d
    logical :: initialized
    initialized = ccpp_d%thrd_no  /= CCPP_DEFAULT_THREAD_NUMBER .or. &
                  ccpp_d%chunk_no /= CCPP_DEFAULT_CHUNK_NUMBER
  end function ccpp_cfg_initialized
  
end module ccpp_config
