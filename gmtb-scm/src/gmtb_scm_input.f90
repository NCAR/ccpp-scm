!> \file gmtb_scm_input.f90
!!  Contains input-related subroutines -- reading in model configuration from file or the command line and reading in the case
!!  initial conditions and forcing; also contains reference profile input (temporarily hard-coded).

module gmtb_scm_input

use gmtb_scm_kinds, only : sp, dp, qp
use netcdf

implicit none

contains

!> \ingroup SCM
!! @{
!! \defgroup input gmtb_scm_input
!! @{
!! Contains input-related subroutines -- reading in model configuration from file or the command line and reading in the case
!! initial conditions and forcing; also contains reference profile input (temporarily hard-coded).

!> Subroutine to get basic model configuration data from a namelist file (whose name is specified as the first argument on the command line) and from
!! data entered on the command line with the format: "var1='string' var2=d.d var3=i". Namelist variables listed on the command line
!! override those that are specified in the external namelist file. The case configuration namelist variables are also written out to
!! a namelist file placed in the output directory. Note: This routine uses GET_COMMAND which is an intrinsic routine in the Fortran 2003 standard. This
!! requires that the compiler supports this standard.
subroutine get_config_nml(scm_state)
  use gmtb_scm_type_defs, only : scm_state_type

  type(scm_state_type), target, intent(inout) :: scm_state

  character(len=80)    :: experiment_name !< name of the experiment configuration file (usually case name)
  character(len=80)    :: model_name !< name of the host model (currently only GFS supported)
  character(len=80)    :: case_name !< name of case initialization and forcing dataset
  real(kind=dp)           :: dt !< time step in seconds
  real(kind=dp)             :: runtime !< total runtime in seconds
  real(kind=dp)             :: output_frequency !< freqency of output writing in seconds
  integer             :: n_levels !< number of model levels (currently only 64 supported)
  integer              :: n_columns !< number of columns to use
  integer              :: n_time_levels
  integer              :: time_scheme !< 1 => forward Euler, 2 => filtered leapfrog
  character(len=80)    :: output_dir !< name of the output directory
  character(len=80)    :: output_file !< name of the output file (without the file extension)
  character(len=80)    :: case_data_dir !< path to the directory containing case initialization and forcing data
  character(len=80)    :: vert_coord_data_dir !< path to the directory containing vertical coordinate data
  character(len=80)    :: case_config_dir !< path to the directory containing case configuration files (relative to build dir)
  integer              :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer              :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  real(kind=dp)              :: relax_time !< relaxation time scale (s)
  logical              :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
  integer              :: sfc_type !< 0: sea surface, 1: land surface, 2: sea-ice surface
  integer              :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere
  integer              :: year, month, day, hour

  character(len=80), allocatable    :: physics_suite(:) !< name of the physics suite name (currently only GFS_operational supported)
  character(len=65), allocatable    :: physics_nml(:)
  character(len=80)                :: physics_suite_dir !< location of the physics suite XML files for the IPD (relative to the executable path)
  real(kind=dp), allocatable       :: column_area(:)
  !integer, allocatable              :: n_phy_fields(:) !< number of fields in the data type sent through the IPD

  integer                           :: i, last_physics_specified
  integer                        :: ioerror(5)

  INTEGER, PARAMETER    :: buflen = 255
  CHARACTER(LEN=buflen) :: buf, opt_namelist_vals
  CHARACTER(1)             :: response

  NAMELIST /case_config/ model_name, n_columns, case_name, dt, time_scheme, runtime, output_frequency, &
    n_levels, output_dir, output_file, case_data_dir, vert_coord_data_dir, thermo_forcing_type, mom_forcing_type, relax_time, &
    sfc_type, sfc_flux_spec, reference_profile_choice, year, month, day, hour

  NAMELIST /physics_config/ physics_suite, physics_suite_dir, physics_nml, column_area

  !>  \section get_config_alg Algorithm
  !!  @{

  !> Define default values for case configuration (to be overridden by external namelist file or command line arguments)
  case_config_dir = '../etc/case_config'

  model_name = 'GFS'
  n_columns = 1
  case_name = 'twpice'
  dt = 600.0
  time_scheme = 2
  runtime = 2138400.0
  output_frequency = 600.0
  n_levels = 64
  output_dir = 'output'
  output_file = 'output'
  case_data_dir = '../data/processed_case_input'
  vert_coord_data_dir = '../data/vert_coord_data'
  thermo_forcing_type = 2
  mom_forcing_type = 3
  relax_time = 7200.0
  sfc_flux_spec = .false.
  sfc_type = 0
  reference_profile_choice = 1
  year = 2006
  month = 1
  day = 19
  hour = 3

  physics_suite_dir = '../../gmtb-ccpp/examples/'

  last_physics_specified = -1

  !> Parse the command line arguments.
  opt_namelist_vals = ''

  !> - Get the entire command used to invoke the SCM (including arguments)
  CALL GET_COMMAND(command=buf, status=ioerror(1))

  if (ioerror(1) /= 0) then
    write(*,*) 'An error was encountered reading the command line: error code: '
  else
    !> - Discard SCM invocation command.
    i = INDEX(buf, ' ') !gets position of first space
    buf = trim(adjustl(buf(i+1:))) !discards invoking command from buf

    !> - Get experiment name from first argument (if not present, use default values in this source file)
    i = INDEX(buf, ' ') !gets position of second space (if present)
    experiment_name = trim(adjustl(buf(:i)))
    buf = trim(adjustl(buf(i+1:))) !discards experiment name from buf

    !> - Put remaining command line arguments into "internal" namelist to read in after the external file.
    opt_namelist_vals = '&case_config '//trim(buf)//' /' !assigns rest of command line into case_config namelist

    !> Attempt to read in external namelist file.
    open(unit=1, file=trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml', status='old', action='read', &
      iostat=ioerror(2))
    if(ioerror(2) /= 0) then
      write(*,*) 'There was an error opening the file '//trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml '//'&
        Error code = ',ioerror(2)
    else
      read(1, NML=case_config, iostat=ioerror(3))
    end if

    if(ioerror(3) /= 0) then
      write(*,*) 'There was an error reading the namelist case_config in the file '&
        //trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml'//'&
        Error code = ',ioerror(3)
    end if

    !Using n_columns, allocate memory for the physics suite names and number of fields needed by each. If there are more physics suites
    !than n_columns, notify the user and stop the program. If there are less physics suites than columns, notify the user and attempt to
    !continue (getting permission from user), filling in the unspecified suites as the same as the last specified suite.
    allocate(physics_suite(n_columns), physics_nml(n_columns), column_area(n_columns))

    do i=1, n_columns
      physics_suite(i) = 'none'
    end do

    if(ioerror(2) == 0) then
      read(1, NML=physics_config, iostat=ioerror(4))
    end if
    close(1)

    if(ioerror(4) /= 0) then
      write(*,*) 'There was an error reading the namelist physics_config in the file '&
        //trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml'//'&
        Error code = ',ioerror(4)
      write(*,*) 'Check to make sure that the number of specified physics suites is not greater than n_columns '&
        //'in the case_config namelist. Stopping...'
      STOP
    else
      !check to see if the number of physics_suite_names matches n_cols; if physics_suite_names is less than n_cols, assume that there are multiple columns using the same physics
      do i=1, n_columns
        if (physics_suite(i) == 'none') then
          if(i == 1 ) then
            write(*,*) 'No physics suites were specified in '//trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml. '&
              //'Please edit this file and start again.'
            STOP
          else
            if(last_physics_specified < 0) last_physics_specified = i-1
            !only ask for response the first time an unspecified physics suite is found for a column
            if(last_physics_specified == i-1) then
              write(*,*) 'Too few physics suites were specified for the number of columns in '&
                //trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml. All columns with unspecified physics are set '&
                //'to the last specified suite. Is this the desired behavior (y/n)?'
              read(*,*) response
            end if
            if (response == 'y' .or. response == 'Y') then
              physics_suite(i) = physics_suite(last_physics_specified)
              physics_nml(i) = physics_nml(last_physics_specified)
              !n_phy_fields(i) = n_phy_fields(last_physics_specified)
            else
              write(*,*) 'Please edit '//trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml to contain the same '&
                //' number of physics suites as columns and start again. Stopping...'
              STOP
            end if
          end if !check on i
        ! else
        !   !check to see if n_phy_fields was initialized for this suites
        !   if (n_phy_fields(i) < 0) then
        !     write(*,*) 'The variable n_phy_fields was not initialized for the physics suite '//trim(physics_suite(i))//'. Please '&
        !       //'edit '//trim(adjustl(case_config_dir))//'/''//trim(experiment_name)//'.nml to contain values of this variable '&
        !       //'for each suite.'
        !   end if
        end if !check on physics_suite
      end do
    end if

    if(ioerror(2) /= 0 .or. ioerror(3) /=0 .or. ioerror(4) /=0) then
      write(*,*) 'Since there was an error reading in the '//trim(adjustl(case_config_dir))//'/'//trim(experiment_name)//'.nml'&
        //' file, the default values of the namelist variables will be used, modified by any values contained on the command line.'
    end if

    !> "Read" in internal namelist variables from command line.
    write(*,*) 'Loading optional namelist variables from command line...'
    read(opt_namelist_vals, NML=case_config, IOSTAT=ioerror(5))

    !> Check for namelist reading errors. If errors are found, output the namelist variables that were loaded to console and ask user for input on whether to proceed.
    if(ioerror(5) /= 0) then
      write(*,*) 'There was an error reading the optional namelist variables from the command line: error code ',ioerror(5)
    end if
    if(ioerror(2) /= 0 .or. ioerror(3) /=0 .or. ioerror(4) /= 0 .or. ioerror(5) /= 0) then
      write(*,*) 'Since there was an error either reading the namelist file or an error reading namelist variables from the &
        command line, the namelist variables that are available to use are some combination of the default values and the values &
        that were able to be read. Please look over the namelist variables below to make sure they are set as intended.'
      write(*,NML=case_config)
      write(*,*) 'Continue with these values? (y/n):'
      read(*,*) response

      if (response /= 'y' .and. response /= 'Y') THEN
        write(*,*) 'Stopping...'
        stop
      end if
    end if
  end if

  !> If there was an error reading from the command line, use the default values of the namelist variables and check with the user whether to continue.
  if(ioerror(1) /= 0) then
    write(*,*) 'There was an error reading from the command line: error code ',ioerror(1)
    write(*,NML=case_config)
    write(*,*) 'Continue with default values? (y/n):'
    read(*,*) response

    if (response /= 'y' .and. response /= 'Y') THEN
      write(*,*) 'Stopping...'
      stop
    end if
  end if

  !> Write the case_config namelist variables to an output file in the output directory for inspection or re-use.
  CALL SYSTEM('mkdir -p '//TRIM(output_dir))
  open(unit=1, file=trim(output_dir)//'/'//trim(experiment_name)//'.nml', status='replace', action='write', iostat=ioerror(1))
  write(1,NML=case_config)
  close(1)

  select case(time_scheme)
    case(1)
      n_time_levels = 1
    case(2)
      n_time_levels = 2
    case default
      n_time_levels = 2
  end select

  call scm_state%create(n_columns, n_levels, n_time_levels)

  scm_state%experiment_name = experiment_name
  scm_state%model_name = model_name
  scm_state%output_dir = output_dir
  scm_state%physics_suite_dir = physics_suite_dir
  scm_state%case_data_dir = case_data_dir
  scm_state%vert_coord_data_dir = vert_coord_data_dir
  scm_state%output_file = output_file
  scm_state%case_name = case_name
  scm_state%physics_suite_name = physics_suite
  scm_state%physics_nml = physics_nml
  scm_state%area(:,1) = column_area

  scm_state%n_cols = n_columns
  scm_state%n_levels = n_levels
  scm_state%n_time_levels = n_time_levels
  scm_state%dt = dt
  scm_state%runtime = runtime
  scm_state%time_scheme = time_scheme
  scm_state%init_year = year
  scm_state%init_month = month
  scm_state%init_day = day
  scm_state%init_hour = hour

  scm_state%output_frequency = output_frequency
  scm_state%thermo_forcing_type = thermo_forcing_type
  scm_state%mom_forcing_type = mom_forcing_type
  scm_state%sfc_flux_spec = sfc_flux_spec
  scm_state%sfc_type = sfc_type
  scm_state%sfc_type_real = DBLE(sfc_type)
  scm_state%reference_profile_choice = reference_profile_choice
  scm_state%relax_time = relax_time

!> @}
end subroutine get_config_nml


!> Subroutine to read the netCDF file containing case initialization and forcing. The forcing files (netCDF4) should be located in the
!! "processed_case_input" directory.
subroutine get_case_init(scm_state, scm_input)
  use gmtb_scm_type_defs, only : scm_state_type, scm_input_type
  type(scm_state_type), intent(in) :: scm_state
  type(scm_input_type), target, intent(inout) :: scm_input

  integer               :: input_nlev !< number of levels in the input file
  integer               :: input_ntimes !< number of times represented in the input file

  ! dimension variables
  real(kind=dp), allocatable  :: input_pres(:) !< input file pressure levels (Pa)
  real(kind=dp), allocatable  :: input_time(:) !< input file times (seconds since the beginning of the case)

  !initial profile variables
  real(kind=dp), allocatable  :: input_height(:) !< height of the pressure levels
  real(kind=dp), allocatable  :: input_thetail(:) !< ice-liquid water potential temperature profile (K)
  real(kind=dp), allocatable  :: input_qt(:) !< total water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_ql(:) !< liquid water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_qi(:) !< ice water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_u(:) !< east-west horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_v(:) !< north-south horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_tke(:) !< TKE profile (m^2 s^-2)
  real(kind=dp), allocatable  :: input_ozone(:) !< ozone profile (kg kg^-1)

  !surface time-series variables
  real(kind=dp), allocatable  :: input_lat(:) !< time-series of column latitude
  real(kind=dp), allocatable  :: input_lon(:) !< time-series of column longitude
  real(kind=dp), allocatable  :: input_pres_surf(:) !< time-series of surface pressure (Pa)
  real(kind=dp), allocatable  :: input_T_surf(:) !< time-series of surface temperature (K)
  real(kind=dp), allocatable  :: input_sh_flux_sfc(:) !< time-series of surface sensible heat flux (K m s^-1)
  real(kind=dp), allocatable  :: input_lh_flux_sfc(:) !< time-series of surface latent heat flux (kg kg^-1 m s^-1)

  !2D (time, pressure) variables
  real(kind=dp), allocatable  :: input_w_ls(:,:) !< 2D vertical velocity (m s^-1)
  real(kind=dp), allocatable  :: input_omega(:,:) !< 2D pressure vertical velocity (Pa s^-1)
  real(kind=dp), allocatable  :: input_u_g(:,:) !< 2D geostrophic east-west wind (m s^-1)
  real(kind=dp), allocatable  :: input_v_g(:,:) !< 2D geostrophic north-south wind (m s^-1)
  real(kind=dp), allocatable  :: input_u_nudge(:,:) !< 2D nudging east-west wind (m s^-1)
  real(kind=dp), allocatable  :: input_v_nudge(:,:) !< 2D nudging north-south wind (m s^-1)
  real(kind=dp), allocatable  :: input_T_nudge(:,:) !< 2D nudging abs. temperature (K)
  real(kind=dp), allocatable  :: input_thil_nudge(:,:) !< 2D nudging liq. pot. temperature (K)
  real(kind=dp), allocatable  :: input_qt_nudge(:,:) !< 2D nudging specific humidity (kg kg^-1)
  real(kind=dp), allocatable  :: input_dT_dt_rad(:,:) !< 2D radiative heating rate (K s^-1)
  real(kind=dp), allocatable  :: input_h_advec_thetail(:,:) !< 2D theta_il tendency due to large-scale horizontal advection (K s^-1)
  real(kind=dp), allocatable  :: input_h_advec_qt(:,:) !< 2D q_t tendency due to large-scale horizontal advection (kg kg^-1 s^-1)
  real(kind=dp), allocatable  :: input_v_advec_thetail(:,:) !< 2D theta_il tendency due to large-scale vertical advection (K s^-1)
  real(kind=dp), allocatable  :: input_v_advec_qt(:,:) !< 2D q_t tendency due to large-scale horizontal vertical (kg kg^-1 s^-1)

  CHARACTER(LEN=nf90_max_name)      :: tmpName
  integer                           :: ncid, varID, grp_ncid, allocate_status

  !>  \section get_case_init_alg Algorithm
  !!  @{

  !> - Open the case input file found in the processed_case_input dir corresponding to the experiment name.
  call check(NF90_OPEN(trim(adjustl(scm_state%case_data_dir))//'/'//trim(adjustl(scm_state%case_name))//'.nc',nf90_nowrite,ncid))

  !> - Get the dimensions (global group).

  call check(NF90_INQ_DIMID(ncid,"levels",varID))
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nlev))
  call check(NF90_INQ_DIMID(ncid,"time",varID))
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_ntimes))

  !> - Allocate the dimension variables.
  allocate(input_pres(input_nlev),input_time(input_ntimes), stat=allocate_status)

  !> - Read in the dimension variables.
  call check(NF90_INQ_VARID(ncid,"levels",varID))
  call check(NF90_GET_VAR(ncid,varID,input_pres))
  call check(NF90_INQ_VARID(ncid,"time",varID))
  call check(NF90_GET_VAR(ncid,varID,input_time))

  !> - Read in the initial conditions.

  !>  - Find group ncid for initial group.
  call check(NF90_INQ_GRP_NCID(ncid,"initial",grp_ncid))

  !>  - Allocate the initial profiles.
  allocate(input_height(input_nlev), input_thetail(input_nlev), input_qt(input_nlev), input_ql(input_nlev), input_qi(input_nlev), &
    input_u(input_nlev), input_v(input_nlev), input_tke(input_nlev), input_ozone(input_nlev), stat=allocate_status)

  !>  - Read in the initial profiles. The variable names in all input files are expected to be identical.
  call check(NF90_INQ_VARID(grp_ncid,"height",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_height))
  call check(NF90_INQ_VARID(grp_ncid,"thetail",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_thetail))
  call check(NF90_INQ_VARID(grp_ncid,"qt",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_qt))
  call check(NF90_INQ_VARID(grp_ncid,"ql",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_ql))
  call check(NF90_INQ_VARID(grp_ncid,"qi",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_qi))
  call check(NF90_INQ_VARID(grp_ncid,"u",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_u))
  call check(NF90_INQ_VARID(grp_ncid,"v",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_v))
  call check(NF90_INQ_VARID(grp_ncid,"tke",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_tke))
  call check(NF90_INQ_VARID(grp_ncid,"ozone",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_ozone))

  !> - Read in the forcing data.

  !>  - Find group ncid for forcing group.
  call check(NF90_INQ_GRP_NCID(ncid,"forcing",grp_ncid))

  !>  - (Recall that multidimensional arrays need to be read in with the order of dimensions reversed from the netCDF file).

  !>  - Allocate the time-series and 2D forcing data.
  allocate(input_lat(input_ntimes), input_lon(input_ntimes), input_pres_surf(input_ntimes), input_T_surf(input_ntimes),            &
    input_sh_flux_sfc(input_ntimes), input_lh_flux_sfc(input_ntimes), input_w_ls(input_ntimes, input_nlev), &
    input_omega(input_ntimes, input_nlev), input_u_g(input_ntimes, input_nlev), input_v_g(input_ntimes, input_nlev), &
    input_dT_dt_rad(input_ntimes, input_nlev), input_h_advec_thetail(input_ntimes, input_nlev), &
    input_h_advec_qt(input_ntimes, input_nlev), input_v_advec_thetail(input_ntimes, input_nlev), &
    input_v_advec_qt(input_ntimes, input_nlev), input_u_nudge(input_ntimes, input_nlev), input_v_nudge(input_ntimes, input_nlev),  &
    input_T_nudge(input_ntimes, input_nlev), input_thil_nudge(input_ntimes, input_nlev), input_qt_nudge(input_ntimes, input_nlev), &
    stat=allocate_status)

  !>  - Read in the time-series and 2D forcing data.
  call check(NF90_INQ_VARID(grp_ncid,"lat",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_lat))
  call check(NF90_INQ_VARID(grp_ncid,"lon",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_lon))
  call check(NF90_INQ_VARID(grp_ncid,"p_surf",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_pres_surf))
  call check(NF90_INQ_VARID(grp_ncid,"T_surf",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_T_surf))
  if(scm_state%sfc_flux_spec) then
    call check(NF90_INQ_VARID(grp_ncid,"sh_flux_sfc",varID))
    call check(NF90_GET_VAR(grp_ncid,varID,input_sh_flux_sfc))
    call check(NF90_INQ_VARID(grp_ncid,"lh_flux_sfc",varID))
    call check(NF90_GET_VAR(grp_ncid,varID,input_lh_flux_sfc))
  else
    input_sh_flux_sfc = 0.0
    input_lh_flux_sfc = 0.0
  end if
  call check(NF90_INQ_VARID(grp_ncid,"w_ls",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_w_ls))
  call check(NF90_INQ_VARID(grp_ncid,"omega",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_omega))
  call check(NF90_INQ_VARID(grp_ncid,"u_g",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_u_g))
  call check(NF90_INQ_VARID(grp_ncid,"v_g",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_v_g))
  call check(NF90_INQ_VARID(grp_ncid,"u_nudge",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_u_nudge))
  call check(NF90_INQ_VARID(grp_ncid,"v_nudge",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_v_nudge))
  call check(NF90_INQ_VARID(grp_ncid,"T_nudge",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_T_nudge))
  call check(NF90_INQ_VARID(grp_ncid,"thil_nudge",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_thil_nudge))
  call check(NF90_INQ_VARID(grp_ncid,"qt_nudge",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_qt_nudge))
  call check(NF90_INQ_VARID(grp_ncid,"dT_dt_rad",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_dT_dt_rad))
  call check(NF90_INQ_VARID(grp_ncid,"h_advec_thetail",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_h_advec_thetail))
  call check(NF90_INQ_VARID(grp_ncid,"h_advec_qt",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_h_advec_qt))
  call check(NF90_INQ_VARID(grp_ncid,"v_advec_thetail",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_v_advec_thetail))
  call check(NF90_INQ_VARID(grp_ncid,"v_advec_qt",varID))
  call check(NF90_GET_VAR(grp_ncid,varID,input_v_advec_qt))
  call check(NF90_CLOSE(NCID=ncid))

  call scm_input%create(input_ntimes, input_nlev)

  scm_input%input_nlev = input_nlev
  scm_input%input_ntimes = input_ntimes

  scm_input%input_pres = input_pres
  scm_input%input_time = input_time
  scm_input%input_height = input_height
  scm_input%input_thetail = input_thetail
  scm_input%input_qt = input_qt
  scm_input%input_ql = input_ql
  scm_input%input_qi = input_qi
  scm_input%input_u = input_u
  scm_input%input_v = input_v
  scm_input%input_tke = input_tke
  scm_input%input_ozone = input_ozone
  scm_input%input_lat = input_lat
  scm_input%input_lon = input_lon
  scm_input%input_pres_surf = input_pres_surf
  scm_input%input_T_surf = input_T_surf
  scm_input%input_sh_flux_sfc = input_sh_flux_sfc
  scm_input%input_lh_flux_sfc = input_lh_flux_sfc
  scm_input%input_w_ls = input_w_ls
  scm_input%input_omega = input_omega
  scm_input%input_u_g = input_u_g
  scm_input%input_v_g = input_v_g
  scm_input%input_dT_dt_rad = input_dT_dt_rad
  scm_input%input_h_advec_thetail = input_h_advec_thetail
  scm_input%input_h_advec_qt = input_h_advec_qt
  scm_input%input_v_advec_thetail = input_v_advec_thetail
  scm_input%input_v_advec_qt = input_v_advec_qt
  scm_input%input_u_nudge = input_u_nudge
  scm_input%input_v_nudge = input_v_nudge
  scm_input%input_T_nudge = input_T_nudge
  scm_input%input_thil_nudge = input_thil_nudge
  scm_input%input_qt_nudge = input_qt_nudge

!> @}
end subroutine get_case_init

!> Subroutine to get reference profile to use above the case data (temporarily hard-coded profile)
subroutine get_reference_profile(scm_state, scm_reference)
  use gmtb_scm_type_defs, only : scm_state_type, scm_reference_type

  type(scm_state_type), target, intent(in) :: scm_state
  type(scm_reference_type), target, intent(inout) :: scm_reference

  integer  :: nlev !< number of pressure levels in the reference profile
  real(kind=dp), allocatable :: pres(:)  !< reference profile pressure (Pa)
  real(kind=dp), allocatable :: T(:) !< reference profile temperature (K)
  real(kind=dp), allocatable :: qv(:) !< reference profile specific humidity (kg kg^-1)
  real(kind=dp), allocatable :: ozone(:) !< reference profile ozone concentration (kg kg^-1)

  integer :: reference_file_choice
  integer :: i, ioerror
  character(len=120)                 :: line
  real :: dummy

  integer                           :: ncid, varID, grp_ncid, allocate_status
  CHARACTER(LEN=nf90_max_name)      :: tmpName

  select case (scm_state%reference_profile_choice)
    case (1)
      open(unit=1, file=trim(adjustl(scm_state%case_data_dir))//'/'//'McCProfiles.dat', status='old', action='read', iostat=ioerror)
      if(ioerror /= 0) then
        write(*,*) 'There was an error opening the file McCprofiles.dat in the processed_case_input directory. &
          Error code = ',ioerror
        stop
      endif

      ! find number of records
      read(1, '(a)', iostat=ioerror) line !first line is header
      nlev = 0
      do
        read(1, '(a)', iostat=ioerror) line
        if (ioerror /= 0) exit
        nlev = nlev + 1
      end do

      allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev))

      rewind(1)
      read(1, '(a)', iostat=ioerror) line !first line is header
      do i=1, nlev
        read(1,'(5ES16.4)', iostat=ioerror) dummy, pres(i), T(i), qv(i), ozone(i)
      END DO
      close(1)
    case (2)
      call check(NF90_OPEN(trim(adjustl(scm_state%case_data_dir))//'/'//'mid_lat_summer_std.nc',nf90_nowrite,ncid))

      call check(NF90_INQ_DIMID(ncid,"height",varID))
      call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, nlev))

      !> - Allocate the dimension variables.
      allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev), stat=allocate_status)

      call check(NF90_INQ_VARID(ncid,"pressure",varID))
      call check(NF90_GET_VAR(ncid,varID,pres))
      call check(NF90_INQ_VARID(ncid,"temperature",varID))
      call check(NF90_GET_VAR(ncid,varID,T))
      call check(NF90_INQ_VARID(ncid,"q_v",varID))
      call check(NF90_GET_VAR(ncid,varID,qv))
      call check(NF90_INQ_VARID(ncid,"o3",varID))
      call check(NF90_GET_VAR(ncid,varID,ozone))

      call check(NF90_CLOSE(NCID=ncid))

  end select

  call scm_reference%create(nlev)

  scm_reference%ref_nlev = nlev
  scm_reference%ref_pres = pres
  scm_reference%ref_T = T
  scm_reference%ref_qv = qv
  scm_reference%ref_ozone = ozone


end subroutine get_reference_profile

!> Subroutine to get reference profile to use above the case data (temporarily hard-coded profile)
subroutine get_reference_profile_old(nlev, pres, T, qv, ozone)
  integer, intent(out)  :: nlev !< number of pressure levels in the reference profile
  real(kind=dp), allocatable, intent(out) :: pres(:)  !< reference profile pressure (Pa)
  real(kind=dp), allocatable, intent(out) :: T(:) !< reference profile temperature (K)
  real(kind=dp), allocatable, intent(out) :: qv(:) !< reference profile specific humidity (kg kg^-1)
  real(kind=dp), allocatable, intent(out) :: ozone(:) !< reference profile ozone concentration (kg kg^-1)

  !> \todo write a more sophisticated reference profile subroutine (can choose between reference profiles)

  !> - For the prototype, the 'McClatchey' sounding used in Jennifer Fletcher's code is hardcoded.
  nlev = 20

  allocate(pres(nlev), T(nlev), qv(nlev), ozone(nlev))

  pres = (/ 1030.0, 902.0, 802.0, 710.0, 628.0, 554.0, 487.0, 426.0, 372.0, 281.0, 209.0, 130.0, 59.5, 27.7, 13.2, 6.52, 3.33, &
    0.951, 0.0671, 0.000300 /)
  pres = pres*100.0

  T = (/ 294., 290., 285., 279., 273., 267., 261., 255., 248., 235., 222., 216., 218., 224., 234., 245., 258., 276., 218., 210. /)

  qv = (/ 11.75, 8.611, 6.047, 3.877, 2.363, 1.387, 0.9388, 0.6364, 0.4019, 0.1546, 0.01976, 0.004002, 0.003999, 0.004011, &
    0.004002, 0.004004, 0.003994, 0.003995, 0.003996, 0.004000 /)
  qv = qv*1.0E-3

  ozone = (/ 6., 6., 6., 6.2, 6.4, 6.6, 6.9, 7.5, 7.9, 9., 12., 19., 34., 30., 20., 9.2, 4.1, 0.43, 0.0086, 0.0000043 /)
  ozone = ozone*1.0E-5

end subroutine get_reference_profile_old

!> Generic subroutine to check for netCDF I/O errors
subroutine check(status)
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "stopped"
  end if
end subroutine check
!> @}
!> @}
end module gmtb_scm_input
