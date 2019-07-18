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
  real(kind=dp)        :: dt !< time step in seconds
  real(kind=dp)        :: runtime !< total runtime in seconds
  real(kind=dp)        :: output_frequency !< freqency of output writing in seconds
  integer              :: n_levels !< number of model levels (currently only 64 supported)
  integer              :: n_soil   !< number of model levels (currently only 4 supported)
  integer              :: n_columns !< number of columns to use
  integer              :: n_time_levels
  integer              :: time_scheme !< 1 => forward Euler, 2 => filtered leapfrog
  character(len=80)    :: output_dir !< name of the output directory
  character(len=80)    :: output_file !< name of the output file (without the file extension)
  character(len=80)    :: case_data_dir !< path to the directory containing case initialization and forcing data
  character(len=80)    :: vert_coord_data_dir !< path to the directory containing vertical coordinate data
  integer              :: thermo_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer              :: mom_forcing_type !< 1: "revealed forcing", 2: "horizontal advective forcing", 3: "relaxation forcing"
  integer              :: C_RES            !< reference "C" resoltiion of FV3 grid (needed for GWD and mountain blocking)
  real(kind=dp)        :: relax_time !< relaxation time scale (s)
  logical              :: sfc_flux_spec !< flag for using specified surface fluxes instead of calling a surface scheme
  real(kind=dp)        :: sfc_roughness_length_cm !< surface roughness length used for calculating surface layer parameters from specified fluxes
  integer              :: sfc_type !< 0: sea surface, 1: land surface, 2: sea-ice surface
  logical              :: model_ics !<  true means have land info too
  integer              :: reference_profile_choice !< 1: McClatchey profile, 2: mid-latitude summer standard atmosphere
  integer              :: year, month, day, hour
  real(kind=dp)        :: column_area

  character(len=80), allocatable  :: physics_suite(:) !< name of the physics suite name (currently only GFS_operational supported)
  character(len=64), allocatable   :: physics_nml(:)

  integer                          :: ioerror

  CHARACTER(LEN=*), parameter :: experiment_namelist = 'input_experiment.nml'

  NAMELIST /case_config/ model_name, n_columns, case_name, dt, time_scheme, runtime, output_frequency, &
    n_levels, output_dir, output_file, case_data_dir, vert_coord_data_dir, thermo_forcing_type, model_ics,C_RES,mom_forcing_type, relax_time, &
    sfc_type, sfc_flux_spec, sfc_roughness_length_cm, reference_profile_choice, year, month, day, hour, column_area
    
  NAMELIST /physics_config/ physics_suite, physics_nml

  !>  \section get_config_alg Algorithm
  !!  @{

  !> Define default values for experiment configuration (to be overridden by external namelist file or command line arguments)
  model_name = 'GFS'
  n_columns = 1
  case_name = 'twpice'
  dt = 600.0
  time_scheme = 2
  runtime = 2138400.0
  output_frequency = 600.0
  n_levels = 64
  n_soil   = 4
  output_dir = 'output'
  output_file = 'output'
  case_data_dir = '../data/processed_case_input'
  vert_coord_data_dir = '../data/vert_coord_data'
  thermo_forcing_type = 2
  mom_forcing_type = 3
  C_RES = 384
  relax_time = 7200.0
  sfc_flux_spec = .false.
  sfc_roughness_length_cm = 1.0
  sfc_type = 0
  model_ics = .false.
  reference_profile_choice = 1
  year = 2006
  month = 1
  day = 19
  hour = 3

  open(unit=10, file=experiment_namelist, status='old', action='read', iostat=ioerror)
  if(ioerror /= 0) then
    write(*,'(a,i0)') 'There was an error opening the file ' // experiment_namelist // &
                      '; error code = ', ioerror
    STOP
  else
    read(10, NML=case_config, iostat=ioerror)
  end if

  if(ioerror /= 0) then
    write(*,'(a,i0)') 'There was an error reading the namelist case_config in the file '&
                      // experiment_namelist // '; error code = ',ioerror
    STOP
  end if

  !The current implementation of GFS physics does not support more than one column, since radiation sub schemes use
  !internal module variables. This means that one cannot specify different ways to treat O3, CO2 etc., and also that
  !the code crashes in GFS_initialize_scm_run and later in radiation_gases.f, because it tries to allocate module
  !variables that are already allocated. For now, throw an error and abort.
  if (n_columns>1) then
    write(*,'(a)') 'The current implementation does not allow to run more than one column at a time.'
    STOP
  end if

  !Using n_columns, allocate memory for the physics suite names and number of fields needed by each. If there are more physics suites
  !than n_columns, notify the user and stop the program. If there are less physics suites than columns, notify the user and attempt to
  !continue (getting permission from user), filling in the unspecified suites as the same as the last specified suite.
  allocate(physics_suite(n_columns), physics_nml(n_columns))

  read(10, NML=physics_config, iostat=ioerror)
  close(10)

  select case(time_scheme)
    case(1)
      n_time_levels = 1
    case(2)
      n_time_levels = 2
    case default
      n_time_levels = 2
  end select

  call scm_state%create(n_columns, n_levels, n_soil, n_time_levels)
  
  scm_state%experiment_name = experiment_name
  scm_state%model_name = model_name
  scm_state%output_dir = output_dir
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
  scm_state%C_RES            = C_RES            
  scm_state%sfc_flux_spec = sfc_flux_spec
  scm_state%sfc_roughness_length_cm(:) = sfc_roughness_length_cm
  scm_state%sfc_type = sfc_type
  scm_state%sfc_type_real = DBLE(sfc_type)
  scm_state%model_ics = model_ics
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
  integer               :: input_nsoil !< number of soil levels in the input file
  integer               :: input_ntimes !< number of times represented in the input file

  ! dimension variables
  real(kind=dp), allocatable  :: input_pres(:) !< input file pressure levels (Pa)
  real(kind=dp), allocatable  :: input_time(:) !< input file times (seconds since the beginning of the case)

  !initial profile variables
  real(kind=dp), allocatable  :: input_thetail(:) !< ice-liquid water potential temperature profile (K)
  real(kind=dp), allocatable  :: input_temp(:) !< temperature profile (K)
  real(kind=dp), allocatable  :: input_qt(:) !< total water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_ql(:) !< liquid water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_qi(:) !< ice water specific humidity profile (kg kg^-1)
  real(kind=dp), allocatable  :: input_u(:) !< east-west horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_v(:) !< north-south horizontal wind profile (m s^-1)
  real(kind=dp), allocatable  :: input_tke(:) !< TKE profile (m^2 s^-2)
  real(kind=dp), allocatable  :: input_ozone(:) !< ozone profile (kg kg^-1)
! additional land info
  real(kind=dp), allocatable  :: input_stc(:) !< soil temperature
  real(kind=dp), allocatable  :: input_smc(:) !< soil moisture    
  real(kind=dp), allocatable  :: input_slc(:) !< soil liquid content
  real(kind=dp), allocatable  :: input_pres_i(:) !< interface pressures
  real(kind=dp), allocatable  :: input_pres_l(:) !< layer pressures
  integer                     :: input_vegsrc !< vegetation source
  integer                     :: input_vegtyp !< vegetation type
  integer                     :: input_soiltyp!< soil type
  integer                     :: input_slopetype !< slope type
  real(kind=dp)               :: input_vegfrac  !< vegetation fraction
  real(kind=dp)               :: input_shdmin  !< minimun vegetation fraction
  real(kind=dp)               :: input_shdmax  !< maximun vegetation fraction
  real(kind=dp)               :: input_zorl    !< surfce roughness length [cm]
  real(kind=dp)               :: input_slmsk   !< sea land ice mask [0,1,2]
  real(kind=dp)               :: input_canopy  !< amount of water stored in canopy
  real(kind=dp)               :: input_hice    !< ice thickness
  real(kind=dp)               :: input_fice    !< ice fraction
  real(kind=dp)               :: input_tisfc   !< ice surface temperature
  real(kind=dp)               :: input_snwdph  !< snow depth
  real(kind=dp)               :: input_snoalb  !< snow albedo
  real(kind=dp)               :: input_sncovr  !< snow cover
  real(kind=dp)               :: input_area    !< surfce area [m^2]
  real(kind=dp)               :: input_tg3     !< maximun vegetation fraction
  real(kind=dp)               :: input_uustar  !< maximun vegetation fraction
  real(kind=dp)               :: input_alvsf !< uv+visible black sky albedo (z=60 degree)"
  real(kind=dp)               :: input_alnsf !< near IR black sky albedo (z=60 degree)"
  real(kind=dp)               :: input_alvwf !< uv+visible white sky albedo"
  real(kind=dp)               :: input_alnwf !< near IR white sky albedo"
  real(kind=dp)               :: input_stddev !<  surface roughness
  real(kind=dp)               :: input_convexity !<  surface roughness
  real(kind=dp)               :: input_ol1 !<  surface roughness
  real(kind=dp)               :: input_ol2 !<  surface roughness
  real(kind=dp)               :: input_ol3 !<  surface roughness
  real(kind=dp)               :: input_ol4 !<  surface roughness
  real(kind=dp)               :: input_oa1 !<  surface roughness
  real(kind=dp)               :: input_oa2 !<  surface roughness
  real(kind=dp)               :: input_oa3 !<  surface roughness
  real(kind=dp)               :: input_oa4 !<  surface roughness
  real(kind=dp)               :: input_sigma !<  surface roughness
  real(kind=dp)               :: input_theta !<  surface roughness
  real(kind=dp)               :: input_gamma !<  surface roughness
  real(kind=dp)               :: input_elvmax!<  surface roughness
  real(kind=dp)               :: input_facsf !< near IR white sky albedo"
  real(kind=dp)               :: input_facwf !< near IR white sky albedo"

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
  integer                           :: ncid, varID, grp_ncid, allocate_status,ierr

  !>  \section get_case_init_alg Algorithm
  !!  @{

  !> - Open the case input file found in the processed_case_input dir corresponding to the experiment name.
  call check(NF90_OPEN(trim(adjustl(scm_state%case_data_dir))//'/'//trim(adjustl(scm_state%case_name))//'.nc',nf90_nowrite,ncid))

  !> - Get the dimensions (global group).

  call check(NF90_INQ_DIMID(ncid,"levels",varID))
  call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nlev))
  if (scm_state%model_ics) then
    call check(NF90_INQ_DIMID(ncid,"nsoil",varID))
    call check(NF90_INQUIRE_DIMENSION(ncid, varID, tmpName, input_nsoil))
  endif
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
  allocate(input_thetail(input_nlev), input_qt(input_nlev), input_ql(input_nlev), input_qi(input_nlev), &
    input_u(input_nlev), input_v(input_nlev), input_tke(input_nlev), input_ozone(input_nlev), stat=allocate_status)
   if (scm_state%model_ics) then
     allocate(input_stc(input_nsoil), input_temp(input_nlev),input_smc(input_nsoil), input_slc(input_nsoil), &
              input_pres_i(input_nlev+1),input_pres_l(input_nlev), stat=allocate_status)
     input_pres_i(:) = -999.9
     input_pres_l(:) = -999.9
  endif

  !>  - Read in the initial profiles. The variable names in all input files are expected to be identical.
  if (.NOT. scm_state%model_ics) then
     call check(NF90_INQ_VARID(grp_ncid,"thetail",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_thetail))
  endif
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
  if (scm_state%model_ics) then
     call check(NF90_INQ_VARID(grp_ncid,"temp",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_temp))
     call check(NF90_INQ_VARID(grp_ncid,"stc",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_stc))
     call check(NF90_INQ_VARID(grp_ncid,"smc",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_smc))
     call check(NF90_INQ_VARID(grp_ncid,"slc",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_slc))
     ierr = NF90_INQ_VARID(grp_ncid,"pres_i",varID)
     if (ierr.EQ.0) then ! input file should have pres_i and pres_l
        call check(NF90_GET_VAR(grp_ncid,varID,input_pres_i))
        call check(NF90_INQ_VARID(grp_ncid,"pres_l",varID))
        call check(NF90_GET_VAR(grp_ncid,varID,input_pres_l))
     endif
     call check(NF90_INQ_GRP_NCID(ncid,"scalars",grp_ncid))
     call check(NF90_INQ_VARID(grp_ncid,"vegsrc",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_vegsrc))
     call check(NF90_INQ_VARID(grp_ncid,"vegtyp",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_vegtyp))
     call check(NF90_INQ_VARID(grp_ncid,"soiltyp",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_soiltyp))
     call check(NF90_INQ_VARID(grp_ncid,"slopetyp",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_slopetype))
     call check(NF90_INQ_VARID(grp_ncid,"vegfrac",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_vegfrac))
     call check(NF90_INQ_VARID(grp_ncid,"shdmin",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_shdmin))
     call check(NF90_INQ_VARID(grp_ncid,"shdmax",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_shdmax))
     call check(NF90_INQ_VARID(grp_ncid,"zorl",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_zorl))
     call check(NF90_INQ_VARID(grp_ncid,"slmsk",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_slmsk))
     call check(NF90_INQ_VARID(grp_ncid,"canopy",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_canopy))
     call check(NF90_INQ_VARID(grp_ncid,"hice",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_hice))
     call check(NF90_INQ_VARID(grp_ncid,"fice",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_fice))
     call check(NF90_INQ_VARID(grp_ncid,"tisfc",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_tisfc))
     call check(NF90_INQ_VARID(grp_ncid,"snwdph",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_snwdph))
     call check(NF90_INQ_VARID(grp_ncid,"snoalb",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_snoalb))
     call check(NF90_INQ_VARID(grp_ncid,"sncovr",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_sncovr))
     call check(NF90_INQ_VARID(grp_ncid,"area",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_area))
     call check(NF90_INQ_VARID(grp_ncid,"tg3",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_tg3))
     call check(NF90_INQ_VARID(grp_ncid,"uustar",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_uustar))
     call check(NF90_INQ_VARID(grp_ncid,"alvsf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_alvsf))
     call check(NF90_INQ_VARID(grp_ncid,"alnsf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_alnsf))
     call check(NF90_INQ_VARID(grp_ncid,"alvwf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_alvwf))
     call check(NF90_INQ_VARID(grp_ncid,"alnwf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_alnwf))
     call check(NF90_INQ_VARID(grp_ncid,"stddev",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_stddev))
     call check(NF90_INQ_VARID(grp_ncid,"convexity",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_convexity))
     call check(NF90_INQ_VARID(grp_ncid,"oa1",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_oa1))
     call check(NF90_INQ_VARID(grp_ncid,"oa2",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_oa2))
     call check(NF90_INQ_VARID(grp_ncid,"oa3",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_oa3))
     call check(NF90_INQ_VARID(grp_ncid,"oa4",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_oa4))
     call check(NF90_INQ_VARID(grp_ncid,"ol1",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_ol1))
     call check(NF90_INQ_VARID(grp_ncid,"ol2",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_ol2))
     call check(NF90_INQ_VARID(grp_ncid,"ol3",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_ol3))
     call check(NF90_INQ_VARID(grp_ncid,"ol4",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_ol4))
     call check(NF90_INQ_VARID(grp_ncid,"sigma",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_sigma))
     call check(NF90_INQ_VARID(grp_ncid,"gamma",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_gamma))
     call check(NF90_INQ_VARID(grp_ncid,"theta",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_theta))
     call check(NF90_INQ_VARID(grp_ncid,"elvmax",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_elvmax))
     call check(NF90_INQ_VARID(grp_ncid,"facsf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_facsf))
     call check(NF90_INQ_VARID(grp_ncid,"facwf",varID))
     call check(NF90_GET_VAR(grp_ncid,varID,input_facwf))
  endif

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
  if (scm_state%model_ics) then
     call scm_input%create_modelics(input_nsoil,input_nlev)
  endif

  scm_input%input_nlev = input_nlev
  scm_input%input_ntimes = input_ntimes

  scm_input%input_pres = input_pres
  scm_input%input_time = input_time
  if (scm_state%model_ics) then
     scm_input%input_temp = input_temp
  else
     scm_input%input_thetail = input_thetail
  endif
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
  if (scm_state%model_ics) then
     scm_input%input_stc   = input_stc  
     scm_input%input_smc   = input_smc  
     scm_input%input_slc   = input_slc  
     scm_input%input_vegsrc   = input_vegsrc
     scm_input%input_vegtyp   = input_vegtyp
     scm_input%input_soiltyp  = input_soiltyp
     scm_input%input_slopetype = input_slopetype
     scm_input%input_vegfrac  = input_vegfrac 
     scm_input%input_shdmin   = input_shdmin  
     scm_input%input_shdmax   = input_shdmax  
     scm_input%input_zorl     = input_zorl    
     scm_input%input_slmsk    = input_slmsk   
     scm_input%input_canopy   = input_canopy  
     scm_input%input_hice     = input_hice    
     scm_input%input_fice     = input_fice    
     scm_input%input_tisfc    = input_tisfc   
     scm_input%input_snwdph   = input_snwdph  
     scm_input%input_snoalb   = input_snoalb  
     scm_input%input_sncovr   = input_sncovr  
     scm_input%input_area     = input_area    
     scm_input%input_tg3      = input_tg3     
     scm_input%input_uustar   = input_uustar  
     scm_input%input_alvsf    = input_alvsf   
     scm_input%input_alnsf    = input_alnsf   
     scm_input%input_alvwf    = input_alvwf   
     scm_input%input_alnwf    = input_alnwf   
     scm_input%input_convexity= input_convexity
     scm_input%input_stddev   = input_stddev  
     scm_input%input_oa1      = input_oa1     
     scm_input%input_oa2      = input_oa2     
     scm_input%input_oa3      = input_oa3     
     scm_input%input_oa4      = input_oa4     
     scm_input%input_ol1      = input_ol1     
     scm_input%input_ol2      = input_ol2    
     scm_input%input_ol3      = input_ol3     
     scm_input%input_ol4      = input_ol4     
     scm_input%input_sigma    = input_sigma   
     scm_input%input_theta    = input_theta   
     scm_input%input_gamma    = input_gamma   
     scm_input%input_elvmax   = input_elvmax  
     scm_input%input_facsf    = input_facsf   
     scm_input%input_facwf    = input_facwf   
     scm_input%input_pres_i   = input_pres_i  
     scm_input%input_pres_l   = input_pres_l  
  endif

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

  integer :: i, ioerror
  character(len=120)                 :: line
  real :: dummy

  integer                           :: ncid, varID, allocate_status
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
