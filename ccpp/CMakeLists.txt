#------------------------------------------------------------------------------
#
# Create list of SCHEME_FILES, HOST_FILES, and SUITE_FILES
# Paths should be relative to CMAKE_SOURCE_DIR (this file's directory)
#
#------------------------------------------------------------------------------

set(SCHEME_FILES_DIR "${CMAKE_SOURCE_DIR}/ccpp/physics/physics")
set(SCHEME_METADATA_FILES
  # time-vary
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/module_ccpp_suite_simulator.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_ccpp_suite_sim_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_time_vary_pre.scm.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_setup.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rad_time_vary.scm.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_phys_time_vary.scm.meta"
  # radiation
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radsw_param.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_rrtmg_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_surface.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rad_sw_pre.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radsw_main.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rrtmg_sw_post.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radlw_main.meta"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/rrtmg_lw_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_radiation_post.meta"
  # physics
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_reset.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_stateout_update.meta"
  "${SCHEME_FILES_DIR}/tools/get_prs_fv3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_1.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_pre.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/dcyc2t3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_inter.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_2.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part1.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diff.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/scm_sfc_flux_spec.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/Noahmp/lnd_iau_mod.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst_pre.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_nst_post.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/Land/Noahmp/noahmpdrv.meta"
  "${SCHEME_FILES_DIR}/SFC_Models/SeaIce/CICE/sfc_sice.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part2.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_post.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diag.meta"
  "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diag_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_pre.meta"
  "${SCHEME_FILES_DIR}/PBL/SATMEDMF/satmedmfvdifq.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_PBL_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_pre.meta"
  "${SCHEME_FILES_DIR}/GWD/unified_ugwp.meta"
  "${SCHEME_FILES_DIR}/GWD/unified_ugwp_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_GWD_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_photochemistry.meta"
  "${SCHEME_FILES_DIR}/tools/get_phi_fv3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_3.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_pre.meta"
  "${SCHEME_FILES_DIR}/CONV/SAMF/samfdeepcnv.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_DCNV_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_pre.meta"
  "${SCHEME_FILES_DIR}/CONV/SAMF/samfshalcnv.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_SCNV_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_suite_interstitial_4.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/cnvc90.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_pre.meta"
  "${SCHEME_FILES_DIR}/MP/Thompson/mp_thompson_pre.meta"
  "${SCHEME_FILES_DIR}/MP/Thompson/mp_thompson.meta"
  "${SCHEME_FILES_DIR}/MP/Thompson/mp_thompson_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_MP_generic_post.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/maximum_hourly_diagnostics.meta"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_physics_post.meta"
)

## DH* TEST
## Run ccpp_validator
#ccpp_validator(VERBOSITY ${CCPP_VERBOSITY}
#               SOURCE_FILES "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diff.f"
#               METADATA_FILES "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diff.meta")
## *DH

# Create list of scheme source files from metadata files
ccpp_source_files(SCHEME_FORTRAN_FILES ${SCHEME_METADATA_FILES})

set(HOST_FILES_DIR "${CMAKE_SOURCE_DIR}/scm/src")
set(HOST_METADATA_FILES
  "${HOST_FILES_DIR}/scm_kinds.meta"
  "${HOST_FILES_DIR}/scm_physical_constants.meta"
  "${HOST_FILES_DIR}/scm_type_defs.meta"
  "${HOST_FILES_DIR}/CCPP_typedefs.meta"
  "${HOST_FILES_DIR}/GFS_typedefs.meta"
  "${HOST_FILES_DIR}/scm_time_integration.meta"
)

ccpp_source_files(HOST_FORTRAN_FILES ${HOST_METADATA_FILES})
message(INFO "${HOST_METADATA_FILES} --> ${HOST_FORTRAN_FILES}")

set(SUITE_FILES_DIR "${CMAKE_SOURCE_DIR}/ccpp/suites")
set(SUITE_FILES
  "${SUITE_FILES_DIR}/suite_SCM_GFS_v17_p8.xml"
  "${SUITE_FILES_DIR}/suite_SCM_GFS_v17_p8_ps.xml"
)

set(HOST "ccpp-scm")

set(CCPP_VERBOSITY "2" CACHE STRING "Verbosity level of output (default: 0)")

# By default, generated caps go in ccpp subdir
set(OUTPUT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/ccpp")

# Run ccpp_validator
ccpp_validator(VERBOSITY ${CCPP_VERBOSITY}
               SOURCE_FILES ${SCHEME_FORTRAN_FILES}
               METADATA_FILES ${SCHEME_METADATA_FILES})
ccpp_validator(VERBOSITY ${CCPP_VERBOSITY}
               SOURCE_FILES ${HOST_FORTRAN_FILES}
               METADATA_FILES ${HOST_METADATA_FILES})

# Run ccpp_capgen_ng
ccpp_capgen(VERBOSITY ${CCPP_VERBOSITY}
            HOSTFILES ${HOST_METADATA_FILES}
            SCHEMEFILES ${SCHEME_METADATA_FILES}
            SUITES ${SUITE_FILES}
            HOST_NAME ${HOST}
            OUTPUT_ROOT "${OUTPUT_ROOT}")

# Retrieve the list of Fortran files required for test host from datatable.xml;
# this includes capgen-generated files and dependencies inferred from metadata
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--dependencies")
set(CAPGEN_DEPENDENCIES ${CCPP_FILES})
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--scheme-files")
set(SCHEME_FORTRAN_FILES ${CCPP_FILES})
ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
              REPORT_NAME "--capgen-files")
set(CAPGEN_FILES ${CCPP_FILES})

message(STATUS "List of capgen dependencies: ${CAPGEN_DEPENDENCIES}")
message(STATUS "List of filtered scheme files: ${SCHEME_FORTRAN_FILES}")
message(STATUS "List of capgen-generated files: ${CAPGEN_FILES}")

# Temporary: manually added dependencies because of working on a subset of schemes only
set(EXTRA_FILES
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radlw_datatb.f"
  "${SCHEME_FILES_DIR}/Radiation/RRTMG/radsw_datatb.f"
  "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/module_ccpp_suite_simulator.F90"
)
# End temporary

add_library(scm-ccpp STATIC
  ${EXTRA_FILES}
  ${CAPGEN_DEPENDENCIES}
  ${SCHEME_FORTRAN_FILES}
  ${HOST_FORTRAN_FILES}
  ${CAPGEN_FILES}
)

#set_target_properties(scm-ccpp PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(scm-ccpp MPI::MPI_Fortran)
if(OPENMP)
  target_link_libraries(scm-ccpp OpenMP::OpenMP_Fortran)
endif()
target_link_libraries(scm-ccpp NetCDF::NetCDF_Fortran)
if(TARGET bacio::bacio)
  set(BACIO_TARGET bacio::bacio)
elseif(TARGET bacio::bacio_4)
  set(BACIO_TARGET bacio::bacio_4)
else()
  message(FATAL_ERROR
    "bacio found but expected target not exported (bacio::bacio or bacio::bacio_4)")
endif()
target_link_libraries(scm-ccpp ${BACIO_TARGET})
if(ip_FOUND)
  target_link_libraries(scm-ccpp ip::ip_d)
else()
  target_link_libraries(scm-ccpp sp::sp_d)
endif()
target_link_libraries(scm-ccpp w3emc::w3emc_d)

set_target_properties(scm-ccpp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
target_include_directories(scm-ccpp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>)

#message(FATAL_ERROR "xxx")
#
#add_executable(test_chunked_data.x
#  ${CAPGEN_DEPENDENCIES}
#  ${SCHEME_FORTRAN_FILES}
#  ${HOST_FORTRAN_FILES}
#  ${CAPGEN_FILES}
#)
#target_link_libraries(test_chunked_data.x PRIVATE MPI::MPI_Fortran)
#if(OPENMP)
#  target_link_libraries(test_chunked_data.x PRIVATE OpenMP::OpenMP_Fortran)
#endif()
#set_target_properties(test_chunked_data.x PROPERTIES LINKER_LANGUAGE Fortran)
#
## Add executable to be called with ctest
#add_test(NAME test_chunked_data
#         COMMAND test_chunked_data.x)
#