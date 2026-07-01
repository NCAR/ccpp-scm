#
# CCPP SCM source code
#

SET(scm_source_files scm.F90
            scm_input.F90
            #scm_utils.F90
            scm_vgrid.F90
            scm_setup.F90
            #scm_forcing.F90
            #scm_time_integration.F90
            scm_output.F90
)

add_executable(scm ${scm_source_files})

target_link_libraries(scm scm-ccpp)

if(OPENMP)
  target_link_libraries(scm OpenMP::OpenMP_Fortran)
endif()
target_link_libraries(scm NetCDF::NetCDF_Fortran)

add_custom_command(
        TARGET scm POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/run_scm.py
                ${CMAKE_BINARY_DIR}/run_scm.py)
