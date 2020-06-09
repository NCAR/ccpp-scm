#!/usr/bin/env python

# CCPP prebuild config for GMTB Single Column Model (SCM)


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "SCM"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics.
VARIABLE_DEFINITION_FILES = [
    'ccpp/physics/physics/machine.F',
    'ccpp/physics/physics/radsw_param.f',
    'ccpp/physics/physics/radlw_param.f',
    'scm/src/GFS_typedefs.F90',
    'scm/src/gmtb_scm_kinds.F90',
    'scm/src/gmtb_scm_type_defs.F90',
    'scm/src/gmtb_scm_physical_constants.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_optical_props.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_source_functions.F90'
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_types' : '',
        'ccpp_t' : 'cdata',
        },
    'machine' : {
        'machine' : '',
        },
    'module_radlw_parameters' : {
        'module_radsw_parameters' : '',
        },
    'module_radlw_parameters' : {
        'module_radlw_parameters' : '',
        },
    'GFS_typedefs' : {
        'GFS_diag_type' : 'physics%Diag(cdata%blk_no)',
        'GFS_control_type' : 'physics%Model(cdata%blk_no)',
        'GFS_cldprop_type' : 'physics%Cldprop(cdata%blk_no)',
        'GFS_tbd_type' : 'physics%Tbd(cdata%blk_no)',
        'GFS_sfcprop_type' : 'physics%Sfcprop(cdata%blk_no)',
        'GFS_coupling_type' : 'physics%Coupling(cdata%blk_no)',
        'GFS_interstitial_type' : 'physics%Interstitial(cdata%blk_no)',
        'GFS_statein_type' : 'physics%Statein(cdata%blk_no)',
        'GFS_radtend_type' : 'physics%Radtend(cdata%blk_no)',
        'GFS_grid_type' : 'physics%Grid(cdata%blk_no)',
        'GFS_stateout_type' : 'physics%Stateout(cdata%blk_no)',
        'GFS_typedefs' : '',
        },
    'gmtb_scm_physical_constants' : {
        'gmtb_scm_physical_constants' : '',
        },
    'gmtb_scm_type_defs' : {
        'gmtb_scm_type_defs' : '',
        'physics_type' : 'physics',
        },
   'mo_gas_concentrations' : {
        'ty_gas_concs' : '',
        'mo_gas_concentrations' : '',
        },
    'mo_gas_optics_rrtmgp' : {
        'ty_gas_optics_rrtmgp' : '',
        'mo_gas_optics_rrtmgp' : '',
        },
    'mo_optical_props' : {
        'ty_optical_props_1scl' : '',
        'ty_optical_props_2str' : '',
        'mo_optical_props' : '',
        },
    'mo_cloud_optics' : {
        'ty_cloud_optics' : '',
        'mo_cloud_optics' : '',
        },
    'mo_source_functions' : {
        'ty_source_func_lw' : '',
        'mo_source_functions' : '',
        },
    }

# Add all physics scheme dependencies relative to basedir - note that the CCPP
# rules stipulate that dependencies are not shared between the schemes!
SCHEME_FILES_DEPENDENCIES = [
    'ccpp/physics/physics/GFDL_parse_tracers.F90',
    'ccpp/physics/physics/aer_cloud.F',
    'ccpp/physics/physics/aerclm_def.F',
    'ccpp/physics/physics/aerinterp.F90',
    'ccpp/physics/physics/calpreciptype.f90',
    'ccpp/physics/physics/cldwat2m_micro.F',
    'ccpp/physics/physics/cldmacro.F',
    'ccpp/physics/physics/date_def.f',
    'ccpp/physics/physics/funcphys.f90',
    'ccpp/physics/physics/gfs_phy_tracer_config.F',
    'ccpp/physics/physics/gocart_tracer_config_stub.f',
    'ccpp/physics/physics/h2o_def.f',
    'ccpp/physics/physics/h2ointerp.f90',
    'ccpp/physics/physics/iccn_def.F',
    'ccpp/physics/physics/iccninterp.F90',
    'ccpp/physics/physics/iounitdef.f',
    'ccpp/physics/physics/machine.F',
    'ccpp/physics/physics/mersenne_twister.f',
    'ccpp/physics/physics/mfpbl.f',
    'ccpp/physics/physics/micro_mg_utils.F90',
    'ccpp/physics/physics/micro_mg2_0.F90',
    'ccpp/physics/physics/micro_mg3_0.F90',
    'ccpp/physics/physics/module_bfmicrophysics.f',
    'ccpp/physics/physics/multi_gases.F90',
    'ccpp/physics/physics/module_gfdl_cloud_microphys.F90',
    'ccpp/physics/physics/module_nst_model.f90',
    'ccpp/physics/physics/module_nst_parameters.f90',
    'ccpp/physics/physics/module_nst_water_prop.f90',
    'ccpp/physics/physics/module_mp_radar.F90',
    'ccpp/physics/physics/module_mp_thompson.F90',
    'ccpp/physics/physics/module_mp_thompson_make_number_concentrations.F90',
    'ccpp/physics/physics/module_MP_FER_HIRES.F90',
#    'ccpp/physics/physics/HWRF_mersenne_twister.F90',
#    'ccpp/physics/physics/HWRF_mcica_random_numbers.F90',
    'ccpp/physics/physics/module_bl_mynn.F90',
    'ccpp/physics/physics/module_sf_mynn.F90',
    'ccpp/physics/physics/module_SF_JSFC.F90',
    'ccpp/physics/physics/module_BL_MYJPBL.F90',
    'ccpp/physics/physics/module_sf_noahmp_glacier.f90',
    'ccpp/physics/physics/module_sf_noahmplsm.f90',
    'ccpp/physics/physics/cires_ugwp_module.F90',
    'ccpp/physics/physics/ugwp_driver_v0.F',
    'ccpp/physics/physics/cires_ugwp_triggers.F90',
    'ccpp/physics/physics/cires_ugwp_initialize.F90',
    'ccpp/physics/physics/cires_ugwp_solvers.F90',
    'ccpp/physics/physics/cires_ugwp_utils.F90',
    'ccpp/physics/physics/cires_orowam2017.f',
    'ccpp/physics/physics/cires_vert_lsatdis.F90',
    'ccpp/physics/physics/cires_vert_orodis.F90',
    'ccpp/physics/physics/cires_vert_wmsdis.F90',
    'ccpp/physics/physics/namelist_soilveg.f',
    'ccpp/physics/physics/mfpblt.f',
    'ccpp/physics/physics/mfpbltq.f',
    'ccpp/physics/physics/mfscu.f',
    'ccpp/physics/physics/mfscuq.f',
    'ccpp/physics/physics/noahmp_tables.f90',
    'ccpp/physics/physics/num_parthds.F',
    'ccpp/physics/physics/ozne_def.f',
    'ccpp/physics/physics/ozinterp.f90',
    'ccpp/physics/physics/physcons.F90',
    'ccpp/physics/physics/physparam.f',
    'ccpp/physics/physics/radcons.f90',
    'ccpp/physics/physics/radiation_aerosols.f',
    'ccpp/physics/physics/radiation_astronomy.f',
    'ccpp/physics/physics/radiation_clouds.f',
    'ccpp/physics/physics/radiation_gases.f',
    'ccpp/physics/physics/radiation_surface.f',
    'ccpp/physics/physics/radlw_datatb.f',
    'ccpp/physics/physics/radlw_param.f',
    'ccpp/physics/physics/radsw_datatb.f',
    'ccpp/physics/physics/radsw_param.f',
    'ccpp/physics/physics/samfaerosols.F',
    'ccpp/physics/physics/sfcsub.F',
    'ccpp/physics/physics/sflx.f',
    'ccpp/physics/physics/set_soilveg.f',
    'ccpp/physics/physics/surface_perturbation.F90',
    'ccpp/physics/physics/cu_gf_deep.F90',
    'ccpp/physics/physics/cu_gf_sh.F90',
    'ccpp/physics/physics/tridi.f',
    'ccpp/physics/physics/wv_saturation.F',
    'ccpp/physics/physics/module_sf_ruclsm.F90',
    'ccpp/physics/physics/namelist_soilveg_ruc.F90',
    'ccpp/physics/physics/set_soilveg_ruc.F90',
    'ccpp/physics/physics/module_soil_pre.F90',
#    'ccpp/physics/physics/module_sf_noahlsm.F90',
#    'ccpp/physics/physics/module_sf_noahlsm_glacial_only.F90',
#    'ccpp/physics/physics/module_sf_exchcoef.f90',
    # RRTMGP
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_optics.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_util_reorder.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_fluxes.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_rte_util_array.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_optical_props.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_rte_kind.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_rte_lw.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_rte_sw.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/mo_source_functions.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_compute_bc.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_fluxes_byband.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_fluxes_byband_kernels.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_fluxes_bygpoint.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_heating_rates.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/mo_rrtmgp_clr_all_sky.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90',
    'ccpp/physics/physics/rte-rrtmgp/extensions/cloud_optics/mo_cloud_sampling.F90',
    # derived data type definitions
    'scm/src/GFS_typedefs.F90',
    'scm/src/gmtb_scm_kinds.F90',
    'scm/src/gmtb_scm_physical_constants.F90',
    'scm/src/gmtb_scm_type_defs.F90',
]

# Add all physics scheme files relative to basedir
SCHEME_FILES = {
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'ccpp/physics/physics/GFS_DCNV_generic.F90'             : ['physics'],
    'ccpp/physics/physics/GFS_GWD_generic.F90'              : ['physics'],
    'ccpp/physics/physics/GFS_MP_generic.F90'               : ['physics'],
    'ccpp/physics/physics/GFS_PBL_generic.F90'              : ['physics'],
    'ccpp/physics/physics/GFS_SCNV_generic.F90'             : ['physics'],
    'ccpp/physics/physics/GFS_phys_time_vary.scm.F90'       : ['physics'],
    'ccpp/physics/physics/GFS_rad_time_vary.scm.F90'        : ['physics'],
    'ccpp/physics/physics/GFS_rrtmg_post.F90'               : ['physics'],
    'ccpp/physics/physics/GFS_rrtmg_pre.F90'                : ['physics'],
    'ccpp/physics/physics/GFS_rrtmg_setup.F90'              : ['physics'],
    'ccpp/physics/physics/GFS_suite_interstitial.F90'       : ['physics'],
    'ccpp/physics/physics/GFS_surface_generic.F90'          : ['physics'],
    'ccpp/physics/physics/GFS_surface_composites.F90'       : ['physics'],
    'ccpp/physics/physics/GFS_surface_loop_control.F90'     : ['physics'],
    'ccpp/physics/physics/GFS_time_vary_pre.scm.F90'        : ['physics'],
    'ccpp/physics/physics/cires_ugwp.F90'                   : ['physics'],
    'ccpp/physics/physics/cires_ugwp_post.F90'              : ['physics'],
    'ccpp/physics/physics/cnvc90.f'                         : ['physics'],
    'ccpp/physics/physics/cs_conv.F90'                      : ['physics'],
    'ccpp/physics/physics/cs_conv_aw_adj.F90'               : ['physics'],
    'ccpp/physics/physics/cu_ntiedtke_pre.F90'              : ['physics'],
    'ccpp/physics/physics/cu_ntiedtke.F90'                  : ['physics'],
    'ccpp/physics/physics/cu_ntiedtke_post.F90'             : ['physics'],
    'ccpp/physics/physics/dcyc2.f'                          : ['physics'],
    'ccpp/physics/physics/drag_suite.F90'                   : ['physics'],
    'ccpp/physics/physics/gcm_shoc.F90'                     : ['physics'],
    'ccpp/physics/physics/get_prs_fv3.F90'                  : ['physics'],
    'ccpp/physics/physics/gfdl_cloud_microphys.F90'         : ['physics'],
#    'ccpp/physics/physics/gfdl_sfc_layer.F90'               : ['physics'],
    'ccpp/physics/physics/gscond.f'                         : ['physics'],
    'ccpp/physics/physics/gwdc.f'                           : ['physics'],
    'ccpp/physics/physics/gwdps.f'                          : ['physics'],
    'ccpp/physics/physics/h2ophys.f'                        : ['physics'],
    'ccpp/physics/physics/samfdeepcnv.f'                    : ['physics'],
    'ccpp/physics/physics/samfshalcnv.f'                    : ['physics'],
    'ccpp/physics/physics/sascnvn.F'                        : ['physics'],
    'ccpp/physics/physics/shalcnv.F'                        : ['physics'],
    'ccpp/physics/physics/maximum_hourly_diagnostics.F90'   : ['physics'],
    'ccpp/physics/physics/m_micro.F90'                      : ['physics'],
    'ccpp/physics/physics/m_micro_interstitial.F90'         : ['physics'],
    'ccpp/physics/physics/cu_gf_driver_pre.F90'             : ['physics'],
    'ccpp/physics/physics/cu_gf_driver.F90'                 : ['physics'],
    'ccpp/physics/physics/cu_gf_driver_post.F90'            : ['physics'],
    'ccpp/physics/physics/moninedmf.f'                      : ['physics'],
    'ccpp/physics/physics/moninshoc.f'                      : ['physics'],
    'ccpp/physics/physics/satmedmfvdif.F'                   : ['physics'],
    'ccpp/physics/physics/satmedmfvdifq.F'                  : ['physics'],
    'ccpp/physics/physics/shinhongvdif.F90'                 : ['physics'],
    'ccpp/physics/physics/ysuvdif.F90'                      : ['physics'],
    'ccpp/physics/physics/module_MYNNPBL_wrapper.F90'       : ['physics'],
    'ccpp/physics/physics/module_MYNNSFC_wrapper.F90'       : ['physics'],
    'ccpp/physics/physics/module_MYNNrad_pre.F90'           : ['physics'],
    'ccpp/physics/physics/module_MYNNrad_post.F90'          : ['physics'],
    'ccpp/physics/physics/module_MYJSFC_wrapper.F90'        : ['physics' ],
    'ccpp/physics/physics/module_MYJPBL_wrapper.F90'        : ['physics' ],
    'ccpp/physics/physics/mp_thompson_pre.F90'              : ['physics'],
    'ccpp/physics/physics/mp_thompson.F90'                  : ['physics'],
    'ccpp/physics/physics/mp_thompson_post.F90'             : ['physics'],
    'ccpp/physics/physics/ozphys.f'                         : ['physics'],
    'ccpp/physics/physics/ozphys_2015.f'                    : ['physics'],
    'ccpp/physics/physics/precpd.f'                         : ['physics'],
#    'ccpp/physics/physics/radlw_main.F90'                   : ['physics'],
#    'ccpp/physics/physics/radsw_main.F90'                   : ['physics'],
    'ccpp/physics/physics/radlw_main.f'                     : ['physics'],
    'ccpp/physics/physics/radsw_main.f'                     : ['physics'],
    'ccpp/physics/physics/rayleigh_damp.f'                  : ['physics'],
    'ccpp/physics/physics/rrtmg_lw_post.F90'                : ['physics'],
    'ccpp/physics/physics/rrtmg_lw_pre.F90'                 : ['physics'],
    'ccpp/physics/physics/rrtmg_sw_post.F90'                : ['physics'],
    'ccpp/physics/physics/rrtmg_sw_pre.F90'                 : ['physics'],
    'ccpp/physics/physics/sfc_diag.f'                       : ['physics'],
    'ccpp/physics/physics/sfc_diag_post.F90'                : ['physics'],
    'ccpp/physics/physics/sfc_drv_ruc.F90'                  : ['physics'],
    'ccpp/physics/physics/lsm_ruc_sfc_sice_interstitial.F90': ['physics'],
    'ccpp/physics/physics/sfc_cice.f'                       : ['physics'],
    'ccpp/physics/physics/sfc_diff.f'                       : ['physics'],
    'ccpp/physics/physics/sfc_drv.f'                        : ['physics'],
#    'ccpp/physics/physics/sfc_noah_wrfv4_interstitial.F90'  : ['physics'],
#    'ccpp/physics/physics/sfc_noah_wrfv4.F90'               : ['physics'],
    'ccpp/physics/physics/sfc_noahmp_drv.f'                 : ['physics'],
    'ccpp/physics/physics/sfc_nst.f'                        : ['physics'],
    'ccpp/physics/physics/sfc_ocean.F'                      : ['physics'],
    'ccpp/physics/physics/sfc_sice.f'                       : ['physics'],
    'ccpp/physics/physics/mp_fer_hires.F90'                 : ['physics'],
    'ccpp/physics/physics/gmtb_scm_sfc_flux_spec.F90'       : ['physics'],
    # RRTMGP
    'ccpp/physics/physics/rrtmg_lw_cloud_optics.F90'        : ['physics'],
    'ccpp/physics/physics/rrtmg_sw_cloud_optics.F90'        : ['physics'],
    'ccpp/physics/physics/rrtmgp_aux.F90'                   : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_gas_optics.F90'         : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_cloud_optics.F90'       : ['physics'],
    'ccpp/physics/physics/rrtmgp_sw_gas_optics.F90'         : ['physics'],
    'ccpp/physics/physics/rrtmgp_sw_cloud_optics.F90'       : ['physics'],
    'ccpp/physics/physics/rrtmgp_sw_aerosol_optics.F90'     : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_rte.F90'                : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_cloud_sampling.F90'     : ['physics'],
    'ccpp/physics/physics/rrtmgp_sw_rte.F90'                : ['physics'],
    'ccpp/physics/physics/rrtmgp_sw_cloud_sampling.F90'     : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_aerosol_optics.F90'     : ['physics'],
    'ccpp/physics/physics/GFS_rrtmgp_setup.F90'             : ['physics'],
    'ccpp/physics/physics/GFS_rrtmgp_pre.F90'               : ['physics'],
    'ccpp/physics/physics/rrtmgp_lw_pre.F90'                : ['physics'],
    'ccpp/physics/physics/GFS_rrtmgp_sw_pre.F90'            : ['physics'],
    'ccpp/physics/physics/GFS_rrtmgp_lw_post.F90'           : ['physics'],
    'ccpp/physics/physics/GFS_rrtmgp_sw_post.F90'           : ['physics'],
    }

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'scm/bin'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = 'ccpp/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = 'ccpp/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = 'ccpp/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE = 'ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE = 'ccpp/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = 'ccpp/physics/CCPP_SCHEMES.sh'

# CCPP host cap in which to insert the ccpp_field_add statements;
# determines the directory to place ccpp_{modules,fields}.inc
TARGET_FILES = [
    'scm/src/gmtb_scm.F90',
    ]

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE = 'ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE = 'ccpp/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = 'ccpp/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = 'ccpp/physics/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'ccpp/suites'

# Optional arguments - only required for schemes that use
# optional arguments. ccpp_prebuild.py will throw an exception
# if it encounters a scheme subroutine with optional arguments
# if no entry is made here. Possible values are: 'all', 'none',
# or a list of standard_names: [ 'var1', 'var3' ].
OPTIONAL_ARGUMENTS = {
    'rrtmg_sw' : {
        'rrtmg_sw_run' : [
            'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels',
            'components_of_surface_downward_shortwave_fluxes',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'rrtmg_lw' : {
        'rrtmg_lw_run' : [
            'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson' : {
        'mp_thompson_init' : [
            'water_friendly_aerosol_number_concentration',
            'ice_friendly_aerosol_number_concentration',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        'mp_thompson_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            'mean_effective_radius_for_liquid_cloud',
            'mean_effective_radius_for_ice_cloud',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson_pre' : {
        'mp_thompson_pre_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        },
    'mp_fer_hires' : {
        'mp_fer_hires_init' : [
            'fraction_of_ice_water_cloud',
            'fraction_of_rain_water_cloud',
            'rime_factor',
            ],
        },
    'rrtmgp_sw_rte' : {
         'rrtmgp_sw_rte_run' : [
             'components_of_surface_downward_shortwave_fluxes',
             'sw_fluxes_sfc',
             'sw_fluxes_toa',
             ],
         },
    'GFS_rrtmgp_sw_post' : {
         'GFS_rrtmgp_sw_post_run' : [
             'components_of_surface_downward_shortwave_fluxes',
             'sw_fluxes_sfc',
             'sw_fluxes_toa',
             ],
         },
    'rrtmgp_lw_rte' : {
         'rrtmgp_lw_rte_run' : [
             'lw_fluxes_sfc',
             'lw_fluxes_toa',
             ],
        },
    'GFS_rrtmgp_lw_post' : {
         'GFS_rrtmgp_lw_post_run' : [
             'lw_fluxes_sfc',
             'lw_fluxes_toa',
             ],
         },
    'GFS_rrtmgp_post' : {
         'GFS_rrtmgp_post_run' : [
             'components_of_surface_downward_shortwave_fluxes',
             ],
         },
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Names of Fortran include files in the host model cap (do not change);
# both files will be written to the directory of each target file
MODULE_INCLUDE_FILE = 'ccpp_modules.inc'
FIELDS_INCLUDE_FILE = 'ccpp_fields.inc'

# Directory where to write static API to
STATIC_API_DIR = 'scm/src/'
STATIC_API_SRCFILE = 'scm/src/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = 'ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = 'ccpp/physics/CCPP_VARIABLES_SCM.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_SCM.tex'


###############################################################################
# Template code to generate include files                                     #
###############################################################################

# Name of the CCPP data structure in the host model cap;
# in the case of SCM, this is a vector with loop index i
CCPP_DATA_STRUCTURE = 'cdata'
