#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Supported suites for CCPP Version 7 release.
#----------------------------------------------------------------------------------------------------------------------------------------------------------
suites_supported_gfortran = [\
             "SCM_GFS_v17_p8_ugwpv1",   "SCM_GFS_v16_RRTMGP",   "SCM_GFS_v16",   "SCM_WoFS_v0",   "SCM_HRRR_gf",                                          \
             "SCM_GFS_v17_p8_ugwpv1_ps","SCM_GFS_v16_RRTMGP_ps","SCM_GFS_v16_ps","SCM_WoFS_v0_ps","SCM_HRRR_gf_ps"]
run_list_supported_gfortran = [\
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "twpice",                "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "twpice",                "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "bomex",                 "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "bomex",                 "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "astex",                 "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "astex",                 "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "COMBLE",                "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "COMBLE",                "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "gabls3",                "suite": "SCM_GFS_v16"}]
#
suites_supported_ifx = suites_supported_gfortran
run_list_supported_ifx = [\
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            #{"case": "twpice",                "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "twpice",                "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            {"case": "bomex",                 "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "bomex",                 "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            #{"case": "astex",                 "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "astex",                 "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            #{"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            #{"case": "COMBLE",                "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "COMBLE",                "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "COMBLE",                "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v17_p8_ugwpv1"},                                                                          \
            #{"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v16_RRTMGP"},                                                                             \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_GFS_v16"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_WoFS_v0"},                                                                                    \
            {"case": "MOSAiC-AMPS",           "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "gabls3",                "suite": "SCM_GFS_v16"}]
#
suites_supported_nvhpc = ["SCM_RAP", "SCM_RAP_ps"]
run_list_supported_nvhpc = [\
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RAP"},                                                                                        \
            {"case": "twpice",                "suite": "SCM_RAP"},                                                                                        \
            {"case": "bomex",                 "suite": "SCM_RAP"},                                                                                        \
            {"case": "astex",                 "suite": "SCM_RAP"},                                                                                        \
            {"case": "LASSO_2016051812",      "suite": "SCM_RAP"}]

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Developmental suites, (w/ supported cases).
#----------------------------------------------------------------------------------------------------------------------------------------------------------
suites_dev_gfortran = [\
            " SCM_GFS_v16_no_nsst",    "SCM_GFS_v17_p8_ugwpv1_no_nsst",    "SCM_RRFS_v1beta_no_nsst",    "SCM_GFS_v17_p8_ugwpv1_tempo",                   \
             "SCM_GFS_v16_no_nsst_ps", "SCM_GFS_v17_p8_ugwpv1_no_nsst_ps", "SCM_RRFS_v1beta_no_nsst_ps", "SCM_GFS_v17_p8_ugwpv1_tempo_ps",                \
             "SCM_GFS_v16_gfdlmpv3",   "SCM_GFS_v15p2_ntiedtke",           "SCM_GFS_v16_debug",       	      	      	      	      	                  \
             "SCM_GFS_v16_gfdlmpv3_ps","SCM_GFS_v15p2_ntiedtke_ps",        "SCM_GFS_v16_debug_ps"]
run_list_dev_gfortran = [\
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_GFS_v16_no_nsst"},                                                                          \
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_GFS_v17_p8_ugwpv1_no_nsst"},                                                                \
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_RRFS_v1beta_no_nsst"},                                                                      \
            {"case": "arm_sgp_summer_1997_A",   "suite": "SCM_GFS_v17_p8_ugwpv1_tempo"},                                                                  \
            {"case": "arm_sgp_summer_1997_A",   "suite": "SCM_GFS_v16_gfdlmpv3"},                                                                         \
            {"case": "twpice",                  "suite": "SCM_GFS_v15p2_ntiedtke"},                                                                       \
            {"case": "bomex",                   "suite": "SCM_GFS_v16_debug"}]
#
suites_dev_ifx = [\
             "SCM_GFS_v16_no_nsst",     "SCM_GFS_v17_p8_ugwpv1_no_nsst",    "SCM_RRFS_v1beta_no_nsst",                                                    \
             "SCM_GFS_v16_no_nsst_ps",  "SCM_GFS_v17_p8_ugwpv1_no_nsst_ps", "SCM_RRFS_v1beta_no_nsst_ps",                                                 \
             "SCM_GFS_v16_gfdlmpv3",    "SCM_GFS_v15p2_ntiedtke",           "SCM_GFS_v16_debug",                                                          \
             "SCM_GFS_v16_gfdlmpv3_ps", "SCM_GFS_v15p2_ntiedtke_ps",        "SCM_GFS_v16_debug_ps"]
run_list_dev_ifx = [\
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_GFS_v16_no_nsst"},                                                                          \
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_GFS_v17_p8_ugwpv1_no_nsst"},                                                                \
            {"case": "atomic_Jan16T22Jan18T06", "suite": "SCM_RRFS_v1beta_no_nsst"},                                                                      \
            #{"case": "arm_sgp_summer_1997_A",   "suite": "SCM_GFS_v17_p8_ugwpv1_tempo"},                                                                  \
            {"case": "arm_sgp_summer_1997_A",   "suite": "SCM_GFS_v16_gfdlmpv3"},                                                                         \
            {"case": "twpice",                  "suite": "SCM_GFS_v15p2_ntiedtke"},                                                                       \
            {"case": "bomex",                   "suite": "SCM_GFS_v16_debug"}]

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Legacy suites, (w/ supported cases).
#----------------------------------------------------------------------------------------------------------------------------------------------------------
suites_legacy_gfortran = [\
             "SCM_GFS_v17_p8",    "SCM_HRRR",    "SCM_RRFS_v1beta",    "SCM_RAP",    "SCM_GFS_v15p2",                                                     \
             "SCM_GFS_v17_p8_ps", "SCM_HRRR_ps", "SCM_RRFS_v1beta_ps", "SCM_RAP_ps", "SCM_GFS_v15p2_ps"]
run_list_legacy_gfortran = [\
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v17_p8"},                                                                                 \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_HRRR"},                                                                                       \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RRFS_v1beta"},                                                                                \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RAP"},                                                                                        \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v15p2"},                                                                                  \
            {"case": "twpice",                "suite": "SCM_GFS_v17_p8"},                                                                                 \
            {"case": "twpice",                "suite": "SCM_HRRR"},                                                                                       \
            {"case": "twpice",                "suite": "SCM_RRFS_v1beta"},                                                                                \
            {"case": "twpice",                "suite": "SCM_RAP"},                                                                                        \
            {"case": "twpice",                "suite": "SCM_GFS_v15p2"},                                                                                  \
            {"case": "bomex",                 "suite": "SCM_GFS_v17_p8"},                                                                                 \
            {"case": "bomex",                 "suite": "SCM_HRRR"},                                                                                       \
            {"case": "bomex",                 "suite": "SCM_RRFS_v1beta"},                                                                                \
            {"case": "bomex",                 "suite": "SCM_RAP"},                                                                                        \
            {"case": "bomex",                 "suite": "SCM_GFS_v15p2"},                                                                                  \
            {"case": "astex",                 "suite": "SCM_GFS_v17_p8"},                                                                                 \
            {"case": "astex",                 "suite": "SCM_HRRR"},                                                                                       \
            {"case": "astex",                 "suite": "SCM_RRFS_v1beta"},                                                                                \
            {"case": "astex",                 "suite": "SCM_RAP"},                                                                                        \
            {"case": "astex",                 "suite": "SCM_GFS_v15p2"},                                                                                  \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v17_p8"},                                                                                 \
            {"case": "LASSO_2016051812",      "suite": "SCM_HRRR"},                                                                                       \
            {"case": "LASSO_2016051812",      "suite": "SCM_RRFS_v1beta"},                                                                                \
            {"case": "LASSO_2016051812",      "suite": "SCM_RAP"},                                                                                        \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v15p2"}]
#
suites_legacy_ifx   = suites_legacy_gfortran
run_list_legacy_ifx = run_list_legacy_gfortran

#----------------------------------------------------------------------------------------------------------------------------------------------------------
# Single precision supported suites, (w/ supported cases).
#----------------------------------------------------------------------------------------------------------------------------------------------------------
suites_sp_gfortran = ["SCM_HRRR_gf", "SCM_HRRR_gf_ps"]
run_list_sp_gfortran = [\
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "twpice",                "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "bomex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "astex",                 "suite": "SCM_HRRR_gf"},                                                                                    \
            {"case": "LASSO_2016051812",      "suite": "SCM_HRRR_gf"}]
#
suites_sp_ifx   = suites_sp_gfortran
run_list_sp_ifx = run_list_sp_gfortran

# make this work with suite_info.py
class gnu_test_cases:
    run_list_supported = run_list_supported_gfortran
    run_list_legacy = run_list_legacy_gfortran
    run_list_dev = run_list_dev_gfortran
    run_list_sp = run_list_sp_gfortran

class intel_test_cases:
    run_list_supported = run_list_supported_ifx
    run_list_legacy = run_list_legacy_ifx
    run_list_dev = run_list_dev_ifx
    run_list_sp = run_list_sp_ifx

class nvhpc_test_cases:
    run_list_supported = run_list_supported_nvhpc
    run_list_legacy = []
    run_list_dev = []
    run_list_sp = []
