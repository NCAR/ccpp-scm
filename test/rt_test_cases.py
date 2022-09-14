run_list = [\
            #----------------------------------------------------------------------------------------------------------------------------------------------
            # CCPP-SCM v5 suites
            #----------------------------------------------------------------------------------------------------------------------------------------------
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v15p2"},                                                                                 \ #NOT WORKING IN SCM
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_csawmg"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GSD_v1"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v15p2"},                                                                                 \ #NOT WORKING IN SCM
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16"},                                                                                   \ #NOT WORKING IN SCM
            {"case": "LASSO_2016051812",      "suite": "SCM_csawmg"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "LASSO_2016051812",      "suite": "SCM_GSD_v1"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "astex",                 "suite": "SCM_GFS_v15p2"},                                                                                 \ #NOT WORKING IN SCM
            {"case": "astex",                 "suite": "SCM_csawmg"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "astex",                 "suite": "SCM_GSD_v1"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "bomex",                 "suite": "SCM_GFS_v15p2"},                                                                                 \ #NOT WORKING IN SCM
            {"case": "bomex",                 "suite": "SCM_csawmg"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "bomex",                 "suite": "SCM_GSD_v1"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "twpice",                "suite": "SCM_GFS_v15p2"},                                                                                 \ #NOT WORKING IN SCM
            {"case": "twpice",                "suite": "SCM_csawmg"},                                                                                    \ #NOT WORKING IN SCM 
            {"case": "twpice",                "suite": "SCM_GSD_v1"},                                                                                    \ #NOT WORKING IN SCM
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RRFS_v1alpha", "namelist": "input_RRFS_v1alpha.nml", "tracer": "tracers_RRFS_v1alpha.txt"},  \ #NOT WORKING IN SCM
            {"case": "LASSO_2016051812",      "suite": "SCM_RRFS_v1alpha", "namelist": "input_RRFS_v1alpha.nml", "tracer": "tracers_RRFS_v1alpha.txt"},  \ #NOT WORKING IN SCM
            {"case": "astex",                 "suite": "SCM_RRFS_v1alpha", "namelist": "input_RRFS_v1alpha.nml", "tracer": "tracers_RRFS_v1alpha.txt"},  \ #NOT WORKING IN SCM
            {"case": "bomex",                 "suite": "SCM_RRFS_v1alpha", "namelist": "input_RRFS_v1alpha.nml", "tracer": "tracers_RRFS_v1alpha.txt"},  \ #NOT WORKING IN SCM
            {"case": "twpice",                "suite": "SCM_RRFS_v1alpha", "namelist": "input_RRFS_v1alpha.nml", "tracer": "tracers_RRFS_v1alpha.txt"},  \ #NOT WORKING IN SCM
            #----------------------------------------------------------------------------------------------------------------------------------------------
            # CCPP-SCM v6 supported suites 
            #----------------------------------------------------------------------------------------------------------------------------------------------
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v16"},                                                                                   \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_GFS_v17_p8"},                                                                                \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_HRRR"},                                                                                      \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RRFS_v1beta"},                                                                               \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_RAP"},                                                                                       \
            {"case": "arm_sgp_summer_1997_A", "suite": "SCM_WoFS_v0"},                                                                                   \
            {"case": "twpice",                "suite": "SCM_GFS_v16"},                                                                                   \
            {"case": "twpice",                "suite": "SCM_GFS_v17_p8"},                                                                                \
            {"case": "twpice",                "suite": "SCM_HRRR"},                                                                                      \
            {"case": "twpice",                "suite": "SCM_RRFS_v1beta"},                                                                               \
            {"case": "twpice",                "suite": "SCM_RAP"},                                                                                       \
            {"case": "twpice",                "suite": "SCM_WoFS_v0"},                                                                                   \
            {"case": "bomex",                 "suite": "SCM_GFS_v16"},                                                                                   \
            {"case": "bomex",                 "suite": "SCM_GFS_v17_p8"},                                                                                \
            {"case": "bomex",                 "suite": "SCM_HRRR"},                                                                                      \
            {"case": "bomex",                 "suite": "SCM_RRFS_v1beta"},                                                                               \
            {"case": "bomex",                 "suite": "SCM_RAP"},                                                                                       \
            {"case": "bomex",                 "suite": "SCM_WoFS_v0"},                                                                                   \
            {"case": "astex",                 "suite": "SCM_GFS_v16"},                                                                                   \
            {"case": "astex",                 "suite": "SCM_GFS_v17_p8"},                                                                                \
            {"case": "astex",                 "suite": "SCM_HRRR"},                                                                                      \
            {"case": "astex",                 "suite": "SCM_RRFS_v1beta"},                                                                               \
            {"case": "astex",                 "suite": "SCM_RAP"},                                                                                       \
            {"case": "astex",                 "suite": "SCM_WoFS_v0"},                                                                                   \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v16"},                                                                                   \
            {"case": "LASSO_2016051812",      "suite": "SCM_GFS_v17_p8"},                                                                                \
            {"case": "LASSO_2016051812",      "suite": "SCM_HRRR"},                                                                                      \
            {"case": "LASSO_2016051812",      "suite": "SCM_RRFS_v1beta"},                                                                               \
            {"case": "LASSO_2016051812",      "suite": "SCM_RAP"},                                                                                       \
            {"case": "LASSO_2016051812",      "suite": "SCM_WoFS_v0"},                                                                                   \
            #---------------------------------------------------------------------------------------------------------------------------------------------------
            # Unsupported suites
            #--------------------------------------------------------------------------------------------------------------------------------------------------- 
            {"case": "gabls3",                "suite": "SCM_GFS_v16"},                                                                                          \
            {"case": "gabls3_noahmp",         "suite": "SCM_GFS_v17_p8"},                                                                                       \
            {"case": "gabls3_ruc",            "suite": "SCM_RAP"},                                                                                              \
#            {"case": "ARMCU_REF",             "suite": "SCM_GFS_v16"},                                                                                          \
            {"case": "fv3_model_point_noah",  "suite": "SCM_GFS_v16"},                                                                                          \
            {"case": "twpice",                "suite": "SCM_GFS_v15p2_RRTMGP",  "namelist": "input_GFS_v15p2_RRTMGP.nml",  "tracer": "tracers_GFS_v15p2.txt"},  \
            {"case": "twpice",                "suite": "SCM_GFS_v16_RRTMGP",    "namelist": "input_GFS_v16_RRTMGP.nml",    "tracer": "tracers_GFS_v16.txt"},    \
            {"case": "twpice",                "suite": "SCM_GFS_v17_p8_RRTMGP", "namelist": "input_GFS_v17_p8_RRTMGP.nml", "tracer": "tracers_GFS_v17_p8.txt"}, \
            {"case": "twpice",                "suite": "SCM_RAP_RRTMGP",        "namelist": "input_RAP_RRTMGP.nml",        "tracer": "tracers_RAP.txt"}]
