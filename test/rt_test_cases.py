run_list = [\
            #----------------------------------------------------------------------------------------------------------------------------------------------
            # Supported suites for CCPP Version 7 release
            #----------------------------------------------------------------------------------------------------------------------------------------------
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
            {"case": "gabls3",                "suite": "SCM_GFS_v16"},                                                                                    \
            #----------------------------------------------------------------------------------------------------------------------------------------------
            # Unsupported suites (w/ supported cases)
            #----------------------------------------------------------------------------------------------------------------------------------------------
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
