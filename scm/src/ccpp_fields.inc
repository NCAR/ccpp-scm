
            call ccpp_field_add(cdata(i), 'FV3-GFS_Cldprop_type', '', c_loc(physics%Cldprop(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Cldprop_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Control_type', '', c_loc(physics%Model(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Control_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Coupling_type', '', c_loc(physics%Coupling(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Coupling_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Diag_type', '', c_loc(physics%Diag(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Diag_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Grid_type', '', c_loc(physics%Grid(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Grid_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Init_type', '', c_loc(physics%Init_parm(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Init_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Interstitial_type', '', c_loc(physics%Interstitial(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Interstitial_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Radtend_type', '', c_loc(physics%Radtend(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Radtend_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Sfccycle_type', '', c_loc(physics%Sfccycle(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Sfccycle_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Sfcprop_type', '', c_loc(physics%Sfcprop(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Sfcprop_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Statein_type', '', c_loc(physics%Statein(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Statein_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Stateout_type', '', c_loc(physics%Stateout(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Stateout_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'FV3-GFS_Tbd_type', '', c_loc(physics%Tbd(i)), ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "FV3-GFS_Tbd_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'Monin-Obukhov_similarity_function_for_heat', physics%Sfcprop(i)%ffhh, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "Monin-Obukhov_similarity_function_for_heat" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'Monin-Obukhov_similarity_function_for_heat_at_2m', physics%Interstitial(i)%fh2, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "Monin-Obukhov_similarity_function_for_heat_at_2m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'Monin-Obukhov_similarity_function_for_momentum', physics%Sfcprop(i)%ffmm, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "Monin-Obukhov_similarity_function_for_momentum" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'Monin-Obukhov_similarity_function_for_momentum_at_10m', physics%Interstitial(i)%fm10, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "Monin-Obukhov_similarity_function_for_momentum_at_10m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90', physics%Tbd(i)%acv, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'accumulated_lwe_thickness_of_precipitation_amount', physics%Diag(i)%totprcp, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "accumulated_lwe_thickness_of_precipitation_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'adjusted_vertical_layer_dimension_for_radiation', physics%Interstitial(i)%lmk, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "adjusted_vertical_layer_dimension_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'adjusted_vertical_level_dimension_for_radiation', physics%Interstitial(i)%lmp, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "adjusted_vertical_level_dimension_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_asymmetry_parameter_for_longwave_bands_01-16', physics%Interstitial(i)%faerlw(:,:,:,3), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_asymmetry_parameter_for_longwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_asymmetry_parameter_for_shortwave_bands_01-16', physics%Interstitial(i)%faersw(:,:,:,3), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_asymmetry_parameter_for_shortwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_optical_depth_for_longwave_bands_01-16', physics%Interstitial(i)%faerlw(:,:,:,1), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_optical_depth_for_longwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_optical_depth_for_shortwave_bands_01-16', physics%Interstitial(i)%faersw(:,:,:,1), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_optical_depth_for_shortwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_single_scattering_albedo_for_longwave_bands_01-16', physics%Interstitial(i)%faerlw(:,:,:,2), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_single_scattering_albedo_for_longwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'aerosol_single_scattering_albedo_for_shortwave_bands_01-16', physics%Interstitial(i)%faersw(:,:,:,2), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "aerosol_single_scattering_albedo_for_shortwave_bands_01-16" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure', physics%Statein(i)%prsl, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure_at_interface', physics%Statein(i)%prsi, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure_at_interface" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure_at_interface_for_radiation_in_hPa', physics%Interstitial(i)%plvl, ierr, 'hPa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure_at_interface_for_radiation_in_hPa" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure_at_layer_for_radiation_in_hPa', physics%Interstitial(i)%plyr, ierr, 'hPa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure_at_layer_for_radiation_in_hPa" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure_at_lowest_model_layer', physics%Statein(i)%prsl(:,1), ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_pressure_difference_between_midlayers', physics%Interstitial(i)%del, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_pressure_difference_between_midlayers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature', physics%Statein(i)%tgrs, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_at_interface_for_radiation', physics%Interstitial(i)%tlvl, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_at_interface_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_at_layer_for_radiation', physics%Interstitial(i)%tlyr, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_at_layer_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_at_lowest_model_layer', physics%Statein(i)%tgrs(:,1), ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_at_previous_time_step', physics%Tbd(i)%phy_f3d(:,:,3), ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_at_previous_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_save', physics%Interstitial(i)%save_t, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_save" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_two_time_steps_back', physics%Tbd(i)%phy_f3d(:,:,1), ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_two_time_steps_back" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'air_temperature_updated_by_physics', physics%Stateout(i)%gt0, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'angle_from_east_of_maximum_subgrid_orographic_variations', physics%Interstitial(i)%theta, ierr, 'degrees')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "angle_from_east_of_maximum_subgrid_orographic_variations" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'anisotropy_of_subgrid_orography', physics%Interstitial(i)%gamma, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "anisotropy_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'array_dimension_of_microphysics', physics%Model(i)%num_p3d, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "array_dimension_of_microphysics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'array_dimension_of_random_number', physics%Model(i)%nrcm, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "array_dimension_of_random_number" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'asymmetry_of_subgrid_orography', physics%Interstitial(i)%oa4, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "asymmetry_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'atmosphere_boundary_layer_thickness', physics%Diag(i)%hpbl, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "atmosphere_boundary_layer_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'atmosphere_heat_diffusivity', physics%Interstitial(i)%dkt, ierr, 'm2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "atmosphere_heat_diffusivity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'atmosphere_heat_diffusivity_background', physics%Model(i)%xkzm_h, ierr, 'm2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "atmosphere_heat_diffusivity_background" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'atmosphere_momentum_diffusivity_background', physics%Model(i)%xkzm_m, ierr, 'm2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "atmosphere_momentum_diffusivity_background" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles', physics%Interstitial(i)%aerodp, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "atmosphere_optical_thickness_due_to_ambient_aerosol_particles" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'bulk_richardson_number_at_lowest_model_level', physics%Interstitial(i)%rb, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "bulk_richardson_number_at_lowest_model_level" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'canopy_upward_latent_heat_flux', physics%Interstitial(i)%evcw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "canopy_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'canopy_water_amount', physics%Sfcprop(i)%canopy, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "canopy_water_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cell_area', physics%Grid(i)%area, ierr, 'm2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cell_area" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cell_size', physics%Grid(i)%dx, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cell_size" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cell_soil_type', physics%Interstitial(i)%soiltype, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cell_soil_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cell_vegetation_type', physics%Interstitial(i)%vegtype, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cell_vegetation_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'change_in_ozone_concentration', physics%Interstitial(i)%dq3dt_loc(:,:,6:6+physics%Interstitial(i)%oz_coeff-1), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "change_in_ozone_concentration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'characteristic_grid_length_scale', physics%Interstitial(i)%dlength, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "characteristic_grid_length_scale" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_area_fraction', physics%Interstitial(i)%cldf, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_area_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_area_fraction_for_radiation', physics%Interstitial(i)%cldsa, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_area_fraction_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_condensed_water_conversion_threshold', physics%Model(i)%wminco, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_condensed_water_conversion_threshold" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_condensed_water_specific_humidity_save', physics%Interstitial(i)%save_qcw, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_condensed_water_specific_humidity_save" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_condensed_water_specific_humidity_updated_by_physics', physics%Stateout(i)%gq0(:,:,scm_state%cloud_water_index), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_condensed_water_specific_humidity_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_ice_specific_humidity', physics%Interstitial(i)%clw(:,:,1), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_ice_specific_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_ice_water_path', physics%Interstitial(i)%clouds(:,:,4), ierr, 'g m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_ice_water_path" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_liquid_water_path', physics%Interstitial(i)%clouds(:,:,2), ierr, 'g m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_liquid_water_path" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_liquid_water_specific_humidity', physics%Interstitial(i)%clw(:,:,2), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_liquid_water_specific_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_rain_water_path', physics%Interstitial(i)%clouds(:,:,6), ierr, 'g m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_rain_water_path" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_snow_water_path', physics%Interstitial(i)%clouds(:,:,8), ierr, 'g m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_snow_water_path" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cloud_work_function', physics%Interstitial(i)%cld1d, ierr, 'm2 s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cloud_work_function" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_c_0', physics%Sfcprop(i)%c_0, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_c_0" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_c_d', physics%Sfcprop(i)%c_d, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_c_d" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_for_evaporation_of_rainfall', physics%Model(i)%evpco, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_for_evaporation_of_rainfall" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_from_cloud_ice_to_snow', physics%Model(i)%psautco, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_from_cloud_ice_to_snow" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_from_cloud_water_to_rain', physics%Model(i)%prautco, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_from_cloud_water_to_rain" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_w_0', physics%Sfcprop(i)%w_0, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_w_0" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'coefficient_w_d', physics%Sfcprop(i)%w_d, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "coefficient_w_d" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'column_precipitable_water', physics%Diag(i)%pwat, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "column_precipitable_water" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'components_of_surface_downward_shortwave_fluxes', '', c_loc(physics%Interstitial(i)%scmpsw), rank=size(shape(physics%Interstitial(i)%scmpsw)), dims=shape(physics%Interstitial(i)%scmpsw), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "components_of_surface_downward_shortwave_fluxes" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'convective_cloud_cover', physics%Interstitial(i)%cnvc, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "convective_cloud_cover" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'convective_cloud_switch', physics%Model(i)%clstp, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "convective_cloud_switch" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'convective_cloud_water_specific_humidity', physics%Interstitial(i)%cnvw, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "convective_cloud_water_specific_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'convective_transportable_tracers', physics%Interstitial(i)%clw, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "convective_transportable_tracers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'convexity_of_subgrid_orography', physics%Interstitial(i)%oc, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "convexity_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cosine_of_latitude', physics%Grid(i)%coslat, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cosine_of_latitude" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cosine_of_solar_declination_angle', physics%Model(i)%cdec, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cosine_of_solar_declination_angle" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cosine_of_zenith_angle', physics%Radtend(i)%coszen, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cosine_of_zenith_angle" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'countergradient_mixing_term_for_temperature', physics%Interstitial(i)%gamt, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "countergradient_mixing_term_for_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'countergradient_mixing_term_for_water_vapor', physics%Interstitial(i)%gamq, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "countergradient_mixing_term_for_water_vapor" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'critical_relative_humidity', physics%Interstitial(i)%rhc, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "critical_relative_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'critical_relative_humidity_at_PBL_top', physics%Interstitial(i)%rhcpbl, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "critical_relative_humidity_at_PBL_top" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'critical_relative_humidity_at_surface', physics%Interstitial(i)%rhcbot, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "critical_relative_humidity_at_surface" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'critical_relative_humidity_at_top_of_atmosphere', physics%Interstitial(i)%rhctop, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "critical_relative_humidity_at_top_of_atmosphere" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cumulative_change_in_temperature_due_to_surface_processes', physics%Diag(i)%dt3dt(:,:,2), ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cumulative_change_in_temperature_due_to_surface_processes" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag', physics%Diag(i)%du3dt(:,:,4), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cumulative_change_in_x_wind_due_to_surface_processes', physics%Diag(i)%du3dt(:,:,2), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cumulative_change_in_x_wind_due_to_surface_processes" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag', physics%Diag(i)%dv3dt(:,:,4), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'cumulative_change_in_y_wind_due_to_surface_processes', physics%Diag(i)%dv3dt(:,:,2), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "cumulative_change_in_y_wind_due_to_surface_processes" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'daytime_points', physics%Interstitial(i)%idxday, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "daytime_points" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'daytime_points_dimension', physics%Interstitial(i)%nday, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "daytime_points_dimension" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'deep_soil_temperature', physics%Sfcprop(i)%tg3, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "deep_soil_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'diffusivity_background_sigma_level', physics%Model(i)%xkzm_s, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "diffusivity_background_sigma_level" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dimensionless_exner_function_at_lowest_model_interface', physics%Statein(i)%prsik(:,1), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dimensionless_exner_function_at_lowest_model_interface" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dimensionless_exner_function_at_lowest_model_layer', physics%Statein(i)%prslk(:,1), ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dimensionless_exner_function_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dimensionless_exner_function_at_model_layers', physics%Statein(i)%prslk, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dimensionless_exner_function_at_model_layers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'diurnal_thermocline_layer_heat_content', physics%Sfcprop(i)%xt, ierr, 'K m')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "diurnal_thermocline_layer_heat_content" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'diurnal_thermocline_layer_thickness', physics%Sfcprop(i)%xz, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "diurnal_thermocline_layer_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'diurnal_thermocline_layer_x_current', physics%Sfcprop(i)%xu, ierr, 'm2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "diurnal_thermocline_layer_x_current" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'diurnal_thermocline_layer_y_current', physics%Sfcprop(i)%xv, ierr, 'm2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "diurnal_thermocline_layer_y_current" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dominant_freezing_rain_type', physics%Interstitial(i)%domzr, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dominant_freezing_rain_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dominant_rain_type', physics%Interstitial(i)%domr, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dominant_rain_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dominant_sleet_type', physics%Interstitial(i)%domip, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dominant_sleet_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dominant_snow_type', physics%Interstitial(i)%doms, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dominant_snow_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'dynamics_to_physics_timestep_ratio', physics%Interstitial(i)%frain, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "dynamics_to_physics_timestep_ratio" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'equation_of_time', physics%Model(i)%slag, ierr, 'radians')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "equation_of_time" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'error_flag', physics%Interstitial(i)%errflg, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "error_flag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'error_message', physics%Interstitial(i)%errmsg, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "error_message" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'extra_top_layer', physics%LTP, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "extra_top_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_TKE_dissipation_heating', physics%Model(i)%dspheat, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_TKE_dissipation_heating" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_deep_convection', physics%Interstitial(i)%kcnv, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_deep_convection" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_diagnostics', physics%Model(i)%lssav, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_diagnostics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_diagnostics_3D', physics%Model(i)%ldiag3d, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_diagnostics_3D" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_guess_run', physics%Interstitial(i)%flag_guess, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_guess_run" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_iteration', physics%Interstitial(i)%flag_iter, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_iteration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_land_surface_scheme', physics%Model(i)%lsm, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_land_surface_scheme" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_mom4_coupling', physics%Model(i)%mom4ice, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_mom4_coupling" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_nsstm_run', physics%Model(i)%nstf_name(1), ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_nsstm_run" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_precipitation_type', physics%Sfcprop(i)%srflag, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_precipitation_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_precipitation_type_algorithm', physics%Model(i)%cal_pre, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_precipitation_type_algorithm" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_for_reduced_drag_coefficient_over_sea', physics%Model(i)%redrag, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_for_reduced_drag_coefficient_over_sea" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_idealized_physics', physics%Model(i)%lsidea, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_idealized_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_print', physics%Model(i)%lprnt, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_print" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_to_calc_lw', physics%Model(i)%lslwr, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_to_calc_lw" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'flag_to_calc_sw', physics%Model(i)%lsswr, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "flag_to_calc_sw" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'forecast_hour', physics%Model(i)%solhr, ierr, 'h')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "forecast_hour" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'forecast_time', physics%Model(i)%fhour, ierr, 'h')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "forecast_time" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'fraction_of_convective_cloud', physics%Cldprop(i)%cv, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "fraction_of_convective_cloud" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height', physics%Interstitial(i)%clx, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'free_convection_layer_thickness', physics%Sfcprop(i)%d_conv, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "free_convection_layer_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'gas_constant_dry_air', con_rd, ierr, 'J kg-1 K-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "gas_constant_dry_air" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'gas_constant_water_vapor', con_rv, ierr, 'J kg-1 K-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "gas_constant_water_vapor" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'geopotential', physics%Statein(i)%phil, ierr, 'm2 s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "geopotential" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'geopotential_at_interface', physics%Statein(i)%phii, ierr, 'm2 s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "geopotential_at_interface" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature', physics%Interstitial(i)%del_gz, ierr, 'm2 s-2 K-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'gravitational_acceleration', con_g, ierr, 'm s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "gravitational_acceleration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'grid_size_related_coefficient_used_in_scale-sensitive_schemes', physics%Interstitial(i)%work1, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "grid_size_related_coefficient_used_in_scale-sensitive_schemes" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement', physics%Interstitial(i)%work2, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'height_above_mean_sea_level_at_lowest_model_layer', physics%Diag(i)%zlvl, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "height_above_mean_sea_level_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'horizontal_dimension', physics%Interstitial(i)%ix, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "horizontal_dimension" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'horizontal_index_of_printed_column', physics%Interstitial(i)%ipr, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "horizontal_index_of_printed_column" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'horizontal_loop_extent', physics%Interstitial(i)%im, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "horizontal_loop_extent" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'index_for_liquid_cloud_condensate', physics%Model(i)%ntcw, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "index_for_liquid_cloud_condensate" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'index_of_dtlm_start', physics%Sfcprop(i)%ifd, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "index_of_dtlm_start" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'index_of_highest_temperature_inversion', physics%Interstitial(i)%kinver, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "index_of_highest_temperature_inversion" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'index_of_time_step', physics%Model(i)%kdt, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "index_of_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_atmosphere_detrainment_convective_mass_flux', physics%Interstitial(i)%dt_mf, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_atmosphere_detrainment_convective_mass_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_atmosphere_downdraft_convective_mass_flux', physics%Interstitial(i)%dd_mf, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_atmosphere_downdraft_convective_mass_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_atmosphere_updraft_convective_mass_flux', physics%Interstitial(i)%ud_mf, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_atmosphere_updraft_convective_mass_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_cosine_of_zenith_angle', physics%Interstitial(i)%xcosz, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_cosine_of_zenith_angle" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_surface_upward_latent_heat_flux', physics%Interstitial(i)%dqsfc1, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_surface_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_surface_upward_sensible_heat_flux', physics%Interstitial(i)%dtsfc1, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_surface_upward_sensible_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_surface_x_momentum_flux', physics%Interstitial(i)%dusfc1, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_surface_x_momentum_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_surface_y_momentum_flux', physics%Interstitial(i)%dvsfc1, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_surface_y_momentum_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_x_stress_due_to_gravity_wave_drag', physics%Interstitial(i)%dusfcg, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_x_stress_due_to_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'instantaneous_y_stress_due_to_gravity_wave_drag', physics%Interstitial(i)%dvsfcg, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "instantaneous_y_stress_due_to_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'iteration_number', physics%Interstitial(i)%iter, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "iteration_number" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'kinematic_surface_upward_latent_heat_flux', physics%Interstitial(i)%evap, ierr, 'kg kg-1 m s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "kinematic_surface_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'kinematic_surface_upward_sensible_heat_flux', physics%Interstitial(i)%hflx, ierr, 'K m s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "kinematic_surface_upward_sensible_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'large_scale_condensate_heating_rate_at_model_layers', physics%Diag(i)%dt3dt(:,:,6), ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "large_scale_condensate_heating_rate_at_model_layers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'large_scale_condensate_moistening_rate_at_model_layers', physics%Diag(i)%dq3dt(:,:,4), ierr, 'kg kg-1 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "large_scale_condensate_moistening_rate_at_model_layers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'largest_cloud_top_vertical_index_encountered_thus_far', physics%Tbd(i)%acvt, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "largest_cloud_top_vertical_index_encountered_thus_far" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'latitude_index_in_debug_printouts', physics%Interstitial(i)%latidxprnt, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "latitude_index_in_debug_printouts" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'latitude_of_ozone_forcing_data_from_host', physics%ozone_lat, ierr, 'degree')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "latitude_of_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'longitude', physics%Grid(i)%xlon, ierr, 'radians')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "longitude" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lw_fluxes_sfc', '', c_loc(physics%Radtend(i)%sfcflw), rank=size(shape(physics%Radtend(i)%sfcflw)), dims=shape(physics%Radtend(i)%sfcflw), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lw_fluxes_sfc" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lw_fluxes_top_atmosphere', '', c_loc(physics%Diag(i)%topflw), rank=size(shape(physics%Diag(i)%topflw)), dims=shape(physics%Diag(i)%topflw), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lw_fluxes_top_atmosphere" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep', physics%Diag(i)%rainc, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lwe_thickness_of_deep_convective_precipitation_amount', physics%Interstitial(i)%raincd, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lwe_thickness_of_deep_convective_precipitation_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lwe_thickness_of_precipitation_amount_on_dynamics_timestep', physics%Diag(i)%rain, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lwe_thickness_of_precipitation_amount_on_dynamics_timestep" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lwe_thickness_of_shallow_convective_precipitation_amount', physics%Interstitial(i)%raincs, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lwe_thickness_of_shallow_convective_precipitation_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'lwe_thickness_of_stratiform_precipitation_amount', physics%Interstitial(i)%rainst, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "lwe_thickness_of_stratiform_precipitation_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'maximum_column_heating_rate', physics%Interstitial(i)%cumabs, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "maximum_column_heating_rate" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'maximum_subgrid_orography', physics%Interstitial(i)%elvmax, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "maximum_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'maximum_vegetation_area_fraction', physics%Sfcprop(i)%shdmax, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "maximum_vegetation_area_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mean_change_over_depth_in_sea_water_temperature', physics%Interstitial(i)%dtzm, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mean_change_over_depth_in_sea_water_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mean_effective_radius_for_ice_cloud', physics%Interstitial(i)%clouds(:,:,5), ierr, 'micron')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mean_effective_radius_for_ice_cloud" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mean_effective_radius_for_liquid_cloud', physics%Interstitial(i)%clouds(:,:,3), ierr, 'micron')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mean_effective_radius_for_liquid_cloud" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mean_effective_radius_for_rain_drop', physics%Interstitial(i)%clouds(:,:,7), ierr, 'micron')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mean_effective_radius_for_rain_drop" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mean_effective_radius_for_snow_flake', physics%Interstitial(i)%clouds(:,:,9), ierr, 'micron')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mean_effective_radius_for_snow_flake" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'minimum_vegetation_area_fraction', physics%Sfcprop(i)%shdmin, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "minimum_vegetation_area_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'model_layer_number_at_cloud_base', physics%Interstitial(i)%mbota, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "model_layer_number_at_cloud_base" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'model_layer_number_at_cloud_top', physics%Interstitial(i)%mtopa, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "model_layer_number_at_cloud_top" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'mpi_rank', physics%Model(i)%me, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "mpi_rank" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'multiplication_factors_for_convective_gravity_wave_drag', physics%Model(i)%cgwf, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "multiplication_factors_for_convective_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag', physics%Model(i)%cdmbgwd, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'natural_log_of_ozone_forcing_data_pressure_levels', physics%Interstitial(i)%oz_pres, ierr, 'log(Pa)')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "natural_log_of_ozone_forcing_data_pressure_levels" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'natural_log_of_ozone_forcing_data_pressure_levels_from_host', physics%ozone_pres, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "natural_log_of_ozone_forcing_data_pressure_levels_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep', physics%Sfcprop(i)%tprcp, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'normalized_soil_wetness', physics%Diag(i)%wet1, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "normalized_soil_wetness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_coefficients_in_ozone_forcing_data', physics%Interstitial(i)%oz_coeff, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_coefficients_in_ozone_forcing_data" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_coefficients_in_ozone_forcing_data_from_host', physics%n_ozone_coefficients, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_coefficients_in_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_equatorial_longitude_points', physics%Model(i)%lonr, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_equatorial_longitude_points" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_hydrometeors', physics%Model(i)%ncld, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_hydrometeors" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_latitutde_points_in_ozone_forcing_data_from_host', physics%n_ozone_lats, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_latitutde_points_in_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_statistical_measures_of_subgrid_orography', physics%Model(i)%nmtvr, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_statistical_measures_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_time_levels_in_ozone_forcing_data_from_host', physics%n_ozone_times, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_time_levels_in_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_vertical_diffusion_tracers', physics%Interstitial(i)%nvdiff, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_vertical_diffusion_tracers" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'number_of_vertical_layers_for_radiation_calculations', physics%Model(i)%levr, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "number_of_vertical_layers_for_radiation_calculations" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ocean_mixed_layer_thickness', physics%Sfcprop(i)%zm, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ocean_mixed_layer_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'omega', physics%Statein(i)%vvl, ierr, 'Pa s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "omega" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'orography', physics%Sfcprop(i)%oro, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'orography_unfiltered', physics%Sfcprop(i)%oro_uf, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "orography_unfiltered" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ozone_concentration_at_layer_for_radiation', physics%Interstitial(i)%olyr, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ozone_concentration_at_layer_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ozone_concentration_updated_by_physics', physics%Stateout(i)%gq0(:,:,scm_state%ozone_index), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ozone_concentration_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ozone_forcing', physics%Tbd(i)%ozpl, ierr, 'various')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ozone_forcing" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ozone_forcing_from_host', physics%ozone_forcing_in, ierr, 'various')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ozone_forcing_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'pi', con_pi, ierr, 'radians')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "pi" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'pressure_at_bottom_of_convective_cloud', physics%Cldprop(i)%cvb, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "pressure_at_bottom_of_convective_cloud" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'pressure_at_top_of_convective_cloud', physics%Cldprop(i)%cvt, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "pressure_at_top_of_convective_cloud" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'pressure_cutoff_for_rayleigh_damping', physics%Model(i)%prslrd0, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "pressure_cutoff_for_rayleigh_damping" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'random_number_array', physics%Tbd(i)%rann, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "random_number_array" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer', physics%Interstitial(i)%work3, ierr, 'ratio')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ratio_of_snowfall_to_rainfall', physics%Diag(i)%sr, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ratio_of_snowfall_to_rainfall" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ratio_of_vapor_to_dry_air_gas_constants_minus_one', con_fvirt, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ratio_of_vapor_to_dry_air_gas_constants_minus_one" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'ratio_of_wind_at_lowest_model_layer_and_wind_at_10m', physics%Sfcprop(i)%f10m, ierr, 'ratio')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "ratio_of_wind_at_lowest_model_layer_and_wind_at_10m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_concentration', physics%Sfcprop(i)%fice, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_concentration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_concentration_for_physics', physics%Interstitial(i)%cice, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_concentration_for_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_temperature', physics%Sfcprop(i)%tisfc, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_temperature_for_physics', physics%Interstitial(i)%tice, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_temperature_for_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_thickness', physics%Sfcprop(i)%hice, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_ice_thickness_for_physics', physics%Interstitial(i)%zice, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_ice_thickness_for_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_land_ice_mask', physics%Interstitial(i)%islmsk, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_land_ice_mask" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_land_ice_mask_real', physics%Sfcprop(i)%slmsk, ierr, 'flag')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_land_ice_mask_real" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_surface_reference_temperature', physics%Sfcprop(i)%tref, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_surface_reference_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sea_water_salinity', physics%Sfcprop(i)%xs, ierr, 'ppt m')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sea_water_salinity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'seed_random_numbers_lw', physics%Tbd(i)%icsdlw, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "seed_random_numbers_lw" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'seed_random_numbers_sw', physics%Tbd(i)%icsdsw, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "seed_random_numbers_sw" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sensible_heat_flux_due_to_rainfall', physics%Sfcprop(i)%qrain, ierr, 'W')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sensible_heat_flux_due_to_rainfall" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sensitivity_of_dtl_heat_content_to_surface_temperature', physics%Sfcprop(i)%xtts, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sensitivity_of_dtl_heat_content_to_surface_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sensitivity_of_dtl_thickness_to_surface_temperature', physics%Sfcprop(i)%xzts, ierr, 'm K-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sensitivity_of_dtl_thickness_to_surface_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sine_of_latitude', physics%Grid(i)%sinlat, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sine_of_latitude" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sine_of_solar_declination_angle', physics%Model(i)%sdec, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sine_of_solar_declination_angle" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'slope_of_subgrid_orography', physics%Interstitial(i)%sigma, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "slope_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'smallest_cloud_base_vertical_index_encountered_thus_far', physics%Tbd(i)%acvb, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "smallest_cloud_base_vertical_index_encountered_thus_far" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'snow_deposition_sublimation_upward_latent_heat_flux', physics%Interstitial(i)%sbsno, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "snow_deposition_sublimation_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'snow_freezing_rain_upward_latent_heat_flux', physics%Interstitial(i)%snohf, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "snow_freezing_rain_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'soil_moisture_content', physics%Diag(i)%soilm, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "soil_moisture_content" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'soil_temperature', physics%Sfcprop(i)%stc, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "soil_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'soil_type', physics%Model(i)%isot, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "soil_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'soil_upward_latent_heat_flux', physics%Interstitial(i)%evbs, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "soil_upward_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'soil_vertical_dimension', physics%Model(i)%lsoil, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "soil_vertical_dimension" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'solar_constant', physics%Model(i)%solcon, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "solar_constant" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'specific_heat_of_dry_air_at_constant_pressure', con_cp, ierr, 'J kg-1 K-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "specific_heat_of_dry_air_at_constant_pressure" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'specific_humidity_at_2m', physics%Sfcprop(i)%q2m, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "specific_humidity_at_2m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'specific_humidity_at_lowest_model_layer', physics%Statein(i)%qgrs(:,1,scm_state%water_vapor_index), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "specific_humidity_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'standard_deviation_of_subgrid_orography', physics%Interstitial(i)%hprime1, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "standard_deviation_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'statistical_measures_of_subgrid_orography', physics%Sfcprop(i)%hprime, ierr, 'various')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "statistical_measures_of_subgrid_orography" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sub-layer_cooling_amount', physics%Sfcprop(i)%dt_cool, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sub-layer_cooling_amount" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sub-layer_cooling_thickness', physics%Sfcprop(i)%z_c, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sub-layer_cooling_thickness" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'subsurface_runoff_flux', physics%Interstitial(i)%drain, ierr, 'g m-2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "subsurface_runoff_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_air_pressure', physics%Statein(i)%pgr, ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_air_pressure" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_air_pressure_at_previous_time_step', physics%Tbd(i)%phy_f2d(:,2), ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_air_pressure_at_previous_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_air_pressure_two_time_steps_back', physics%Tbd(i)%phy_f2d(:,1), ierr, 'Pa')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_air_pressure_two_time_steps_back" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_air_temperature_for_radiation', physics%Interstitial(i)%tsfa, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_air_temperature_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_albedo_due_to_UV_and_VIS_diffused', physics%Interstitial(i)%sfcalb(:,4), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_albedo_due_to_UV_and_VIS_diffused" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_albedo_due_to_UV_and_VIS_direct', physics%Interstitial(i)%sfcalb(:,3), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_albedo_due_to_UV_and_VIS_direct" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_albedo_due_to_near_IR_diffused', physics%Interstitial(i)%sfcalb(:,2), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_albedo_due_to_near_IR_diffused" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_albedo_due_to_near_IR_direct', physics%Interstitial(i)%sfcalb(:,1), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_albedo_due_to_near_IR_direct" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_diffused_shortwave_albedo', physics%Radtend(i)%sfalb, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_diffused_shortwave_albedo" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_diffuse_near_infrared_shortwave_flux', physics%Interstitial(i)%adjnirdfd, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_diffuse_near_infrared_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%nirdfdi, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux', physics%Interstitial(i)%adjvisdfd, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%visdfdi, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_direct_near_infrared_shortwave_flux', physics%Interstitial(i)%adjnirbmd, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_direct_near_infrared_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%nirbmdi, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux', physics%Interstitial(i)%adjvisbmd, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%visbmdi, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_longwave_flux', physics%Interstitial(i)%adjsfcdlw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_longwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_longwave_flux_absorbed_by_ground', physics%Interstitial(i)%gabsbdlw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_longwave_flux_absorbed_by_ground" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_longwave_flux_on_radiation_time_step', physics%Coupling(i)%sfcdlw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_longwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_shortwave_flux', physics%Interstitial(i)%adjsfcdsw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_downwelling_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%sfcdsw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_downwelling_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_drag_coefficient_for_heat_and_moisture_in_air', physics%Interstitial(i)%cdq, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_drag_coefficient_for_heat_and_moisture_in_air" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_drag_coefficient_for_momentum_in_air', physics%Interstitial(i)%cd, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_drag_coefficient_for_momentum_in_air" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_drag_mass_flux_for_heat_and_moisture_in_air', physics%Diag(i)%chh, ierr, 'kg m-2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_drag_mass_flux_for_heat_and_moisture_in_air" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_drag_wind_speed_for_momentum_in_air', physics%Diag(i)%cmm, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_drag_wind_speed_for_momentum_in_air" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_friction_velocity', physics%Sfcprop(i)%uustar, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_friction_velocity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_ground_temperature_for_radiation', physics%Interstitial(i)%tsfg, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_ground_temperature_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_longwave_emissivity', physics%Radtend(i)%semis, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_longwave_emissivity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_midlayer_air_temperature_in_longwave_radiation', physics%Radtend(i)%tsflw, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_midlayer_air_temperature_in_longwave_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_net_downwelling_shortwave_flux', physics%Interstitial(i)%adjsfcnsw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_net_downwelling_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_net_downwelling_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%sfcnsw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_net_downwelling_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_roughness_length', physics%Sfcprop(i)%zorl, ierr, 'cm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_roughness_length" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_runoff', physics%Diag(i)%srunoff, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_runoff" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_runoff_flux', physics%Interstitial(i)%runoff, ierr, 'g m-2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_runoff_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_skin_temperature', physics%Sfcprop(i)%tsfc, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_skin_temperature" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_skin_temperature_after_iteration', physics%Interstitial(i)%tsurf, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_skin_temperature_after_iteration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_skin_temperature_for_nsst', physics%Interstitial(i)%tseal, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_skin_temperature_for_nsst" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_slope_classification', physics%Interstitial(i)%slopetype, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_slope_classification" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_snow_area_fraction', physics%Interstitial(i)%snowc, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_snow_area_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_snow_area_fraction_for_diagnostics', physics%Sfcprop(i)%sncovr, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_snow_area_fraction_for_diagnostics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_snow_melt', physics%Interstitial(i)%snowmt, ierr, 'm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_snow_melt" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_snow_thickness_water_equivalent', physics%Sfcprop(i)%snowd, ierr, 'mm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_snow_thickness_water_equivalent" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_specific_humidity', physics%Interstitial(i)%qss, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_specific_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upward_potential_latent_heat_flux', physics%Interstitial(i)%ep1d, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upward_potential_latent_heat_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_diffuse_near_infrared_shortwave_flux', physics%Interstitial(i)%adjnirdfu, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_diffuse_near_infrared_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%nirdfui, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux', physics%Interstitial(i)%adjvisdfu, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%visdfui, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_direct_near_infrared_shortwave_flux', physics%Interstitial(i)%adjnirbmu, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_direct_near_infrared_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%nirbmui, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux', physics%Interstitial(i)%adjvisbmu, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', physics%Coupling(i)%visbmui, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_upwelling_longwave_flux', physics%Interstitial(i)%adjsfculw, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_upwelling_longwave_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_wind_enhancement_due_to_convection', physics%Tbd(i)%phy_f2d(:,physics%Model(i)%num_p2d), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_wind_enhancement_due_to_convection" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'surface_wind_stress', physics%Interstitial(i)%stress, ierr, 'm2 s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "surface_wind_stress" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sw_fluxes_sfc', '', c_loc(physics%Radtend(i)%sfcfsw), rank=size(shape(physics%Radtend(i)%sfcfsw)), dims=shape(physics%Radtend(i)%sfcfsw), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sw_fluxes_sfc" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'sw_fluxes_top_atmosphere', '', c_loc(physics%Diag(i)%topfsw), rank=size(shape(physics%Diag(i)%topfsw)), dims=shape(physics%Diag(i)%topfsw), ierr=ierr)
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "sw_fluxes_top_atmosphere" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'temperature_at_2m', physics%Sfcprop(i)%t2m, ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "temperature_at_2m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step', physics%Tbd(i)%htlw0, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step', physics%Tbd(i)%htlwc, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_model_physics', physics%Interstitial(i)%dtdt, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_model_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky', physics%Interstitial(i)%dtdtc, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step', physics%Tbd(i)%htsw0, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step', physics%Tbd(i)%htswc, ierr, 'K s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_rain_water_mixing_ratio_due_to_model_physics', physics%Interstitial(i)%rainp, ierr, 'kg kg-1 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_rain_water_mixing_ratio_due_to_model_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_tracers_due_to_model_physics', physics%Interstitial(i)%dqdt, ierr, 'kg kg-1 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_tracers_due_to_model_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_x_wind_due_to_convective_gravity_wave_drag', physics%Interstitial(i)%gwdcu, ierr, 'm s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_x_wind_due_to_convective_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_x_wind_due_to_model_physics', physics%Interstitial(i)%dudt, ierr, 'm s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_x_wind_due_to_model_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_y_wind_due_to_convective_gravity_wave_drag', physics%Interstitial(i)%gwdcv, ierr, 'm s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_y_wind_due_to_convective_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tendency_of_y_wind_due_to_model_physics', physics%Interstitial(i)%dvdt, ierr, 'm s-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tendency_of_y_wind_due_to_model_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'threshold_volume_fraction_of_condensed_water_in_soil', physics%Diag(i)%smcref2, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "threshold_volume_fraction_of_condensed_water_in_soil" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_integral_of_x_stress_due_to_gravity_wave_drag', physics%Diag(i)%dugwd, ierr, 'Pa s')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_integral_of_x_stress_due_to_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_integral_of_y_stress_due_to_gravity_wave_drag', physics%Diag(i)%dvgwd, ierr, 'Pa s')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_integral_of_y_stress_due_to_gravity_wave_drag" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_levels_in_ozone_forcing_data_from_host', physics%ozone_time, ierr, 'day')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_levels_in_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_scale_for_rayleigh_damping', physics%Model(i)%ral_ts, ierr, 'd')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_scale_for_rayleigh_damping" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_step_for_dynamics', physics%Model(i)%dtf, ierr, 's')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_step_for_dynamics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_step_for_physics', physics%Model(i)%dtp, ierr, 's')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_step_for_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'time_step_for_radiation', physics%Interstitial(i)%raddt, ierr, 's')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "time_step_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'total_cloud_fraction', physics%Interstitial(i)%clouds(:,:,1), ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "total_cloud_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'total_runoff', physics%Diag(i)%runoff, ierr, 'kg m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "total_runoff" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'tracer_concentration', physics%Statein(i)%qgrs, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "tracer_concentration" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'transpiration_flux', physics%Interstitial(i)%trans, ierr, 'kg m-2 s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "transpiration_flux" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'upper_bound_on_max_albedo_over_deep_snow', physics%Sfcprop(i)%snoalb, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "upper_bound_on_max_albedo_over_deep_snow" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'upward_heat_flux_in_soil', physics%Interstitial(i)%gflx, ierr, 'W m-2')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "upward_heat_flux_in_soil" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vegetation_area_fraction', physics%Interstitial(i)%sigmaf, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vegetation_area_fraction" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vegetation_type', physics%Model(i)%ivegsrc, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vegetation_type" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_dimension', physics%Model(i)%levs, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_dimension" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_dimension_of_ozone_forcing_data', physics%Interstitial(i)%levozp, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_dimension_of_ozone_forcing_data" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_dimension_of_ozone_forcing_data_from_host', physics%n_ozone_layers, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_dimension_of_ozone_forcing_data_from_host" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_at_cloud_base', physics%Interstitial(i)%kbot, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_at_cloud_base" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_at_cloud_top', physics%Interstitial(i)%ktop, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_at_cloud_top" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_at_top_of_atmosphere_boundary_layer', physics%Interstitial(i)%kpbl, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_at_top_of_atmosphere_boundary_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_difference_between_inout_and_local', physics%Interstitial(i)%kd, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_difference_between_inout_and_local" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_difference_between_layer_and_lower_bound', physics%Interstitial(i)%kb, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_difference_between_layer_and_lower_bound" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_index_difference_between_layer_and_upper_bound', physics%Interstitial(i)%kt, ierr, 'index')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_index_difference_between_layer_and_upper_bound" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_interface_dimension', physics%Interstitial(i)%levi, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_interface_dimension" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_layer_dimension_for_radiation', physics%Interstitial(i)%lm, ierr, 'count')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_layer_dimension_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_temperature_average_range_lower_bound', physics%Model(i)%nstf_name(4), ierr, 'mm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_temperature_average_range_lower_bound" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'vertical_temperature_average_range_upper_bound', physics%Model(i)%nstf_name(5), ierr, 'mm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "vertical_temperature_average_range_upper_bound" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_fraction_of_condensed_water_in_soil_at_wilting_point', physics%Diag(i)%smcwlt2, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_fraction_of_condensed_water_in_soil_at_wilting_point" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_fraction_of_soil_moisture', physics%Sfcprop(i)%smc, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_fraction_of_soil_moisture" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_fraction_of_unfrozen_soil_moisture', physics%Sfcprop(i)%slc, ierr, 'frac')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_fraction_of_unfrozen_soil_moisture" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_ccl4', physics%Interstitial(i)%gasvmr(:,:,9), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_ccl4" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_cfc11', physics%Interstitial(i)%gasvmr(:,:,6), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_cfc11" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_cfc113', physics%Interstitial(i)%gasvmr(:,:,10), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_cfc113" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_cfc12', physics%Interstitial(i)%gasvmr(:,:,7), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_cfc12" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_cfc22', physics%Interstitial(i)%gasvmr(:,:,8), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_cfc22" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_ch4', physics%Interstitial(i)%gasvmr(:,:,3), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_ch4" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_co', physics%Interstitial(i)%gasvmr(:,:,5), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_co" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_co2', physics%Interstitial(i)%gasvmr(:,:,1), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_co2" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_n2o', physics%Interstitial(i)%gasvmr(:,:,2), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_n2o" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'volume_mixing_ratio_o2', physics%Interstitial(i)%gasvmr(:,:,4), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "volume_mixing_ratio_o2" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_equivalent_accumulated_snow_depth', physics%Sfcprop(i)%weasd, ierr, 'mm')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_equivalent_accumulated_snow_depth" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity', physics%Statein(i)%qgrs(:,:,scm_state%water_vapor_index), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity_at_layer_for_radiation', physics%Interstitial(i)%qlyr, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity_at_layer_for_radiation" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity_at_previous_time_step', physics%Tbd(i)%phy_f3d(:,:,4), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity_at_previous_time_step" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity_save', physics%Interstitial(i)%save_qv, ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity_save" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity_two_time_steps_back', physics%Tbd(i)%phy_f3d(:,:,2), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity_two_time_steps_back" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'water_vapor_specific_humidity_updated_by_physics', physics%Stateout(i)%gq0(:,:,scm_state%water_vapor_index), ierr, 'kg kg-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "water_vapor_specific_humidity_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'wind_speed_at_lowest_model_layer', physics%Interstitial(i)%wind, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "wind_speed_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'x_wind', physics%Statein(i)%ugrs, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'x_wind_at_10m', physics%Diag(i)%u10m, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind_at_10m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'x_wind_at_lowest_model_layer', physics%Statein(i)%ugrs(:,1), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'x_wind_save', physics%Interstitial(i)%save_u, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind_save" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'x_wind_updated_by_physics', physics%Stateout(i)%gu0, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'y_wind', physics%Statein(i)%vgrs, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'y_wind_at_10m', physics%Diag(i)%v10m, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind_at_10m" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'y_wind_at_lowest_model_layer', physics%Statein(i)%vgrs(:,1), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind_at_lowest_model_layer" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'y_wind_save', physics%Interstitial(i)%save_v, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind_save" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'y_wind_updated_by_physics', physics%Stateout(i)%gv0, ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind_updated_by_physics" to CCPP data structure')
            end if
            call ccpp_field_add(cdata(i), 'zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes', physics%Interstitial(i)%xmu, ierr, 'none')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes" to CCPP data structure')
            end if