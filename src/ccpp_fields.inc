
            call ccpp_fields_add(cdata(i), 'air_temperature', scm_state%state_T(i,:,cdata_time_index), ierr, 'K')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "air_temperature" to CCPP data structure')
                return
            end if
            call ccpp_fields_add(cdata(i), 'x_wind', scm_state%state_u(i,:,cdata_time_index), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "x_wind" to CCPP data structure')
                return
            end if
            call ccpp_fields_add(cdata(i), 'y_wind', scm_state%state_v(i,:,cdata_time_index), ierr, 'm s-1')
            if (ierr /= 0) then
                call ccpp_error('Unable to add field "y_wind" to CCPP data structure')
                return
            end if