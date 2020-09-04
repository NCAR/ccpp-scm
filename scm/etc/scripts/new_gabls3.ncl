load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;Define constants

  g     = 9.80665 ;gravity (m s-2)
  R_dry = 287.    ;ideal gas dry air constant (J kg-1 K-1)
  R_vap = 461.5   ;gas constant for water vapor
  c_p   = 1004.   ;specific heat at const pressure (J kg-1 K-1)
  l_v   = 2.5e6   ;latent heat of vaporization (J kg-1)
  l_s   = 2.836e6 ;latent heat of sublimation (J kg-1)
  P0    = 100000. ;intial pressure (Pa)
  kappa = R_dry/c_p

; GABLS initial vertical profiles and forcings. For levels with no entries,
; values are linearly interpolated based on nearby points
; Data source: http://projects.knmi.nl/gabls/setup.html

; Time, 01Jul2006 12:00:00 UTC - 02Jul2006 12:00:00 UTC
; Must be in seconds
  time           = fspan(0,86400,145)
  time@long_name = "elapsed time since the beginning of the simulation" 
  time@units     = "s"

;  timb = fspan(0,86400,145)
;  print(timb)

; Lat and Lon of the column
  lat = 51.9711
  lat@long_name = "latitude of column"
  lat@units = "degrees N"
  lon = 4.9267
  lon@long_name = "longitude of column"
  lon@units = "degrees E"

; Height levels for initial values of temperature and specific humidity
; Other variables will need to be interpolated
  height = (/ 0., 2., 10., 20., 40., 80., 140., 205., 1800., 2200., \
              5000., 12000., 14000., 50000. /)
  height@long_name = "physical height at pressure levels"
  height@units     = "m"

; Total water specific humidity (q) at heights listed above
  qt = (/ 9.3e-3, 9.3e-3, 8.5e-3, 8.4e-3, 8.3e-3, 8.2e-3, 8.1e-3, \
          8.0e-3, 7.5e-3, 2.0e-3, 0.3e-3, 0.01e-3, 0.003e-3, 0. /) 
  qt@long_name = "initial profile of total water specific humidity"
  qt@units     = "kg kg^-1"

; Converting levels (m) to pressures (Pa) using Exner eqn
; Calculate Exner and pressure at surface
  p_surf_hpa = 1024.4 ;(hPa, given on GABLS setup link)
  p_surf     = new((/ dimsizes(time) /),float)
  p_surf(0)  = p_surf_hpa * 100. ;now in Pa, rather than hPa

  do i=1,dimsizes(time)-1,1
    p_surf(i) = p_surf(0)
  end do

  p_surf@long_name = "surface pressure"
  p_surf@units     = "Pa"

  levels = new((/ dimsizes(height) /),float)
  levels@long_name = "pressure levels"
  levels@units     = "Pa"

; Air temperature profile (not needed in netCDF input file, but necessary
; for calculating other variables that gmtb-scm will find necessary)
  T_C = (/ 27.0, 27.0, 26.4, 26.2, 25.9, 25.5, 24.8, 24.3, 9.0, 9.0, \
           -6.4, -61., -54., -50. /)
  T_K = T_C + 273.15 ;convert from degC to Kelvin

  T_surf = new((/ dimsizes(time) /),float)
  T_surf = 23.4 + 273.15

;  do i=0,dimsizes(time)-1,1
;    T_surf(i) = T_K(0)
;  end do
  T_surf@long_name = "surface absolute temperature"
  T_surf@units     = "K"

; Calculate approximate virtual temperature to solve hypsometric equation
; for the model required pressure levels
  T_V = new ((/ dimsizes(levels) /),float)
  do i=0,dimsizes(levels)-1,1
    T_V(i) = T_K(i)*((1 + (qt(i)/0.622))/(1+qt(i)))
  end do

; Solve the hypsometric equation for p2 to attain each new pressure level
; and fill the array (the variable to fill is levels, created above.) 
; The value for levels(0) = p_surf(0)
  levels(0) = p_surf(0)
  do i=1,dimsizes(levels)-1,1
    levels(i) = levels(i-1) / (exp((g*(height(i)-height(i-1)))/(R_dry \
                * T_V(i))))
  end do

; Calculate potential temperature (theta) in order to calculate the
; necessary thetail (ice-liquid water potential temperature)
  theta = new((/ dimsizes(levels) /),float)
  do i=0,dimsizes(levels)-1,1
    theta(i) = T_K(i)*((P0/levels(i))^(kappa))
  end do

; Calculate thetail. Because ql and qi are not given, this is not possible,
; so the initial profile of thetail will be assumed theta.
; The remainder of the equation used in the gmtb-scm user's guide would be 
; zero, leaving you with only theta.
  thetail           = theta
  thetail@long_name = "initial profile of ice-liquid water potential "+ \
                      "temperature"
  thetail@units     = "K"

; Horizontal motion (u and v winds), needs to be interpolated
  u_specified = (/ 0.0, -4.0, -5.5, -5.5, -2.0, -2.0 /)
  u_h         = (/ 0., 10., 353., 1238., 2000., 50000. /)

  u = linint1(u_h,u_specified,False,height,0)
  u@long_name = "initial profile of E-W horizontal wind"
  u@units     = "m s^-1"

  v_specified = (/ 0.0, -0.4, -0.5, -0.5, 2., 2. /)
  v_h         = u_h

  v = linint1(v_h,v_specified,False,height,0)
  v@long_name = "initial profile of N-S horizontal wind"
  v@units     = "m s^-1"

; Use GLDAS to derive the sensible and latent heat flux for the Cabauw 
; location, then you can compare the GLDAS H/Le to the given intial Bowen 
; ratio. Begin by reading in the data from each file.
  gl_dir = "/glade/scratch/damico/GLDAS_GABLS3/"
  gl_00  = addfile(gl_dir+"GLDAS_NOAH025_3H.A20060701.0000.021.nc4.SUB.nc4", \
                   "r")

  gl_lat = gl_00->lat
  gl_lon = gl_00->lon

  glny = dimsizes(gl_lat)
  glnx = dimsizes(gl_lon)

  delete(gl_00)

  start_date = 0.
  end_date   = 0.

  start_date@units = "hours since 2006-07-01 12:00:00"
  end_date@units   = "hours since 2006-07-02 12:00:00"
  gldt             = 3

  end_date_on_start = cd_convert(end_date,start_date@units)

  gl_time_hr       = ispan(0,toint(end_date_on_start),gldt)
  gl_time_hr@units = start_date@units

  gl_time       = tofloat(gl_time_hr) * 60 * 60
  gl_time@units = "seconds since 2006-07-01 12:00:00"

  time_yyyymmddhh = toint(cd_calendar(gl_time,0))

  glnt = dimsizes(gl_time_hr)

  gl_shf = new((/ glnt,glny,glnx /),float)
  gl_lhf = new((/ glnt,glny,glnx /),float)
  do i=0,glnt-1,1
    fname = gl_dir + "GLDAS_NOAH025_3H.A" + \
            sprinti("%0.4i",time_yyyymmddhh(i,0)) + \
            sprinti("%0.2i",time_yyyymmddhh(i,1)) + \
            sprinti("%0.2i",time_yyyymmddhh(i,2)) +"."+ \
            sprinti("%0.2i",time_yyyymmddhh(i,3))+ \
            sprinti("%0.2i",time_yyyymmddhh(i,4))+ \
            ".021.nc4.SUB.nc4"
    if (isfilepresent(fname)) then
      gl_ff         = addfile(fname,"r")
      gl_shf(i,:,:) = gl_ff->Qh_tavg(0,:,:)
      gl_lhf(i,:,:) = gl_ff->Qle_tavg(0,:,:)
      print("Read->"+fname)
    else
      print(fname+" MISSING")
    end if
  end do

; Interpolate GLDAS shf and lhf to the single Cabauw point.
  shf_pt = linint2_points_Wrap(gl_lon,gl_lat,gl_shf,False,lon,lat,0)
  lhf_pt = linint2_points_Wrap(gl_lon,gl_lat,gl_lhf,False,lon,lat,0)

; Interpolate with time
  shf = linint1(gl_time,shf_pt(:,0),False,time,0)
  lhf = linint1(gl_time,lhf_pt(:,0),False,time,0)

; Convert from W m^-2 to K m s^-1 and kg kg^-1 m s^-1
  sh_flux_sfc = shf * R_dry * T_surf / (c_p * p_surf)
  lh_flux_sfc = lhf * R_dry * T_surf / (l_v * p_surf)

  sh_flux_sfc@long_name = "surface sensible heat flux"
  sh_flux_sfc@units     = "K m s^-1"
  lh_flux_sfc@long_name = "surface latent heat flux"
  lh_flux_sfc@units     = "kg kg^-1 m s^-1"

; Calculate vertical motion, needed in both w and omega, given in omega. Omega
; begins at 0.12 Pa s^-1 at 1500 m at the beginning of the study period. After
; five hours, omega starts to decrease to 0, ending at 0 after seven hours.
; Omega is measured at 1500 m, and decreases from that value to 0 m. It also
; decreases from 0.12 to 5000 m, above which the value is zero.
  om_specified = (/ 0.12, 0.12,   0.,     0. /)
  omega_t      = (/    0,18000,25200, 86400  /)
  om           = linint1(omega_t,om_specified,False,time,0)

  omega     = new((/ dimsizes(levels),dimsizes(time) /),float)
  omega_h   = (/ 0.,1500.,5000.,50000. /)
  omega_int = new((/ dimsizes(omega_h),dimsizes(time) /),float)
  do i=0,dimsizes(time)-1,1
    omega_int(:,i) = (/ 0.,om(i),0.,0. /)
    omega(:,i)     = linint1(omega_h,omega_int(:,i),False,height,0)
  end do

  omega@long_name = "large scale pressure vertical velocity"
  omega@units     = "Pa s^-1"

; Convert omega into w, also required input
; Begin by solving for rho
  rho = new((/ dimsizes(levels) /),float)
  do i=0,dimsizes(levels)-1,1
    rho(i) = levels(i)/(R_dry*T_K(i))
  end do

; Solve for w (equation commented below)
; w = -om / (rho*g)
  w_ls = new((/ dimsizes(levels),dimsizes(time) /),float)
  do i=0,dimsizes(levels)-1,1
    do j=0,dimsizes(time)-1,1
      w_ls(i,j) = -omega(i,j) / (rho(i)*g)
    end do
  end do
  w_ls@long_name = "large scale vertical velocity"
  w_ls@units     = "m s^-1"

; Geostrophic winds are given at certain times but only at the surface.
; Instructions from GABLS suggest interpolation to surface winds up to 
; 2000 m, and then constant from that point. This should be done for each
; time step. Starting with time.
  u_geo_spec = (/ -7.8, -7.8, -6.5, -5.0, -5.0, -6.5 /)
  u_geo_t    = (/ 0, 21600, 39600, 54000, 64800, 86400 /)
  u_geo      = linint1(u_geo_t,u_geo_spec,False,time,0)

  v_geo_spec = (/ 0.0, 0.0, 4.5, 4.5, 4.5, 2.5 /)
  v_geo_t    = (/ 0, 21600, 39600, 54000, 64800, 86400 /)
  v_geo      = linint1(v_geo_t,v_geo_spec,False,time,0)

; Interpolate with height, linearly, from the surface values above to -2.0 
; m s^-1 for u and 2.0 m s^-1 for v. Above 2000 m, geostrophic wind is
; constant with height (according to GABLS guide).
  u_g       = new((/ dimsizes(levels),dimsizes(time) /),float)
  u_geo_h   = (/ 0., 2000., 50000. /)
  u_geo_int = new((/ dimsizes(u_geo_h),dimsizes(time) /),float)

  do i=0,dimsizes(time)-1,1
    u_geo_int(:,i) = (/ u_geo(i), -2., -2. /)
    u_g(:,i)       = linint1(u_geo_h,u_geo_int(:,i),False,height,0)
  end do

  v_g       = new((/ dimsizes(levels),dimsizes(time) /),float)
  v_geo_int = new((/ dimsizes(u_geo_h),dimsizes(time) /),float)

  v_geo_h   = u_geo_h

  do i=0,dimsizes(time)-1,1
    v_geo_int(:,i) = (/ v_geo(i), 2., 2. /)
    v_g(:,i)       = linint1(v_geo_h,v_geo_int(:,i),False,height,0)
  end do

  u_g@long_name = "large scale geostrophic E-W wind"
  v_g@long_name = "large scale geostrophic N-S wind"
  u_g@units     = "m s^-1"
  v_g@unit      = "m s^-1"

; The last remaining variable given in the data that the gmtb-scm can 
; use is horizontal specific humidity tendency. The values are given 
; between 200-1000m, decreasing to zero down towards the surface and
; decreasing to zero from 1000m to 1500, with zero entirely above. Once 
; again, there needs to be interpolation in time and with height.
  h_advec_spec = (/ 0., 8.e-8, 0., -8.e-8, 0., 0. /)
  h_advec_time = (/ 0, 32400, 43200, 50400, 61200, 86400 /)
  h_advec      = linint1(h_advec_time,h_advec_spec,False,time,0)

; Interpolate with height
  h_advec_qt  = new((/ dimsizes(levels),dimsizes(time) /),float)
  h_advec_h   = (/ 0., 200., 1000., 1500., 50000. /)
  h_advec_int = new((/ dimsizes(h_advec_h),dimsizes(time) /),float)

  do i=0,dimsizes(time)-1,1
    h_advec_int(:,i) = (/ 0., h_advec(i), h_advec(i), 0., 0. /)
    h_advec_qt(:,i)  = linint1(h_advec_h,h_advec_int(:,i),False,height,0)
  end do

  h_advec_qt@long_name = "prescribed q_t tendency due to horizontal " + \
                         "advection"
  h_advec_qt@units     = "kg kg^-1 s^-1"

; The following variables are not in the GABLS data available or cannot
; be derived with the information given. They are still required by the
; model.
; Initial profiles
  ql           = new((/ dimsizes(levels) /),float)
  ql           = 0.
  ql@long_name = "initial profile of liquid water specific humidity"
  ql@units     = "kg kg^-1"

  qi           = new((/ dimsizes(levels) /),float)
  qi           = 0.
  qi@long_name = "initial profile of ice water specific humidity"
  qi@units     = "kg kg^-1"

  tke           = new((/ dimsizes(levels) /),float)
  tke           = 0.
  tke@long_name = "initial profile of turbulence kinetic energy"
  tke@units     = "m^2 s^-2"

  ozone           = new((/ dimsizes(levels) /),float)
  ozone           = 0.
  ozone@long_name = "initial profile of ozone mass mixing ratio"
  ozone@units     = "kg kg^-1"

; Forcing

  u_nudge           = new((/ dimsizes(levels),dimsizes(time) /),float)
  u_nudge           = 0.
  u_nudge@long_name = "E-W wind to nudge toward"
  u_nudge@units     = "m s^-1"

  v_nudge           = new((/ dimsizes(levels),dimsizes(time) /),float)
  v_nudge           = 0.
  v_nudge@long_name = "N-S wind to nudge toward"
  v_nudge@units     = "m s^-1"

  T_nudge           = new((/ dimsizes(levels),dimsizes(time) /),float)
  T_nudge           = 0.
  T_nudge@long_name = "absolute temperature to nudge toward"
  T_nudge@units     = "K"

  thil_nudge           = new((/ dimsizes(levels),dimsizes(time) /),float)
  thil_nudge           = 0.
  thil_nudge@long_name = "potential temperature to nudge toward"
  thil_nudge@units     = "K"

  qt_nudge           = new((/ dimsizes(levels),dimsizes(time) /),float)
  qt_nudge           = 0.
  qt_nudge@long_name = "q_t to nudge toward"
  qt_nudge@units     = "kg kg^-1"

  dT_dt_rad           = new((/ dimsizes(levels),dimsizes(time) /),float)
  dT_dt_rad           = 0.
  dT_dt_rad@long_name = "prescribed radiative heating rate"
  dT_dt_rad@units     = "K s^-1"

  h_advec_thetail           = new((/ dimsizes(levels),dimsizes(time) /), \
                              float)
  h_advec_thetail           = 0.
  h_advec_thetail@long_name = "prescribed theta_il tendency due to " + \
                              "horizontal advection"
  h_advec_thetail@units     = "K s^-1"

  v_advec_thetail           = new((/ dimsizes(levels),dimsizes(time) /), \
                              float)
  v_advec_thetail           = 0.
  v_advec_thetail@long_name = "prescribed theta_il tendency due to " + \
                              "vertical advection"
  v_advec_thetail@units     = "K s^-1"                            

  v_advec_qt           = new((/ dimsizes(levels),dimsizes(time) /),float)
  v_advec_qt           = 0.
  v_advec_qt@long_name = "prescribe q_t tendency due to vertical "+\
                         "advection"
  v_advec_qt@units     = "kg kg^-1 s^-1"

; Define dimension sizes for file creation
  ntim = dimsizes(time)
  klev = dimsizes(levels)

; Write the netCDF output file, titled "gabls3.nc"
  setfileoption("nc","FileStructure","Advanced")
  setfileoption("nc","Format","NetCDF4")
  fout = "gabls3.nc"           ;File name of output file
  system("/bin/rm -fv "+fout) ;If a version already exists, delete it.
  fo = addfile(fout,"c")      ;Create the netCDF file

  setfileoption(fo,"DefineMode",True) ;Enter define mode
  fAtt = True                         ;Set file attributes to True
  fAtt@title         = "GMTB SCM forcing file for GABLS3 case"
  fAtt@creation_date = systemfunc("date")
  fileattdef(fo,fAtt)                 ;Add attibutes to output file

  dimNames = (/ "time","levels" /)
  dimSizes = (/ ntim,klev /)
  dimUnlim = (/ False,True /)
;  dimNames = (/ "levels","time" /)
;  dimSizes = (/ klev,ntim /)
;  dimUnlim = (/ True,False /)

  filedimdef(fo,dimNames,dimSizes,dimUnlim)

;  chunkSizes = (/ dimsizes(time),dimsizes(levels)/2 /)
;  filechunkdimdef(fo,dimNames,chunkSizes,dimUnlim)

  filevardef(fo,"time",typeof(time),"time")
  filevarattdef(fo,"time",time)
  fo->time = (/ time /)

  filevardef(fo,"levels",typeof(levels),"levels")
  filevarattdef(fo,"levels",levels)
  fo->levels = (/ levels /)

; Define group names for netCDF4 file
  grpnames = (/ "scalars","initial","forcing" /)
  filegrpdef(fo,grpnames)

; Scalars  --> Lat and Lon

  g1 = fo=>/scalars
  g2 = fo=>/initial
  g3 = fo=>/forcing

;  filedimdef(g1,dimNames,dimSizes,dimUnlim)

; Initial

  filedimdef(g2,dimNames,dimSizes,dimUnlim)
  filedimdef(g3,dimNames,dimSizes,dimUnlim)
  
  filevardef(g1,"lat",typeof(lat),"ncl_scalar")
  filevarattdef(g1,"lat",lat)
  g1->lat = lat

  filevardef(g1,"lon",typeof(lon),"ncl_scalar")
  filevarattdef(g1,"lon",lon)
  g1->lon = lon
  
  filevardef(g2,"height",typeof(height),"levels")
  filevarattdef(g2,"height",height)
  g2->height = (/height/)

  filevardef(g2,"thetail",typeof(thetail),"levels")
  filevarattdef(g2,"thetail",thetail)
  g2->thetail = (/thetail/)

  filevardef(g2,"qt",typeof(qt),"levels")
  filevarattdef(g2,"qt",qt)
  g2->qt = (/qt/)

  filevardef(g2,"ql",typeof(ql),"levels")
  filevarattdef(g2,"ql",ql)
  g2->ql = (/ql/)

  filevardef(g2,"qi",typeof(qi),"levels")
  filevarattdef(g2,"qi",qi)
  g2->qi = (/qi/)

  filevardef(g2,"u",typeof(u),"levels")
  filevarattdef(g2,"u",u)
  g2->u = (/u/)

  filevardef(g2,"v",typeof(v),"levels")
  filevarattdef(g2,"v",v)
  g2->v = (/v/)

  filevardef(g2,"tke",typeof(tke),"levels")
  filevarattdef(g2,"tke",tke)
  g2->tke = (/tke/)

  filevardef(g2,"ozone",typeof(ozone),"levels")
  filevarattdef(g2,"ozone",ozone)
  g2->ozone = (/ozone/)

; Forcing

  filevardef(g3,"p_surf",typeof(p_surf),"time")
  filevarattdef(g3,"p_surf",p_surf)
  g3->p_surf = (/ p_surf /)

  filevardef(g3,"T_surf",typeof(T_surf),"time")
  filevarattdef(g3,"T_surf",T_surf)
  g3->T_surf = (/ T_surf /)

  filevardef(g3,"sh_flux_sfc",typeof(sh_flux_sfc),"time")
  filevarattdef(g3,"sh_flux_sfc",sh_flux_sfc)
  g3->sh_flux_sfc = (/ sh_flux_sfc /)

  filevardef(g3,"lh_flux_sfc",typeof(lh_flux_sfc),"time")
  filevarattdef(g3,"lh_flux_sfc",lh_flux_sfc)
  g3->lh_flux_sfc = (/ lh_flux_sfc /)

  filevardef(g3,"w_ls",typeof(w_ls),(/"levels","time"/))
  filevarattdef(g3,"w_ls",w_ls)
  g3->w_ls = (/ w_ls /)

  filevardef(g3,"omega",typeof(omega),(/"levels","time"/))
  filevarattdef(g3,"omega",omega)
  g3->omega = (/ omega /)

  filevardef(g3,"u_g",typeof(u_g),(/"levels","time"/))
  filevarattdef(g3,"u_g",u_g)
  g3->u_g = (/ u_g /)

  filevardef(g3,"v_g",typeof(v_g),(/"levels","time"/))
  filevarattdef(g3,"v_g",v_g)
  g3->v_g = (/ v_g /)

  filevardef(g3,"u_nudge",typeof(u_nudge),(/"levels","time"/))
  filevarattdef(g3,"u_nudge",u_nudge)
  g3->u_nudge = (/ u_nudge /)

  filevardef(g3,"v_nudge",typeof(v_nudge),(/"levels","time"/))
  filevarattdef(g3,"v_nudge",v_nudge)
  g3->v_nudge = (/ v_nudge /)

  filevardef(g3,"T_nudge",typeof(T_nudge),(/"levels","time"/))
  filevarattdef(g3,"T_nudge",T_nudge)
  g3->T_nudge = (/ T_nudge /)

  filevardef(g3,"thil_nudge",typeof(thil_nudge),(/"levels","time"/))
  filevarattdef(g3,"thil_nudge",thil_nudge)
  g3->thil_nudge = (/ thil_nudge /)

  filevardef(g3,"qt_nudge",typeof(qt_nudge),(/"levels","time"/))
  filevarattdef(g3,"qt_nudge",qt_nudge)
  g3->qt_nudge = (/ qt_nudge /)

  filevardef(g3,"dT_dt_rad",typeof(dT_dt_rad),(/"levels","time"/))
  filevarattdef(g3,"dT_dt_rad",dT_dt_rad)
  g3->dT_dt_rad = (/ dT_dt_rad /)

  filevardef(g3,"h_advec_thetail",typeof(h_advec_thetail), \
            (/"levels","time"/))
  filevarattdef(g3,"h_advec_thetail",h_advec_thetail)
  g3->h_advec_thetail = (/ h_advec_thetail /)

  filevardef(g3,"v_advec_thetail",typeof(v_advec_thetail), \
            (/"levels","time"/))
  filevarattdef(g3,"v_advec_thetail",v_advec_thetail)
  g3->v_advec_thetail = (/ v_advec_thetail /)

  filevardef(g3,"h_advec_qt",typeof(h_advec_qt),(/"levels","time"/))
  filevarattdef(g3,"h_advec_qt",h_advec_qt)
  g3->h_advec_qt = (/ h_advec_qt /)

  filevardef(g3,"v_advec_qt",typeof(v_advec_qt),(/"levels","time"/))
  filevarattdef(g3,"v_advec_qt",v_advec_qt)
  g3->v_advec_qt = (/ v_advec_qt /)

  print("wrote new file: "+fout)

end
