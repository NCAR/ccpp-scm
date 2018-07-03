;************************************************
; load required libraries and begin
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Constants
g     = 9.80665  ; gravity (m s^-1)
R_dry = 287.0    ; dry air constant (J deg^-1 kg^-1)
c_p   = 1004.0   ; specific heat at constant pressure (J deg^-1 kg^-1)
P0    = 100000.0 ; initial pressure (Pa)
kappa = R_dry/c_p

;;; BOMEX initial mean vertical profiles and forcings. At levels where no entries are
;;; found, values are to be linearly interpolated based on values at surrounding points.
;;; This information can be found in Siebesma et al. 2003 (JAS).

;; time as specified in Siebesma et al. 2003
time = fspan(0, 21600, 37) 
time@long_name = "elapsed time since the beginning of the simulation"
time@units = "s"

;; latitude and longitude of column
lat = 15.0
lat@long_name = "latitude of column"
lat@units = "degrees N"
lon = -56.0
lon@long_name = "longitude of column"
lon@units = "degrees E"

;; height levels as specified in Siebesma et al. 2003
height = (/0.,300.,500.,520.,700.,1480.,1500.,2000.,2100.,3000./)
height@long_name = "physical height at pressure levels"
height@units = "m"

;; total water specific humidity
qt_specified = (/17.0,16.3,10.7,4.2,3.0/)
qt_specified_kgkg = qt_specified/1000 ; convert g/kg to kg/kg
qt_h = (/0,520,1480,2000,3000/)
qt = linint1(qt_h,qt_specified_kgkg,False,height,0)
qt@long_name = "initial profile of total water specific humidity"
qt@units = "kg kg^-1"

;; liquid water potential temperature
thetail_specified = (/298.7,298.7,302.4,308.2,311.85/)
thetail_h = (/0,520,1480,2000,3000/)
thetail = linint1(thetail_h,thetail_specified,False,height,0)
thetail@long_name = "initial profile of ice-liquid water potential temperature"
thetail@units = "K"

;; horizontal velocity components 
u_specified = (/-8.75,-8.75,-4.61/)
u_h = (/0,700,3000/)
u = linint1(u_h,u_specified,False,height,0)
u@long_name = "initial profile of E-W horizontal wind"
u@units = "m s^-1"
v_specified = (/0,0/)
v_h = (/0,3000/)
v = linint1(v_h,v_specified,False,height,0)
v@long_name = "initial profile of N-S horizontal wind"
v@units = "m s^-1"

;; vertical velocity component
w_specified = (/0,-0.65,0/)
w_specified_ms = w_specified/100 ; convert cm/s to m/s
w_h = (/0,1500,2100/)
w = linint1(w_h,w_specified_ms,False,height,0)
w = where(ismissing(w),0,w)

w_ls = new((/dimsizes(w),dimsizes(time)/),"float")
w_ls@long_name = "large scale vertical velocity"
w_ls@units = "m s^-1"
do i=0,dimsizes(time)-1,1
  w_ls(:,i) = w
end do

;; forcing by radiative heating
rad_heat_specified = (/-2.0,-2.0,0/)
rad_heat_specified_ks = rad_heat_specified/86400 ; convert from K/day to K/sec
rad_heat_h = (/0,1500,3000/)
rad_heat = linint1(rad_heat_h,rad_heat_specified_ks,False,height,0)

dT_dt_rad = new((/dimsizes(rad_heat),dimsizes(time)/),"float")
dT_dt_rad@long_name = "prescribed radiative heating rate"
dT_dt_rad@units = "K/day" ; check units
do i=0,dimsizes(time)-1,1
  dT_dt_rad(:,i) = rad_heat
end do

;; horizontal advection of qt
hor_advec_qt_specified = (/-1.2e-8,-1.2e-8,0/)
hor_advec_qt_specified_kgkg = hor_advec_qt_specified/1000 ; convert g/kg to kg/kg
hor_advec_qt_h = (/0,300,500/)
hor_advec_qt = linint1(hor_advec_qt_h,hor_advec_qt_specified_kgkg,False,height,0)
hor_advec_qt = where(ismissing(hor_advec_qt),0,hor_advec_qt)

h_advec_qt = new((/dimsizes(hor_advec_qt),dimsizes(time)/),"float")
h_advec_qt@long_name = "prescribed q_t tendency due to horizontal advection"
h_advec_qt@units = "kg kg^-1 s^-1" ; check units
do i=0,dimsizes(time)-1,1
  h_advec_qt(:,i) = hor_advec_qt
end do

;; sensible and latent heat fluxes
sh_flux_sfc = 8e-3
sh_flux_sfc@long_name = "surface sensible heat flux"
sh_flux_sfc@units     = "K m s^-1"

lh_flux_sfc = 5.2e-5
lh_flux_sfc@long_name = "surface latent heat flux"
lh_flux_sfc@units     = "m s^-1"  ; check units

;; need to convert level (m) to pressures (Pa) using Exner eqn
; calculate Exner and pressure at surface
p_surf = 101500.0 ; from Siebesma et al. 2003
p_surf@long_name = "surface pressure"
p_surf@units     = "Pa"

exner = new((/dimsizes(height)/),float)
levels = new((/dimsizes(height)/),float)
levels@long_name = "Pa"
levels@units = "pressure levels"

; calulcate the Exner value at the surface
exner_surf = (p_surf/P0)^kappa

; calculate the Exner value and pressure at first level (which is also the surface in this case)
exner(0) = exner_surf - (height(0) - height(0))*(g/c_p)/((1.0 + 0.61*qt(0))*thetail(0))
levels(0) = P0*exner(0)^(1.0/kappa)

; calculate exner at every z-level above
do i = 1, dimsizes(height)-1
  exner(i) = exner(i-1) - (height(i)-height(i-1))*(g/c_p)/((1.0 + 0.61*0.5*(qt(i-1)+qt(i)))*0.5*(thetail(i-1)+thetail(i)))
  levels(i) = P0*exner(i)^(1.0/kappa)
end do

;; covert w to omega
; solve for T first
T = thetail*(levels/P0)^kappa
; solve for rho
rho = levels/(R_dry*T)
; solve for omega
om = -w*rho*g

omega = new((/dimsizes(om),dimsizes(time)/),"float")
omega@long_name = "large scale pressure vertical velocity"
omega@units = "Pa s^-1"
do i=0,dimsizes(time)-1,1
  omega(:,i) = om
end do

; using T from omega calculation set surface temperature
T_surf                 = T(0)
T_surf@long_name       = "surface absolute temperature"
T_surf@units           = "K"

;;; The following variables are not in Siebesma et al. 2003, but are necessary
;; Initial profiles
ql                = new((/dimsizes(levels)/),float)
ql                = 0.0
ql@long_name      = "initial profile of liquid water specific humidity"
ql@units          = "kg kg^-1"

qi                = new((/dimsizes(levels)/),float)
qi                = 0.0
qi@long_name      = "initial profile of ice water specific humidity"
qi@units          = "kg kg^-1"

tke               = new((/dimsizes(levels)/),float)
tke               = 0.0
tke@long_name     = "initial profile of turbulence kinetic energy"
tke@units         = "m^2 s^-2"

ozone             = new((/dimsizes(levels)/),float)
ozone             = 0.0
ozone@long_name   = "initial profile of ozone mass mixing ratio"
ozone@units       = "kg kg^-1"

;; Forcings
h_advec_thetail        = new((/dimsizes(levels),dimsizes(time)/),float)
h_advec_thetail(:,:)   = 0.0
h_advec_thetail@long_name = "prescribed theta_il tendency due to horizontal advection"
h_advec_thetail@units     = "K s^-1"

v_advec_thetail           = new((/dimsizes(levels),dimsizes(time)/),float)
v_advec_thetail(:,:)      = 0.0
v_advec_thetail@long_name = "prescribed theta_il tendency due to vertical advection"
v_advec_thetail@units     = "K s^-1"

v_advec_qt             = new((/dimsizes(levels),dimsizes(time)/),float)
v_advec_qt(:,:)        = 0.0
v_advec_qt@long_name   = "prescribed q_t tendency due to vertical advection"
v_advec_qt@units       = "kg kg^-1 s^-1"

u_g                    = new((/dimsizes(levels),dimsizes(time)/),float)
u_g(:,:)               = 0.0
u_g@long_name          = "large scale geostrophic E-W wind"
u_g@units              = "m s^-1"

v_g                    = new((/dimsizes(levels),dimsizes(time)/),float)
v_g(:,:)               = 0.0
v_g@long_name          = "large scale geostrophic N-S wind"
v_g@units              = "m s^-1"

u_nudge                = new((/dimsizes(levels),dimsizes(time)/),float)
u_nudge(:,:)           = 0.0
u_nudge@long_name      = "E-W wind to nudge toward"
u_nudge@units          = "m s^-1"

v_nudge                = new((/dimsizes(levels),dimsizes(time)/),float)
v_nudge(:,:)           = 0.0
v_nudge@long_name      = "N-S wind to nudge toward"
v_nudge@units          = "m s^-1"

T_nudge                = new((/dimsizes(levels),dimsizes(time)/),float)
T_nudge(:,:)           = 0.0
T_nudge@long_name      = "absolute temperature to nudge toward"
T_nudge@units          = "K"

thil_nudge             = new((/dimsizes(levels),dimsizes(time)/),float)
thil_nudge(:,:)        = 0.0
thil_nudge@long_name   = "potential temperature to nudge toward"
thil_nudge@units       = "K"

qt_nudge               = new((/dimsizes(levels),dimsizes(time)/),float)
qt_nudge(:,:)          = 0.0
qt_nudge@long_name     = "q_t to nudge toward"
qt_nudge@units         = "kg kg^-1"

;;; Write output to BOMEX forcing file (named: bomex_input.nc)
setfileoption("nc","Format","NetCDF4")                 ; Se to output NetCDF4
fout = "bomex.nc"                                ; Name of outut file
system("/bin/rm -f "+fout)                             ; Remove any existing file with the same output name
fo = addfile(fout,"c")                                 ; Create the netcdf file

setfileoption(fo,"DefineMode",True)                    ; Enter define mode
fAtt = True                                            ; Set file attributes to True
fAtt@title = "GMTB SCM forcing file for BOMEX case"    ; Set the title of the file
fileattdef(fo,fAtt)                                    ; Add the attributes to the output file

; Output file dimensions
;dimNames = (/"levels","time"/)
;dimSizes = (/dimsizes(levels),dimsizes(time)/)
;dimUnlim = (/True,False/)                        ; Do the variables have unlimited dimension?
dimNames = (/"time","levels"/)
dimSizes = (/dimsizes(time),dimsizes(levels)/)
dimUnlim = (/False,True/)                        ; Do the variables have unlimited dimension?
filedimdef(fo,dimNames,dimSizes,dimUnlim)

filevardef(fo,"time",typeof(time),"time")
filevarattdef(fo,"time",time)
fo->time = (/time/)

filevardef(fo,"levels",typeof(levels),"levels")
filevarattdef(fo,"levels",levels)
fo->levels = (/levels/)

; Group names for NetCDF4 file
grpnames = (/"scalars", "initial", "forcing"/)
filegrpdef(fo, grpnames)

; Scalars -- none

; Initial
g2 = fo=>/initial

filedimdef(g2,dimNames,dimSizes,dimUnlim)

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
g3 = fo=>/forcing

filedimdef(g3,dimNames,dimSizes,dimUnlim)

filevardef(g3,"lat",typeof(lat),"time")
filevarattdef(g3,"lat",lat)
g3->lat = (/lat/)

filevardef(g3,"lon",typeof(lon),"time")
filevarattdef(g3,"lon",lon)
g3->lon = (/lon/)

filevardef(g3,"p_surf",typeof(p_surf),"time")
filevarattdef(g3,"p_surf",p_surf)
g3->p_surf = (/p_surf/)

filevardef(g3,"T_surf",typeof(T_surf),"time")
filevarattdef(g3,"T_surf",T_surf)
g3->T_surf = (/T_surf/)

filevardef(g3,"sh_flux_sfc",typeof(sh_flux_sfc),"time")
filevarattdef(g3,"sh_flux_sfc",sh_flux_sfc)
g3->sh_flux_sfc = (/sh_flux_sfc/)

filevardef(g3,"lh_flux_sfc",typeof(lh_flux_sfc),"time")
filevarattdef(g3,"lh_flux_sfc",lh_flux_sfc)
g3->lh_flux_sfc = (/lh_flux_sfc/)

filevardef(g3,"w_ls",typeof(w_ls),(/"levels","time"/))
filevarattdef(g3,"w_ls",w_ls)
g3->w_ls = (/w_ls/)

filevardef(g3,"omega",typeof(omega),(/"levels","time"/))
filevarattdef(g3,"omega",omega)
g3->omega = (/omega/)

filevardef(g3,"u_g",typeof(u_g),(/"levels","time"/))
filevarattdef(g3,"u_g",u_g)
g3->u_g = (/u_g/)

filevardef(g3,"v_g",typeof(v_g),(/"levels","time"/))
filevarattdef(g3,"v_g",v_g)
g3->v_g = (/v_g/)

filevardef(g3,"u_nudge",typeof(u_nudge),(/"levels","time"/))
filevarattdef(g3,"u_nudge",u_nudge)
g3->u_nudge = (/u_nudge/)

filevardef(g3,"v_nudge",typeof(v_nudge),(/"levels","time"/))
filevarattdef(g3,"v_nudge",v_nudge)
g3->v_nudge = (/v_nudge/)

filevardef(g3,"T_nudge",typeof(T_nudge),(/"levels","time"/))
filevarattdef(g3,"T_nudge",T_nudge)
g3->T_nudge = (/T_nudge/)

filevardef(g3,"thil_nudge",typeof(thil_nudge),(/"levels","time"/))
filevarattdef(g3,"thil_nudge",thil_nudge)
g3->thil_nudge = (/thil_nudge/)

filevardef(g3,"qt_nudge",typeof(qt_nudge),(/"levels","time"/))
filevarattdef(g3,"qt_nudge",qt_nudge)
g3->qt_nudge = (/qt_nudge/)

filevardef(g3,"dT_dt_rad",typeof(dT_dt_rad),(/"levels","time"/))
filevarattdef(g3,"dT_dt_rad",dT_dt_rad)
g3->dT_dt_rad = (/dT_dt_rad/)

filevardef(g3,"h_advec_thetail",typeof(h_advec_thetail),(/"levels","time"/))
filevarattdef(g3,"h_advec_thetail",h_advec_thetail)
g3->h_advec_thetail = (/h_advec_thetail/)

filevardef(g3,"v_advec_thetail",typeof(v_advec_thetail),(/"levels","time"/))
filevarattdef(g3,"v_advec_thetail",v_advec_thetail)
g3->v_advec_thetail = (/v_advec_thetail/)

filevardef(g3,"h_advec_qt",typeof(h_advec_qt),(/"levels","time"/))
filevarattdef(g3,"h_advec_qt",h_advec_qt)
g3->h_advec_qt = (/h_advec_qt/)

filevardef(g3,"v_advec_qt",typeof(v_advec_qt),(/"levels","time"/))
filevarattdef(g3,"v_advec_qt",v_advec_qt)
g3->v_advec_qt = (/v_advec_qt/)


