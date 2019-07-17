#!/usr/bin/env python
from netCDF4 import Dataset
import numpy as np

# define date
YYYY=2014
MM=8
DD=1
HH=0
MI=0
SC=0

# define path to FV3 GFS initial and boundary conditions
icpath='../../data/raw_case_input/UFS_test_ics'
fixpath='../../data/raw_case_input/UFS_test_ics'
print (icpath)
#  defint i,j and tile to extract colmn
ipt=16
jpt=41
tilenum=2

ipt2=ipt*2+1
jpt2=jpt*2+1

#  I have surface cycle on the sfc_data initial conditions to get the proper surface fields
infile1='%s/gfs_data.tile%i.nc' % (icpath,tilenum)
infile2='%s/sfc_data.tile%i.nc' %(icpath,tilenum)
infile3='%s/C96_oro_data.tile%i.nc' %(fixpath,tilenum)
infile4='%s/C96_grid.tile%i.nc' %(fixpath,tilenum)
infile5='%s/gfs_ctrl.nc' %icpath
ncin1=Dataset(infile1)
ncin2=Dataset(infile2)
ncin3=Dataset(infile3)
ncin4=Dataset(infile4)
ncin5=Dataset(infile5)
# assume model contains one less level than the cold start spectral GFS initial conditions
nlevs=len(ncin1.dimensions['lev'])-1
# pick off lat and lon from i and j point defined above
lon0=ncin4['x'][jpt2,ipt2]
lat0=ncin4['y'][jpt2,ipt2]

# extract out area of grid cell
area_in=ncin4['area'][jpt2-1:jpt2+1,ipt2-1:ipt2+1]
print (lat0,lon0)

# upper air fields from initial conditions
zh=ncin1['zh'][::-1,jpt,ipt]
uw1=ncin1['u_w'][::-1,jpt,ipt]
uw2=ncin1['u_w'][::-1,jpt,ipt+1]
us1=ncin1['u_s'][::-1,jpt,ipt]
us2=ncin1['u_s'][::-1,jpt+1,ipt]
vw1=ncin1['v_w'][::-1,jpt,ipt]
vw2=ncin1['v_w'][::-1,jpt,ipt+1]
vs1=ncin1['v_s'][::-1,jpt,ipt]
vs2=ncin1['v_s'][::-1,jpt+1,ipt]
ucomp=0.25*(uw1+uw2+us1+us2)  # estimate u winds on the a grid
vcomp=0.25*(vw1+vw2+vs1+vs2)  # estimate v winds on the a grid
sphum=ncin1['sphum'][::-1,jpt,ipt]
# o3 and qv are taken from ics. 
o3=ncin1['o3mr'][::-1,jpt,ipt]
liqwat=ncin1['liq_wat'][:-1,jpt,ipt]

# surface pressure and skin temperature
ps=ncin1['ps'][jpt,ipt]
ts=ncin2['tsea'][jpt,ipt]

# land state
stc_in=ncin2['stc'][:,jpt,ipt]
smc_in=ncin2['smc'][:,jpt,ipt]
slc_in=ncin2['slc'][:,jpt,ipt]
tg3_in=ncin2['tg3'][jpt,ipt]

# surface properties
uustar_in=ncin2['uustar'][jpt,ipt]
alvsf=ncin2['alvsf'][jpt,ipt]
alvwf=ncin2['alvwf'][jpt,ipt]
alnsf=ncin2['alnsf'][jpt,ipt]
alnwf=ncin2['alnwf'][jpt,ipt]
facsf_in=ncin2['facsf'][jpt,ipt]
facwf_in=ncin2['facwf'][jpt,ipt]
styp_in=ncin2['stype'][jpt,ipt]
slope_in=ncin2['slope'][jpt,ipt]
vtyp_in=ncin2['vtype'][jpt,ipt]
vfrac_in=ncin2['vfrac'][jpt,ipt]
shdmin_in=ncin2['shdmin'][jpt,ipt]
shdmax_in=ncin2['shdmax'][jpt,ipt]
zorl_in=ncin2['zorl'][jpt,ipt]
slmsk_in=ncin2['slmsk'][jpt,ipt]
canopy_in=ncin2['canopy'][jpt,ipt]
hice_in=ncin2['hice'][jpt,ipt]
fice_in=ncin2['fice'][jpt,ipt]
tisfc_in=ncin2['tisfc'][jpt,ipt]
snwdph_in=ncin2['snwdph'][jpt,ipt]
snoalb_in=ncin2['snoalb'][jpt,ipt]

# orographyic properties
stddev_in=ncin3['stddev'][jpt,ipt]
convexity_in=ncin3['convexity'][jpt,ipt]
oa1_in=ncin3['oa1'][jpt,ipt]
oa2_in=ncin3['oa2'][jpt,ipt]
oa3_in=ncin3['oa3'][jpt,ipt]
oa4_in=ncin3['oa4'][jpt,ipt]
ol1_in=ncin3['ol1'][jpt,ipt]
ol2_in=ncin3['ol2'][jpt,ipt]
ol3_in=ncin3['ol3'][jpt,ipt]
ol4_in=ncin3['ol4'][jpt,ipt]
theta_in=ncin3['theta'][jpt,ipt]
gamma_in=ncin3['gamma'][jpt,ipt]
sigma_in=ncin3['sigma'][jpt,ipt]
elvmax_in=ncin3['elvmax'][jpt,ipt]

# vertical coordinate definition
ak=ncin5['vcoord'][0,::-1]
bk=ncin5['vcoord'][1,::-1]

#calculate temperature
rdgas  = 287.05
rvgas  = 461.50
zvir = rvgas/rdgas - 1.
grav=9.80665
gz=zh*grav
pn1=np.zeros([nlevs+1])
temp=np.zeros([nlevs])
for k in range(nlevs+1):
  pn1[k]=np.log(ak[k]+ps*bk[k])
for k in range(nlevs):
  temp[k] = (gz[k]-gz[k+1])/( rdgas*(pn1[k+1]-pn1[k])*(1.+zvir*sphum[k]) )

# open output file
nc = Dataset('../../data/processed_case_input/fv3_model_point.nc', mode='w')
nc.description = "FV3GFS model profile input (no forcing)"

time = nc.createDimension('time',None)
levels = nc.createDimension('levels',None)
nsoil  = nc.createDimension('nsoil',None)
t = nc.createVariable('time',np.float64,('time',))
t.units = "s" 
t.description = "elapsed time since the beginning of the simulation" 
z = nc.createVariable('levels',np.float64,('levels',))
z.units = "Pa"
z.description = "pressure levels"
#scalars
iyr = nc.createVariable('scalars/init_year',np.int32)
imo = nc.createVariable('scalars/init_month',np.int32)
idy = nc.createVariable('scalars/init_day',np.int32)
ihr = nc.createVariable('scalars/init_hour',np.int32)
imi = nc.createVariable('scalars/init_minute',np.int32)
isc = nc.createVariable('scalars/init_second',np.int32)
ivegsrc  = nc.createVariable('scalars/vegsrc',np.int32)
ivegtyp  = nc.createVariable('scalars/vegtyp',np.int32)
isoiltyp = nc.createVariable('scalars/soiltyp',np.int32)
islopetyp = nc.createVariable('scalars/slopetyp',np.int32)
ivegfrac = nc.createVariable('scalars/vegfrac',np.float)
ishdmin  = nc.createVariable('scalars/shdmin',np.float)
ishdmax  = nc.createVariable('scalars/shdmax',np.float)
izorl    = nc.createVariable('scalars/zorl',np.float)
islmsk   = nc.createVariable('scalars/slmsk',np.float)
icanopy    = nc.createVariable('scalars/canopy',np.float)
ihice    = nc.createVariable('scalars/hice',np.float)
ifice    = nc.createVariable('scalars/fice',np.float)
itisfc    = nc.createVariable('scalars/tisfc',np.float)
isnwdph    = nc.createVariable('scalars/snwdph',np.float)
isnoalb    = nc.createVariable('scalars/snoalb',np.float)
isncovr    = nc.createVariable('scalars/sncovr',np.float)
itg3     = nc.createVariable('scalars/tg3',np.float)
iuustar  = nc.createVariable('scalars/uustar',np.float)

iyr.units = "years"  
iyr.description = "year at time of initial values"  
imo.units = "months"  
imo.description = "month at time of initial values"  
idy.units = "days"  
idy.description = "day at time of initial values"  
ihr.units = "hours"  
ihr.description = "hour at time of initial values"  
imi.units = "minutes"  
imi.description = "minute at time of initial values"  
isc.units = "seconds"  
isc.description = "second at time of initial values"  
ivegsrc.description = "vegetation soure (1-2)"
ivegtyp.description = "vegetation type (1-12)"
isoiltyp.description = "soil type (1-12)"
islopetyp.description = "slope type (1-9)"
ivegfrac.description = "vegetation fraction"
ishdmin.description = "minimum vegetation fraction"
ishdmax.description = "maximum vegetation fraction"
izorl.description = "surface roughness length"
islmsk.description = "land-sea-ice mask"
icanopy.description = "canopy moisture"
ihice.description = "ice thickness"
ifice.description = "ice fraction"
itisfc.description = "ice temperature"
isnwdph.description = "snow depth"
isnoalb.description = "snow albedo"
isncovr.description = "snow cover"
itg3.description = "deep soil temperature"
itg3.units = "K"  
iuustar.description = "frication velocity"
iuustar.units = "m2s-2?"  


#initial
ic_t = nc.createVariable('initial/temp',np.float64,('levels',))
ic_qt = nc.createVariable('initial/qt',np.float64,('levels',))
ic_ql = nc.createVariable('initial/ql',np.float64,('levels',))
ic_qi = nc.createVariable('initial/qi',np.float64,('levels',))
ic_u = nc.createVariable('initial/u',np.float64,('levels',))
ic_v = nc.createVariable('initial/v',np.float64,('levels',))
ic_tke = nc.createVariable('initial/tke',np.float64,('levels',))
ic_o3  = nc.createVariable('initial/ozone',np.float64,('levels',))
ic_stc = nc.createVariable('initial/stc',np.float64,('nsoil',))
ic_smc = nc.createVariable('initial/smc',np.float64,('nsoil',))
ic_slc = nc.createVariable('initial/slc',np.float64,('nsoil',))
ic_t.units = "K" 
#ic_t.description = "initial profile of ice-liquid water potential temperature" 
ic_t.description = "initial profile of temperature" 
ic_qt.units = "kg kg^-1" 
ic_qt.description = "initial profile of total water specific humidity" 
ic_ql.units = "kg kg^-1" 
ic_ql.description = "initial profile of liquid water specific humidity" 
ic_qi.units = "kg kg^-1" 
ic_qi.description = "initial profile of ice water specific humidity" 
ic_u.units = "m s^-1" 
ic_u.description = "initial profile of E-W horizontal wind" 
ic_v.units = "m s^-1" 
ic_v.description = "initial profile of N-S horizontal wind" 
ic_tke.units = "m^2 s^-2" 
ic_tke.description = "initial profile of turbulence kinetic energy" 
ic_o3.units = "kg kg^-1" 
ic_o3.description = "initial profile of ozone mass mixing ratio" 
ic_stc.units = "K"
ic_stc.description = "initial profile of soil temperature"
ic_smc.units = "kg"
ic_smc.description = "initial profile of soil moisture"
ic_slc.units = "kg"
ic_slc.description = "initial profile of soil liquid moisture"

lat=nc.createVariable('forcing/lat',np.float64,('time',))
lat.units = "degrees N" 
lat.description = "latitude of column" 
lon=nc.createVariable('forcing/lon',np.float64,('time',))
lon.units = "degrees E" 
lon.description = "longitude of column" 
p_surf=nc.createVariable('forcing/p_surf',np.float64,('time',))
p_surf.units = "Pa" 
p_surf.description = "surface pressure" 
T_surf=nc.createVariable('forcing/T_surf',np.float64,('time',))
T_surf.units = "K" 
T_surf.description = "surface absolute temperature" 
area1=nc.createVariable('scalars/area',np.float64,('time',))
alb1=nc.createVariable('scalars/alvsf',np.float64,('time',))
alb2=nc.createVariable('scalars/alnsf',np.float64,('time',))
alb3=nc.createVariable('scalars/alvwf',np.float64,('time',))
alb4=nc.createVariable('scalars/alnwf',np.float64,('time',))
stddev=nc.createVariable('scalars/stddev',np.float64,('time',))
convexity=nc.createVariable('scalars/convexity',np.float64,('time',))
oa1=nc.createVariable('scalars/oa1',np.float64,('time',))
oa2=nc.createVariable('scalars/oa2',np.float64,('time',))
oa3=nc.createVariable('scalars/oa3',np.float64,('time',))
oa4=nc.createVariable('scalars/oa4',np.float64,('time',))
ol1=nc.createVariable('scalars/ol1',np.float64,('time',))
ol2=nc.createVariable('scalars/ol2',np.float64,('time',))
ol3=nc.createVariable('scalars/ol3',np.float64,('time',))
ol4=nc.createVariable('scalars/ol4',np.float64,('time',))
theta=nc.createVariable('scalars/theta',np.float64,('time',))
gamma=nc.createVariable('scalars/gamma',np.float64,('time',))
sigma=nc.createVariable('scalars/sigma',np.float64,('time',))
elvmax=nc.createVariable('scalars/elvmax',np.float64,('time',))
facsf=nc.createVariable('scalars/facsf',np.float64,('time',))
facwf=nc.createVariable('scalars/facwf',np.float64,('time',))
area1.units = "m^2" 
alb1.units = "None" 
alb2.units = "None" 
alb3.units = "None" 
alb4.units = "None" 
facsf.units = "None" 
facwf.units = "None" 
area1.description = "grid cell area"
alb1.description = "uv+visible black sky albedo (z=60 degree)"
alb2.description = "near IR black sky albedo (z=60 degree)"
alb3.description = "uv+visible white sky albedo"
alb4.description = "near IR white sky albedo"
stddev.description = "surface orography standard deviation"
facsf.description = "fraction of grid cell with strong sun angle albedo dependence"
facwf.description = "fraction of grid cell with weak sun angle albedo dependence"
w_ls=nc.createVariable('forcing/w_ls',np.float64,('levels','time',))
w_ls.units = "m s^-1" 
w_ls.description = "large scale vertical velocity" 
omega=nc.createVariable('forcing/omega',np.float64,('levels','time',))
omega.units = "Pa s^-1" 
omega.description = "large scale pressure vertical velocity" 
u_g=nc.createVariable('forcing/u_g',np.float64,('levels','time',))
u_g.units = "m s^-1" 
u_g.description = "large scale geostrophic E-W wind" 
v_g=nc.createVariable('forcing/v_g',np.float64,('levels','time',))
v_g.units = "m s^-1" 
v_g.description = "large scale geostrophic N-S wind" 
u_nudge=nc.createVariable('forcing/u_nudge',np.float64,('levels','time',))
u_nudge.units = "m s^-1" 
u_nudge.description = "E-W wind to nudge toward" 
v_nudge=nc.createVariable('forcing/v_nudge',np.float64,('levels','time',))
v_nudge.units = "m s^-1" 
v_nudge.description = "N-S wind to nudge toward" 
T_nudge=nc.createVariable('forcing/T_nudge',np.float64,('levels','time',))
T_nudge.units = "K" 
T_nudge.description = "absolute temperature to nudge toward" 
thil_nudge=nc.createVariable('forcing/thil_nudge',np.float64,('levels','time',))
thil_nudge.units = "K" 
thil_nudge.description = "potential temperature to nudge toward" 
qt_nudge=nc.createVariable('forcing/qt_nudge',np.float64,('levels','time',))
qt_nudge.units = "kg kg^-1" 
qt_nudge.description = "q_t to nudge toward" 
dT_dt_rad=nc.createVariable('forcing/dT_dt_rad',np.float64,('levels','time',))
dT_dt_rad.units = "K s^-1" 
dT_dt_rad.description = "prescribed radiative heating rate" 
h_advec_thetail=nc.createVariable('forcing/h_advec_thetail',np.float64,('levels','time',))
h_advec_thetail.units = "K s^-1" 
h_advec_thetail.description = "prescribed theta_il tendency due to horizontal advection" 
v_advec_thetail=nc.createVariable('forcing/v_advec_thetail',np.float64,('levels','time',))
v_advec_thetail.units = "K s^-1" 
v_advec_thetail.description = "prescribed theta_il tendency due to vertical advection" 
h_advec_qt=nc.createVariable('forcing/h_advec_qt',np.float64,('levels','time',))
h_advec_qt.units = "kg kg^-1 s^-1" 
h_advec_qt.description = "prescribed q_t tendency due to horizontal advection" 
v_advec_qt=nc.createVariable('forcing/v_advec_qt',np.float64,('levels','time',))
v_advec_qt.units = "kg kg^-1 s^-1" 
v_advec_qt.description = "prescribed q_t tendency due to vertical advection" 


# date
iyr[:]=YYYY
imo[:]=MM
idy[:]=DD
ihr[:]=HH
imi[:]=MI
isc[:]=SC
# axes
t[0]=0

#ics
ic_t[:]   = temp[0:nlevs]
ic_qt[:]  = sphum[0:nlevs]
ic_ql[:]  = liqwat[0:nlevs]
ic_qi[:]  = 0.0
ic_u[:]   = ucomp[0:nlevs]
ic_v[:]   = vcomp[0:nlevs]
ic_tke[:] = 0.0
ic_o3[:]  = o3[0:nlevs]
ic_stc[:] = stc_in
ic_smc[:] = smc_in
ic_slc[:] = slc_in

lat[:]=lat0
lon[:]=lon0

p_surf[:]=ps
T_surf[:]=ts
area1[:]=area_in.sum()
alb1[:]=alvsf
alb2[:]=alnsf
alb3[:]=alvwf
alb4[:]=alnwf
stddev[:]=stddev_in
convexity[:]=convexity_in
oa1[:]=oa1_in
oa2[:]=oa2_in
oa3[:]=oa3_in
oa4[:]=oa4_in
ol1[:]=ol1_in
ol2[:]=ol2_in
ol3[:]=ol3_in
ol4[:]=ol4_in
theta[:]=theta_in
gamma[:]=gamma_in
sigma[:]=sigma_in
elvmax[:]=elvmax_in
facsf[:]=facsf_in
facwf[:]=facwf_in
z[:]=np.exp(pn1[0:nlevs])
w_ls[:]=0.0
omega[:]=0.0
u_g[:]=0.0
v_g[:]=0.0
u_nudge[:]=0.0
v_nudge[:]=0.0
T_nudge[:]=0.0
thil_nudge[:]=0.0
qt_nudge[:]=0.0
dT_dt_rad[:]=0.0
h_advec_thetail[:]=0.0
v_advec_thetail[:]=0.0
h_advec_qt[:]=0.0
v_advec_qt[:]=0.0

ivegsrc[:] = 1
ivegtyp[:] = vtyp_in
isoiltyp[:] = styp_in
islopetyp[:] = slope_in
ishdmin[:] = shdmin_in
ishdmax[:] = shdmax_in
izorl[:] = zorl_in
islmsk[:] = slmsk_in
icanopy[:] = canopy_in
ihice[:] = hice_in
ifice[:] = fice_in
itisfc[:] = tisfc_in
isnwdph[:] = snwdph_in
isnoalb[:] = snoalb_in
isncovr[:] = 0.0
itg3[:] = tg3_in
iuustar[:] = uustar_in
ivegfrac[:]=vfrac_in

nc.close()

