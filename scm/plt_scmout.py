##########################################################################################
##########################################################################################
import os,netCDF4,numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from pylab import figure

to = 0
tf = 6

exp    = 'twpice'
dirG   = 'output_'+exp+'_SCM_GFS_v15'
dirGP  = 'output_'+exp+'_SCM_GFS_v15_RRTMGP_input_GFS_v15_RRTMGP'
fileG  = 'bin/'+dirG+'/output.nc'
fileGP = 'bin/'+dirGP+'/output.nc'

fileOUT  = 'scm_GvGP_'+exp+'.eps'
fileOUTd = 'scm_GmGP_'+exp+'.eps'

# Read in data
# RRTMG
dataG                   = netCDF4.Dataset(fileG,'r')
time                    = dataG.variables['time'][:]
pres                    = dataG.variables['pres'][:,:,:]
pres_i                  = dataG.variables['pres_i'][:,:,:]
sigma                   = dataG.variables['sigma'][:,:,:]
sigma_i                 = dataG.variables['sigma_i'][:,:,:]
cldcovR                 = dataG.variables['cldcov'][:,:,:]
lw_rad_heating_rate_REF = dataG.variables['lw_rad_heating_rate'][:,:,:]
sw_rad_heating_rate_REF = dataG.variables['sw_rad_heating_rate'][:,:,:]
dT_dt_lwrad_REF         = dataG.variables['dT_dt_lwrad'][:,:,:]
dT_dt_swrad_REF         = dataG.variables['dT_dt_swrad'][:,:,:]
sw_up_TOA_totR          = dataG.variables['sw_up_TOA_tot'][:,:]
sw_dn_TOA_totR          = dataG.variables['sw_dn_TOA_tot'][:,:]
sw_up_TOA_clrR          = dataG.variables['sw_up_TOA_clr'][:,:]
sw_up_SFC_totR          = dataG.variables['sw_up_sfc_tot'][:,:]
sw_dn_SFC_totR          = dataG.variables['sw_dn_sfc_tot'][:,:]
sw_up_SFC_clrR          = dataG.variables['sw_up_sfc_clr'][:,:]
sw_dn_SFC_clrR          = dataG.variables['sw_dn_sfc_clr'][:,:]
lw_up_TOA_totR          = dataG.variables['lw_up_TOA_tot'][:,:]
lw_up_TOA_clrR          = dataG.variables['lw_up_TOA_clr'][:,:]
lw_up_SFC_totR          = dataG.variables['lw_up_sfc_tot'][:,:]
lw_dn_SFC_totR          = dataG.variables['lw_dn_sfc_tot'][:,:]
lw_up_SFC_clrR          = dataG.variables['lw_up_sfc_clr'][:,:]
lw_dn_SFC_clrR          = dataG.variables['lw_dn_sfc_clr'][:,:]

# RRTMGP
dataGP              = netCDF4.Dataset(fileGP,'r')
lw_rad_heating_rate = dataGP.variables['lw_rad_heating_rate'][:,:,:]
sw_rad_heating_rate = dataGP.variables['sw_rad_heating_rate'][:,:,:]
dT_dt_lwrad         = dataGP.variables['dT_dt_lwrad'][:,:,:]
dT_dt_swrad         = dataGP.variables['dT_dt_swrad'][:,:,:]
cldcov              = dataGP.variables['cldcov'][:,:,:]
sw_up_TOA_tot       = dataGP.variables['sw_up_TOA_tot'][:,:]
sw_dn_TOA_tot       = dataGP.variables['sw_dn_TOA_tot'][:,:]
sw_up_TOA_clr       = dataGP.variables['sw_up_TOA_clr'][:,:]
sw_up_SFC_tot       = dataGP.variables['sw_up_sfc_tot'][:,:]
sw_dn_SFC_tot       = dataGP.variables['sw_dn_sfc_tot'][:,:]
sw_up_SFC_clr       = dataGP.variables['sw_up_sfc_clr'][:,:]
sw_dn_SFC_clr       = dataGP.variables['sw_dn_sfc_clr'][:,:]
lw_up_TOA_clr       = dataGP.variables['lw_up_TOA_clr'][:,:]
lw_up_TOA_tot       = dataGP.variables['lw_up_TOA_tot'][:,:]
lw_up_SFC_tot       = dataGP.variables['lw_up_sfc_tot'][:,:]
lw_dn_SFC_tot       = dataGP.variables['lw_dn_sfc_tot'][:,:]
lw_up_SFC_clr       = dataGP.variables['lw_up_sfc_clr'][:,:]
lw_dn_SFC_clr       = dataGP.variables['lw_dn_sfc_clr'][:,:]
ntime = len(time)
nlay = len(pres[0,:])

# Compute mean pressure (for contouring)
mean_pres = np.sum(pres[:,:],axis=0)/len(time)

nrows = 6
ncols = 4

# Make plots...


##########################################################################
# 1) Overlapping fields
##########################################################################
fig  = plt.figure(0, figsize=(12,14), dpi=80, facecolor='w', edgecolor='k')

# Plot heating-rate profiles
# 1) SW
plt.subplot(nrows,ncols,1)
plt.plot(time/3600.,sw_dn_TOA_totR[:,0],'-b')
plt.plot(time/3600.,sw_dn_TOA_tot[:,0],'-r')
plt.axis([to,tf,0,1400])
plt.xticks(rotation=45)
plt.ylabel('Flux (W/m2)')
plt.xlabel('Forecast time (hours)')
plt.title('SW down @ TOA (all)')
#
plt.subplot(nrows,ncols,2)
plt.plot(time/3600.,sw_dn_SFC_totR[:,0],'-b')
plt.plot(time/3600.,sw_dn_SFC_tot[:,0],'-r')
plt.axis([to,tf,0,1200])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW down @ surface (all)')
#
plt.subplot(nrows,ncols,3)
plt.plot(time/3600.,sw_up_SFC_totR[:,0],'-b')
plt.plot(time/3600.,sw_up_SFC_tot[:,0],'-r')
plt.axis([to,tf,0,100])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ surface (all)')
#
plt.subplot(nrows,ncols,4)
plt.plot(time/3600.,sw_up_TOA_totR[:,0],'-b')
plt.plot(time/3600.,sw_up_TOA_tot[:,0],'-r')
plt.axis([to,tf,0,600])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ TOA (all)')
#
plt.subplot(nrows,ncols,6)
plt.plot(time/3600.,sw_dn_SFC_clrR[:,0],'-b')
plt.plot(time/3600.,sw_dn_SFC_clr[:,0],'-r')
plt.axis([to,tf,0,1100])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.ylabel('Flux (W/m2)')
plt.title('SW down @ surface (clr)')
#
plt.subplot(nrows,ncols,7)
plt.plot(time/3600.,sw_up_SFC_clrR[:,0],'-b')
plt.plot(time/3600.,sw_up_SFC_clr[:,0],'-r')
plt.axis([to,tf,0,100])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ surface (clr)')
#
plt.subplot(nrows,ncols,8)
plt.plot(time/3600.,sw_up_TOA_clrR[:,0],'-b')
plt.plot(time/3600.,sw_up_TOA_clr[:,0],'-r')
plt.axis([to,tf,0,200])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ TOA (clr)')
#LW
plt.subplot(nrows,ncols,10)
plt.plot(time/3600.,lw_dn_SFC_totR[:,0],'-b')
plt.plot(time/3600.,lw_dn_SFC_tot[:,0],'-r')
plt.axis([to,tf,0,500])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.ylabel('Flux (W/m2)')
plt.title('LW down @ surface (all)')
#
plt.subplot(nrows,ncols,11)
plt.plot(time/3600.,lw_up_SFC_totR[:,0],'-b')
plt.plot(time/3600.,lw_up_SFC_tot[:,0],'-r')
plt.axis([to,tf,0,600])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('LW up @ surface (all)')
#
plt.subplot(nrows,ncols,12)
plt.plot(time/3600.,lw_up_TOA_totR[:,0],'-b')
plt.plot(time/3600.,lw_up_TOA_tot[:,0],'-r')
plt.axis([to,tf,0,500])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('LW up @ TOA (all)')
#
plt.subplot(nrows,ncols,14)
plt.plot(time/3600.,lw_dn_SFC_clrR[:,0],'-b')
plt.plot(time/3600.,lw_dn_SFC_clr[:,0],'-r')
plt.axis([to,tf,0,500])
plt.ylabel('Flux (W/m2)')
plt.xlabel('Forecast time (hours)')
plt.title('LW down @ surface (clr)')
#
plt.subplot(nrows,ncols,15)
plt.plot(time/3600.,lw_up_SFC_clrR[:,0],'-b')
plt.plot(time/3600.,lw_up_SFC_clr[:,0],'-r')
plt.axis([to,tf,0,600])
plt.xticks(rotation=45)
plt.xlabel('Forecast time (hours)')
plt.title('LW up @ surface (clr)')
#
plt.subplot(nrows,ncols,16)
plt.plot(time/3600.,lw_up_TOA_clrR[:,0],'-b')
plt.plot(time/3600.,lw_up_TOA_clr[:,0],'-r')
plt.axis([to,tf,0,400])
plt.xticks(rotation=45)
plt.xlabel('Forecast time (hours)')
plt.title('LW up @ TOA (clr)')

#
levels = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
plt.subplot(nrows,ncols,21)
plt.contourf(time/3600.,mean_pres[:,0]/100.,np.transpose(cldcovR[:,:,0])*100.,levels, cmap='YlGnBu')
plt.axis([to,tf,1000,0])
plt.ylabel('P (hPa)')
plt.xlabel('Forecast time (hours)')
plt.title('Cloud-cover (G)')
plt.colorbar()
#
plt.subplot(nrows,ncols,22)
plt.contourf(time/3600.,mean_pres[:,0]/100.,np.transpose(cldcov[:,:,0])*100.,levels, cmap='YlGnBu')
plt.axis([to,tf,1000,0])
plt.xlabel('Forecast time (hours)')
plt.title('Cloud-cover (GP)')
plt.colorbar()
#
levelsD = np.arange(-50, 50, 5)
plt.subplot(nrows,ncols,23)
plt.contourf(time/3600.,mean_pres[:,0]/100.,np.transpose(cldcovR[:,:,0]-cldcov[:,:,0])*100.,levelsD, cmap='seismic')
plt.axis([to,tf,1000,0])
plt.xlabel('Forecast time (hours)')
plt.title('Change Cloud-cover (G-GP)')
plt.colorbar()

#
plt.show()
fig.savefig(fileOUT, bbox_inches='tight')
plt.close()

##########################################################################
# 2) Differences
##########################################################################
fig  = plt.figure(1, figsize=(12,14), dpi=80, facecolor='w', edgecolor='k')

yo = -10
yf = 10
# Plot heating-rate profiles
# 1) SW
plt.subplot(nrows,ncols,1)
plt.plot(time/3600.,sw_dn_TOA_totR[:,0]- sw_dn_TOA_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
plt.xticks(rotation=45)
plt.ylabel('Flux (W/m2)')
plt.xlabel('Forecast time (hours)')
plt.title('SW down @ TOA (all)')
#
plt.subplot(nrows,ncols,2)
plt.plot(time/3600.,sw_dn_SFC_totR[:,0] - sw_dn_SFC_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW down @ surface (all)')
#
plt.subplot(nrows,ncols,3)
plt.plot(time/3600.,sw_up_SFC_totR[:,0] - sw_up_SFC_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ surface (all)')
#
plt.subplot(nrows,ncols,4)
plt.plot(time/3600.,sw_up_TOA_totR[:,0] - sw_up_TOA_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ TOA (all)')
#
plt.subplot(nrows,ncols,6)
plt.plot(time/3600.,sw_dn_SFC_clrR[:,0] - sw_dn_SFC_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.ylabel('Flux (W/m2)')
plt.title('SW down @ surface (clr)')
#
plt.subplot(nrows,ncols,7)
plt.plot(time/3600.,sw_up_SFC_clrR[:,0] - sw_up_SFC_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ surface (clr)')
#
plt.subplot(nrows,ncols,8)
plt.plot(time/3600.,sw_up_TOA_clrR[:,0] - sw_up_TOA_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('SW up @ TOA (clr)')
#LW
plt.subplot(nrows,ncols,10)
plt.plot(time/3600.,lw_dn_SFC_totR[:,0] - lw_dn_SFC_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.ylabel('Flux (W/m2)')
plt.title('LW down @ surface (all)')
#
plt.subplot(nrows,ncols,11)
plt.plot(time/3600.,lw_up_SFC_totR[:,0] - lw_up_SFC_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('LW up @ surface (all)')
#
plt.subplot(nrows,ncols,12)
plt.plot(time/3600.,lw_up_TOA_totR[:,0] - lw_up_TOA_tot[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
ax = plt.gca()
ax.tick_params(axis='x',labelbottom='off')
plt.title('LW up @ TOA (all)')
#
plt.subplot(nrows,ncols,14)
plt.plot(time/3600.,lw_dn_SFC_clrR[:,0] - lw_dn_SFC_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
plt.ylabel('Flux (W/m2)')
plt.xlabel('Forecast time (hours)')
plt.title('LW down @ surface (clr)')
#
plt.subplot(nrows,ncols,15)
plt.plot(time/3600.,lw_up_SFC_clrR[:,0] - lw_up_SFC_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
plt.xticks(rotation=45)
plt.xlabel('Forecast time (hours)')
plt.title('LW up @ surface (clr)')
#
plt.subplot(nrows,ncols,16)
plt.plot(time/3600.,lw_up_TOA_clrR[:,0] - lw_up_TOA_clr[:,0],'-b')
plt.plot([0,100],[0,0],'g--')
plt.axis([to,tf,yo,yf])
plt.xticks(rotation=45)
plt.xlabel('Forecast time (hours)')
plt.title('LW up @ TOA (clr)')

#
plt.show()
fig.savefig(fileOUTd, bbox_inches='tight')
plt.close()


##########################################################################################
# END PROGRAM
##########################################################################################
