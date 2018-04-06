from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
from statsmodels.stats.weightstats import DescrStatsW
from scipy.stats import ttest_ind
from matplotlib.markers import TICKDOWN
import matplotlib.patches as mpatches

area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['maizeData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
latnc=N.flipud(latnc)


ncvar_maizef= N.zeros((2160, 4320))
ncvar_maizef=ncvar_maize[0,1,:,:]
ncvar_maize1=ncvar_maize[0,4,:,:]
ncvar_mask= N.zeros((2160, 4320))
ncvar_mask=ncvar_maize[0,0,:,:]


ncvar_maizef[N.isnan(ncvar_maizef)] = -9999
ncvar_mask[N.isnan(ncvar_mask)] = -9999
ncvar_maize1[N.isnan(ncvar_maize1)] = -9999

ncvar_maizef = ma.masked_where(ncvar_maizef<=0,ncvar_maizef)
#ncvar_maizef= ma.masked_where(ncvar_mask<0.01,ncvar_maizef)
ncvar_maizef = ma.masked_where(ncvar_maize1<=0,ncvar_maizef)

ncvar_maize1=N.flipud(ncvar_maize1)
ncvar_maizef=N.flipud(ncvar_maizef)
ncvar_mask=N.flipud(ncvar_mask)


lon2,lat2 = N.meshgrid(gridlon,gridlat)

iyield = interp(ncvar_maizef,lonnc,latnc,lon2,lat2  , order=1)
iarea=interp(ncvar_mask,lonnc,latnc,lon2,lat2  , order=1)


region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
#maitrop = region1.variables['soy_trop'][99,:,:]
#maitemp = region1.variables['soy_temp'][99,:,:]
#maitropi=region1.variables['soy_trop_irrig'][99,:,:]
#maitempi=region1.variables['soy_temp_irrig'][99,:,:]

maitrop = N.average(region1.variables['maize_trop'][95:105,:,:],axis=0)
maitemp = N.average(region1.variables['maize_temp'][95:105,:,:],axis=0)
maitropi=N.average(region1.variables['maize_trop_irrig'][95:105,:,:],axis=0)
maitempi=N.average(region1.variables['maize_temp_irrig'][95:105,:,:],axis=0)

gridarea = region1.variables['area'][:,:]
maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)

maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)
maizetor=maitrop+maitemp
maizetoi=maitropi+maitempi
maizeto = maitrop+maitemp+maitropi+maitempi

 




dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1y = N.average(dat.variables['totalyield'][95:105,:,:],axis=0)
iyield1yn = dat.variables['totalyield'][95:105,:,:]

latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2y = N.average(dat2.variables['totalyield'][95:105,:,:],axis=0)
iyield2yn = dat2.variables['totalyield'][95:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3y = N.average(dat3.variables['totalyield'][95:105,:,:],axis=0)
iyield3yn = dat3.variables['totalyield'][95:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constco2_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield4y = N.average(dat4.variables['totalyield'][95:105,:,:],axis=0)
iyield4yn = dat4.variables['totalyield'][95:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/maizetemp_historical_constclim_irrig_fert_0.5x0.5.nc_0.5x0.5.nc','r')
iyield5y = N.average(dat5.variables['totalyield'][95:105,:,:],axis=0)
iyield5yn = dat5.variables['totalyield'][95:105,:,:]

dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield6y = N.average(dat6.variables['totalyield'][95:105,:,:],axis=0)



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/maizetemp_historical_co2_irrig_fert_0.5x0.5_evp.nc','r')
iyield1ye = N.average(dat.variables['g_ET'][95:105,0,:,:],axis=0)
iyield1yne = dat.variables['g_ET'][95:105,0,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/maizetemp_historical_co2_rf_fert_0.5x0.5_evp.nc','r')
iyield2ye = N.average(dat2.variables['g_ET'][95:105,0,:,:],axis=0)
iyield2yne = dat2.variables['g_ET'][95:105,0,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/maizetemp_historical_co2_irrig_nofert_0.5x0.5_evp.nc','r')
iyield3ye = N.average(dat3.variables['g_ET'][95:105,0,:,:],axis=0)
iyield3yne = dat3.variables['g_ET'][95:105,0,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/maizetemp_historical_constco2_irrig_fert_0.5x0.5_evp.nc','r')
iyield4ye = N.average(dat4.variables['g_ET'][95:105,0,:,:],axis=0)
iyield4yne = dat4.variables['g_ET'][95:105,0,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/maizetemp_historical_constclim_irrig_fert_0.5x0.5_evp.nc','r')
iyield5ye = N.average(dat5.variables['g_ET'][95:105,0,:,:],axis=0)
iyield5yne = dat5.variables['g_ET'][95:105,0,:,:]


iyield1y= ma.masked_where(iyield1y<=0.,iyield1y)
iyield2y= ma.masked_where(iyield2y<=0.,iyield2y)
iyield3y= ma.masked_where(iyield3y<=0.,iyield3y)
iyield4y= ma.masked_where(iyield4y<=0.,iyield4y)
iyield5y= ma.masked_where(iyield5y<=0.,iyield5y)
iyield6y= ma.masked_where(iyield6y<=0.,iyield6y)

iyield1,lonisam1 = shiftgrid(180.5,iyield1y,lonisam,start=False)
iyield2,lonisam1 = shiftgrid(180.5,iyield2y,lonisam,start=False)
iyield3,lonisam1 = shiftgrid(180.5,iyield3y,lonisam,start=False)
iyield4,lonisam1 = shiftgrid(180.5,iyield4y,lonisam,start=False)
iyield5,lonisam1 = shiftgrid(180.5,iyield5y,lonisam,start=False)
iyield6,lonisam1 = shiftgrid(180.5,iyield6y,lonisam,start=False)

iyield1= ma.masked_where(maizeto<=0.,iyield1)
iyield2= ma.masked_where(maizeto<=0.,iyield2)
iyield3= ma.masked_where(maizeto<=0.,iyield3)
iyield4= ma.masked_where(maizeto<=0.,iyield4)
iyield5= ma.masked_where(maizeto<=0.,iyield5)
iyield6= ma.masked_where(maizeto<=0.,iyield6)

iyield1=ma.filled(iyield1, fill_value=0.)
iyield2=ma.filled(iyield2, fill_value=0.)
iyield3=ma.filled(iyield3, fill_value=0.)
iyield4=ma.filled(iyield4, fill_value=0.)
iyield5=ma.filled(iyield5, fill_value=0.)
iyield6=ma.filled(iyield6, fill_value=0.)

iyield1= ma.masked_where(iyield1<=0.,iyield1)
iyield2= ma.masked_where(iyield2<=0.,iyield2)
iyield3= ma.masked_where(iyield3<=0.,iyield3)
iyield4= ma.masked_where(iyield4<=0.,iyield4)
iyield5= ma.masked_where(iyield5<=0.,iyield5)
iyield6= ma.masked_where(iyield6<=0.,iyield6)








iyield1ye= ma.masked_where(iyield1ye<=0.,iyield1ye)
iyield2ye= ma.masked_where(iyield2ye<=0.,iyield2ye)
iyield3ye= ma.masked_where(iyield3ye<=0.,iyield3ye)
iyield4ye= ma.masked_where(iyield4ye<=0.,iyield4ye)
iyield5ye= ma.masked_where(iyield5ye<=0.,iyield5ye)

iyield1e,lonisam1 = shiftgrid(180.5,iyield1ye,lonisam,start=False)
iyield2e,lonisam1 = shiftgrid(180.5,iyield2ye,lonisam,start=False)
iyield3e,lonisam1 = shiftgrid(180.5,iyield3ye,lonisam,start=False)
iyield4e,lonisam1 = shiftgrid(180.5,iyield4ye,lonisam,start=False)
iyield5e,lonisam1 = shiftgrid(180.5,iyield5ye,lonisam,start=False)

iyield1e= ma.masked_where(maizeto<=0.,iyield1e)
iyield2e= ma.masked_where(maizeto<=0.,iyield2e)
iyield3e= ma.masked_where(maizeto<=0.,iyield3e)
iyield4e= ma.masked_where(maizeto<=0.,iyield4e)
iyield5e= ma.masked_where(maizeto<=0.,iyield5e)

iyield1e=ma.filled(iyield1e, fill_value=0.)
iyield2e=ma.filled(iyield2e, fill_value=0.)
iyield3e=ma.filled(iyield3e, fill_value=0.)
iyield4e=ma.filled(iyield4e, fill_value=0.)
iyield5e=ma.filled(iyield5e, fill_value=0.)

iyield1e= ma.masked_where(iyield1e<=0.,iyield1e)
iyield2e= ma.masked_where(iyield2e<=0.,iyield2e)
iyield3e= ma.masked_where(iyield3e<=0.,iyield3e)
iyield4e= ma.masked_where(iyield4e<=0.,iyield4e)
iyield5e= ma.masked_where(iyield5e<=0.,iyield5e)


wue1=iyield1/iyield1e
wue2=iyield2/iyield2e
wue3=iyield3/iyield3e
wue4=iyield4/iyield4e
wue5=iyield5/iyield5e


lon2,lat2 = N.meshgrid(lonisam1,latisam)
cmap = plt.cm.bwr
bounds=[-12,-10,-8,-6,-4,-2,0,8,16,24,32,40,48]
norm = colors.BoundaryNorm(bounds, cmap.N)



fig = plt.figure(figsize=(10,10))
ax2 = fig.add_subplot(411)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

x,y = map(lon2,lat2)


map.drawcoastlines(color='gray')

#aa=(iyield1-iyield2)/iyield1*100
aa=(wue1-wue2)/wue1*100
aa= ma.masked_where(aa==0,aa)

cs1 = map.pcolormesh(x,y,aa,cmap=plt.cm.bwr,norm=norm)

#cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",extend='both')

#cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(412)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines(color='grey')
map.drawmapboundary()
#bb=(iyield1-iyield3)/iyield1*100
#bb= ma.masked_where(bb==0,bb)
bb=(wue1-wue3)/wue1*100
bb= ma.masked_where(bb==0,bb)

cs1 = map.pcolormesh(x,y,bb,cmap=plt.cm.bwr,norm=norm)

#cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",extend='both')
#cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(413)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines(color='grey')
map.drawmapboundary()
#cc=(iyield1-iyield4)/iyield1*100
cc=(wue1-wue4)/wue1*100
cc= ma.masked_where(cc==0,cc)


cs1 = map.pcolormesh(x,y,cc,cmap=plt.cm.bwr,norm=norm)
#cbar = map.colorbar(cs1,location='right',size="5%",pad="2%",extend='both')
#cbar.ax.tick_params(labelsize=15)

plt.axis('off')

ax2 = fig.add_subplot(414)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')

map.drawcoastlines(color='grey')
#map.drawcountries()
map.drawmapboundary()

#dd=(iyield1-iyield5)/iyield1*100
#dd= ma.masked_where(dd==0,dd)
dd=(wue1-wue5)/wue5*100
dd= ma.masked_where(dd==0,dd)


cs1 = map.pcolormesh(x,y,dd,cmap=plt.cm.bwr,norm=norm)

cbar = map.colorbar(cs1,location='bottom',size="8%",pad="4%",extend='both')
cbar.ax.tick_params(labelsize=16)




plt.axis('off')
#fig.colorbar(cs1)

aa= ma.masked_where(maizeto<=0,aa)
bb= ma.masked_where(maizeto<=0,bb)
cc= ma.masked_where(maizeto<=0,cc)
dd= ma.masked_where(maizeto<=0,dd)
aa=ma.filled(aa, fill_value=0.)
bb=ma.filled(bb, fill_value=0.)
cc=ma.filled(cc, fill_value=0.)
dd=ma.filled(dd, fill_value=0.)

ax=0
alla=0
allb=0
allc=0
alld=0
for i in range(0,360):
	for j in range(0,720):
		if (maizeto[i,j]>0):
			alla=aa[i,j]+alla
			allb=bb[i,j]+allb
                        allc=cc[i,j]+allc
                        alld=dd[i,j]+alld

			ax=ax+1
da=alla/ax
db=allb/ax
dc=allc/ax
dd1=alld/ax
print da,db,dc,dd1
#print aa.shape,maizeto.shape
print N.average(aa),N.std(aa),N.average(bb),N.std(bb),N.average(cc),N.std(cc),N.average(dd),N.std(dd)
print N.average(aa,weights=maizeto),N.average(bb,weights=maizeto),N.average(cc,weights=maizeto),N.average(dd,weights=maizeto)
#weighted_stats = DescrStatsW(aa, weights=maizeto, ddof=0)
#print weighted_stats.mean, weighted_stats.std,weighted_stats.std_mean


plt.savefig('maidiff_1996_2005_wue.jpg',dpi=300,bbox_inches='tight')
plt.show()

