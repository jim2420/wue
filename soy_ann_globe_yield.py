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
import datetime
from matplotlib.dates import DateFormatter


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf_isam.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]
gridlat=area.variables['lat'][:]
gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)
#print gridlon
nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['soybeanData'][:]
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
maitrop = region1.variables['soy_trop'][0:105,:,:]
maitemp = region1.variables['soy_temp'][0:105,:,:]
maitropi=region1.variables['soy_trop_irrig'][0:105,:,:]
maitempi=region1.variables['soy_temp_irrig'][0:105,:,:]


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
maizetropo=maitrop+maitropi
maizetempo=maitemp+maitempi

 



dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
iyield1y = N.average(dat.variables['totalyield'][95:105,:,:],axis=0)
iyield1ynew = dat.variables['totalyield'][0:105,:,:]
#print iyield1ynew.shape
latisam=dat.variables['lat'][:]
lonisam=dat.variables['lon'][:]
dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
iyield2y = N.average(dat2.variables['totalyield'][95:105,:,:],axis=0)
iyield2ynew = dat2.variables['totalyield'][0:105,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
iyield3y = N.average(dat3.variables['totalyield'][95:105,:,:],axis=0)
iyield3ynew = dat3.variables['totalyield'][0:105,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constco2_irrig_fert_0.5x0.5.nc','r')
iyield4y = N.average(dat4.variables['totalyield'][95:105,:,:],axis=0)
iyield4ynew = dat4.variables['totalyield'][0:105,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_constclim_irrig_fert_0.5x0.5.nc','r')
iyield5y = N.average(dat5.variables['totalyield'][95:105,:,:],axis=0)
iyield5ynew = dat5.variables['totalyield'][0:105,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
iyield6y = N.average(dat6.variables['totalyield'][95:105,:,:],axis=0)
iyield6ynew = dat6.variables['totalyield'][0:105,:,:]






dat=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_co2_irrig_fert_0.5x0.5_evp.nc','r')
iyield1ye = N.average(dat.variables['g_ET'][95:105,1,:,:],axis=0)
iyield1ynewe = dat.variables['g_ET'][0:105,1,:,:]

dat2=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_co2_rf_fert_0.5x0.5_evp.nc','r')
iyield2ye = N.average(dat2.variables['g_ET'][95:105,1,:,:],axis=0)
iyield2ynewe = dat2.variables['g_ET'][0:105,1,:,:]

dat3=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_co2_irrig_nofert_0.5x0.5_evp.nc','r')
iyield3ye = N.average(dat3.variables['g_ET'][95:105,1,:,:],axis=0)
iyield3ynewe = dat3.variables['g_ET'][0:105,1,:,:]

dat4=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_constco2_irrig_fert_0.5x0.5_evp.nc','r')
iyield4ye = N.average(dat4.variables['g_ET'][95:105,1,:,:],axis=0)
iyield4ynewe = dat4.variables['g_ET'][0:105,1,:,:]

dat5=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_constclim_irrig_fert_0.5x0.5_evp.nc','r')
iyield5ye = N.average(dat5.variables['g_ET'][95:105,1,:,:],axis=0)
iyield5ynewe = dat5.variables['g_ET'][0:105,1,:,:]


dat6=NetCDFFile('/scratch2/scratchdirs/tslin2/isam/cheyenne/yieldout/isamhistorical_cru/heat/fertfao/fixedirr_new1/evp/soytemp_historical_co2_rf_nofert_0.5x0.5_evp.nc','r')
iyield6ye = N.average(dat6.variables['g_ET'][95:105,1,:,:],axis=0)
iyield6ynewe = dat6.variables['g_ET'][0:105,1,:,:]


iyield1y= ma.masked_where(iyield1y<=0.,iyield1y)
iyield2y= ma.masked_where(iyield2y<=0.,iyield2y)
iyield3y= ma.masked_where(iyield3y<=0.,iyield3y)
iyield4y= ma.masked_where(iyield4y<=0.,iyield4y)
iyield5y= ma.masked_where(iyield5y<=0.,iyield5y)
iyield6y= ma.masked_where(iyield6y<=0.,iyield6y)


iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)

iyield1,lonisam1 = shiftgrid(180.5,iyield1y,lonisam,start=False)
iyield2,lonisam1 = shiftgrid(180.5,iyield2y,lonisam,start=False)
iyield3,lonisam1 = shiftgrid(180.5,iyield3y,lonisam,start=False)
iyield4,lonisam1 = shiftgrid(180.5,iyield4y,lonisam,start=False)
iyield5,lonisam1 = shiftgrid(180.5,iyield5y,lonisam,start=False)
iyield6,lonisam1 = shiftgrid(180.5,iyield6y,lonisam,start=False)



iyield1ye= ma.masked_where(iyield1ye<=0.,iyield1ye)
iyield2ye= ma.masked_where(iyield2ye<=0.,iyield2ye)
iyield3ye= ma.masked_where(iyield3ye<=0.,iyield3ye)
iyield4ye= ma.masked_where(iyield4ye<=0.,iyield4ye)
iyield5ye= ma.masked_where(iyield5ye<=0.,iyield5ye)
iyield6ye= ma.masked_where(iyield6ye<=0.,iyield6ye)


iyield1ynewe= ma.masked_where(iyield1ynewe<=0.,iyield1ynewe)
iyield2ynewe= ma.masked_where(iyield2ynewe<=0.,iyield2ynewe)
iyield3ynewe= ma.masked_where(iyield3ynewe<=0.,iyield3ynewe)
iyield4ynewe= ma.masked_where(iyield4ynewe<=0.,iyield4ynewe)
iyield5ynewe= ma.masked_where(iyield5ynewe<=0.,iyield5ynewe)
iyield6ynewe= ma.masked_where(iyield6ynewe<=0.,iyield6ynewe)

iyield1e,lonisam1 = shiftgrid(180.5,iyield1ye,lonisam,start=False)
iyield2e,lonisam1 = shiftgrid(180.5,iyield2ye,lonisam,start=False)
iyield3e,lonisam1 = shiftgrid(180.5,iyield3ye,lonisam,start=False)
iyield4e,lonisam1 = shiftgrid(180.5,iyield4ye,lonisam,start=False)
iyield5e,lonisam1 = shiftgrid(180.5,iyield5ye,lonisam,start=False)
iyield6e,lonisam1 = shiftgrid(180.5,iyield6ye,lonisam,start=False)


#print lonisam
#print lonisam1
maizeto1,lonisam2=shiftgrid(0.5,maizeto,lonisam1,start=True)
maizeto1i,lonisam2=shiftgrid(0.5,maitropi,lonisam1,start=True)
maizeto1r,lonisam2=shiftgrid(0.5,maitrop,lonisam1,start=True)
maizete1i,lonisam2=shiftgrid(0.5,maitempi,lonisam1,start=True)
maizete1r,lonisam2=shiftgrid(0.5,maitemp,lonisam1,start=True)



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



iyield1ynew=ma.filled(iyield1ynew, fill_value=0.)
iyield2ynew=ma.filled(iyield2ynew, fill_value=0.)
iyield3ynew=ma.filled(iyield3ynew, fill_value=0.)
iyield4ynew=ma.filled(iyield4ynew, fill_value=0.)
iyield5ynew=ma.filled(iyield5ynew, fill_value=0.)
iyield6ynew=ma.filled(iyield6ynew, fill_value=0.)

iyield1ynew= ma.masked_where(iyield1ynew<=0.,iyield1ynew)
iyield2ynew= ma.masked_where(iyield2ynew<=0.,iyield2ynew)
iyield3ynew= ma.masked_where(iyield3ynew<=0.,iyield3ynew)
iyield4ynew= ma.masked_where(iyield4ynew<=0.,iyield4ynew)
iyield5ynew= ma.masked_where(iyield5ynew<=0.,iyield5ynew)
iyield6ynew= ma.masked_where(iyield6ynew<=0.,iyield6ynew)



iyield1e=ma.filled(iyield1e, fill_value=0.)
iyield2e=ma.filled(iyield2e, fill_value=0.)
iyield3e=ma.filled(iyield3e, fill_value=0.)
iyield4e=ma.filled(iyield4e, fill_value=0.)
iyield5e=ma.filled(iyield5e, fill_value=0.)
iyield6e=ma.filled(iyield6e, fill_value=0.)

iyield1e= ma.masked_where(iyield1e<=0.,iyield1e)
iyield2e= ma.masked_where(iyield2e<=0.,iyield2e)
iyield3e= ma.masked_where(iyield3e<=0.,iyield3e)
iyield4e= ma.masked_where(iyield4e<=0.,iyield4e)
iyield5e= ma.masked_where(iyield5e<=0.,iyield5e)
iyield6e= ma.masked_where(iyield6e<=0.,iyield6e)

iyield1ynewe=ma.filled(iyield1ynewe, fill_value=0.)
iyield2ynewe=ma.filled(iyield2ynewe, fill_value=0.)
iyield3ynewe=ma.filled(iyield3ynewe, fill_value=0.)
iyield4ynewe=ma.filled(iyield4ynewe, fill_value=0.)
iyield5ynewe=ma.filled(iyield5ynewe, fill_value=0.)
iyield6ynewe=ma.filled(iyield6ynewe, fill_value=0.)

iyield1ynewe= ma.masked_where(iyield1ynewe<=0.,iyield1ynewe)
iyield2ynewe= ma.masked_where(iyield2ynewe<=0.,iyield2ynewe)
iyield3ynewe= ma.masked_where(iyield3ynewe<=0.,iyield3ynewe)
iyield4ynewe= ma.masked_where(iyield4ynewe<=0.,iyield4ynewe)
iyield5ynewe= ma.masked_where(iyield5ynewe<=0.,iyield5ynewe)
iyield6ynewe= ma.masked_where(iyield6ynewe<=0.,iyield6ynewe)




maizeto2=N.zeros((105,360,720))
maizeto2r=N.zeros((105,360,720))
maizeto2i=N.zeros((105,360,720))
maizete2r=N.zeros((105,360,720))
maizete2i=N.zeros((105,360,720))

for i in range(0,105):
	for x in range(0,360):
		for y in range(0,720):
			maizeto2[i,x,y]=maizeto1[i,x,y]
                        maizeto2r[i,x,y]=maizeto1r[i,x,y]
                        maizeto2i[i,x,y]=maizeto1i[i,x,y]
                        maizete2r[i,x,y]=maizete1r[i,x,y]
                        maizete2i[i,x,y]=maizete1i[i,x,y]
##yield (t/ha)
allynew=N.average(iyield1ynew,weights=maizeto2,axis=(1,2))
allyinew=N.average(iyield2ynew,weights=maizeto2,axis=(1,2))
allynnew=N.average(iyield3ynew,weights=maizeto2,axis=(1,2))
allycnew=N.average(iyield4ynew,weights=maizeto2,axis=(1,2))
allyclinew=N.average(iyield5ynew,weights=maizeto2,axis=(1,2))
allyinnew=N.average(iyield6ynew,weights=maizeto2,axis=(1,2))

##et (mm)
allynewe=N.average(iyield1ynewe,weights=maizeto2,axis=(1,2))
allyinewe=N.average(iyield2ynewe,weights=maizeto2,axis=(1,2))
allynnewe=N.average(iyield3ynewe,weights=maizeto2,axis=(1,2))
allycnewe=N.average(iyield4ynewe,weights=maizeto2,axis=(1,2))
allyclinewe=N.average(iyield5ynewe,weights=maizeto2,axis=(1,2))
allyinnewe=N.average(iyield6ynewe,weights=maizeto2,axis=(1,2))


##wue (kg/ha/mm)
allynewew=N.average(iyield1ynew*1000/iyield1ynewe,weights=maizeto2,axis=(1,2))
allyinewew=N.average(iyield2ynew*1000/iyield2ynewe,weights=maizeto2,axis=(1,2))
allynnewew=N.average(iyield3ynew*1000/iyield3ynewe,weights=maizeto2,axis=(1,2))
allycnewew=N.average(iyield4ynew*1000/iyield4ynewe,weights=maizeto2,axis=(1,2))
allyclinewew=N.average(iyield5ynew*1000/iyield5ynewe,weights=maizeto2,axis=(1,2))
allyinnewew=N.average(iyield6ynew*1000/iyield6ynewe,weights=maizeto2,axis=(1,2))

#finalnew=((N.average(allynew)-N.average(allyinnew))/N.average(allyinnew)*100,(N.average(allynnew)-N.average(allyinnew))/N.average(allyinnew)*100,(N.average(allyinew)-N.average(allyinnew))/N.average(allyinnew)*100)

##yield relative change 
finalnew=(allynew-allyinnew)/allynew*100
finalnnew=(allynew-allynnew)/allynew*100
finalinew=(allynew-allyinew)/allynew*100
finalclinew=(allynew-allyclinew)/allynew*100
finalcnew=(allynew-allycnew)/allynew*100

##water use efficiency relative change 
wfinalnew=(allynewe-allyinnewe)/allynewe*100
wfinalnnew=(allynewe-allynnewe)/allynewe*100
wfinalinew=(allynewe-allyinewe)/allynewe*100
wfinalclinew=(allynewe-allyclinewe)/allynewe*100
wfinalcnew=(allynewe-allycnewe)/allynewe*100




#finalstd=(N.std(allyinew)/N.std(allynew),N.std(allynnew)/N.std(allynew),N.std(allycnew)/N.std(allynew),N.std(allyclinew)/N.std(allynew))

fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(212)
xx=range(1901,2006)
xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
ax.plot_date(xdates,finalcnew,"k-",label=r"$\Delta$ CO$_{2}$",linewidth='3')
ax.plot_date(xdates,finalclinew,"r-",label=r"$\Delta$ Climate",linewidth='3')
ax.plot_date(xdates,finalnnew,"g-",label=r"$\Delta$ NF",linewidth='3')
ax.plot_date(xdates,finalinew,"b-",label=r"$\Delta$ Irrigation",linewidth='3')
plt.xlim(xdates[1],xdates[104])

leg = plt.legend(loc=4,fancybox=True, fontsize=18)
leg.get_frame().set_alpha(0.5)

#ax.set_ylim([0,14])
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
plt.tick_params(axis='both',labelsize=18)

plt.xlabel("Year",fontsize=18)

plt.ylabel('Effect on soybean yield (%)',fontsize=18)


ax = fig.add_subplot(211)

params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
ax.plot_date(xdates,allynew,"y-",label=" $S_{all}$",linewidth='3')
ax.plot_date(xdates,allycnew,"k-",label=" $S_{CO2}$",linewidth='3')
ax.plot_date(xdates,allyclinew,"r-",label="$S_{Climate}$",linewidth='3')
ax.plot_date(xdates,allynnew,"g-",label=" $S_{NF}$",linewidth='3')
ax.plot_date(xdates,allyinew,"b-",label=" $S_{Irrigation}$",linewidth='3')
plt.xlim(xdates[1],xdates[104])

leg = plt.legend(loc=2,fancybox=True, fontsize=18)
leg.get_frame().set_alpha(0.5)

#ax.set_ylim([0,14])
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
plt.tick_params(axis='both',labelsize=18)

plt.xlabel("Year",fontsize=18)

plt.ylabel('Soybean yield (t/ha)',fontsize=18)




plt.tight_layout()

plt.savefig('soybeanhis_yield_1902_2005.png')
plt.show()
