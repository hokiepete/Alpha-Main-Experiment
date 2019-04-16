import h5py as hp
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import numpy.ma as ma
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import time
import calendar
plt.close('all')
tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
timestep=0
lat_max = 41.4548797564
lat_min = 41.0250244184
lon_max = -70.2909317079
lon_min = -70.916465753
        
tstart = calendar.timegm(time.strptime('Aug 7, 2018 @ 12:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

wthresh = 4
nthresh = 4
lthresh = 4

#wthresh = 0
#nthresh = 0

with hp.File('windageLCS.hdf5','r') as loadfile:
        nftle = loadfile['nftle'][:]
        nconcav = loadfile['nconcavity'][:]
        ndirdiv = loadfile['ndirectionalderivative'][:]
        nlon = loadfile['nlon'][:]
        nlat = loadfile['nlat'][:]
        wftle = loadfile['wftle'][:]
        wconcav = loadfile['wconcavity'][:]
        wdirdiv = loadfile['wdirectionalderivative'][:]
        wlon = loadfile['wlon'][:]
        wlat = loadfile['wlat'][:]
        lftle = loadfile['lftle'][:]
        lconcav = loadfile['lconcavity'][:]
        ldirdiv = loadfile['ldirectionalderivative'][:]
        llon = loadfile['llon'][:]
        llat = loadfile['llat'][:]
        t = loadfile['time'][:]
        loadfile.close()
print(time.gmtime(t[-1]+tstart))

nftle = ma.masked_where(nftle==nftle.max(axis=None),nftle)
wftle = ma.masked_where(wftle==wftle.max(axis=None),wftle)
lftle = ma.masked_where(lftle==lftle.max(axis=None),lftle)

ndirdiv = np.ma.masked_where(nconcav>0,ndirdiv)
ndirdiv = np.ma.masked_where(nftle<=nthresh,ndirdiv)
wdirdiv = np.ma.masked_where(wconcav>0,wdirdiv)
wdirdiv = np.ma.masked_where(wftle<=wthresh,wdirdiv)
ldirdiv = np.ma.masked_where(lconcav>0,ldirdiv)
ldirdiv = np.ma.masked_where(lftle<=lthresh,ldirdiv)

nftle = ma.masked_where(nftle<nthresh,nftle)
wftle = ma.masked_where(wftle<wthresh,wftle)
lftle = ma.masked_where(lftle<lthresh,lftle)


nlon,nlat = np.meshgrid(nlon,nlat)
wlon,wlat = np.meshgrid(wlon,wlat)
llon,llat = np.meshgrid(llon,llat)
lon_min = wlon.min()
lon_max = wlon.max()
lat_min = wlat.min()
lat_max = wlat.max()
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'h',
            area_thresh=0.,
            )

width=7
height=width
pltsize = [width,height]
#pltsize = [5,4]
tup = (0,2,5,12,14,17)
#tup = (6,8,10,12,14,16)
#tup = (7,8,9,13,14,15)
tup = (0,1,2,12,13,14) 
fig = plt.figure(1,figsize=pltsize, dpi=150)
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],colors='red',latlon=True)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],colors='blue',latlon=True)
lridge = m.contour(llon,llat,ldirdiv,levels =[0],colors='green',latlon=True)
h1,_ = nridge.legend_elements()
h2,_ = wridge.legend_elements()
h3,_ = lridge.legend_elements()
plt.legend([h1[0], h2[0], h3[0]], ['No Windage','0.019 Windage', '0.009 Windage'],loc='lower right')
#plt.legend((nridge,wridge,lridge),('No Windage','0.019 Windage', '0.009 Windage'))
plt.title('Windage LCS at 11am EDT')
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('WindageLCS.png', transparent=False, bbox_inches='tight')

fig = plt.figure(2,figsize=pltsize, dpi=150)
m = Basemap(llcrnrlon=-70.6,
            llcrnrlat=41.2,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],colors='red',latlon=True)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],colors='blue',latlon=True)
lridge = m.contour(llon,llat,ldirdiv,levels =[0],colors='green',latlon=True)
h1,_ = nridge.legend_elements()
h2,_ = wridge.legend_elements()
h3,_ = lridge.legend_elements()
plt.legend([h1[0], h2[0], h3[0]], ['No Windage','0.019 Windage', '0.009 Windage'],loc='lower right')
#plt.legend((nridge,wridge,lridge),('No Windage','0.019 Windage', '0.009 Windage'))
plt.title('Windage LCS at 11am EDT')
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('close_WindageLCS.png', transparent=False, bbox_inches='tight')
