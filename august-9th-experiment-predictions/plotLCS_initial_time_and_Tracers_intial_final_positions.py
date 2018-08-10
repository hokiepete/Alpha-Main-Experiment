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
root = Dataset('windageTracers.nc','r')
loadfile = root.variables
wtlat = loadfile['lat'][:]
wtlon = loadfile['lon'][:]
t = loadfile['time'][:]
root.close()

root = Dataset('nowindageTracers.nc','r')
loadfile = root.variables
ntlat = loadfile['lat'][:]
ntlon = loadfile['lon'][:]
#t = loadfile['time'][:]
root.close()
        
tstart = calendar.timegm(time.strptime('Aug 7, 2018 @ 12:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))

wthresh = 4
nthresh = 4

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
        t = loadfile['time'][:]
        loadfile.close()
print(time.gmtime(t[-1]+tstart))

nftle = ma.masked_where(nftle==nftle.max(axis=None),nftle)
wftle = ma.masked_where(wftle==wftle.max(axis=None),wftle)
ndirdiv = np.ma.masked_where(nconcav>0,ndirdiv)
ndirdiv = np.ma.masked_where(nftle<=nthresh,ndirdiv)
wdirdiv = np.ma.masked_where(wconcav>0,wdirdiv)
wdirdiv = np.ma.masked_where(wftle<=wthresh,wdirdiv)
nftle = ma.masked_where(nftle<nthresh,nftle)
wftle = ma.masked_where(wftle<wthresh,wftle)

nlon,nlat = np.meshgrid(nlon,nlat)
wlon,wlat = np.meshgrid(wlon,wlat)
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

pltsize = [15,12]
#pltsize = [5,4]
tup = (0,2,5,12,14,17)
#tup = (6,8,10,12,14,16)
#tup = (7,8,9,13,14,15)
tup = (0,1,2,12,13,14) 
fig = plt.figure(1,figsize=pltsize, dpi=150)
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],colors='red',latlon=True)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],colors='blue',latlon=True)

m.scatter(ntlon[tup,-1],ntlat[tup,-1], latlon=True,color='black',label="No Windage Tracers")
m.scatter(wtlon[tup,-1],wtlat[tup,-1], latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[tup,0],wtlat[tup,0], latlon=True,color='red',label="Initial Position")
'''
m.scatter(ntlon[:,-1],ntlat[:,-1], latlon=True,color='black',label="No Windage Tracers")
m.scatter(wtlon[:,-1],wtlat[:,-1], latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[:,0],wtlat[:,0], latlon=True,color='red',label="Initial Position")
#'''
plt.title('Windage LCS at 2pm EDT')
plt.legend()
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#plt.savefig('WindageLCS.png', transparent=False, bbox_inches='tight')

"""
fig = plt.figure(2,figsize=pltsize, dpi=150)
m = Basemap(llcrnrlon=-70.9,
            llcrnrlat=41.25,
            urcrnrlon=-70.7,
            urcrnrlat=41.4,
            #lat_0=(lat_max - lat_min)/2,
            #lon_0=(lon_max-lon_min)/2,
            projection='merc',
            resolution = 'f',
            area_thresh=0.,
            )

m.scatter(ntlon[:,:9],ntlat[:,:9], latlon=True,color='orange',label="No Windage Tracers")
m.scatter(wtlon[:,:9],wtlat[:,:9], latlon=True,color='cyan',label="Windage Tracers")
plt.legend()
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('WindageLCS_Closeup.png', transparent=False, bbox_inches='tight')


fig = plt.figure(3,figsize=pltsize, dpi=150)
'''
m.plot(nlon[:,0],nlat[:,0], latlon=True,color='r',label="No Windage LCS Initial")
m.plot(nlon[:,-1],nlat[:,-1], latlon=True,color='y',label="No Windage LCS Final")
m.plot(wlon[:,0],wlat[:,0], latlon=True,color='b',label="Windage LCS Initial")
m.plot(wlon[:,-1],wlat[:,-1], latlon=True,color='g',label="Windage LCS Final")
'''
m.scatter(nlon[:,0],nlat[:,0],latlon=True,color='r',label="No Windage LCS 2pm")
m.scatter(nlon[:,2],nlat[:,2],latlon=True,color='y',label="No Windage LCS 4pm")
#m.scatter(nlon[:,-1],nlat[:,-1],latlon=True,color='k',label="No Windage LCS 6pm")
m.scatter(wlon[:,0],wlat[:,0],latlon=True,color='b',label="Windage LCS 2pm")
m.scatter(wlon[:,2],wlat[:,2],latlon=True,color='g',label="Windage LCS 4pm")
#m.scatter(wlon[:,-1],wlat[:,-1],latlon=True,color='m',label="Windage LCS 6pm")
m.scatter(wtlon[:,:9],wtlat[:,:9], latlon=True,color='cyan',label="Windage Tracers")
plt.legend()
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('WindageTracers.png', transparent=False, bbox_inches='tight')


fig = plt.figure(4,figsize=pltsize, dpi=150)
'''
m.plot(nlon[:,0],nlat[:,0], latlon=True,color='r',label="No Windage LCS Initial")
m.plot(nlon[:,-1],nlat[:,-1], latlon=True,color='y',label="No Windage LCS Final")
m.plot(wlon[:,0],wlat[:,0], latlon=True,color='b',label="Windage LCS Initial")
m.plot(wlon[:,-1],wlat[:,-1], latlon=True,color='g',label="Windage LCS Final")
'''
m.scatter(nlon[:,0],nlat[:,0],latlon=True,color='r',label="No Windage LCS 2pm")
m.scatter(nlon[:,2],nlat[:,2],latlon=True,color='y',label="No Windage LCS 4pm")
#m.scatter(nlon[:,-1],nlat[:,-1],latlon=True,color='k',label="No Windage LCS 6pm")
m.scatter(wlon[:,0],wlat[:,0],latlon=True,color='b',label="Windage LCS 2pm")
m.scatter(wlon[:,2],wlat[:,2],latlon=True,color='g',label="Windage LCS 4pm")
#m.scatter(wlon[:,-1],wlat[:,-1],latlon=True,color='m',label="Windage LCS 6pm")
m.scatter(ntlon[:,:9],ntlat[:,:9], latlon=True,color='orange',label="No Windage Tracers")
plt.legend()
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('NoWindageTracers.png', transparent=False, bbox_inches='tight')



"""