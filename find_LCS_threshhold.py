# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 22:14:16 2018

@author: pnola
"""

import h5py as hp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import numpy.ma as ma
import time
import calendar
plt.close('all')
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
print('begin plot')
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
#"""
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],latlon=True)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],latlon=True)

tracers = np.genfromtxt('m_drifters.txt', delimiter=',')

plt.figure(2)
plt.subplot(121)
pc = m.contourf(wlon,wlat,wftle,levels=np.linspace(wftle.min(),wftle.max(),301),latlon=True,vmin=0, vmax=wftle.max(),cmap='Blues')#,alpha=0.4)
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.title("Blue: windage = 0.019; Red: No windage, 2pm EDT")
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(pc, cax=cax)

plt.subplot(122)
pc = m.contourf(nlon,nlat,nftle,levels=np.linspace(nftle.min(),nftle.max(),301),latlon=True,vmin=0, vmax=nftle.max(),cmap='Reds')#,alpha=0.4)
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(pc, cax=cax)
#plt.title(time.gmtime(t[-1]+tstart))
plt.figure(1)
plt.subplot(111)

m.drawcoastlines()
nridge = m.contour(nlon,nlat,ndirdiv,levels =[0],colors='red',latlon=True,alpha=0.6)
wridge = m.contour(wlon,wlat,wdirdiv,levels =[0],colors='blue',latlon=True,alpha=0.6)
m.scatter(tracers[:,0],tracers[:,1],latlon=True,color='black')
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
plt.title('Blue = Windage 0.019, Red = No Windage, 11am EDT')
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ax = plt.gca()
def format_coord(x, y):
    #return 'x={0[0]:.4f}, y={0[1]:.4f}'.format(m(x, y, inverse = True))

    decimal_degrees = m(x, y, inverse = True)
    degrees_lon,decimals_lon = np.divmod(decimal_degrees[0],-1)
    degrees_lat,decimals_lat = np.divmod(decimal_degrees[1],1)
    return 'lon={0}deg,{2:.4f}min ; lat={1}deg,{3:.4f}min'.format(int(degrees_lon),int(degrees_lat),decimals_lon*-60,decimals_lat*60)

ax.format_coord = format_coord
plt.show()
#"""

'''
npp = nridge.collections[0].get_paths()
wpp = wridge.collections[0].get_paths()
nindex = [104,75,31,32,26] #no windage
windex = [78,54,16,17,10,19] # windage = 0.019
nv = np.empty([0,2])
wv = np.empty([0,2])
plt.close('all')
for i in range(len(nindex)):
    if i<3:
        nv = np.concatenate((nv,list(reversed(npp[nindex[i]].vertices))))
    else:
        nv = np.concatenate((nv,npp[nindex[i]].vertices))

for i in range(len(windex)):
    if i<3:
        wv = np.concatenate((wv,list(reversed(wpp[windex[i]].vertices))))
    else:
        wv = np.concatenate((wv,wpp[windex[i]].vertices))
        
nx,ny = m(nv[:,0],nv[:,1],inverse=True)
wx,wy = m(wv[:,0],wv[:,1],inverse=True)
m.plot(nx,ny, latlon=True,color='r')
m.plot(wx,wy, latlon=True,color='b')
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

ax = plt.gca()
def format_coord(x, y):
    return 'x=%.4f, y=%.4f'%(m(x, y, inverse = True))
ax.format_coord = format_coord
plt.show()
#'''