# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 23:34:03 2018

@author: pnola
"""

import h5py as hp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma
plt.close('all')

wthresh = 4
nthresh = 5

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


#ridgelines = m.contour(nlon,nlat,ndirdiv,levels =[0], latlon=True)
ridgelines = m.contour(wlon,wlat,wdirdiv,levels =[0], latlon=True)
pp = ridgelines.collections[0].get_paths()
###View contours
plt.close('all')
for p in range(len(pp)):
    v = pp[p].vertices
    x = v[:,0]
    y = v[:,1]
    if x.size > 10:
        m.plot(x,y)#, latlon=True)
        m.drawcoastlines()
        plt.title('{:03d}'.format(p))
        plt.savefig('{:03d}.png'.format(p))
        plt.close('all')
        
plt.show()        

