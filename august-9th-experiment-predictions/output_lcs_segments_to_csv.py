#import h5py as hp
import h5py as hp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma
import csv

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
            resolution = 'c',
            area_thresh=0.,
            )

fig=plt.figure()
ridgelines = m.contour(wlon,wlat,wdirdiv,levels =[0], latlon=True)
plt.close(fig)
pp = ridgelines.collections[0].get_paths()

##OUTPUT desired contours
#from operator import itemgetter 
#index = [104,75,31,32,26] #no windage
index = [29] # windage = 0.019
v = np.empty([0,2])
for i in range(len(index)):
    if i<0:
        v = np.concatenate((v,list(reversed(pp[index[i]].vertices))))
    else:
        v = np.concatenate((v,pp[index[i]].vertices))
    #v = nvpp[index[i]].vertices

plt.figure()
plt.plot(v[:,0],v[:,1])


x,y = m(v[:,0],v[:,1],inverse=True)
del v

v = np.stack((x,y),axis=1)
with open("windage0,019.txt", "w") as f:
    writer = csv.writer(f)
    for row in v:
        writer.writerow(row)
    f.close()

fig = plt.figure()
ridgelines = m.contour(nlon,nlat,ndirdiv,levels =[0], latlon=True)
plt.close(fig)
pp = ridgelines.collections[0].get_paths()

##OUTPUT desired contours
#from operator import itemgetter 
#index = [104,75,31,32,26] #no windage
index = [81,105] # windage = 0.019
v = np.empty([0,2])
for i in range(len(index)):
    if i<0:
        v = np.concatenate((v,list(reversed(pp[index[i]].vertices))))
    else:
        v = np.concatenate((v,pp[index[i]].vertices))
    #v = nvpp[index[i]].vertices

#plt.figure()
plt.plot(v[:,0],v[:,1])

x,y = m(v[:,0],v[:,1],inverse=True)
del v

v = np.stack((x,y),axis=1)
with open("nowindage.txt", "w") as f:
    writer = csv.writer(f)
    for row in v:
        writer.writerow(row)
    f.close()


