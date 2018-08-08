#import h5py as hp
from netCDF4 import Dataset
import numpy.ma as ma
import numpy as np
import time
import calendar
tstart = calendar.timegm(time.strptime('Jul 13, 2018 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
import h5py as hp
timestep=0
dx = 200
dy = 200
root = Dataset('no_windage.nc','r')
londfile = root.variables
nlat = londfile['lat'][:]
nlon = londfile['lon'][:]
nftle = londfile['FTLE'][timestep,:,:,0]
ncg = londfile['CauchyGreen'][timestep,:,:,:,:]
t = londfile['time'][:]
print(time.gmtime(t[0]+tstart))
ydim , xdim = nftle.shape
#print(time.gmtime(t[0]*24*60*60+tstart))
nftle = ma.masked_where(nftle==999,nftle)

root = Dataset('windage=0,019.nc','r')
londfile = root.variables
wlat = londfile['lat'][:]
wlon = londfile['lon'][:]
wftle = londfile['FTLE'][timestep,:,:,0]
wcg = londfile['CauchyGreen'][timestep,:,:,:,:]
wftle = ma.masked_where(wftle==999,wftle)


ndfdy,ndfdx = np.gradient(nftle,dy,dx,edge_order=2)
ndfdydy,ndfdydx = np.gradient(ndfdy,dy,dx,edge_order=2)
ndfdxdy,ndfdxdx = np.gradient(ndfdx,dy,dx,edge_order=2)
wdfdy,wdfdx = np.gradient(wftle,dy,dx,edge_order=2)
wdfdydy,wdfdydx = np.gradient(wdfdy,dy,dx,edge_order=2)
wdfdxdy,wdfdxdx = np.gradient(wdfdx,dy,dx,edge_order=2)

ndirdiv = np.ma.empty([ydim,xdim])
wdirdiv = np.ma.empty([ydim,xdim])
nconcav = np.ma.empty([ydim,xdim])
wconcav = np.ma.empty([ydim,xdim])
for i in range(ydim):
    print(str(i/ydim*100)+" % done")
    for j in range(xdim):
        if (ndfdx[i,j] and ndfdy[i,j] and ndfdxdy[i,j] and ndfdydy[i,j] and ndfdxdx[i,j] and ndfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(ncg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            ndirdiv[i,j] = np.dot([ndfdx[i,j],ndfdy[i,j]],eigenVectors[:,0])
            nconcav[i,j] = np.dot(np.dot([[ndfdxdx[i,j],ndfdxdy[i,j]],[ndfdydx[i,j],ndfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        else:
            ndirdiv[i,j] = np.ma.masked
            nconcav[i,j] = np.ma.masked

        if (wdfdx[i,j] and wdfdy[i,j] and wdfdxdy[i,j] and wdfdydy[i,j] and wdfdxdx[i,j] and wdfdydx[i,j]) is not np.ma.masked:    
            eigenValues, eigenVectors = np.linalg.eig(wcg[i,j,:,:])
            idx = eigenValues.argsort()[::-1]   
            eigenVectors = eigenVectors[:,idx]
            wdirdiv[i,j] = np.dot([wdfdx[i,j],wdfdy[i,j]],eigenVectors[:,0])
            wconcav[i,j] = np.dot(np.dot([[wdfdxdx[i,j],wdfdxdy[i,j]],[wdfdydx[i,j],wdfdydy[i,j]]],eigenVectors[:,0]),eigenVectors[:,0])
        else:
           wdirdiv[i,j] = np.ma.masked
           wconcav[i,j] = np.ma.masked

with hp.File('windageLCS.hdf5','w') as savefile:
        savefile.create_dataset('nftle',shape=nftle.shape,data=nftle)
        savefile.create_dataset('nconcavity',shape=nconcav.shape,data=nconcav)
        savefile.create_dataset('ndirectionalderivative',shape=ndirdiv.shape,data=ndirdiv)
        savefile.create_dataset('nlon',shape=nlon.shape,data=nlon)
        savefile.create_dataset('nlat',shape=nlat.shape,data=nlat)
        savefile.create_dataset('wftle',shape=wftle.shape,data=wftle)
        savefile.create_dataset('wconcavity',shape=wconcav.shape,data=wconcav)
        savefile.create_dataset('wdirectionalderivative',shape=wdirdiv.shape,data=wdirdiv)
        savefile.create_dataset('wlon',shape=wlon.shape,data=wlon)
        savefile.create_dataset('wlat',shape=wlat.shape,data=wlat)
        savefile.create_dataset('time',shape=t.shape,data=t)
        savefile.close()
