# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:17:52 2018

@author: pnola
"""
from netCDF4 import Dataset
import numpy as np
import time
import calendar
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
plt.close('all')
tstart = calendar.timegm(time.strptime('Aug 7, 2018 @ 12:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
wc = 1#0.019 #Windage Coefficent

ncfile="MIT_nsf_alpha200m_2018080700_2018080800_2018081000_01h-depth-00m-wind.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t=15
dx=200 #m
dy=200 #m
lat = vars["lat"][:]
lon = vars["lon"][:]
data_time = vars["time"][:]+tstart
print(time.gmtime(data_time[t]))
u_sea = vars["East_vel"][t,:,:]
v_sea = vars["North_vel"][t,:,:]
u_wind = vars["East_wind"][t,:,:]
v_wind = vars["North_wind"][t,:,:]
u = u_sea + wc*u_wind
v = v_sea + wc*v_wind
#u = u_wind
#v = v_wind

root.close()
ydim = lat.shape[0]
xdim = lon.shape[0]
lon_min = lon.min()
lon_max = lon.max()
lat_min = lat.min()
lat_max = lat.max()
lon, lat = np.meshgrid(lon,lat)
dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)
s1 = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = np.min(np.linalg.eig(S)[0])

        else:
            s1[i,j] = np.ma.masked
s1_windage = 3600*s1 #convert to hr^-1
del u,v,s1
u = u_sea
v = v_sea
dudy,dudx = np.gradient(u,dy,dx,edge_order=2)
dvdy,dvdx = np.gradient(v,dy,dx,edge_order=2)
s1 = np.ma.empty([ydim,xdim])
J = np.array([[0, 1], [-1, 0]])
for i in range(ydim):
    for j in range(xdim):
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j] and u[i,j] and v[i,j]) is not np.ma.masked:    
            Utemp = np.array([u[i, j], v[i, j]])
            Grad = np.array([[dudx[i, j], dudy[i, j]], [dvdx[i, j], dvdy[i, j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            s1[i,j] = np.min(np.linalg.eig(S)[0])

        else:
            s1[i,j] = np.ma.masked
s1_nowindage = 3600*s1 #convert to hr^-1
del s1
print('begin plot')
#'''
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'h',
            area_thresh=0.,
            )
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
'''
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
#''
root = Dataset('m_drifters_windage.nc','r')
loadfile = root.variables
wtlat = loadfile['lat'][:]
wtlon = loadfile['lon'][:]
t = loadfile['time'][:]
root.close()

root = Dataset('m_drifters_nowindage.nc','r')
loadfile = root.variables
ntlat = loadfile['lat'][:]
ntlon = loadfile['lon'][:]
#t = loadfile['time'][:]
root.close()
'''
tracers = np.genfromtxt('m_drifters.txt', delimiter=',')
date_time = time.gmtime(data_time[t])
#timestamp = str(date_time.tm_mon)+'/'+str(date_time.tm_mday)+'/'+str(date_time.tm_year)+' @ {0:02d}'.format(date_time.tm_hour)+'
timestamp = "{0:02d}/{1:02d}/{2:02d} @ {3:02d}:{4:02d} UTC".format(date_time.tm_mon,date_time.tm_mday,date_time.tm_year,date_time.tm_hour,date_time.tm_min)
print(timestamp)
pltsize = [15,12]
fig = plt.figure(1,figsize=pltsize, dpi=150)
plt.subplot(1,2,1)
sc = m.contourf(lon,lat,s1_nowindage,levels=np.linspace(s1_nowindage.min(axis=None),s1_nowindage.max(axis=None),301),latlon=True)
m.scatter(tracers[:,0],tracers[:,1],latlon=True,color='orange')
'''
m.scatter(ntlon,ntlat, latlon=True,color='orange',label="No Windage Tracers")
m.scatter(wtlon,wtlat, latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[:,0],wtlat[:,0], latlon=True,color='red',label="Initial Position")
'''
plt.legend()
plt.title("sea s$_{1}$ "+timestamp)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(sc, cax=cax)
#plt.savefig('No_Windage_S1_Closeup.png', transparent=False, bbox_inches='tight')


#fig = plt.figure(2,figsize=pltsize, dpi=150)
plt.subplot(1,2,2)
sc = m.contourf(lon,lat,s1_windage,levels=np.linspace(s1_windage.min(axis=None),s1_windage.max(axis=None),301),latlon=True)
m.scatter(tracers[:,0],tracers[:,1],latlon=True,color='orange')

'''
m.scatter(ntlon,ntlat, latlon=True,color='orange',label="No Windage Tracers")
m.scatter(wtlon,wtlat, latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[:,0],wtlat[:,0], latlon=True,color='red',label="Initial Position")
'''
plt.legend()
plt.title("wind s$_{1}$ "+timestamp)
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(sc, cax=cax)
#plt.savefig('Windage_S1_{0:02d}.png'.format(t), transparent=False, bbox_inches='tight')