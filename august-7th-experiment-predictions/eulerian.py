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
tstart = calendar.timegm(time.strptime('Aug 5, 2018 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
wc = 0.019 #Windage Coefficent

ncfile="MIT_nsf_alpha200m_2018080500_2018080600_2018080800_01h-depth-00m-wind.nc"
root = Dataset(ncfile,'r') #read the data
vars = root.variables #dictionary, all variables in dataset\
print(vars.keys())
t=18
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
s1_windage = s1
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
s1_nowindage = s1
del s1
print('begin plot')
'''
m = Basemap(llcrnrlon=lon_min,
            llcrnrlat=lat_min,
            urcrnrlon=lon_max,
            urcrnrlat=lat_max,
            projection='merc',
            resolution = 'i',
            area_thresh=0.,
            )
#'''

parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
pltsize = [15,12]
#'''
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
'''
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

fig = plt.figure(1,figsize=pltsize, dpi=150)
sc = m.contourf(lon,lat,s1_nowindage,levels=np.linspace(s1_nowindage.min(axis=None),s1_nowindage.max(axis=None),301),latlon=True)
m.scatter(ntlon,ntlat, latlon=True,color='orange',label="No Windage Tracers")
m.scatter(wtlon,wtlat, latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[:,0],wtlat[:,0], latlon=True,color='red',label="Initial Position")
plt.legend()
plt.title("s$_{1}$ at 2pm")
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(sc, cax=cax)
plt.savefig('No_Windage_S1_Closeup.png', transparent=False, bbox_inches='tight')


fig = plt.figure(2,figsize=pltsize, dpi=150)
sc = m.contourf(lon,lat,s1_windage,levels=np.linspace(s1_windage.min(axis=None),s1_windage.max(axis=None),301),latlon=True)
m.scatter(ntlon,ntlat, latlon=True,color='orange',label="No Windage Tracers")
m.scatter(wtlon,wtlat, latlon=True,color='cyan',label="Windage Tracers")
m.scatter(wtlon[:,0],wtlat[:,0], latlon=True,color='red',label="Initial Position")
plt.legend()
plt.title("windage s$_{1}$ at 2pm")
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(sc, cax=cax)
plt.savefig('Windage_S1_Closeup.png', transparent=False, bbox_inches='tight')
'''
u = u_sea + wc*u_wind
v = v_sea + wc*v_wind
fig = plt.figure(3,figsize=pltsize, dpi=150)
m.quiver(lon[::10,::10],lat[::10,::10],u_wind[::10,::10],v_wind[::10,::10],latlon=True,color='blue',label='Wind Velocity',scale=100)
m.quiver(lon[::10,::10],lat[::10,::10],u[::10,::10],v[::10,::10],latlon=True,color='black',label='Hybrid Velocity',scale=10)
m.quiver(lon[::10,::10],lat[::10,::10],u_sea[::10,::10],v_sea[::10,::10],latlon=True,color='orange',label='Sea Velocity',scale=10)
plt.legend()
plt.title("Velocity fields at 2pm")
m.drawcoastlines()
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.show()
plt.savefig('quiver.png', transparent=False, bbox_inches='tight')
