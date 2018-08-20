#import h5py as hp
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import Basemap
import time
import calendar
plt.close('all')
#tstart = calendar.timegm(time.strptime('Jun 1, 2017 @ 00:00:00 UTC', '%b %d, %Y @ %H:%M:%S UTC'))
timestep=0
lat_max = 41.4548797564
lat_min = 41.0250244184
lon_max = -70.2909317079
lon_min = -70.916465753
n=[48,49,50,51,52]
initial = np.genfromtxt('manikan_initial_pos.txt', delimiter=',')
final = np.genfromtxt('manikan_final_pos.txt', delimiter=',')

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
fig = plt.figure(1,figsize=pltsize, dpi=150)
ax=plt.subplot(111)
m.scatter(initial[:,1],initial[:,0],latlon=True,color='y',label="10:40 position")
m.scatter(final[:,1],final[:,0],latlon=True,color='g',label="2:30 position")
for i in range(initial.shape[0]):
    m.plot([initial[i,1],final[i,1]],[initial[i,0],final[i,0]],latlon=True)
plt.legend()
m.drawcoastlines()
parallels = np.arange(round(lat_min,1),lat_max+0.1,0.1)
meridians = np.arange(round(lon_max,1),lon_min-0.1,-0.1)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.savefig('WindageLCS.png', transparent=False, bbox_inches='tight')


