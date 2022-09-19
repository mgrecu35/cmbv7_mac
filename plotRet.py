from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

fname1="/Users/mgrecu/cmbv7/Data/2B-CS-CONUS.GPM.DPRGMI.CORRA2022.20190108-S083514-E084324.027624.V07A.HDF5"
fname2="outDir/2B.GPM.DPRGMI.20190108.027624.HDF5"

fh1=Dataset(fname1)
sfcRain=fh1["KuKaGMI/nearSurfPrecipTotRate"][:]
lon=fh1["KuKaGMI/Longitude"][:]
lat=fh1["KuKaGMI/Latitude"][:]
binNodes=fh1["KuKaGMI/phaseBinNodes"][:]

fh2=Dataset(fname2)
sfcRain1=fh2["KuKaGMI/nearSurfPrecipTotRate"][:]
lon1=fh2["KuKaGMI/Longitude"][:]
lat1=fh2["KuKaGMI/Latitude"][:]
binNodes2=fh2["KuKaGMI/phaseBinNodes"][:]

a=np.nonzero(sfcRain>0)
b=np.nonzero(binNodes[:,:,2][a]>=binNodes[:,:,4][a])
plt.figure(figsize=(5,5))
ax=plt.subplot(111)
plt.scatter(sfcRain[a][b],sfcRain1[a][b])
plt.xlabel('V7 rate [mm/h]')
plt.ylabel('New rate [mm/h]')
plt.title('Near surface snow')
plt.plot(range(0,6),range(0,6))
plt.plot(np.arange(0,6),np.arange(0,6)*1.5)
plt.xlim(0,5)
plt.ylim(0,5)
plt.legend(['Data','1:1 line','1:1.5 line'])
ax.set_aspect('equal')
plt.savefig('nearSurfPRate.scatter.027624.png')

import cartopy
import cartopy.crs as ccrs
fig=plt.figure(figsize=(6,6))
ax=fig.add_subplot(211,projection=ccrs.PlateCarree())
plt.pcolormesh(lon[300:599,:],lat[300:599,:],sfcRain[300:599,:],\
               cmap='jet',norm=matplotlib.colors.LogNorm(vmin=0.1,vmax=10))
ax.coastlines()
ax.gridlines()
plt.title('V7')
#plt.colorbar()
ax=fig.add_subplot(212,projection=ccrs.PlateCarree())
plt.pcolormesh(lon[300:599,:],lat[300:599,:],sfcRain1[300:599,:],\
               cmap='jet',norm=matplotlib.colors.LogNorm(vmin=0.1,vmax=10))
plt.title('New')
ax.coastlines()
ax.gridlines()
plt.tight_layout()
axc = fig.add_axes([0.8, 0.3, 0.025, 0.4])
cbar=plt.colorbar(orientation='vertical',cax=axc)
cbar.ax.set_title('mm/h')
plt.savefig('nearSurfPRate.map.027624.png')
