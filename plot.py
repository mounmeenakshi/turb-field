import numpy as np
import pickle
import matplotlib.pyplot as plt

dbfile = open('mag3d_z.pkl','rb') 
mag = pickle.load(dbfile)

dbfile.close()

print (mag.shape)
N=128
#mag=np.reshape(mag.real,(N,N,N))

#print (np.min(mag),np.max(mag))
mag_r=mag[20,:,:].real
#mag_r=mag_r.T

arr1=np.arange(0,N)
arr2=np.arange(0,N)
arr3=np.arange(0,N)


fig=plt.figure()

cmap = plt.get_cmap('jet')

c=plt.pcolormesh(arr1,arr2,np.log10(mag_r), cmap=cmap)#,vmin=-6,vmax=-2)
cbar=plt.colorbar(c)
#cbar.ax.set_ylabel(r'$\log(\mathrm{n[cm^{-3}]})$',size=14)

#plt.title(r'Density plot (R-$\phi(\theta\approx 90)$)', fontweight ="bold")
plt.tight_layout()
#plt.grid()
#fig.savefig('/mnt/home/student/cmeenakshi/public/folders/images_236/temp_t_90_shock_100.png',dpi=300)
plt.show()