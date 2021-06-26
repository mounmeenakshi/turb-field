import numpy as np
import pickle
import matplotlib.pyplot as plt

dbfile = open('mag3d_z.pkl','rb') 
mag = pickle.load(dbfile)

dbfile.close()


mag1=np.reshape(mag.real,(100,100,100))
print (np.min(mag),np.max(mag))
mag_r=mag1[:,70,:]

Nx=100
arr1=np.arange(-50,50)
arr2=np.arange(-50,50)
arr3=np.arange(-50,50)


fig=plt.figure()

cmap = plt.get_cmap('jet')

c=plt.pcolormesh(arr1,arr2,np.log10(mag_r), cmap=cmap)#,vmin=-3,vmax=-1)#,vmin=-4.177)# vmin=np.min(den_profile),vmax=np.max(den_profile)
cbar=plt.colorbar(c)
cbar.ax.set_ylabel(r'$\log(\mathrm{n[cm^{-3}]})$',size=14)

#plt.title(r'Density plot (R-$\phi(\theta\approx 90)$)', fontweight ="bold")
plt.tight_layout()
#plt.grid()
#fig.savefig('/mnt/home/student/cmeenakshi/public/folders/images_236/temp_t_90_shock_100.png',dpi=300)
plt.show()
