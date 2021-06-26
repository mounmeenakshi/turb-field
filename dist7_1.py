'''
This code computed the A rayleigh distribution for 3-D real fields and its power spectrum.

'''


import numpy as np
from matplotlib import pyplot as plt
import math
import cmath
from scipy import integrate
from scipy.optimize import curve_fit
import pickle

                          #**********Function  for Rayleigh distribution ******************************************************

def my_dist(A,ks):
	
	try:
		sig2=ks**(-zeta)    #sigma**2
		probab=(A/(sig2))*np.exp(-A**2/(2.*sig2))

	except: probab=0.0
	return probab
	



def func(x,m,c):
	return m*x+c


 

zeta=3.
Nx=Ny=Nz=128

L=10.

y1=np.zeros(Ny)  

x1=np.zeros(Nx)
z1=np.zeros(Nz)

arr=np.arange(0,Nx)
print (len(arr))
#print (arr)
    #*********The k-space array *************************************
y1[0:int(Ny/2)]=arr[0:int(Ny/2)]/L
y1[int(Ny/2):,]=-(Ny-arr[int(Ny/2):,])/L

#print (np.fft.fftfreq(Nx))
print (y1)

x1[0:int(Nx/2)]=arr[0:int(Nx/2)]/L
x1[int(Nx/2):,]=-(Nx-arr[int(Nx/2):,])/L

z1[0:int(Nz/2)]=arr[0:int(Nz/2)]/L
z1[int(Nz/2):,]=-(Nz-arr[int(Nz/2):,])/L



#k1=np.zeros((Nz*Ny*Nx)) 

kx,ky,kz=np.meshgrid(x1,y1,z1)
k1=np.sqrt(ky**2+kx**2+kz**2)
#k1=np.sqrt(z1[np.newaxis,np.newaxis,:]**2+y1[np.newaxis,:,np.newaxis]**2+x1[:,np.newaxis,np.newaxis]**2)

print (k1.shape)

#print (k1[126:200])
print (k1.min(),k1.max())


#print (y1,x1)
#print (k)



'''
 # To plot the Rayleigh distribution for random A values for the K-space.  To check how Rayleigh distribution looks.
A_k1=np.random.uniform(low=0,high=10,size=2000)
zipped_lists = zip(A_k1,k) #sorting k baised on A 

sorted_lists = sorted(zipped_lists)

k_sorted = [x for k, x in sorted_lists]
A_k1.sort()
k_sorted=np.asarray(k_sorted)
print (A_k1[0:2],k_sorted[0:2])
A_k1=np.asarray(A_k1)

#plt.plot(A_k1,my_dist(A_k1,abs(k_sorted)))


#plt.show()

'''
sigma=np.zeros((Nx,Ny,Nz))
A_arr=np.zeros((Nx,Ny,Nz))

A_kx=np.zeros((Nx,Ny,Nz))
A_ky=np.zeros((Nx,Ny,Nz))
A_kz=np.zeros((Nx,Ny,Nz))

 # We assign the A_k for each K-value, based on the sigma for its distribution, we take a range here (of 1 \sigma)


k_min=1
mask=k1>k_min
sigma[mask]=np.sqrt(k1[mask]**(-zeta)) 
#sigma[mask]=2.0
delta=1.  #deviation around the max probability point
A_kx[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+3*sigma[mask])
A_ky[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+3*sigma[mask])
A_kz[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+3*sigma[mask])

A_arr=np.sqrt(A_kz**2+A_ky**2+A_kx**2) 



'''
# saving the vector potential in real space after ifftn
print (np.min(A_arr),np.max(A_arr))  
#Ax[k1<0.1]=0;Ay[k1<0.1]=0;Az[k1<0.1]=0.0;
#A_arr[k<0.1]=0

A_kx=np.fft.ifftn(A_kx)
A_ky=np.fft.ifftn(A_ky)
A_kz=np.fft.ifftn(A_kz)

  
dbfile=open('mag3d_x.pkl','wb')
pickle.dump(A_kx,dbfile)
dbfile.close()

dbfile=open('mag3d_y.pkl','wb')
pickle.dump(A_ky,dbfile)
dbfile.close()

dbfile=open('mag3d_z.pkl','wb')
pickle.dump(A_kz,dbfile)
dbfile.close()


# to check sigma 


#using A_k and k values to find the probability distribution

A_arr=np.ndarray.flatten(A_arr)
k1=np.ndarray.flatten(k1)
p = my_dist(A_arr,k1)

p=np.asarray(p)[A_arr!=0].tolist()
A_arr=A_arr[A_arr!=0].tolist()

print (np.min(p),np.max(p))

zipped_lists = zip(A_arr,p)

sorted_lists = sorted(zipped_lists)


p_sorted = [x for p, x in sorted_lists]   # sorting out the Probability in order of A increasing for plotting
A_arr.sort()

#plt.xlim(0,5)
#print (y1t[o:5])
plt.plot(A_arr,p_sorted)     #plotting the probability distribution and one can check if this is similar to the trial above

print (p_sorted[0:10])

plt.title(r"$\mathrm{Ak\, vs\, Prob.\,}$")
plt.show()

A_arr=np.ndarray.flatten(A_arr)
k1=np.ndarray.flatten(k1)
k1=np.asarray(k1)[A_arr!=0].tolist()
A_arr=A_arr[A_arr!=0].tolist()


zipped_lists = zip(k1,A_arr)  #To check the power-spectrum of A_k vs k   (comment the sorted A_arr)

sorted_lists = sorted(zipped_lists)

arr_sorted = [x for A, x in sorted_lists]
k1.sort()
k1=np.asarray(k1)
arr_sorted=np.asarray(arr_sorted)

#arr_sorted[y1<1]=0.0
popt, pcov = curve_fit(func,np.log10(k1),np.log10(arr_sorted**2))

plt.plot(np.log10(k1),np.log10(arr_sorted**2))
plt.plot(np.log10(k1),func(np.log10(k1),*popt), 'r-',label='fit: y= %5.3f x + %5.3f' % tuple(popt))
plt.legend()

plt.title(r"$\mathrm{k\, vs\, |Ak|^2\, (log-log\,plot)}$")
plt.show()

#************************ Magnetic field *************************************8



B=np.zeros(Nx*Ny*Nz)
Bx=np.zeros((Nx,Ny,Nz))
By=np.zeros((Nx,Ny,Nz))
Bz=np.zeros((Nx,Ny,Nz))
l=0
for k in range(Nz):   # y is rows, x is columns in the z=constant plane
	for i in range(Ny):
		for j in range(Nx):
			Bx[j,i,k]=(ky[j,i,k]*A_kz[j,i,k]-kz[j,i,k]*A_ky[j,i,k])
			By[j,i,k]=-(kx[j,i,k]*A_kz[j,i,k]-kz[j,i,k]*A_kx[j,i,k])
			Bz[j,i,k]=(kx[j,i,k]*A_ky[j,i,k]-ky[j,i,k]*A_kx[j,i,k])
			l=l+1
	



'''


Bx = ky*A_kz - kz*A_ky
By = - kx*A_kz + kz*A_kx
Bz = kx*A_ky - ky*A_kx


#print (np.min(B),np.max(B))
Bx_fft=np.fft.ifftn(Bx*1.j)
By_fft=np.fft.ifftn(By*1.j)
Bz_fft=np.fft.ifftn(Bz*1.j)



B=np.sqrt(Bx**2+By**2+Bz**2)
#mag1=np.sqrt(Bx_fft.imag**2+By_fft.imag**2+Bz_fft.imag**2+Bx_fft.real**2+By_fft.real**2+Bz_fft.real**2)

#print (np.min(mag1),np.max(mag1))

#mag=np.reshape(mag1,(Nz,Ny,Nx))

#plt.plot(k1,B)
#plt.show()


#plt.xlim(0.1,1)
#********************* plot Mag field power spectra ***************************************************


B=np.ndarray.flatten(B)
k1=np.ndarray.flatten(k1)


k1=np.asarray(k1)[B!=0].tolist()
B=B[B!=0].tolist()

zipped_lists = zip(k1,B)  #To check the power-spectrum of A_k vs k   (comment the sorted A_arr)

sorted_lists = sorted(zipped_lists)

b_sorted = [x for B, x in sorted_lists]
k1.sort()
k1=np.asarray(k1)
b_sorted=np.asarray(b_sorted)
#plt.loglog(k1,b_sorted)
#arr_sorted[y1<1]=0.0
plt.plot(np.log10(k1),np.log10(b_sorted**2))
popt, pcov = curve_fit(func,np.log10(k1),np.log10(b_sorted**2))

plt.plot(np.log10(k1),np.log10(b_sorted**2))
plt.plot(np.log10(k1),func(np.log10(k1),*popt), 'r-',label='fit: y= %5.3f x + %5.3f' % tuple(popt))
plt.legend()

plt.title(r"$\mathrm{k\, vs\, |Bk|^2\, (log-log\,plot)}$")

plt.show()



#Bx_fft=np.fft.ifftn(Bx*1.j)
#By_fft=np.fft.ifftn(By*1.j)
#Bz_fft=np.fft.ifftn(Bz*1.j)


dbfile=open('mag3d_x.pkl','wb')
pickle.dump(Bx_fft,dbfile)
dbfile.close()

dbfile=open('mag3d_y.pkl','wb')
pickle.dump(By_fft,dbfile)
dbfile.close()

dbfile=open('mag3d_z.pkl','wb')
pickle.dump(Bz_fft,dbfile)
dbfile.close()

dbfile=open('mag3d.pkl','wb')
pickle.dump(B,dbfile)
dbfile.close()



del A_kz
del A_ky
del A_kx
del sigma
del k1

del Bx
del By
del Bz
del B


