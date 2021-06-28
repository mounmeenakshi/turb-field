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
		#sig2=4.
		probab=(A/(sig2))*np.exp(-A**2/(2.*sig2))

	except: probab=0.0
	return probab
	



def func(x,m,c):
	return m*x+c


 

zeta=5.
Nx=128
Ny=128
Nz=128

L=2.

y1=np.zeros(Ny)  

x1=np.zeros(Nx)
z1=np.zeros(Nz)

arrx=np.arange(0,Nx)
arry=np.arange(0,Ny)
arrz=np.arange(0,Nz)

#print (len(arr))
#print (arr)
    #*********The k-space array *************************************


'''
y1[0:int(Ny/2)]=arry[0:int(Ny/2)]/L
y1[int(Ny/2):,]=-(Ny-arry[int(Ny/2):,])/L

print (np.fft.fftfreq(Nx))
print (y1)

x1[0:int(Nx/2)]=arrx[0:int(Nx/2)]/L
x1[int(Nx/2):,]=-(Nx-arrx[int(Nx/2):,])/L

z1[0:int(Nz/2)]=arrz[0:int(Nz/2)]/L
z1[int(Nz/2):,]=-(Nz-arrz[int(Nz/2):,])/L

'''

#k1=np.zeros((Nz*Ny*Nx)) 

x1 = np.fft.fftfreq(Nx,d=0.1)  # L=10 now
y1 = np.fft.fftfreq(Ny,d=0.1) #L =12 here    min 1/(120*0.1)
z1 = np.fft.fftfreq(Nz,d=0.1)  # L=12.8 



ky,kx,kz=np.meshgrid(y1,x1,z1)

#print ('ky shape',ky.shape)


k1=np.sqrt(ky**2+kx**2+kz**2)


#print (k1.shape)


print ('wavenumber k min and max',np.min(k1[k1!=0]),k1.max())



'''
k1=np.ndarray.flatten(k1)

 # To plot the Rayleigh distribution for random A values for the K-space.  To check how Rayleigh distribution looks.
A_k1=np.random.uniform(low=0,high=10,size=2000)
zipped_lists = zip(A_k1,k1) #sorting k baised on A 

sorted_lists = sorted(zipped_lists)

k_sorted = [x for k1, x in sorted_lists]
A_k1.sort()
k_sorted=np.asarray(k_sorted)
print (A_k1[0:2],k_sorted[0:2])
A_k1=np.asarray(A_k1)

plt.plot(A_k1,my_dist(A_k1,abs(k_sorted)))


plt.show()
'''

sigma=np.zeros((Nx,Ny,Nz))


A_arr=np.zeros((Nx,Ny,Nz))

A_kx_r=np.zeros((Nx,Ny,Nz))
A_ky_r=np.zeros((Nx,Ny,Nz))
A_kz_r=np.zeros((Nx,Ny,Nz))
A_kx_i=np.zeros((Nx,Ny,Nz))
A_ky_i=np.zeros((Nx,Ny,Nz))
A_kz_i=np.zeros((Nx,Ny,Nz))

 # We assign the A_k for each K-value, based on the sigma for its distribution, we take a range here (of 1 \sigma)


#k_min=0.1

#print ('k_max,k_min',k_max,k_min)

k_max=max(np.max(abs(x1)),np.max(abs(y1)),np.max(abs(z1)))

k_min=max(np.min(abs(x1[1:])),np.min(abs(y1[1:])),np.min(abs(z1[1:])))

print ('k_min and k_max is',k_min,k_max)
mask=k1>k_min
mask&=k1<k_max

#print (mask)
#mask &=k1<6.
sigma[mask]=np.sqrt(k1[mask]**(-zeta)) 
#sigma[mask]=2.0
delta=0.5  #deviation around the max probability point
phi=2*np.pi*np.random.random((Nx,Ny,Nz))
A_mag=np.zeros((Nx,Ny,Nz))

A_mag[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+delta*sigma[mask])

#print (A_mag)
#print (phi)
A_kx_r[mask]=A_mag[mask]*np.cos(phi[mask])
A_kx_i[mask]=A_mag[mask]*np.sin(phi[mask])


del A_mag

A_mag=np.zeros((Nx,Ny,Nz))
phi=2*np.pi*np.random.random((Nx,Ny,Nz))
A_mag[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+delta*sigma[mask])

A_ky_r[mask]=A_mag[mask]*np.cos(phi[mask])
A_ky_i[mask]=A_mag[mask]*np.sin(phi[mask])


del A_mag

A_mag=np.zeros((Nx,Ny,Nz))
phi=2*np.pi*np.random.random((Nx,Ny,Nz))
A_mag[mask]=np.random.uniform(sigma[mask]-delta*sigma[mask],sigma[mask]+delta*sigma[mask])


A_kz_r[mask]=A_mag[mask]*np.cos(phi[mask])
A_kz_i[mask]=A_mag[mask]*np.sin(phi[mask])

del A_mag
del phi

#print ('max_Aky',np.max(A_ky_i[mask]))

A_arr=np.sqrt(A_kz_i**2+A_kz_r**2+A_ky_i**2+A_ky_r**2+A_kx_i**2+A_kx_r**2) 

#A_arr=np.sqrt(A_kz**2+A_ky**2+A_kx**2) 
#print (A_arr[mask])
#print (np.exp(2*np.pi*0.234*1j))



# saving the vector potential in real space after ifftn
print (np.min(A_arr),np.max(A_arr))  
#Ax[k1<0.1]=0;Ay[k1<0.1]=0;Az[k1<0.1]=0.0;
#A_arr[k<0.1]=0


A_kx=np.fft.ifftn(A_kx_r+A_kx_i*1.j)
A_ky=np.fft.ifftn(A_ky_r+A_ky_i*1.j)
A_kz=np.fft.ifftn(A_kz_r+A_kz_i*1.j)


dbfile=open('ak3d_x.pkl','wb')
pickle.dump(A_kx,dbfile)
dbfile.close()

dbfile=open('ak3d_y.pkl','wb')
pickle.dump(A_ky,dbfile)
dbfile.close()

dbfile=open('ak3d_z.pkl','wb')
pickle.dump(A_kz,dbfile)
dbfile.close()



# to check sigma 



#using A_k and k values to find the probability distribution


'''
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

#Bx=np.zeros((Nx,Ny,Nz))
#By=np.zeros((Nx,Ny,Nz))
#Bz=np.zeros((Nx,Ny,Nz))

#print ([mask])

#print (np.max(A_kz))
Bx= (ky*(A_kz_r*1.j-A_kz_i) - kz*(A_ky_r*1.j-A_ky_i))
By =(- kx*(A_kz_r*1.j-A_kz_i) + kz*(A_kx_r*1.j-A_kx_i))     # i**2= -1.
Bz = (kx*(A_ky_r*1.j-A_ky_i) - ky*(A_kx_r*1.j-A_kx_i))


print ('sum',np.sum(Bx+By+Bz))

print ('Bx_real',np.max(Bx.imag),np.max(Bx.real))
#print (np.min(B),np.max(B))
Bx_fft=np.fft.ifftn(Bx)
By_fft=np.fft.ifftn(By)
Bz_fft=np.fft.ifftn(Bz)

#B=np.zeros(Nx*Ny*Nz)
'''
#########Check if the divergence of magnetic field is zeros in real space, for check in k-space, just multiply Bi by ki

Bx_fft=np.ndarray.flatten(Bx_fft)
By_fft=np.ndarray.flatten(By_fft)
Bz_fft=np.ndarray.flatten(Bz_fft)

Bx_check=(Bx_fft[1:]-Bx_fft[-1:])/5.
By_check=(By_fft[1:]-By_fft[-1:])/5.
Bz_check=(Bz_fft[1:]-Bz_fft[-1:])/5.




Bx_check=(Bx_fft[1:,1:,1:]-Bx_fft[-1:,-1:,-1:])/0.1
By_check=(By_fft[1:,1:,1:]-By_fft[-1:,-1:,-1:])/0.1
Bz_check=(Bz_fft[1:,1:,1:]-Bz_fft[-1:,-1:,-1:])/0.1

print (np.sum(Bx_check+By_check+Bz_check))


'''
B=np.sqrt(Bx.real**2+Bx.imag**2+By.real**2+By.imag**2+Bz.real**2+Bz.imag**2)
#mag1=np.sqrt(Bx_fft.imag**2+By_fft.imag**2+Bz_fft.imag**2+Bx_fft.real**2+By_fft.real**2+Bz_fft.real**2)

#print (np.min(mag1),np.max(mag1))

#mag=np.reshape(mag1,(Nz,Ny,Nx))

#plt.plot(k1,B)
#plt.show()


#plt.xlim(0.1,1)
#********************* plot Mag field power spectra ***************************************************

del A_kz_r,A_kz_i
del A_ky_r,A_ky_i
del A_kx_r,A_kx_i
del sigma


B=np.ndarray.flatten(B)
k1=np.ndarray.flatten(k1)

print (B)
'''
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

'''
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

del k1


del Bx
del By
del Bz
del B

