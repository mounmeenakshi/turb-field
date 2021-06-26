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

def my_dist(A,k1):
	zeta=3.
	sig2=k1**(-zeta)    #sigma**2
	return (A/(sig2))*np.exp(-A**2/(2.*sig2))   #((A/(2.*np.pi*sig**2))*np.exp(-A**2/(2.*sig**2)))
	



def func(x,m,c):
	return m*x+c

N = 2000
 

Nx=Ny=Nz=100

L=10

y1=np.zeros(Ny)  

x1=np.zeros(Nx)
z1=np.zeros(Nz)

arr=np.arange(0,Nx)
print (len(arr))
    #*********The k-space array *************************************
y1[0:int(Ny/2)]=arr[0:int(Ny/2)]/L
y1[int(Ny/2):,]=-(Nx-arr[int(Ny/2):,])/L


x1[0:int(Nx/2)]=arr[0:int(Nx/2)]/L
x1[int(Nx/2):,]=-(Nx-arr[int(Nx/2):,])/L

z1[0:int(Nz/2)]=arr[0:int(Nz/2)]/L
z1[int(Nz/2):,]=-(Nz-arr[int(Nz/2):,])/L

y1[0]=0.000001
x1[0]=0.00001
z1[0]=0.00001
#print (y1)

#print (y1[0],y1[1])


#print (np.min(y1t),np.max(y1t))
k1=np.zeros((Nz*Ny*Nx)) 
#l=0
#for i in range(Ny):
#	k[i,:]=np.sqrt(y1[i]**2+x1[:]**2)
	
l=0
for k in range(Nz):   # y is rows, x is columns in the z=constant plane
	for i in range(Ny):
		for j in range(Nx):
			k1[l]=np.sqrt(z1[k]**2+y1[i]**2+x1[j]**2)
			l=l+1



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

A_arr=np.zeros(Nz*Ny*Nx)
Ax=np.zeros(Nz*Ny*Nx)
Ay=np.zeros(Nz*Ny*Nx)
Az=np.zeros(Nz*Ny*Nx)

zeta=3.
 # We assign the A_k for each K-value, based on the sigma for its distribution, we take a range here (of 1 \sigma)
l=0
for k in range(Nz):
	for i in range(Ny):
		for j in range(Nx):
			if (k1[l])>0.1:
				sigma=np.sqrt(k1[l]**(-zeta))  #sigma-- it determined the max of probability A_k =sigma
				A_ky=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)  # choosing randomly the A_x and A_y components of A(K)
				A_kx=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)
				A_kz=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)
				Ax[l]=A_kx;Ay[l]=A_ky;Az[l]=A_kz;
				A_arr[l]=np.sqrt(A_kz**2+A_ky**2+A_kx**2)   # taking the magnitude
			l=l+1
				#print ('arr',A_arr[i],sigma)


print (np.min(A_arr),np.max(A_arr))  
#Ax[k1<0.1]=0;Ay[k1<0.1]=0;Az[k1<0.1]=0.0;
#A_arr[k<0.1]=0

'''
p = my_dist(A_arr,k1)  #using A_k and k values to find the probability distribution

print (np.min(p),np.max(p))


# to check sigma 
zipped_lists = zip(A_arr,p)

sorted_lists = sorted(zipped_lists)

p_sorted = [x for p, x in sorted_lists]   # sorting out the Probability in order of A increasing for plotting
A_arr.sort()

plt.xlim(0,5)
#print (y1t[o:5])
plt.plot(A_arr,p_sorted)     #plotting the probability distribution and one can check if this is similar to the trial above



plt.title(r"$\mathrm{Ak\, vs\, Prob.\,}$")
plt.show()





zipped_lists = zip(k1,A_arr)  #To check the power-spectrum of A_k vs k   (comment the sorted A_arr)

sorted_lists = sorted(zipped_lists)

arr_sorted = [x for A_arr, x in sorted_lists]
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

'''
#************************ Magnetic field *************************************8


Bx=np.zeros(Nz*Ny*Nx)
By=np.zeros(Nz*Ny*Nx)
Bz=np.zeros(Nz*Ny*Nx)
l=0
for k in range(Nz):   # y is rows, x is columns in the z=constant plane
	for i in range(Ny):
		for j in range(Nx):
			Bx[l]=(y1[i]*Az[l]-z1[k]*Ay[l])
			By[l]=(-x1[j]*Az[l]+z1[k]*Ax[l])
			Bz[l]=(x1[j]*Ay[l]-y1[i]*Ax[l])
			l=l+1
			


B=np.sqrt(Bx**2+By**2+Bz**2)	    #Check divergence of B is zeros, 

#print (np.min(B),np.max(B))
Bx_fft=np.fft.ifft(Bx*1j)
By_fft=np.fft.ifft(By*1j)
Bz_fft=np.fft.ifft(Bz*1j)
mag1=np.sqrt(Bx_fft.imag**2+By_fft.imag**2+Bz_fft.imag**2+Bx_fft.real**2+By_fft.real**2+Bz_fft.real**2)

print (np.min(mag1),np.max(mag1))

mag=np.reshape(mag1,(Nz,Ny,Nx))

#plt.plot(k1,B)
#plt.show()

B[B==0]=0.1
plt.xlim(0.1,1)
#*********** plot Mag field power spectra ********************************************************************
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


'''


