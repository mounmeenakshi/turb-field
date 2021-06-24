'''
This code computed the A rayleigh distribution for 2-D and its power spectrum.

'''


import numpy as np
from matplotlib import pyplot as plt

from scipy import integrate



                          #**********Function  for Rayleigh distribution ******************************************************

def my_dist(A,k1):
	zeta=1
	sig2=k1**(-zeta)    #sigma**2
	return (A/(2.*np.pi*sig2))*np.exp(-A**2/(2.*sig2))   #((A/(2.*np.pi*sig**2))*np.exp(-A**2/(2.*sig**2)))
	





N = 2000
 

Nx=Ny=500

L=10


y1=np.zeros(Ny)   # y is rows, x is columns

x1=np.zeros(Nx)
arr=np.arange(0,Nx)
print (len(arr))
    #*********The k-space array *************************************
y1[0:int(Ny/2)]=arr[0:int(Ny/2)]/L
y1[int(Ny/2):,]=-(Nx-arr[int(Ny/2):,])/L


x1[0:int(Nx/2)]=arr[0:int(Nx/2)]/L
x1[int(Nx/2):,]=-(Nx-arr[int(Nx/2):,])/L

y1[0]=0.000001
x1[0]=0.00001
#print (y1)

print (y1[0],y1[1])


#print (np.min(y1t),np.max(y1t))
k=np.zeros((Ny*Nx)) # y is rows, x is columns
#l=0
#for i in range(Ny):
#	k[i,:]=np.sqrt(y1[i]**2+x1[:]**2)
	
l=0
for i in range(Ny):
	for j in range(Nx):
		k[l]=np.sqrt(y1[i]**2+x1[j]**2)
		l=l+1

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

A_arr=np.zeros(Ny*Nx)

zeta=1.
 # We assign the A_k for each K-value, based on the sigma for its distribution, we take a range here (of 1 \sigma)
l=0
for i in range(Ny):
	for j in range(Nx):
		sigma=np.sqrt(abs(k[l])**(-zeta))  #sigma-- it determined the max of probability A_k =sigma
		A_ky=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)  # choosing randomly the A_x and A_y components of A(K)
		A_kx=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)
		A_arr[l]=np.sqrt(A_ky**2+A_kx**2)   # taking the magnitude
		l=l+1
		#print ('arr',A_arr[i],sigma)
	
	

print (np.min(A_arr),np.max(A_arr))  


p = my_dist(A_arr,k)  #using A_k and k values to find the probability distribution

print (np.min(p),np.max(p))


# to check sigma 
zipped_lists = zip(A_arr,p)

sorted_lists = sorted(zipped_lists)

p_sorted = [x for p, x in sorted_lists]   # sorting out the Probability in order of A increasing for plotting
A_arr.sort()

plt.xlim(0,1)
#print (y1t[o:5])
plt.plot(A_arr,p_sorted)     #plotting the probability distribution and one can check if this is similar to the trial above

plt.show() 
'''


#k=np.abs(k)
zipped_lists = zip(k,A_arr)  #To check the power-spectrum of A_k vs k   (comment the sorted A_arr)

sorted_lists = sorted(zipped_lists)

arr_sorted = [x for A_arr, x in sorted_lists]
k.sort()
k=np.asarray(k)
arr_sorted=np.asarray(arr_sorted)

#arr_sorted[y1<1]=0.0
plt.plot(np.log10(k),np.log10(arr_sorted**2))

plt.show()
'''

