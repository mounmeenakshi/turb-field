'''
This code computed the A rayleigh distribution for 1-D and its power spectrum.

'''


import numpy as np
from matplotlib import pyplot as plt

from scipy import integrate



                          #**********Function  for Rayleigh distribution ******************************************************

def my_dist(A,y12):
	zeta=1
	sig2=y12**(-zeta)    #sigma**2
	return (A/(2.*np.pi*sig2))*np.exp(-A**2/(2.*sig2))   #((A/(2.*np.pi*sig**2))*np.exp(-A**2/(2.*sig**2)))
	





N = 2000
 

Ny=2000

L=10


y1=np.zeros(Ny)
arr=np.arange(0,Ny)
print (len(arr))
    #*********The k-space array *************************************
y1[0:int(Ny/2)]=arr[0:int(Ny/2)]/L
y1[int(Ny/2):,]=-(Ny-arr[int(Ny/2):,])/L


y1[0]=0.00001

#print (y1)

print (y1[0],y1[1])


#print (np.min(y1t),np.max(y1t))







 # To plot the Rayleigh distribution for random A values for the K-space.  To check how Rayleigh distribution looks.
A_k1=np.random.uniform(low=0,high=10,size=2000)
zipped_lists = zip(A_k1,y1) #sorting k baised on A 

sorted_lists = sorted(zipped_lists)

y1_sorted = [x for y1, x in sorted_lists]
A_k1.sort()
y1_sorted=np.asarray(y1_sorted)
print (A_k1[0:2],y1_sorted[0:2])
A_k1=np.asarray(A_k1)

#plt.plot(A_k1,my_dist(A_k1,abs(y1_sorted)))


A_arr=np.zeros(len(y1))

 # We assign the A_k for each K-value, based on the sigma for its distribution, we take a range here (of 1 \sigma)
for i in range(len(y1)):
	zeta=1.
	sigma=np.sqrt(abs(y1[i])**(-zeta))  #sigma-- it determined the max of probability A_k =sigma
	A_k2=np.random.uniform(sigma-0.5*sigma,sigma+0.5*sigma)
	A_arr[i]=A_k2
	#print ('arr',A_arr[i],sigma)
	

	

print (np.min(A_arr),np.max(A_arr))  
p = my_dist(A_arr,abs(y1))  #using A_k and k values to find the probability distribution

print (np.min(p),np.max(p))


# to check sigma 
zipped_lists = zip(A_arr,p)

sorted_lists = sorted(zipped_lists)

p_sorted = [x for p, x in sorted_lists]   # sorting out the Probability in order of A increasing for plotting
#A_arr.sort()

#plt.xlim(0,1)
#print (y1t[o:5])
#plt.plot(A_arr,p_sorted)     #plotting the probability distribution and one can check if this is similar to the trial above

##plt.show() 



y1=np.abs(y1)
zipped_lists = zip(y1,A_arr)  #To check the power-spectrum of A_k vs k   (comment the sorted A_arr)

sorted_lists = sorted(zipped_lists)

arr_sorted = [x for A_arr, x in sorted_lists]
y1.sort()
y1=np.asarray(y1)
arr_sorted=np.asarray(arr_sorted)

#arr_sorted[y1<1]=0.0
plt.plot(np.log10(y1),np.log10(arr_sorted**2))

plt.show()


