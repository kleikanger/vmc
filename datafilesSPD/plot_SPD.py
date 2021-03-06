#call:
#
# python plotSPD <filename> <number of bins (resolution)>

import matplotlib.pyplot as plt
import numpy as np
import linecache as lc
import sys
from pylab import *


#dataf = sys.argv[2]
n_bins = int(sys.argv[1])
	
datafileVMC = 'spd_2012.5.3.181449.dat'
datafileNoJas = 'spd_2012.5.3.18432.dat'
#datafileDMC = 'spd_2012.5.6.163232.dat'

#datafileVMC = 'spd_2012.5.6.17513.dat'
#datafileDMC = 'spd_2012.5.6.18357.dat'

def sortData(datafile):
	#import sample points
	data=np.genfromtxt(fname=datafile, skip_footer=1)
	#make array with all radial distances
	data=data**2
	r=sqrt(np.add(data[:,0],data[:,1]))
	#sort the data in n bins
	if (np.size(data[0,:])==3): #if a third row exists: dmc-data: weight all values
		n, bins = np.histogram(r, n_bins, weights=data[:,2])#, normed='true')
	else:
		n, bins = np.histogram(r, n_bins)#, normed='true')
	#the mean value of the bins
	radial_bins=np.add(bins[0:n_bins],bins[1:n_bins+1])/2
	#the radial density (density as a func of r, density central symmetrical)
	radial_values=n/(radial_bins*2*np.pi)
	#normalize
	radial_values*=n_bins/(np.sum(radial_values)*np.max(radial_bins))
	del data
	#plot
	return radial_bins, radial_values

b1,v1 = sortData(datafileVMC)
b2,v2 = sortData(datafileNoJas)
#b3,v3 = sortData(datafileDMC)

#fig = plt.figure(figsize=plt.figaspect(2.))
#ax = fig.add_subplot(2, 1, 1)
#ax.set_ylabel('E',size=20)
#ax.set_xlabel('iterations',size=20)
#l = ax.plot(arange(0,len(dat)),dat[:,1],'r')
#l = ax.plot(arange(0,len(dat)),dat[:,0],'k--')
#ax.grid()
#ax = fig.add_subplot(2, 1, 2)
#ax.set_ylabel('number of walkers',size=20)      
#ax.set_xlabel('iterations',size=20)
#l = ax.plot(arange(0,len(dat)),dat[:,2],'r')
#ax.grid()
#plt.show()

fig = plt.figure()
l1=plt.plot(b1[1::],v1[1::],'b--')
#l2=plt.plot(b2[1::],v2[1::],'k')
l2=plt.plot(b2[1::],v2[1::],'k')
plt.ylabel(r'$|\Psi(|\mathbf{r}|)|^2$',fontsize=20)
plt.xlabel(r'$|\mathbf{r}|$',fontsize=20)
plt.legend((l1,l2),('with Jastrow','without Jastrow'))
plt.show()





