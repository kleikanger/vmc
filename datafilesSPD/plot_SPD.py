#call:
#
# python plotSPD <filename> <number of bins (resolution)>

import matplotlib.pyplot as plt
import numpy as np
import linecache as lc
import sys
from pylab import *


datafile = sys.argv[1]
n_bins = int(sys.argv[2])

#import sample points
data=np.genfromtxt(fname=datafile, skip_footer=1)
#make array with all radial distances
data=data**2
r=sqrt(np.add(data[:,0],data[:,1]))
#sort the data in n bins
n, bins = np.histogram(r, n_bins)#, normed='true')
#the mean value of the bins
radial_bins=np.add(bins[0:n_bins],bins[1:n_bins+1])/2
#the radial density (density as a func of r, density central symmetrical)
radial_values=n/(radial_bins*2*np.pi)
#normalize
radial_values*=n_bins/(np.sum(radial_values)*np.max(radial_bins))
#plot
plt.plot(radial_bins[1::],radial_values[1::],'k')
plt.ylabel(r'$|\Psi(|\mathbf{r}|)|^2$',fontsize=20)
plt.xlabel(r'$|\mathbf{r}|$',fontsize=20)
plt.show()
