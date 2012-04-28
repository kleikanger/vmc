import numpy as np
import matplotlib as plt
import os
import sys
from pylab import *

def plot(data):
	fig = plt.figure(figsize=plt.figaspect(2.))
	ax = fig.add_subplot(2, 1, 1)
	ax.set_ylabel(r'$\alpha$',size=20)
	ax.set_xlabel('$N_C$',size=20)
	l = ax.plot(data[:,0],data[:,3],'r')
	l = ax.plot(data[:,0],data[:,4],'k--')
	ax = fig.add_subplot(2, 1, 2)
	ax.set_ylabel(r'$\beta$',size=20)      
	ax.set_xlabel('$N_C$',size=20)
	l = ax.plot(data[:,0],data[:,1],'r')
	l = ax.plot(data[:,0],data[:,2],'k--')
	plt.show()

data = np.genfromtxt(fname='test.txt')
plot(data)
