#arguments infile(a0.dat=a) nomber_of_files(filename<number>.dat) outfile(plot)
#runs blocking.cpp and generates plot

import sys
import os
import matplotlib as plt
import numpy as np
from pylab import *

datafile=""
numprocs=0
outfile=""

if (len(sys.argv)==4):
	datafile=sys.argv[1]
	numprocs=np.int(sys.argv[2])
	outfile=sys.argv[3]
if (len(sys.argv)!=4):
	print('error: wrong number of args')

os.system("g++ -O3 blocking.cpp -o blocking.out")
os.system("./blocking.out 100 3000 2 %i %s  %s"%(numprocs,datafile,datafile))

#datafile=datafile[0:len(datafile)-4]
datafile=datafile+".txt" 
print("writing to datafile %s"%datafile)
data = np.genfromtxt(fname=datafile)

fig=plt.figure()
plt.plot(data[:,0],data[:,2],'k+')
plt.xlabel(r'$\tau_{trial}$', size=20)
plt.ylabel(r'$\epsilon$', size=20)
plt.xlim(np.min(data[:,0]),np.max(data[:,0]))
plt.ylim(np.min(data[:,2]),np.max(data[:,2]))
fig.savefig(outfile+'.eps',format='eps')
plt.show()
