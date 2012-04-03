#arguments infileA+ infileB+ <-p> ...
# 
#<-p> plot results
#
#ex: python blocking.py a+ b+ : plotting all files starting with a or b
#                             : automatically finding number of procs
#NB: filename must be on the form XXXXr0.dat, XXXXr1.dat, ...
#
#runs blocking.cpp and generates plot (.eps) with the 
#same name as the infile

import sys
import os
import matplotlib as plt
import numpy as np
from pylab import *

plot_res=False

def blca(datafile,numprocs): 
	os.system("g++ -O3 blocking.cpp -o blocking.out")
	print('%s%s'%(numprocs,datafile))
	os.system("./blocking.out 100 3000 2 %i %s"%(numprocs,datafile))

	#datafile=datafile[0:len(datafile)-4]
	print('')
	print("writing to datafile %s"%datafile+'.txt')

	data = np.genfromtxt(fname=datafile+'.txt')

	fig=plt.figure()
	plt.plot(data[:,0],data[:,2],'k+')
	plt.xlabel(r'$\tau_{trial}$', size=20)
	plt.ylabel(r'$\epsilon$', size=20)
	plt.xlim(np.min(data[:,0]),np.max(data[:,0]))
	plt.ylim(np.min(data[:,2]),np.max(data[:,2]))
	fig.savefig(datafile+'.eps',format='eps')
	if plot_res:
		os.system('evince %s%s '%(datafile+'.eps','&'))
	print("plot saved : %s"%(datafile+'.eps'))

#if sys.argv[1][-1]=='+':
def blc_plus(ifname):
	print('processing all datafiles %s*'%ifname[1][:-1])
	s=os.listdir('.')
	#find names starting with <sys.argv[1][:-1]>
	for nam in s:
		if (len(nam)>len(ifname)+3):
			if (nam[:len(ifname)-1]==ifname[:-1]) & (nam[-5]=='0') & (nam[-4:]=='.dat'):
				#find numprocs, searching through all files in directory
				nprc=0
				for nal in s:
					if ((nam[:-5]+nam[-4:])==(nal[:-5]+nal[-4:])):
						if (int(nal[-5])>nprc):
							nprc=int(nal[-5])
				#run
				blca(nam[:-5],nprc+1)

#loop through all imput arguments
for imp in sys.argv[1::]:
	if imp[0:2]=='-p':
		plot_res=True
for imp in sys.argv[1::]:
	if imp[-1]=='+':
		print(imp)
		blc_plus(imp)


