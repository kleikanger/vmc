#call:
#
# python plotSPD <filename> <number of bins (resolution)>

import matplotlib.pyplot as plt
import numpy as np
import linecache as lc
import sys
import os
from pylab import *


os.system("c++ -O3 radialSPD.cpp -lcblas -o radialSPD.out")
def genFromFile(df, rs):
	print df
	os.system("./radialSPD.out %s %i"%(df,rs))
	return np.genfromtxt(fname=(df+".txt"))
	
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

#resolution = 150
#DMC6ptw028Data = genFromFile('spd_2012.5.7.11438r',resolution)
#VMC6ptw028Data = genFromFile('spd_2012.5.7.133836r0d',resolution)
#VMC2ptw1Data= genFromFile('spd_2012.5.7.135746.datr0i',resolution)

#fig = plt.figure()
#l1=plt.plot(DMC6ptw028Data[:,1],DMC6ptw028Data[:,0],'b--')
#l1=plt.plot(VMC6ptw028Data[:,1],VMC6ptw028Data[:,0],'b--')
#l2=plt.plot(b2[1::],v2[1::],'k')
#plt.ylabel(r'$|\Psi(|\mathbf{r}|)|^2$',fontsize=20)
#plt.xlabel(r'$|\mathbf{r}|$',fontsize=20)
#plt.legend((l1,l2),('with Jastrow','without Jastrow'))
#plt.show()

##################################################################
######Two particle results DMC,VMC,VMC NO JAST for w=1 and w=0.01#
##################################################################

#DMC2ptw1 = genFromFile('spd_2012.5.7.194530r0i',resolution)
#DMC2ptw01 = genFromFile('spd_2012.5.7.203223r0i',resolution)
#JAS2ptw1 = genFromFile('spd_2012.5.7.205436r0i',resolution)
#VMC2ptw1 = genFromFile('spd_2012.5.7.195219r0i',resolution)
#VMC2ptw01 = genFromFile('spd_2012.5.7.203545r0i',resolution)
#JAS2ptw01 = genFromFile('spd_2012.5.7.21322r0i',resolution)

#fig = plt.figure(figsize=plt.figaspect(2.))
#ax = fig.add_subplot(2, 1, 1)
#ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
#ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
#l1=plt.plot(DMC2ptw1[:,1],DMC2ptw1[:,0],'k')
#l2=plt.plot(VMC2ptw1[:,1],VMC2ptw1[:,0],'r--')
#l3=plt.plot(JAS2ptw1[:,1],JAS2ptw1[:,0],'b-')
#plt.legend((l1,l2,l3),('DMC','VMC','VMC no Jastrow'),numpoints='2')
#ax.grid()
#ax = fig.add_subplot(2, 1, 2)
#ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
#ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
#l1=plt.plot(DMC2ptw01[:,1],DMC2ptw01[:,0],'k')
#l2=plt.plot(VMC2ptw01[:,1],VMC2ptw01[:,0],'r--')
#l3=plt.plot(JAS2ptw01[:,1],JAS2ptw01[:,0],'b-')
#plt.legend((l1,l2,l3),('DMC','VMC','VMC no Jastrow'), numpoints='2')
#ax.grid()
#plt.show()


def plot(DMC1,VMC1,DMC01,VMC01):
	fig = plt.figure(figsize=plt.figaspect(2.))
	ax = fig.add_subplot(2, 1, 1)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(DMC1[:,1],DMC1[:,0],'k--')
	l2=plt.plot(VMC1[:,1],VMC1[:,0],'r')
	#l3=plt.plot(JAS2ptw1[:,1],JAS6ptw1[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'),numpoints='2')
	#ax.grid()
	ax = fig.add_subplot(2, 1, 2)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(DMC01[:,1],DMC01[:,0],'k--')
	l2=plt.plot(VMC01[:,1],VMC01[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	#ax.grid()
	plt.show()
def plot6(ADMC1,AVMC1,ADMC01,AVMC01,BDMC1,BVMC1,BDMC01,BVMC01,CDMC1,CVMC1,CDMC01,CVMC01):
	fig = plt.figure(figsize=plt.figaspect(1.))
	ax = fig.add_subplot(3, 2, 1)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(ADMC1[:,1],ADMC1[:,0],'k--')
	l2=plt.plot(AVMC1[:,1],AVMC1[:,0],'r')
	#l3=plt.plot(JAS2ptw1[:,1],JAS6ptw1[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'),numpoints='2')
	plt.title('A)', fontsize='large', horizontalalignment='left')
	#ax.grid()
	ax = fig.add_subplot(3, 2, 2)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(ADMC01[:,1],ADMC01[:,0],'k--')
	l2=plt.plot(AVMC01[:,1],AVMC01[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	plt.title('B)', fontsize='large', horizontalalignment='left')
	#ax.grid()
	ax = fig.add_subplot(3, 2, 3)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(BDMC1[:,1],BDMC1[:,0],'k--')
	l2=plt.plot(BVMC1[:,1],BVMC1[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	plt.title('C)')
	#ax.grid()
	ax = fig.add_subplot(3, 2, 4)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(BDMC01[:,1],BDMC01[:,0],'k--')
	l2=plt.plot(BVMC01[:,1],BVMC01[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	plt.title('D)', fontsize='large', horizontalalignment='left')
	#ax.grid()
	ax = fig.add_subplot(3, 2, 5)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(CDMC1[:,1],CDMC1[:,0],'k--')
	l2=plt.plot(CVMC1[:,1],CVMC1[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	plt.title('E)')
	#ax.grid()
	ax = fig.add_subplot(3, 2, 6)
	ax.set_ylabel(r'$|\Psi(\mathbf{r})|^2$',fontsize=20)
	ax.set_xlabel(r'$|\mathbf{r}|$',fontsize=20)
	l1=plt.plot(CDMC01[:,1],CDMC01[:,0],'k--')
	l2=plt.plot(CVMC01[:,1],CVMC01[:,0],'r')
	#l3=plt.plot(JAS6ptw01[:,1],JAS6ptw01[:,0],'b-')
	plt.legend((l1,l2),('DMC','VMC'), numpoints='2')
	plt.title('F)', fontsize='large', horizontalalignment='left')
	#ax.grid()
	plt.show()

######################################################
######Six particle DMC,VMC results for w=1 and w=0.01#
######################################################
resolution = 150
DMC6ptw1 = genFromFile('spd_2012.5.7.224444r0i',resolution)
DMC6ptw01 = genFromFile('spd_2012.5.8.077r0i',resolution)
VMC6ptw1 = genFromFile('spd_2012.5.7.232527r0i',resolution)
VMC6ptw01 = genFromFile('spd_2012.5.7.233933r0i',resolution)
#plot(DMC6ptw1,VMC6ptw1,DMC6ptw01,DMC6ptw01)

#########################################################
######Twelve particle DMC,VMC results for w=1 and w=0.01#
#########################################################

#reading from two different files. Taking average.
#resolution = 150
DMC12ptw1 = genFromFile('spd_2012.5.8.172630r0i',resolution)
DMC12ptw1 += genFromFile('spd_2012.5.8.172630r1i',resolution)
DMC12ptw1 /= 2.

DMC12ptw01 = genFromFile('spd_2012.5.8.17568r0i',resolution)
DMC12ptw01 += genFromFile('spd_2012.5.8.17568r1i',resolution)
DMC12ptw01 /= 2.

VMC12ptw1 = genFromFile('spd_2012.5.8.182628r0i',resolution)
VMC12ptw1 += genFromFile('spd_2012.5.8.182628r1i',resolution)
VMC12ptw1 /= 2.

VMC12ptw01 = genFromFile('spd_2012.5.8.185737r0i',resolution)
VMC12ptw01 += genFromFile('spd_2012.5.8.185737r1i',resolution)
VMC12ptw01 /= 2.

#plot(DMC12ptw1,VMC12ptw1,DMC12ptw01,VMC12ptw01)

plot6(DMC6ptw1,VMC6ptw1,DMC6ptw01,VMC6ptw01,DMC12ptw1,VMC12ptw1,DMC12ptw01,VMC12ptw01\
		,DMC12ptw1,VMC12ptw1,DMC12ptw01,VMC12ptw01)
#plot4(VMC12ptw01,VMC12ptw01,VMC12ptw01,VMC12ptw01,VMC12ptw01,VMC12ptw01,VMC12ptw01,VMC12ptw01,)




