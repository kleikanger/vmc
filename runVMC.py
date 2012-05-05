#####################################################
# this program writes h - files with precompiler    #
# definitions, calls Makefile and runs vmc program  #
#####################################################

import sys
import os
import time

#################
#sampling method#
#################

#####################
#optimization method#
#####################

#choose only one should be true for a fast code
conjugate_gradient 	= False
sample_on_grid 		= False
use_dmc_sampler 	= True

###########
#variables#
###########

#vmc variables
#(for cgm and , min_alpha, min_beta is the starting point)
omega 				= .28
delta_t				= .005
min_alpha 	 		= .899#8718 # 0.93 # 0.98 #init value cgm-method and DMC
max_alpha 		 	= 0.5
alpha_variations 	= 1 	#min 1
min_beta 	 		= .414#.6677 #0.56 # 0.4 	#init value cgm-method and DMC
max_beta 		 	= 0.9
beta_variations 	= 1 	#min 1
number_of_particles = 2
sampling_cycles 	= 1e7 	#total number on all procs

thermal_cycles 		= 5e5 	#also used in initialization of DMC 

#dmc variables : note : delta_t:.01 |2:.98,.4,3.0004|6:.93,.56|12:87,68,dt=0.005-0.001 (0.001 converging to slowly?) |omg .5,2pt,a95-96,b.32 E1.66978 
number_of_walkers 	= 1000 #total number on all procs
num_cycles_main_loop= 200
num_c_ET_upd_loop 	= 200 #O(100)-O(1000)
num_c_equilibri_loop= 3000
initial_e_trial 	= 3.0004

####################
#running parameters#
####################

#mpirun flags
number_of_processors= 1
#write running parameters to log (then all data will be traceable)
log_run				= False
#running mode #NOT ACTIVE, find out how to change CC in Makefile
debug 				= False
profile 			= False

#############
#definitions#
#############

#write blockingdata to file blocking/E_<running parameters>.dat
write_blocking_data = 'false'
filepath_blck 		= 'blocking/blc_'
#write variational data to file
write_var_result 	= 'false'
filepath_var 		= 'datafilesVAR/var_'
#write one particle density to file 
write_opd 			= 'false'
filepath_opd 		= 'datafilesSPD/spd_'

#set random number generator
RAN_NORM 			= 'DRanNormalZig32'
RAN_NORM_SET 		= 'RanNormalSetSeedZig32'
RAN_UNI 			= 'DRan_MWC8222'
RAN_UNI_SET			= 'RanSetSeed_MWC8222'

####################################
#writing definitions to headerfiles#
####################################i

#generate string with system time
def time_now():
	return ('%s.%s.%s.%s%s%s'%(time.localtime()[0:6]))
#random number generators (include in walker.cpp)
os.system("echo '#define RAN_NORM %s'  > definitions/randomNumberGenerators.h"%RAN_NORM)
os.system("echo '#define RAN_NORM_SET %s' >> definitions/randomNumberGenerators.h"%RAN_NORM_SET)
os.system("echo '#define RAN_UNI %s' >> definitions/randomNumberGenerators.h"%RAN_UNI)
os.system("echo '#define RAN_UNI_SET %s' >> definitions/randomNumberGenerators.h"%RAN_UNI_SET)
#write blocking data and oneparticle density (include in sampler.cpp)
use_cgm_minimization = ''
def gen_sampl_h():
	#use conjugate gradient method minimization
	if conjugate_gradient:
		use_cgm_minimization = 'false'
	else: 
		use_cgm_minimization = 'false'
	
	f_name_B='"%s%s"'%(filepath_blck,time_now())
	f_name_C='"%s%s.dat"'%(filepath_opd,time_now())
	os.system("echo '#define WRITEOFB %s'  > definitions/sampler_Def.h"%write_blocking_data)
	os.system("echo '#define OFPATHB  %s' >> definitions/sampler_Def.h"%f_name_B)
	os.system("echo '#define WRITEOFC %s' >> definitions/sampler_Def.h"%write_opd)
	os.system("echo '#define OFPATHC %s' >> definitions/sampler_Def.h"%f_name_C)
	os.system("echo '#define CONJGRAD %s' >> definitions/sampler_Def.h"%use_cgm_minimization)
#write variational data to file (include in vmcmain.cpp)
def gen_vmcmain_h():
	f_name_A='"%s%s.dat"'%(filepath_var,time_now())
	os.system("echo '#define WRITEOFA %s'  > definitions/mcongrid_Def.h"%write_var_result)
	os.system("echo '#define OFPATHA  %s' >> definitions/mcongrid_Def.h"%f_name_A)

###############
#preparing run#
###############
#cyc = sampling_cycles / number_of_processors
#a_inc=abs(max_alpha-min_alpha)/alpha_variations 
#a_ini=min_alpha-a_inc
#b_inc=abs(max_beta-min_beta)/beta_variations 
#b_ini=min_beta-b_inc

def write_to_log():
	os.system(("echo\
	'time:%s, num_particles:%s, a_min:%s, a_max:%s,\
 a_variations:%s, b_min:%s, b_max:%s, b_variations:%s,\
 omega:%s, delta_t:%s, sampl_cyc:%s, therm_cyc:%s'\
	>> logfile_runs.txt")%(\
	time_now(),number_of_particles,\
	min_alpha,max_alpha,alpha_variations,\
	min_beta,max_beta,beta_variations,\
	omega,delta_t,sampling_cycles,\
	thermal_cycles))

def all_param():
	return (omega,delta_t,\
		min_alpha-abs(max_alpha-min_alpha)/alpha_variations,#a_ini
		abs(max_alpha-min_alpha)/alpha_variations,alpha_variations,#a_inc
		min_beta-abs(max_beta-min_beta)/beta_variations,#b_ini
		abs(max_beta-min_beta)/beta_variations,#b_inc
		beta_variations,\
		number_of_particles,\
		sampling_cycles / number_of_processors,#cyc
		thermal_cycles,\
		int(conjugate_gradient),\
		int(sample_on_grid),\
		initial_e_trial,
		num_cycles_main_loop,
		num_c_ET_upd_loop,
		num_c_equilibri_loop,
		number_of_walkers/number_of_processors,
		int(use_dmc_sampler))

#rewrite h-files, recompile and run with current parameters
def run():
	if log_run:
		write_to_log()
	gen_sampl_h()
	gen_vmcmain_h()
	running_arguments = (('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s')%(all_param()))
#if (os.system('make profile --silent CC=mpic++ DEBUG=')==0): #--silent
	if (os.system('make --silent CC=mpic++ DEBUG=')==0): #--silent
		os.system("echo 'Compilation successful.'")
		os.system('mpirun -n %i runVMC.out %s'%(number_of_processors,running_arguments))
	else:
		exit(1)
		#os.system('make --debug')

#############
#run program#
#############

#generate one particle densities
conjugate_gradient 	= True
sample_on_grid 		= False
use_dmc_sampler 	= False
number_of_particles = 2
num_cycles_main_loop= 80000
number_of_walkers 	= 30 #total number on all procs
num_c_ET_upd_loop 	= 1 #O(100)-O(1000)

conjugate_gradient 	= True
sample_on_grid 		= True
use_dmc_sampler 	= False
write_opd 			= 'true'
number_of_particles = 2
delta_t 			= 0.05
sampling_cycles     = 1e7
thermal_cycles 		= 10000
min_alpha 			= .988
min_beta 			= .399
omega 				= .01
#run()
write_opd 			= 'false'

conjugate_gradient 	= True
sample_on_grid 		= True
use_dmc_sampler 	= False
number_of_processors= 2
import numpy as np
omega 				= 1
number_of_particles = 6
run()
for omm in np.linspace(0.01,0.075,3):
	omega 				= omm
#run()
for omm in np.linspace(0.1,1,10):
	omega 				= omm
#run()

number_of_particles = 6
omega 				= 1
#run()

conjugate_gradient 	= False
sample_on_grid 		= False
use_dmc_sampler 	= True
number_of_processors= 2
number_of_walkers 	= 1000 #total number on all procs
num_cycles_main_loop= 200
num_c_ET_upd_loop 	= 200 #O(100)-O(1000)
num_c_equilibri_loop= 3000
initial_e_trial 	= 3.0004
min_alpha 			= .924
min_beta 			= .556
omega 				= 1
number_of_particles = 6
omega 				= 1
#reproduces results ref to in LE's master
delta_t=0.005
#run()
delta_t=0.001
#run()

omega 				= .5
delta_t 			= .05
number_of_particles = 12
thermal_cycles 	 	= 4e5


number_of_walkers 	= 5000 #total number on all procs
num_cycles_main_loop= 4000
num_c_ET_upd_loop 	= 200 #O(100)-O(1000)
num_c_equilibri_loop= 3000


number_of_processors= 2
number_of_particles = 2
conjugate_gradient 	= False
sample_on_grid 		= True
use_dmc_sampler 	= False
sampling_cycles 	= 1e8
min_alpha 		 	= .988
min_beta 	 		= .399
#run()

#MINIMIZATION
omega 				= .1
conjugate_gradient 	= True
sample_on_grid 		= False
use_dmc_sampler 	= False
number_of_particles = 2
num_cycles_main_loop= 500000 #30000
number_of_walkers 	= 30 #total number on all procs
num_c_ET_upd_loop 	= 1 #O(100)-O(1000)
number_of_processors= 2
#run()
#run()
#run()
#run()

omega = .1
number_of_processors= 2
delta_t=0.05
number_of_particles = 2
#run()
#run()
#run()
#run()
number_of_particles = 6
#run()
#run()
#run()
#run()
number_of_particles = 12
#run()
#run()
#run()
thermal_cycles 		= 10000
number_of_particles = 20
#run()
#run()


conjugate_gradient 	= False
sample_on_grid 		= True
use_dmc_sampler 	= False
number_of_particles = 2
sampling_cycles     = 1e7
thermal_cycles 		= 0#10000
omega 				= 1

#import numpy as np
#for a in np.linspace(0.001,.4,40):
#	delta_t=a
#run()


number_of_processors= 1
delta_t = 0.05
number_of_particles = 2
min_beta 			= 0.399
min_alpha   		= 0.988
#run()
#os.system("gprof ./runVMC.out")
min_beta 			= 0.557
min_alpha   		= 0.924
number_of_particles = 6
#run()
#os.system("gprof ./runVMC.out")
min_beta 			= 0.548
min_alpha   		= 0.877
number_of_particles = 12
#run()
#os.system("gprof ./runVMC.out")
min_beta 			= 0.68
min_alpha   		= 0.87
number_of_particles = 20
#run()
#os.system("gprof ./runVMC.out")


omega 				= .1
conjugate_gradient 	= False
sample_on_grid 		= False
use_dmc_sampler 	= True
number_of_particles = 2
number_of_processors= 2

number_of_walkers 	= 1000 #total number on all procs
num_cycles_main_loop= 3000
num_c_ET_upd_loop 	= 100 #O(100)-O(1000)
num_c_equilibri_loop= 5000#3000
initial_e_trial 	= 3.0004
number_of_particles = 2
min_beta 			= 0.177
min_alpha   		= 0.949
#delta_t = 0.007
run()
#delta_t = 0.005
#run()
delta_t = 0.003
#run()
delta_t = 0.001
#run()

min_beta 			= 0.219
min_alpha   		= 0.810
number_of_particles = 6
delta_t = 0.007
#run()
write_blocking_data = 'true'
delta_t = 0.005
#run()
write_blocking_data = 'false'
delta_t = 0.003
#run()
delta_t = 0.001
#run()
number_of_particles = 12
min_beta 			= 0.253
min_alpha   		= 0.727
delta_t = 0.007
#run()
write_blocking_data = 'true'
delta_t = 0.005
#run()
write_blocking_data = 'false'
delta_t = 0.003
#run()
delta_t = 0.001
#run()
number_of_particles = 20
min_beta 			= 0.289
min_alpha   		= 0.650
delta_t = 0.007
#run()
#write_blocking_data = 'true'
delta_t = 0.005
#run()
#write_blocking_data = 'false'
delta_t = 0.003
#run()
delta_t = 0.001
#run()

conjugate_gradient 	= False
sample_on_grid 		= True
use_dmc_sampler 	= False
delta_t = 0.05
#run()


#CGM:
#
#
# PLOT CURVES OF E FOR ALL VALUES
#
#
#Test: Do a number of thermalization cycles with a high Delta t eg:0.01	
#
#Do runs with eg:delta_t = .0005,.001,.0015,.002,.004,.006,.008,.01
#an do an extrapolation of the energy to delta_t = 0
#Get mean, variance, energies (6pt.) (+ blocking data for some runs?)
#
#Do testruns with both fixed node approx:kill and deny. check variance and <e>
#Do testruns without the fixed node approx
#
#Plot single particle densities. Plot/find the covariance of the particles <r1-r2> and the tail of the func for the variational system
#and the DMC system. 
#
#Compare the running time of VMC and DMC
#Compare the kinetical, potential, exchange and correlation energy of VMC and DMC
#
#Compare results obtained by using the mixed estimator and the generational generator
#
#


#
# MINIMIZATION RUNS OBTAINING THE ENERGY MINIMA USED IN PROJECT 1
#

delta_t             = .05
conjugate_gradient  = True
sample_on_grid      = False
number_of_particles = 2
sampling_cycles = 1e7
min_alpha = 0.97
min_beta  = 0.25 
omega=.28
#run()
min_alpha = 0.989
min_beta  = 0.30 
omega=.5
#run()
min_alpha = 0.99
min_beta  = 0.40 
omega=1
#run()

delta_t             = .05
conjugate_gradient  = True
sample_on_grid      = False
number_of_particles = 6
sampling_cycles = 1e7
min_alpha = 0.88 
min_beta  = 0.33
omega=.28
#run()
min_alpha = 0.90 
min_beta  = 0.42
omega=.5
#run()
min_alpha = 0.87 
min_beta  = 0.68
omega=1
#run()

delta_t             = .05
conjugate_gradient  = True
sample_on_grid      = False
number_of_particles = 12
sampling_cycles = 1e7
min_alpha = 0.81 
min_beta  = 0.38 
omega=.28
#run()
min_alpha = 0.90 
min_beta  = 0.42
omega=.5
#run()
min_alpha = 0.88 
min_beta  = 0.66
omega=1
#run()

###
#
# Production runs. 
#
###
### Two particles, (write blocking data to file)
conjugate_gradient  = False
sample_on_grid      = True
number_of_particles = 2
sampling_cycles = 1e8
###
omega               = .28
min_alpha = 0.971
min_beta  = 0.252 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = .5
min_alpha = 0.981
min_beta  = 0.309 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = 1.
min_alpha = 0.988
min_beta  = 0.399 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
### Six particles. (write blocking data to file)
conjugate_gradient  = False
sample_on_grid      = True
number_of_particles = 6
sampling_cycles = 1e8
##
omega               = .28
min_alpha = 0.873 
min_beta  = 0.326 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = .5
min_alpha = 0.900
min_beta  = 0.413 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = 1.
min_alpha = 0.924
min_beta  = 0.557 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
### Twelve particles. (write blocking data to file)
conjugate_gradient  = False
sample_on_grid      = True
number_of_particles = 12
sampling_cycles = 1e8
##
omega               = .28
min_alpha = 0.809 
min_beta  = 0.378 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = .5
min_alpha = 0.845
min_beta  = 0.482 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()
###
omega               = 1.
min_alpha = 0.877
min_beta  = 0.658 #my minima 1e7 CGM - samples, gtol 1e-8
delta_t             = .01
#run()
delta_t             = .025
#run()
delta_t             = .05
#run()

