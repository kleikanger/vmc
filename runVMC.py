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
omega 				= 1.0
delta_t				= .005
min_alpha 	 		= 0.93 # 0.98 	#init value cgm-method and DMC
max_alpha 		 	= 0.9
alpha_variations 	= 1 	#min 1
min_beta 	 		= 0.56 # 0.4 	#init value cgm-method and DMC
max_beta 		 	= 0.9
beta_variations 	= 1 	#min 1
number_of_particles = 6
sampling_cycles 	= 1e6 	#total number on all procs

thermal_cycles 		= 4e5 	#also used in initialization of DMC 

#dmc variables : note : delta_t:.01 |2:.98,.4,3.0004|6:.93,.56|12:87,68,dt=0.005-0.001 (0.001 converging to slowly?)
number_of_walkers 	= 4000 #total number on all procs
num_cycles_main_loop= 100
num_c_ET_upd_loop 	= 400 #O(100)-O(1000)
num_c_equilibri_loop= 2000
initial_e_trial 	= 3.190004

####################
#running parameters#
####################

#mpirun flags
number_of_processors= 2
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
#use conjugate gradient method minimization
use_cgm_minimization = ''
if conjugate_gradient:
	use_cgm_minimization = 'true'
else: 
	use_cgm_minimization = 'false'

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
def gen_sampl_h():
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

cyc = sampling_cycles / number_of_processors
a_inc=abs(max_alpha-min_alpha)/alpha_variations 
a_ini=min_alpha-a_inc
b_inc=abs(max_beta-min_beta)/beta_variations 
b_ini=min_beta-b_inc
conjugate_gradient 	= int(conjugate_gradient) 
sample_on_grid 		= int(sample_on_grid)
use_dmc_sampler 	= int(use_dmc_sampler)
number_of_walkers/=number_of_processors;

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
		a_ini,a_inc,alpha_variations,\
		b_ini,b_inc,beta_variations,\
		number_of_particles,\
		cyc,thermal_cycles,\
		conjugate_gradient,\
		sample_on_grid,\
		initial_e_trial,
		num_cycles_main_loop,
		num_c_ET_upd_loop,
		num_c_equilibri_loop,
		number_of_walkers,
		use_dmc_sampler)

#rewrite h-files, recompile and run with current parameters
def run():
	if log_run:
		write_to_log()
	gen_sampl_h()
	gen_vmcmain_h()
	running_arguments = (('%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s')%(all_param()))
	if (os.system('make --silent CC=mpic++ DEBUG=')==0): #--silent
		os.system("echo 'Compilation successful.'")
		os.system('mpirun -n %i runVMC.out %s'%(number_of_processors,running_arguments))
	else:
		exit(1)
		#os.system('make --debug')

#############
#run program#
#############

#do simulations : change values between runs
#delta_t = 0.01
run()
#delta_t = 0.025
#run()
#delta_t = 0.05
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