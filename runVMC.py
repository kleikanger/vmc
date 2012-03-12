#####################################################
# this program writes h - files with precompiler    #
# definitions, calls Makefile and runs vmc program  #
#####################################################

import sys
import os
import time

###########
#variables#
###########

#vmc variables
omega 				= 1.0
delta_t				=.05
min_alpha 	 		= 0.987
max_alpha 		 	= 0.987
alpha_variations 	= 1 #min 1
min_beta 	 		= 0.398
max_beta 		 	= 0.398
beta_variations 	= 1 #min 1
number_of_particles = 2
sampling_cycles 	= 2e6 #total number on all procs
thermal_cycles 		= 4e5

####################
#running parameters#
####################

#write running parameters to log (then all data will be traceable)
log_run				= False 
#mpirun flags
number_of_processors= 1
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
def gen_sampl_h():
	f_name_B='"%s%s"'%(filepath_blck,time_now())
	f_name_C='"%s%s.dat"'%(filepath_opd,time_now())
	os.system("echo '#define WRITEOFB %s'  > definitions/sampler_Def.h"%write_blocking_data)
	os.system("echo '#define OFPATHB  %s' >> definitions/sampler_Def.h"%f_name_B)
	os.system("echo '#define WRITEOFC %s' >> definitions/sampler_Def.h"%write_opd)
	os.system("echo '#define OFPATHC %s' >> definitions/sampler_Def.h"%f_name_C)
#write variational data to file (include in vmcmain.cpp)
def gen_vmcmain_h():
	f_name_A='"%s%s.dat"'%(filepath_var,time_now())
	os.system("echo '#define WRITEOFA %s'  > definitions/vmcmain_Def.h"%write_var_result)
	os.system("echo '#define OFPATHA  %s' >> definitions/vmcmain_Def.h"%f_name_A)

###############
#preparing run#
###############

cyc = sampling_cycles / number_of_processors
a_inc=abs(max_alpha-min_alpha)/alpha_variations 
a_ini=min_alpha-a_inc
b_inc=abs(max_beta-min_beta)/beta_variations 
b_ini=min_beta-b_inc

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
		number_of_processors)

#rewrite h-files, recompile and run with current parameters
def run():
	if log_run:
		write_to_log()
	gen_sampl_h()
	gen_vmcmain_h()
	running_arguments = (('%s %s %s %s %s %s %s %s %s %s %s')%(all_param()[0:-1]))
	if (os.system('make ')==0): #--silent
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




