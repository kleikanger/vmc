#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "sga.h"
#include "../newmatrix/newmatrix.h"

//ref: Jurgen A. Doornik 2005
#include "../ziggurat/zigrandom.h"
#include "../ziggurat/zignor.h"
//Defining random number generators RAN_NORM,RAN_NORM_SET,RAN_UNI,RAN_UNI_SET
//File automatically generated in python script and looks approx like this:
//
// #define RAN_NORM DRanNormalZig32
// #define RAN_NORM_SET RanNormalSetSeedZig32
// #define RAN_UNI DRan_MWC8222
// #define RAN_UNI_SET RanSetSeed_MWC8222
#include "../definitions/randomNumberGenerators.h"

using std::cout;
using std::cerr;
//using std::ofile;
using std::ofstream;
using std::setprecision;
using std::setiosflags;
using std::setw;
using std::ios;
using std::ostringstream;

#include "../definitions/sampler_Def.h"

sga::sga(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank, int nprocs, MPI_Status status)
{/*//startvimfold*/
	this->myrank=myrank;
	this->nprocs=nprocs;
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;
	this->status=status;
}/*//endvimfold*/

sga::~sga()
{/*//startvimfold*/
	delete [] quantum_dot;
}//End function /*//endvimfold*/

//remove e_trial, num_c_dmc_equilibriation_loop
void sga::SGAMin(
		int num_cycles, 
		int thermalization, 
		double* var_par, 
		double delta_t, 
		int num_c_dmc_main_loop,
		int	num_c_dmc_inner_loop,
		int initial_number_of_walkers
		)
{/*//startvimfold*/
	
	//init random number generators
	int idum_d=time(NULL)*(myrank+1);
	RAN_UNI_SET(&idum_d,5);	
	
	int corr_length = 3000; //number of steps between init of particles
	int sga_out_upd = 100; //Start collecting data after 100 cycles (therm..)
	
	int i, n_sga=0;
	int sga_upd=initial_number_of_walkers*num_c_dmc_inner_loop;
	double e_local_temp, e_mean, len_energy_gradient;
	double **e_grad_cum = (double**)matrix(2,2,sizeof(double));
	double energy_gradient[2], energy_gradient_old[2];
	double var_par_cum[2], var_par_cum2[2], m_sga_inc[2];
	double e_local_sga=0., len_energy_gradient_cum[2];
	double m_v_sga[2], e_grad_temp[2];;  

	//array of walkers	
	quantum_dot = new walker*[initial_number_of_walkers]; 
	
	e_grad_cum[0][0] = e_grad_cum[1][0] 
	= e_grad_cum[0][1] = e_grad_cum[1][1] = 0.0;
	len_energy_gradient_cum[0]=len_energy_gradient_cum[1]=0.;
	m_sga_inc[0]=m_sga_inc[1]=1.0;
	m_v_sga[0]=m_v_sga[1]=1.0;

	//Open file stream
	ofstream ofilecga;
	ofilecga.open("cgadata/test.txt");

	//write to screen
	if (myrank==0)
	cout<<"\nDMC simulation - \n"
		<<"initial_number_of_walkers/procs: "<<initial_number_of_walkers<<"\n"
		<<"nprocs:                          "<<nprocs<<"\n"
		<<"cycles main loop:                "<<num_c_dmc_main_loop<<"\n"
		<<"cycles inner loop:               "<<num_c_dmc_inner_loop<<"\n"
		<<"thermalization cycles:           "<<thermalization<<"\n"
		<<"delta_t:                         "<<delta_t<<"\n"
		<<"omega:                           "<<var_par[2]<<"\n"
		<<"number of particles:             "<<num_part<<"\n\n";
	
	//********************************
	// initializing walker objects
	//********************************
	//init one walker at a random position
	quantum_dot[0] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
	quantum_dot[0]->initWalker(var_par, delta_t);
	//init empty walkers
	for (i=1;i<initial_number_of_walkers;i++)
	{
		quantum_dot[i] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
		quantum_dot[i]->initEmptyWalker(var_par, delta_t);
	}
	//Initialize walkers, using VMC routine
	initWalkers(initial_number_of_walkers, thermalization, corr_length, var_par);

	//**********************************
	// minimize trial function
	//**********************************
	for (int loop_main=0; loop_main<num_c_dmc_main_loop; loop_main++) //XXX TODO While some convergence criteria || loop<max_iter
	{
		//set some variables to 0		
		e_grad_cum[0][0]
			=e_grad_cum[0][1]
			=e_grad_cum[1][0]
			=e_grad_cum[1][1]=0.0;
		e_mean=0.0;
		//loop over walkers
		for (int loop_p=0;loop_p<initial_number_of_walkers;loop_p++)
		{
			//energy and gradient accumulating cycles 
			for (int loop_c=0;loop_c<num_c_dmc_inner_loop;loop_c++)
			{
				//move all particles, one at a time
				for (i=0; i<num_part; i++)
					//metropolis-hastings test
					if (quantum_dot[loop_p]->tryRandomStep(i)) 
						quantum_dot[loop_p]->acceptStep(i);
					else  
						quantum_dot[loop_p]->rejectStep(i);
				//find local energy
				e_local_temp = quantum_dot[loop_p]->calcLocalEnergy(var_par);
				//collecting gradient
				quantum_dot[loop_p]->getVarParGrad(e_grad_temp);
				e_grad_cum[0][0]+=e_grad_temp[0];
				e_grad_cum[0][1]+=e_grad_temp[1];
				e_grad_cum[1][0]+=e_grad_temp[0]*e_local_temp;
				e_grad_cum[1][1]+=e_grad_temp[1]*e_local_temp;
				//collecting energy
				e_mean+=e_local_temp;
			}
		}

		//****collect and broadcast. MPI only optimized for 1 computer, too much communication!!
		MPI_Allreduce(MPI_IN_PLACE, e_grad_cum[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, e_grad_cum[1], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &e_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		e_grad_cum[0][0]/=nprocs;
		e_grad_cum[0][1]/=nprocs;
		e_grad_cum[1][0]/=nprocs;
		e_grad_cum[1][1]/=nprocs;
		e_mean/=nprocs;
		MPI_Barrier(MPI_COMM_WORLD);//XXX necc? (feels safe)

		//energy minimization
		e_mean /= (double)sga_upd;
		energy_gradient[0] = 2*(e_grad_cum[1][0]-e_grad_cum[0][0]*e_mean)/(double)sga_upd;
		energy_gradient[1] = 2*(e_grad_cum[1][1]-e_grad_cum[0][1]*e_mean)/(double)sga_upd;
			
		if (n_sga>30&&n_sga<70)
		{
			//len_energy_gradient_cum+=len_energy_gradient;
			len_energy_gradient_cum[1]+=fabs(energy_gradient[1]);
			len_energy_gradient_cum[0]+=fabs(energy_gradient[0]);
		}
		else if (n_sga==70)
		{
			//XXX .8, 300, 20 are tunable parameters XXX
			m_v_sga[0] = pow(300.*(len_energy_gradient_cum[0]/40.),.8); 
			m_v_sga[1] = pow(300.*(len_energy_gradient_cum[1]/40.),.8);
			m_sga_inc[0]  = m_v_sga[0]/20.;
			m_sga_inc[1]  = m_v_sga[1]/20.;
		}
		//setting max length to move to .01. Vector still has the correct direction. 
		len_energy_gradient = sqrt(cblas_ddot(2,energy_gradient,1,energy_gradient,1));	
		if ( (len_energy_gradient*pow((double)m_v_sga[0],-.8)>.01)
		|| (len_energy_gradient*pow((double)m_v_sga[1],-.8)>.01) )
		{
			energy_gradient[0]*=sqrt(0.01)/len_energy_gradient;
			energy_gradient[1]*=sqrt(0.01)/len_energy_gradient;
		}
		
		//can be proven to converge mathematically
		for (i=0;i<2;i++)
		{
			//m++ if the gradient changed its sign
			if (energy_gradient[i]/energy_gradient_old[i]<0) m_v_sga[i]+=m_sga_inc[i];
			//store old gradient
			energy_gradient_old[i]=energy_gradient[i];
			//new value of the variational parameters
			var_par[i]-=energy_gradient_old[i]*pow((double)m_v_sga[i],-.8);
			//prevent unphysical negative values of alpha,beta
			if (var_par[i]<0.01) var_par[i]=0.01;
		}
		n_sga++;
		//reset variational parameter 
		for (i=0;i<initial_number_of_walkers;i++)
			quantum_dot[i]->wSetVarPar(var_par);

		if (n_sga>sga_out_upd) //make better test
		{
			//collect data
			var_par_cum[0]+=var_par[0];
			var_par_cum[1]+=var_par[1];
			var_par_cum2[0]+=var_par[0]*var_par[0];
			var_par_cum2[1]+=var_par[1]*var_par[1];
			
			//calculate error+std.dev
			if (myrank==0)
			if (n_sga%100==0)
			{
				cout<<"\r"<<"                                                                                                ";
				cout<<"\r"<<myrank
					<<setprecision(8)//<<setw(10)
					<<"b_C: "<<var_par_cum[0]/(double)(n_sga-sga_out_upd)
					<<" a_C: "<<var_par_cum[1]/(double)(n_sga-sga_out_upd)
					<<" b: "<<var_par[0]
					<<" a: "<<var_par[1]
					<<" b_V: "<<sqrt(-pow(var_par_cum[0]/(double)(n_sga-sga_out_upd),2)+var_par_cum2[0]/(double)(n_sga-sga_out_upd))
					<<" a_V: "<<sqrt(-pow(var_par_cum[1]/(double)(n_sga-sga_out_upd),2)+var_par_cum2[1]/(double)(n_sga-sga_out_upd))
					<<" c: "<<n_sga*sga_upd*nprocs
					<<" m_B: "<<m_v_sga[0]
					<<" m_A: "<<m_v_sga[1]
					<<"---------";
					fflush(stdout);
			}
		//Write to file	
		if (myrank==0)
			{
				ofilecga << setiosflags(ios::showpoint | ios::uppercase);
				for (i=0;i<dimension;i++)
				{
					ofilecga << setw(16) << setprecision(16) 
						<< n_sga*sga_upd*nprocs << " " 							//loops
						<< var_par[0] << " " 								//alpha
						<< var_par_cum[0]/(double)(n_sga-sga_out_upd) <<" " //alpha mean
						<< var_par[1] << " " 								//beta
						<< var_par_cum[1]/(double)(n_sga-sga_out_upd) <<" ";//beta mean
				}
				ofilecga<<"\n";
			}
		}
	}

	//Collect results
	MPI_Allreduce(MPI_IN_PLACE, var_par_cum, 2, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
	//Write results to screen
	if (myrank==0)
	{	
				cout<<"\n"
					<<"                                                         "
					<<"                                                         "
					<<"                                                         ";
				cout<<"\rFinal result:" 
					<<"alpha="<<var_par_cum[1]/((double)(n_sga-sga_out_upd)*nprocs)
					<<", beta="<<var_par_cum[0]/((double)(n_sga-sga_out_upd)*nprocs)
					<<"\n\n";
	}
}/*//endvimfold*/
/*
	Initializing <initial_number_of_walkers> walkers. <corr_length> number of MC-steps between the walkers. 
   */
void sga::initWalkers(int initial_number_of_walkers, int thermalization, int corr_length, double *var_par)
{	/*//startvimfold*/
	//factory class
	popCtr = new popControl(num_part, spin_up_cutoff, dimension, status);
	int i, loop_c=0, num_init=1; 
	if (myrank==0)
		{
			cout<<"\rStarting initialization of walkers... ";
			fflush(stdout);
		}
	while (num_init<initial_number_of_walkers)
	{
		loop_c++;
		//move all particles, one at a time
		for (i=0; i<num_part; i++)
		{
			//metropolis-hastings test
			if ( quantum_dot[0]->tryRandomStep(i) ) 
				quantum_dot[0]->acceptStep(i);
			else  
				quantum_dot[0]->rejectStep(i);
		}//All particles moved
		//if thermalization finished, copy walker every corr_length step
		if ((loop_c>thermalization) && (loop_c%corr_length==0)) 
		{ 	
			popCtr->cloneWalker(quantum_dot[0],quantum_dot[num_init]);
			num_init++;	
			//print to screen
			if (myrank==0&&num_init%((initial_number_of_walkers)/200+1)==0)
			{
				cout<<"\r"<<"                                            ";
				cout<<"\rinitializing walkers: "<<setprecision(3)
					<<(num_init)/(double)initial_number_of_walkers*100<<"%";
				fflush(stdout);
			}
		}
	}
	//INITIALIZING quantum_dot[0]
	for (loop_c=0;loop_c<corr_length;loop_c++)
		for (int active_part=0; active_part<num_part; active_part++)
		{
			//metropolis-hastings test
			if ( quantum_dot[0]->tryRandomStep(active_part) ) 
				quantum_dot[0]->acceptStep(active_part);
			else  
				quantum_dot[0]->rejectStep(active_part);
		}//All particles moved
	
	//print to screen
	if (myrank==0)
	{
		cout<<"\r"<<"                                            ";
		cout<<"\rinitialization of walkers finished.\n";
		fflush(stdout);
	}
}/*//endvimfold*/

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
