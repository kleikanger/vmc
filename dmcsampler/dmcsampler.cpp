#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include "dmcsampler.h"
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
//using namespace std;
//name and path of ofile
//blocking data
//#ifndef WRITEOFB
//#define WRITEOFB true
//#endif
//single particle density
//#ifndef WRITEOFC
//#define WRITEOFC false
//#endif
//#ifndef OFPATHB
//#define OFPATHB "/home/karleik/masterProgging/vmc/blocking/E_"//.dat
//#endif
//#ifndef OFPATHC
//#define OFPATHC "/home/karleik/masterProgging/vmc/datafiles/xXXspd_2partdt05eopt.dat"
//#endif
//#define CONJGRAD true

dmcsampler::dmcsampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank, int nprocs)
{/*//startvimfold*/
	this->myrank=myrank;
	this->nprocs=nprocs;
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;
	//factory class for manipulating walkers
	popCtr = new popControl(num_part, spin_up_cutoff, dimension);
}/*//endvimfold*/

dmcsampler::~dmcsampler()
{/*//startvimfold*/
	delete [] quantum_dot;
	delete popCtr;
}//End function /*//endvimfold*/

void dmcsampler::sampleDMC(
		int num_cycles, 
		int thermalization, 
		double* var_par, 
		double delta_t, 
		double e_trial, 
		int num_c_dmc_main_loop,
		int	num_c_dmc_inner_loop,
		int	num_c_dmc_equilibriation_loop,
		int initial_number_of_walkers
		)
{/*//startvimfold*/
	

#if WRITEOFB/*//startvimfold*/
	//Only for rank 0 proc
	//if (rank=0) ?
	//if (num_cycles>1.1e7)
    //ofile.open(OFPATHB);
	//double* all_energies = new double[num_cycles + 1];
	//cout<<writing to file ((std::string)OFPATHB)((std::string)OFNAMEB);
#endif
#if WRITEOFC
	ofstream ofilec;
	ofilec.open((OFPATHC));
#endif/*//endvimfold*/

	int idum_d=time(NULL)*(myrank+1);
	RAN_UNI_SET(&idum_d,5);

	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;

	int i,j,k,active_part;
	int accepted=0;
	
	//array of walkers	
	// init 6x number_of_walkers. when 2x number_of_walkers, system renormalize to number_of_walkers. 
	quantum_dot = new walker*[initial_number_of_walkers*6]; 
	bool* occupancy = new bool[initial_number_of_walkers*6];
	double* e_local_old = new double[initial_number_of_walkers*6];

	if (myrank==0)
	cout<<"\nDMC simulation - \n"
		<<" alpha:                          "<<var_par[1]<<", beta: "<<var_par[0]<<"\n"
		<<"initial trial energy:            "<<e_trial<<"\n"
		<<"initial_number_of_walkers/procs: "<<initial_number_of_walkers<<"\n"
		<<"nprocs:                          "<<nprocs<<"\n"
		<<"equilibriation loops:            "<<num_c_dmc_equilibriation_loop<<"\n"
		<<"cycles main loop:                "<<num_c_dmc_main_loop<<"\n"
		<<"cycles inner loop:               "<<num_c_dmc_inner_loop<<"\n"
		<<"delta_t:                         "<<delta_t<<"\n"
		<<"number of particles:             "<<num_part<<"\n";
	
	//WITH MPI: SPLIT INTO nprocs initialization precesses. 
	//initializing walker objects
	quantum_dot[0] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
	quantum_dot[0]->initWalker(var_par, delta_t);
	for (i=1;i<initial_number_of_walkers*6;i++)
	{
		quantum_dot[i] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
		quantum_dot[i]->initEmptyWalker(var_par, delta_t);
		//quantum_dot[i]->initWalker(var_par, delta_t);
	}
	//setting occupancy to 'true' for first initial_number_of_walkers.
	for (i=0;i<initial_number_of_walkers;i++)
		occupancy[i] = true;
	for (i=initial_number_of_walkers;i<initial_number_of_walkers*6;i++)
		occupancy[i] = false;
	
	//TODO input variable corr_length
	int loop_c=0, num_init=1, corr_length = 3000; //number of steps between init of particles

	cout<<"RANK "<<myrank<<" process: starting initialization of walkers...\n";

	//******** Initialize walkers ***********
	/*//startvimfold*/

	//initializing quantum_dot[1] to quantum_dot[num_init-1]
	while (num_init<initial_number_of_walkers)
	{
		loop_c++;
		//move all particles, one at a time
		for (int active_part=0; active_part<num_part; active_part++)
		{
			//metropolis-hastings test
			if ( quantum_dot[0]->tryRandomStep(active_part) ) 
				quantum_dot[0]->acceptStep(active_part);
			else  
				quantum_dot[0]->rejectStep(active_part);
		}//All particles moved
		//if thermalization finished, copy walker every corr_length step
		if ((loop_c>thermalization) && (loop_c%corr_length==0)) 
		{ 	
			popCtr->cloneWalker(quantum_dot[0],quantum_dot[num_init]);
			e_local_old[num_init] = quantum_dot[0]->calcLocalEnergy(var_par);
			num_init++;	
			//obs walker 0 and num_init-1 in same position
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
	e_local_old[0] = quantum_dot[0]->calcLocalEnergy(var_par);
	
	/*//endvimfold*/
	//************************** END OF initialization phase *****************

	cout<<"RANK "<<myrank<<" process: initialization of walkers successful. starting equilibriation of walkers...\n";

	int killsd=0; //
	int num_alive = initial_number_of_walkers;
	double branching_factor;
	int num_new_walkers;
	int num_resurrected = 0;
	double e_cumulative=0.0;
	double e2_cumulative=0.0;
	long int loop_c_cumulative;
	long int total_loop_c_cumulative=0;
	double e_ref=e_trial;
	int loop_main;

	//************************** START equilibriation phase ******************

	/*//startvimfold*/

	//TODO while not equilibriated :: find some criteria
	for (loop_main=0; loop_main<num_c_dmc_equilibriation_loop; loop_main++)
	{
		e_local=0.0;
		num_resurrected=0;
		killsd=0;
		loop_c_cumulative=0;;
		//loop over walkers
		for (int loop_p=0;loop_p<num_alive;loop_p++)
		{
			//energy accumulating cycles 
			for (loop_c=0;loop_c<1;loop_c++) //num_c_dmc_inner_loop/nproc*myrank to ...
			{
				//move all particles, one at a time
				for (int active_part=0; active_part<num_part; active_part++)
				{
					//metropolis-hastings test
					bool mh_test = (quantum_dot[loop_p]->tryRandomStep(active_part)); 
					//reject walkers that has crossed node or that fails the mh_test
					if (mh_test && !quantum_dot[loop_p]->nodeCrossed())
						quantum_dot[loop_p]->acceptStep(active_part);
					else  
						quantum_dot[loop_p]->rejectStep(active_part);
				}//All particles moved

				//find local energy
				e_local_temp = quantum_dot[loop_p]->calcLocalEnergy(var_par);
				//calculate new branching factor
				branching_factor = exp(-delta_t*(0.5*(e_local_temp+e_local_old[loop_p])-e_ref));
				//collect energy weighted by the branching factor
				e_local+=branching_factor*e_local_temp;
				//kill and resurrect walkers
				num_new_walkers=(int)(branching_factor+RAN_UNI());
				if (num_new_walkers==0)
				{
					occupancy[loop_p]=false;
					killsd++;
					loop_c++;//to accumulate the correct nr of cycles!
					break;
				}
				else if (num_new_walkers>=2)
				{
					j=num_alive+num_resurrected;
					while (num_new_walkers>1)
					{
						popCtr->cloneWalker(quantum_dot[loop_p],quantum_dot[j]);
						e_local_old[j]=e_local_temp; //TODO Maybe better be part of walker class. not so much copying, only pointer swapping
						occupancy[j]=true;
						num_resurrected++;
						j++;
						num_new_walkers--;	//LIMIT NUMBER OF NEW PARTICLES TO AVOID NUMERICAL INSTABILITIES. EG: MAX 2 NEW PARTICLES ?	
					}
				}
				//update energy if the walker survived
				e_local_old[loop_p] = e_local_temp;
			}
		//keeping track of the total number of samples
		loop_c_cumulative+=loop_c;
		}
		//sorting walkers. alive first
		int sorted=0; i=0; 
		j=num_alive+num_resurrected-1;//starting on the last one
		while (sorted<killsd)
		{
			walker* swich;
			if (occupancy[i]==false)
			{
				//check that j i occupied
				while (!occupancy[j]) j--;
				if (i>=j) break;//sorting finished

				//pointer swapping
				swich=*(quantum_dot+i);
				*(quantum_dot+i)=*(quantum_dot+j);
				*(quantum_dot+j)=swich;

				occupancy[i]=true;
				occupancy[j]=false;
				e_local_old[i]=e_local_old[j];
				sorted++;
				j--;
			}		
			i++;
			if (i>=j) break;
		}
		num_alive=num_alive-killsd+num_resurrected;
		
		//New reference energy
		e_cumulative+=e_local;
		total_loop_c_cumulative+=loop_c_cumulative;
		//0.01 test other values
		double e_c_temp; long int t_l_c_temp;
		MPI_Allreduce(&e_cumulative, &e_c_temp, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
		MPI_Allreduce(&total_loop_c_cumulative, &t_l_c_temp, 1, MPI_LONG, MPI_SUM,  MPI_COMM_WORLD);
		e_ref = e_c_temp/(double)t_l_c_temp-log((double)num_alive
				/ (double)initial_number_of_walkers)/delta_t*0.01;
		
		//XXX 
		//Not so many send and recieves (5000 here). Seems to work well.
		//e_cumulative+=e_local/(double)loop_c_cumulative;	
		//e_ref = e_cumulative/(double)(loop_main+1)-.01*log((double)num_alive
		// 		/ (double)initial_number_of_walkers)/delta_t;
		
		//renormalize (kill some of the walkers if they get to many) TODO Test algo!
		if (num_alive>initial_number_of_walkers*2)
		{
			int number_to_remove = num_alive-initial_number_of_walkers;
			for (int k=0;k<number_to_remove;k++, --num_alive)
			{
				i=(int)(RAN_UNI()*num_alive+1); //obs ran_uni in (0,1) i in [0,num_alive-1] XXX XXX
				if (i!=(num_alive))
				{
					//removing particle i, and moving the last particle to its position.
					walker* swich;
					swich=*(quantum_dot+i);
					*(quantum_dot+i)=*(quantum_dot+j); //num_alive+num_resurrected-sorted
					*(quantum_dot+num_alive)=swich;
					e_local_old[i]=e_local_old[num_alive];
					occupancy[num_alive]=false;
				}
				else
				{
					occupancy[i]=false;
				}
			}
		}
		//RENORMALIZE IF TO FEW ?
		/*
		if (myrank==0)
		{
			cout<<"\nRANK 0: \nalive/resurrected/killsd<"<<num_alive<<"/"<<num_resurrected<<"/"<<killsd<<"\n";
			cout<<"accumulated energy/number of loops<"
				<<setprecision(16)<<e_cumulative/(double)(loop_main+1)<<"\n";
			cout<<"cycles in main loop: "<<loop_main+1<<" of "<<num_c_dmc_main_loop<<"\n";
			cout<<"\n e_ref<"<<e_ref<<"\n";
		}*/
	}	
	
	MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_loop_c_cumulative, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	e_cumulative/=(double)total_loop_c_cumulative;

	//using different routine for equilibriation with fewer MPI_Calls
	//MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//e_cumulative/=(double)nprocs*loop_main;
	
	cout<<"RANK "<<myrank<<" process: equilibriation of walkers successful.\n";
	if (myrank==0)
	{
		cout<<"current mean E<"<<setprecision(16)<<e_cumulative<<"\n";
		cout<<"e_ref<"<<e_ref<<"\n";
	}
	/*//endvimfold*/
	
	//************************** END OF equilibriation phase *****************
	//************************** START sampling phase ************************
	
	//reset variables/*//startvimfold*/
	//e_trial = initial guess ?
	e_trial=e_cumulative;
	total_loop_c_cumulative=0;
	e_cumulative=0.0;
	for (loop_main=0; loop_main<num_c_dmc_main_loop; loop_main++) //XXX TODO While some convergence criteria || loop<max_iter
	{
		e_local=0.0;
		num_resurrected=0;
		killsd=0;
		loop_c_cumulative=0;
		//loop over walkers
		for (int loop_p=0;loop_p<num_alive;loop_p++)
		{
			//energy accumulating cycles 
			for (loop_c=0;loop_c<num_c_dmc_inner_loop;loop_c++) //num_c_dmc_inner_loop/nproc*myrank to ...
			{
				//move all particles, one at a time
				for (int active_part=0; active_part<num_part; active_part++)
				{/*//startvimfold*/
					//metropolis-hastings test
					bool mh_test = (quantum_dot[loop_p]->tryRandomStep(active_part)); 
					//reject walkers that has crossed node or that fails the mh_test
					if (mh_test && !quantum_dot[loop_p]->nodeCrossed())
					//if (quantum_dot[loop_p]->tryRandomStep(active_part))
					{
					/*
				    //The other method: killing walkers that has crossed a node.	   
					//kill walkers that has crossed node (error linear in delta_t)
					if (quantum_dot[loop_p]->nodeCrossed())
						{
							occupancy[loop_p]=false;
							killsd++;
							break;
							//not adding to loop_c. loop_c has value from previous loop
							//which is correct since the last local_energy was collected then.
						}*/	
					quantum_dot[loop_p]->acceptStep(active_part);
					//if (loop_c>thermalization) { accepted++; } 
					}
					else  
					{
						quantum_dot[loop_p]->rejectStep(active_part);
					}/*//endvimfold*/
				}//All particles moved

				//find local energy
				e_local_temp = quantum_dot[loop_p]->calcLocalEnergy(var_par);
				//calculate new branching factor
				branching_factor = exp(-delta_t*(0.5*(e_local_temp+e_local_old[loop_p])-e_trial)); //e_trial is the (vmc) calculated mean
				
				//collect energy weighted by the branching factor
				e_local+=branching_factor*e_local_temp; //TODO use intermediate variable to calculate b_f*e_l_t
				e2_cumulative+=pow(branching_factor*e_local_temp,2);
				//kill and resurrect walkers
				num_new_walkers=(int)(branching_factor+RAN_UNI());
				if (num_new_walkers==0)
				{
					occupancy[loop_p]=false;
					killsd++;
					loop_c++;//to accumulate the correct nr of cycles! not adding to loop_c if breaking
					break;
				}
				else if (num_new_walkers>=2)//only one test
				{
					j=num_alive+num_resurrected;
					while (num_new_walkers>1)
					{
						popCtr->cloneWalker(quantum_dot[loop_p],quantum_dot[j]);
						e_local_old[j]=e_local_temp; //TODO Maybe better be part of walker class. not so much copying, only pointer swapping
						occupancy[j]=true;
						num_resurrected++;
						j++;
						num_new_walkers--;	//LIMIT NUMBER OF NEW PARTICLES TO AVOID NUMERICAL INSTABILITIES. EG: MAX 2 NEW PARTICLES ?	
					}
				}
				//update energy if the walker survived
				e_local_old[loop_p] = e_local_temp;
			}
		//keeping track of the total number of samples
		loop_c_cumulative+=loop_c;
		}
			
		//collecting energy
		e_cumulative+=e_local;
		//e2_cumulative+=e_local*e_local; //TODO MAYBE collect variance in e_trial
		//total number of moves must be collected for final statistics
		total_loop_c_cumulative+=loop_c_cumulative;
		
		//****collect and broadcast
		//e_trial = e_local/(double)loop_c_cumulative; //GIVING -nan first loop after all walkers dead
		//MPI_Allreduce(MPI_IN_PLACE, &e_trial, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//e_trial/=(double)nprocs;
		
		//****collect and broadcast Find new e_trial
		MPI_Allreduce(MPI_IN_PLACE, &e_local, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &loop_c_cumulative, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		e_trial=e_local/(double)loop_c_cumulative;

		//for (i=0;i<200;i++)
		//	cout<<i<<" "<<occupancy[i]<<" \n";
		//make more intellegent sorting algo
		int sorted=0; i=0; 
		j=num_alive+num_resurrected-1;//starting on the last one
		while (sorted<killsd)
		{
			walker* swich;
			if (occupancy[i]==false)
			{
				//check that j i occupied
				while (!occupancy[j]) j--;
				if (i>=j) break;//sorting finished

				//pointer swapping
				swich=*(quantum_dot+i);
				*(quantum_dot+i)=*(quantum_dot+j);; //num_alive+num_resurrected-sorted
				*(quantum_dot+j)=swich;

				occupancy[i]=true;
				occupancy[j]=false;
				e_local_old[i]=e_local_old[j];
				sorted++;
				j--;
			}		
			i++;
			if (i>=j) break;
		}
		num_alive=num_alive-killsd+num_resurrected;

		//renormalize (kill some of the walkers if they get to many) TODO Test algo!
		if (num_alive>initial_number_of_walkers*2)
		{
			int number_to_remove = num_alive-initial_number_of_walkers;
			for (int k=0;k<number_to_remove;k++, --num_alive)
			{
				i=(int)(RAN_UNI()*num_alive+1); //obs ran_uni in (0,1) i in [0,num_alive-1] XXX XXX
				if (i!=(num_alive))
				{
					//removing particle i, and moving the last particle to its position.
					walker* swich;
					swich=*(quantum_dot+i);
					*(quantum_dot+i)=*(quantum_dot+j); //num_alive+num_resurrected-sorted
					*(quantum_dot+num_alive)=swich;
					e_local_old[i]=e_local_old[num_alive];
					occupancy[num_alive]=false;
				}
				else
				{
					occupancy[i]=false;
				}
			}
		}
	//	if (num_alive<initial_number_of_walkers/2.) //randomly duplicate walkers??
	//	{
	//		//randomly create walkers
	//		j = num_alive-1;
	//		quantum_dot[j]->initialize;
	//		num_alive++;
	//	}
		
		cout<<setprecision(10)<<"RANK "<<myrank<<" process: (alive,resur,killsd)=(<"
			<<num_alive<<","<<num_resurrected<<","<<killsd<<")"
			<<"current mean E"<<e_cumulative/(double)total_loop_c_cumulative<<"\n";
		if (myrank==0)
		{
			cout<<"cycles in main loop: "<<loop_main+1<<" of "
				<<num_c_dmc_main_loop<<" e_trial<"<<e_trial<<"\n";
		}/*//endvimfold*/

	} //************************** END OF DMC sampling **************************

	//collecting results	

	//USE MPI_Reduce(MPI_IN_PLACE,...)	
	MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &e2_cumulative, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_loop_c_cumulative, 1, MPI_LONG, MPI_SUM,  MPI_COMM_WORLD);
	e_cumulative/=(double)total_loop_c_cumulative;
	e2_cumulative/=(double)total_loop_c_cumulative;

	if (myrank==0)
	{	
		cout<<"\nFinal energy<"<<setprecision(16)<<e_cumulative<<"\n";
		cout<<"variance"<<e2_cumulative-pow(e_cumulative,2);
	}

//###################################################################################
//         
// 			NOTES:  total number of accumulated cycles should be accumulated
// 					final statistics summed wrongly if num part not const!
//
// 			MEASURE: variance. Should be exactly 0
// 			RENORMALIZE: more often, then the workload will be approx the same on the different procs
// 					   : should also find a way to add particles to the system
// 					   : some thermalization loops in main loop before starting to collect data
// 					   : implement blocking (e,e2->disk)
//
// 					How do we know that the system is converged?? 
// 					1: the number of walkers stable.
// 					2: the total energy changing slowly
// 					3: chech statistical error with blocking
// 					
//###################################################################################

	//variance = (e_local_squared-e_local*e_local/(double)num_cycles)/(double)num_cycles;
#if WRITEOFB/*//startvimfold*/
	ofstream blockofile;
  	//char *blockoutfilename;
	ostringstream ost;
	ost <<OFPATHB<<
		//"p"<<num_part<<"a"<<var_par[1]<<"b"<<var_par[0]
		//<<"w"<<var_par[2]<<
		"r"<<myrank<<".dat";
	blockofile.open(ost.str().c_str(),ios::out|ios::binary);
	blockofile.write((char*)(all_energies+1),num_cycles * sizeof (double));
	blockofile.close();
	delete [] all_energies;
#endif
#if WRITEOFC
	ofilec.close();
#endif/*//endvimfold*/

}/*//endvimfold*/

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
