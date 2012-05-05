#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
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

dmcsampler::dmcsampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank, int nprocs, MPI_Status status)
{/*//startvimfold*/
	this->myrank=myrank;
	this->nprocs=nprocs;
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;
	this->status=status;
	//factory class for manipulating walkers
	popCtr = new popControl(num_part, spin_up_cutoff, dimension, status);
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
	ofstream blockofile;
  	//char *blockoutfilename;
	ostringstream ost;
	ost <<OFPATHB<<
		//"p"<<num_part<<"a"<<var_par[1]<<"b"<<var_par[0]
		//<<"w"<<var_par[2]<<
		"r"<<myrank<<".dat";
	blockofile.open(ost.str().c_str(), ios::out | ios::binary );
	//blockofile.write((char*)(all_energies+1),n_all_e * sizeof (double));
	//blockofile.close();
	//delete [] all_energies;
	//Only for rank 0 proc
	//if (rank=0) ?
	//if (num_cycles>1.1e7)
    //ofile.open(OFPATHB);
	//double* all_energies = new double[num_c_dmc_inner_loop*num_c_dmc_main_loop*num_part*initial_number_of_walkers*2 + 1]; //resize evt write directly to file
	//int n_all_e=0;
	if (myrank==0)
		cout<<"writing blockingdata to file"<<(OFPATHB)<<"\n";
#endif
#if WRITEOFC
	ofstream ofilec;
	ofilec.open((OFPATHC));
#endif/*//endvimfold*/

	int idum_d=time(NULL)*(myrank+1);
	RAN_UNI_SET(&idum_d,5);
	double e_local=0.0;
	double redist_threshh=.01; //if particles deviates mean*redist_threshh from mean, redistribute particles
	double e_local_squared=0.0;
	double e_local_temp=0.0;
	int corr_length = 3000; //number of steps between init of particles
	int i,j,k,active_part;
	int accepted=0;
	int num_alive = initial_number_of_walkers;
	double branching_factor;
	int num_new_walkers;
	int num_resurrected = 0;
	double e_cumulative=0.0;
	double e2_cumulative=0.0;
	long int loop_c_cumulative;
	long int total_loop_c_cumulative=0;
	double e_ref=e_trial;
	int loop_main, loop_c, killsd=0; 
	
	if (myrank==0)
	cout<<"\nDMC simulation - \n"
		<<"alpha:                           "<<var_par[1]<<"\n"
		<<"beta:                            "<<var_par[0]<<"\n"
		<<"initial trial energy:            "<<e_trial<<"\n"
		<<"initial_number_of_walkers/procs: "<<initial_number_of_walkers<<"\n"
		<<"nprocs:                          "<<nprocs<<"\n"
		<<"equilibriation loops:            "<<num_c_dmc_equilibriation_loop<<"\n"
		<<"cycles main loop:                "<<num_c_dmc_main_loop<<"\n"
		<<"cycles inner loop:               "<<num_c_dmc_inner_loop<<"\n"
		<<"thermalization cycles:           "<<thermalization<<"\n"
		<<"delta_t:                         "<<delta_t<<"\n"
		<<"number of particles:             "<<num_part<<"\n\n";
	//array of walkers	
	// initializing 6x number_of_walkers. when 2x number_of_walkers, system renormalize to number_of_walkers. 
	quantum_dot = new walker*[initial_number_of_walkers*6]; 
	bool* occupancy = new bool[initial_number_of_walkers*6];
	double* e_local_old = new double[initial_number_of_walkers*6];
	
	//********************************
	//initializing walker objects
	//********************************
	quantum_dot[0] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
	quantum_dot[0]->initWalker(var_par, delta_t);
	for (i=1;i<initial_number_of_walkers*6;i++)
	{
		quantum_dot[i] = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
		quantum_dot[i]->initEmptyWalker(var_par, delta_t);//Faster algo
		//quantum_dot[i]->initWalker(var_par, delta_t);
	}
	//setting occupancy to 'true' for first initial_number_of_walkers.
	for (i=0;i<initial_number_of_walkers;i++)
		occupancy[i] = true;
	for (i=initial_number_of_walkers;i<initial_number_of_walkers*6;i++)
		occupancy[i] = false;
	//Initialize walkers, using VMC routine
	initializeSys(initial_number_of_walkers, thermalization, corr_length, var_par, e_local_old);

	
	//************************** START equilibriation phase ******************
	/*//startvimfold*/
	loop_c=0; 
	//TODO while not equilibriated :: find some criteria
	for (loop_main=0; loop_main<num_c_dmc_equilibriation_loop; loop_main++)
	{
		if (myrank==0&&loop_main%(((int)num_c_dmc_equilibriation_loop)/(200)+1)==0)
		{
			fflush(stdout);
			cout<<"\r"<<"                                    ";
			cout<<"\requilibriating walkers: "<<setprecision(3)<<(loop_main+1)/(double)(num_c_dmc_equilibriation_loop)*100<<"%";
		}
		//reset summation variables
		e_local=0.0; num_resurrected=0; killsd=0; loop_c_cumulative=0;;
		//loop over walkers
		for (int loop_p=0;loop_p<num_alive;loop_p++)
		{
			//energy accumulating cycles XXX changing E_T each cycle, maybe not neccesary? 
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
				e_local_temp = quantum_dot[loop_p]->calcLocalEnergy();
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
						e_local_old[j]=e_local_temp;
						occupancy[j]=true;
						num_resurrected++;
						j++;
						num_new_walkers--;	
					}
				}
				//update energy if the walker survived
				e_local_old[loop_p] = e_local_temp;
			}
		//keeping track of the total number of samples
		loop_c_cumulative+=loop_c;
		}
		//sorting walkers, num_alive called by reference and reset
		sortWalkers(num_alive, killsd, num_resurrected, occupancy, e_local_old);	
		
		//New reference energy
		e_cumulative+=e_local;
		total_loop_c_cumulative+=loop_c_cumulative;
		double e_c_temp; long int t_l_c_temp;
		MPI_Allreduce(&e_cumulative, &e_c_temp, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
		MPI_Allreduce(&total_loop_c_cumulative, &t_l_c_temp, 1, MPI_LONG, MPI_SUM,  MPI_COMM_WORLD);
		//0.01 test other values
		e_ref = e_c_temp/(double)t_l_c_temp-log((double)num_alive /
				(double)initial_number_of_walkers)/delta_t*0.01;
		
		//XXX 
		//using different routine for equilibriation with fewer MPI_Calls
		//Not so many send and recieves (5000 here). Seems to work well.
		//e_cumulative+=e_local/(double)loop_c_cumulative;	
		//e_ref = e_cumulative/(double)(loop_main+1)-.01*log((double)num_alive
		// 		/ (double)initial_number_of_walkers)/delta_t;
	}	
	if (myrank==0)
	{
		fflush(stdout);
		cout<<"\r"<<"                                    ";
		cout<<"\requilibriatinon of walkers finished";
		cout<<"\n";
	}
	
	MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_loop_c_cumulative, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	e_cumulative/=(double)total_loop_c_cumulative;

	//using different routine for equilibriation with fewer MPI_Calls
	//MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//e_cumulative/=(double)nprocs*loop_main;
	
	if (myrank==0)
	{
		cout<<"\nequilibriation energy: "
		    <<setprecision(16)<<e_cumulative<<"\n\n";
	}
	/*//endvimfold*/
	
	//************************** END OF equilibriation phase *****************
	//************************** START sampling phase ************************

	ofstream ofiletest;
	if (myrank==0)	
	{
		ofiletest.open("dmctest.dat");
	}
	//XXX XXX test XXX XXX
	double e_tot=0.0;

	//reset variables/*//startvimfold*/
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
					//kill walkers that has crossed node (error linear inst of quadratic? in delta_t)
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
				e_local_temp = quantum_dot[loop_p]->calcLocalEnergy();
				//calculate new branching factor
				branching_factor = exp(-delta_t*(0.5*(e_local_temp+e_local_old[loop_p])-e_trial)); //e_trial is the (vmc) calculated mean
				
				//collect energy weighted by the branching factor
				e_local+=branching_factor*e_local_temp; //TODO use intermediate variable to calculate b_f*e_l_t
				e2_cumulative+=pow(branching_factor*e_local_temp,2);
#if WRITEOFB
		//blocking analysis data to file
		//blockofile<<setiosflags(ios::binary);//showpoint | ios::uppercase );
		//blockofile<<(char)branching_factor*e_local_temp<<"\n";
		double eb_temp = branching_factor*e_local_temp;
		blockofile.write((char*)&eb_temp, sizeof (double));
		//all_energies[n_all_e]=branching_factor*e_local_temp;
		//n_all_e++;
#endif
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
		//total number of moves must be collected for final statistics
		total_loop_c_cumulative+=loop_c_cumulative;
		
		//****collect and broadcast Find new e_trial
		MPI_Allreduce(MPI_IN_PLACE, &e_local, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &loop_c_cumulative, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		e_trial=e_local/(double)loop_c_cumulative;
		//sorting walkers, num_alive called by reference and reset
		sortWalkers(num_alive, killsd, num_resurrected, occupancy, e_local_old);	
				
		//*******************************************************************
		//redistribute using mpi
		//*******************************************************************
		//only redistribute if too many/too few alive
		//if (abs(num_alive-mean)<threshh)
		// 	MPI_Bcast(redistribute==true);
		//MPI_Barrier(MPI_COMM_WORLD);
		//run redistr algo and update mean & threshh
		//if (redistribute) 
		redistributeWalkers(myrank, nprocs, num_alive, occupancy, e_local_old, redist_threshh);
		//redistribute=false;

		e_tot+=e_trial;
		if (myrank==0)
		{
			cout<<"\rcycles in main loop: "<<loop_main+1<<" of "
				<<num_c_dmc_main_loop<<" e_trial<"<<e_trial<<" e_tot approx:"<< e_tot/(double)(loop_main+1);
			fflush(stdout);
			//write to file
			ofiletest<<setprecision(16)<<e_tot/(double)(loop_main+1)<<" "<<e_trial<<" "<<num_alive*2.<<"\n";
		}/*//endvimfold*/

	} //************************** END OF DMC sampling **************************

	//collecting results	

	//USE MPI_Reduce(MPI_IN_PLACE,...)	
	MPI_Allreduce(MPI_IN_PLACE, &e_cumulative, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &e2_cumulative, 1, MPI_DOUBLE, MPI_SUM,  MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &total_loop_c_cumulative, 1, MPI_LONG, MPI_SUM,  MPI_COMM_WORLD);
	e_cumulative/=(double)total_loop_c_cumulative;
	e2_cumulative/=(double)total_loop_c_cumulative;

	//Write results to screen
	if (myrank==0)
	{	
		cout<<"\nFinal energy<"<<setprecision(16)<<e_cumulative<<"\n";
		cout<<"variance"<<e2_cumulative-pow(e_cumulative,2)<<"\n";
	}

#if WRITEOFB/*//startvimfold*/
	//ofstream blockofile;
	//ostringstream ost;
	//ost <<OFPATHB<<
	//	//"p"<<num_part<<"a"<<var_par[1]<<"b"<<var_par[0]
	//	//<<"w"<<var_par[2]<<
	//	"r"<<myrank<<".dat";
	//blockofile.open(ost.str().c_str(),ios::out|ios::binary);
	//blockofile.write((char*)(all_energies+1),n_all_e * sizeof (double));
	blockofile.close();
	//delete [] all_energies;
#endif
#if WRITEOFC
	ofilec.close();
#endif/*//endvimfold*/
	if (myrank==0) ofiletest.close();

}/*//endvimfold*/
void dmcsampler::sortWalkers(int &num_alive, int killsd, 
		int num_resurrected, bool *occupancy, double *e_local_old)	
{/*//startvimfold*/
		int sorted=0, i=0; 
		int j=num_alive+num_resurrected-1;//starting on the last one
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
}/*//endvimfold*/
void dmcsampler::redistributeWalkers(int myrank, int nprocs, int &num_alive, bool *occupancy, double *e_local_old, double redist_threshh)
{/*//startvimfold*/

	//input: redist threshh,num_alive,myrank, walker *quantum_dot, popControl *popCtr
	int i, j; 
	int threshh, nalive_temp, lefto, mean=0, ntransfer;
	int rank[nprocs];
	int nalive[nprocs];
	
	for (i=0;i<nprocs;i++)
		rank[i]=i;
	nalive[myrank]=num_alive;

	for (i=0;i<nprocs;i++)
		MPI_Bcast(&nalive[i], 1, MPI_INT, i, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);//XXX necc? (feels safe)

	//mean of num_alive and the modulus lefto
	for (i=0;i<nprocs;i++)
		mean+=nalive[i];
	lefto=mean%nprocs;
	mean/=nprocs;
	threshh=(int)((double)redist_threshh*mean);
	//evt check if max element is larger than abs(mean-max_elem)>thresh+lefto before sorting

	//The philosophy behind this algo is that no walker should be transfered twice.
	//if we transfer fron node 0 to node 1, there should never be any transfers from node 1.
	//Transfering data is time consuming!
	//All nodes runs to this algo simultaniously to keep communication on a minimum level.
	//Only the data transfer is private to the nodes.
	//communication between the different nodes will happen simultaniously.
	bool cont_iter=true;
	while (cont_iter)
	{
		//continiue iterating while
		cont_iter=false;

		int rank_temp;
		//picksort:Num.Rec. inspired algo.
		//Sorting nalive and rank, Smallest nalive and corresponding rank first
		for (i=nprocs-2;i>=0;i--)
		{
			rank_temp=rank[i];
			nalive_temp=nalive[i];
			//rank_temp=i;//rank[i];
			for (j=i;j<=nprocs-2;j++)
			{
				nalive[j]=nalive[j+1];
				rank[j]=rank[j+1];	
				if (nalive[j+1]>=nalive_temp) break;
			}
			nalive[j]=nalive_temp;
			rank[j]=rank_temp;
		}
		//Transfering between nodes	
		for (i=0;i<nprocs/2;i++) //obs: note integer division
		{
			//simultaneously transfering betw all procs
			//if all ifs fail, coun_iter remains 'false' and the algo is finished
			if  ( (abs(nalive[i]-mean)>threshh+lefto) || 
					(abs(nalive[nprocs-i-1]-mean)>threshh+lefto) )
				//if both nodes has a higher/lower number of particles 
				//than the mean, dont transfer, wait till next loop when sorted	
				if ( (nalive[i]<mean) != (nalive[nprocs-i-1]<mean) )	
			{  
				//transfer walkers till one of the nodes has a number eq to the mean
				//if one is eq to mean: ntransfer=0, resorting and new proposal
				ntransfer = ( (mean-nalive[i]) >= (nalive[nprocs-i-1]-mean)) ? 
					nalive[nprocs-i-1]-mean : mean-nalive[i];
				if  ((myrank==rank[i]) || (myrank==rank[nprocs-i-1]))
				{
					//redistribute send ntransfer particles from nprocs-i-1->i usinc popControl method	
					for (j=0;j<ntransfer;j++)
					{
						if (myrank==rank[nprocs-i-1])
						{
							num_alive--;
							popCtr->transmitWalker(quantum_dot[num_alive-1], myrank, rank[i], myrank);
							occupancy[num_alive]=false;
							MPI_Send(&e_local_old[num_alive], 1, MPI_DOUBLE, rank[i], 499, MPI_COMM_WORLD);	
						}
						if (myrank==rank[i])
						{
							popCtr->transmitWalker(quantum_dot[num_alive], rank[nprocs-i-1], myrank, myrank);
							occupancy[num_alive]=true;
							MPI_Recv(&e_local_old[num_alive], 1, MPI_DOUBLE, rank[nprocs-i-1], 499, MPI_COMM_WORLD, &status);	
							num_alive++;
						}
					}
					//XXX TESTING XXX
//					if (myrank==rank[i])
//						cout<<rank[nprocs-i-1]<<"to"<<rank[i]<<" number  "<<ntransfer<<"\n"; 
				}
				nalive[i]+=ntransfer;
				nalive[nprocs-i-1]-=ntransfer;
			
				//last run if all 
				if (abs(nalive[nprocs-i-1]-mean)>=lefto) cont_iter=true;
				if (abs(nalive[i]-mean)>=lefto) cont_iter=true;
			}
		}
	}
}/*//endvimfold*/

void dmcsampler::initializeSys(int initial_number_of_walkers, int thermalization, int corr_length, double *var_par, double *e_local_old)
{	/*//startvimfold*/
	int loop_c=0, num_init=1; 
	if (myrank==0)
		{
			cout<<"\rStarting initialization of walkers... ";
			fflush(stdout);
		}
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
			e_local_old[num_init] = quantum_dot[0]->calcLocalEnergy();
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
	e_local_old[0] = quantum_dot[0]->calcLocalEnergy();
	
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
