#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include "sampler.h"

using std::cout;
using std::cerr;
//using std::ofile;
using std::ofstream;
using std::setprecision;
using std::setiosflags;
using std::setw;
using std::ios;

//using namespace std;

//name and path of ofile
#ifndef WRITEOFB
#define WRITEOFB false
#endif
#ifndef OFPATHB
#define OFPATHB "/home/karleik/masterProgging/vmc/datafiles/zerotermalization.dat"
#endif

sampler::sampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank)
{/*//startvimfold*/
	this->myrank=myrank;
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;
	quantum_dot = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
}/*//endvimfold*/

sampler::~sampler()
{/*//startvimfold*/
	delete quantum_dot;
}//End function /*//endvimfold*/

void sampler::sample(int num_cycles, int thermalization, double* var_par, double delta_t, double* result)
{/*//startvimfold*/

#if WRITEOFB
	//Only for rank 0 proc
	//if (rank=0) ?
	if (num_cycles>1.1e7)
	{
		//
		cerr<<"Warning: sampler::sample(). Max size for num_cycles (=" 
			<<num_cycles<<")is set to 1e6 when WRITEOFB == true. Terminating";
		exit(1);
	}
	ofstream ofile;
	ofile.open((OFPATHB));
	ofile<<"num_cycles:     "<<num_cycles<<"\n";
	ofile<<"thermalization: "<<thermalization<<"\n";
	ofile<<"delta_t:        "<<delta_t<<"\n";
	ofile<<"num_part:       "<<num_part<<"\n";
	ofile<<"omega           "<<OMG<<"\n"; 
	ofile<<"alpha           "<<var_par[1]<<"\n"; 
	ofile<<"beta            "<<var_par[0]<<"\n"; 
	ofile<<"------------------------------------------------------\n";
	//cout<<writing to file ((std::string)OFPATHB)((std::string)OFNAMEB);
#endif

	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;

	int i,j,k,active_part;
	int accepted=0;
	
	//diffusion const D * delta_t = 0.5 * delta_t
	double dt_x_D = 0.5 * delta_t;

	//INITIATE ATOM
	quantum_dot->initWalker(var_par, delta_t);

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
		//move all particles, one at a time
		for (int active_part=0; active_part<num_part; active_part++)
		{
			//possible to move two part. at a time, one in each SD?
			bool metropolis_test = quantum_dot->tryRandomStep(active_part);

			//metropolis-hastings test
			if ( metropolis_test ) 
			{
				quantum_dot->acceptStep(active_part);
				if (loop_c>thermalization) { accepted++; } 
			}
			else  
			{
				quantum_dot->rejectStep(active_part);
			}
		}//All particles moved
		//if thermalization finished, start collecting data.
		if (loop_c<thermalization) { continue; }
		e_local_temp = quantum_dot->calcLocalEnergy(var_par);
		e_local += e_local_temp;
		e_local_squared += e_local_temp*e_local_temp; 	
#if WRITEOFB
		//evt:
		/*
		e_chunk+=e_local;
		//chunksize = 10
		if (!(loop_c%10)) {
		ofile.write(e_chunk)
		e_chunk=0.0;
		}
		   */
		if (myrank==0)
		{
			ofile << setiosflags(ios::showpoint | ios::uppercase);
			ofile << setw(16) << setprecision(16) << e_local_temp <<"\n";
		}
#endif
	} //************************** END OF MC sampling **************************

	cout<<setprecision(4)<<"alpha: "<<var_par[1];
	cout<<setprecision(4)<<" beta: "<<var_par[0]<<"\t";
	cout<<setprecision(10)<< "Local energy: "<<(e_local/(double)num_cycles);	
	cout<<setprecision(5)<<" variance^2: "<<(e_local_squared-e_local*e_local/num_cycles)/num_cycles;
	cout<<setprecision(5)<<" Acc.rate: "<<accepted/(double)(num_cycles*num_part)<<"\n";
	//cout<<"beta :\t\t"<<var_par[0]<<"\n";

	//return values
	result[0]=(e_local/(double)num_cycles);
	result[1]=(e_local_squared-e_local*e_local/(double)num_cycles)/(double)num_cycles;

#if WRITEOFB
	ofile.close();
#endif

}/*//endvimfold*/


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
