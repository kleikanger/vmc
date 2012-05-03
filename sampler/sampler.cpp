#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>
#include <sstream>
#include "sampler.h"
#include "../newmatrix/newmatrix.h"

using std::cout;
using std::cerr;
//using std::ofile;
using std::ofstream;
using std::setprecision;
using std::setiosflags;
using std::setw;
using std::ios;
using std::ostringstream;

//include definitions
#include "../definitions/sampler_Def.h"

sampler::sampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank)
{/*//startvimfold*/
	this->myrank=myrank;
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;
	quantum_dot = new walker(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
	
	//OBS : OMEGA PART OF VARPAR TODO CHANGE 
	// NUM OF VAR PAR = 3 : TODO set to 2
//#if CONJGRAD
	energy_gradient = new double[num_of_var_par];
//#endif
}/*//endvimfold*/

sampler::~sampler()
{/*//startvimfold*/
	delete quantum_dot;
}//End function /*//endvimfold*/

//REMOVE XXX NOT NECC WHEN CGM REMOVED
double sampler::getEnergyGrad(int i)
{/*//startvimfold*/
	return energy_gradient[i];
}/*//endvimfold*/

void sampler::sample(int num_cycles, int thermalization, double* var_par, double delta_t, double* result)
{/*//startvimfold*/
#if WRITEOFB/*//startvimfold*/
	double* all_energies = new double[num_cycles + 1];
#endif
#if WRITEOFC
	ofstream ofilec;
	ofilec.open((OFPATHC));
#endif/*//endvimfold*/
	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;
	//kin,pot and ocs energies cum
	double energies[3];
	for (int i=0;i<3;i++) energies[i] = 0.0;

	int i,j,k,active_part;
	int accepted=0;

	//init walker object
	quantum_dot->initWalker(var_par, delta_t);

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
		//move all particles, one at a time
		for (int active_part=0; active_part<num_part; active_part++)
		{
			//metropolis-hastings test
			if ( quantum_dot->tryRandomStep(active_part) ) 
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
		e_local_temp = quantum_dot->calcLocalEnergy(energies);
		//e_local_temp = quantum_dot->calcLocalEnergy();
		e_local += e_local_temp;
		e_local_squared += e_local_temp*e_local_temp;	
#if WRITEOFB //write to file blocking data
		all_energies[loop_c-thermalization]=e_local_temp;
#endif
#if WRITEOFC // write to file : single particle density
		if (myrank==0)
		{
			double* x_v = new double[dimension];
			quantum_dot->getRi(0,x_v);	
			{
				ofilec << setiosflags(ios::showpoint | ios::uppercase);
				for (i=0;i<dimension;i++)
				{
					ofilec << setw(16) << setprecision(16) << x_v[i] <<"\t";
				}
				ofilec<<"\n";
			}
			delete [] x_v;
		}
#endif
	} //************************** END OF MC sampling **************************
	//write to screen
	cout<<setprecision(4)<< "(alpha,beta) = ("<<var_par[1]<<","<<var_par[0]<<"),";
	cout<<setprecision(10)<< " kinE: "<<(energies[0]/(double)num_cycles);	
	cout<<setprecision(10)<< ", pot(osc)E: "<<(energies[1]/(double)num_cycles);	
	cout<<setprecision(10)<< ", pot(ee)E: "<<(energies[2]/(double)num_cycles);	
	cout<<setprecision(5)<< ", Acc.rate: "<<accepted/(double)(num_cycles*num_part)<<"\n";
	//return values
	result[0] = (e_local/(double)num_cycles);
	result[1] = (e_local_squared-e_local*e_local/(double)num_cycles)/(double)num_cycles;
#if WRITEOFB
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
#endif

}/*//endvimfold*/

// For vim users: Defining vimfolds.
/// vim:fdm=marker:fmr=//startvimfold,//endvimfold
