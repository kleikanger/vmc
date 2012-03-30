#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <fstream>

#include "mcongrid.h"
#include "../sampler/sampler.h"

//include defines (OFPATHA, OFPATHB, WRITEOFA, WRITEOFA)
#include "../definitions/mcongrid_Def.h"
//#define WRITEOFB true
//#define OFPATHB  "blocking/B_2012.3.9.95041"
//#define WRITEOFA true
//#define OFPATHA "datafiles/opd_2012.3.9.95041.dat"

using std::cout;
using std::cerr;
//using std::ofile;
using std::ofstream;
using std::setprecision;
using std::setiosflags;
using std::setw;
using std::ios;

void writeToFile(double** variational_result_e, double** variational_result_e2,
		double* var_par,double* var_par_inc,int* var_par_cyc, int num_cycles, 
		int thermalization, double delta_t, int num_part, int nprocs);

void writeToScreen(double* var_par, double* var_par_inc, int* var_par_cyc, 
		double** variational_result_e, double** variational_result_e2);

mcongrid::mcongrid(){}
mcongrid::~mcongrid(){}

void mcongrid::runMCongrid(
			int dimension,
			int num_cycles,
			int thermalization,
			int num_part,
			int spin_up_cutoff,
			double delta_t,
			double* var_par,
			double* var_par_inc,
			int* var_par_cyc,
			int num_of_var_par,
			int myrank,
			int nprocs  )
{

	//construct sampler object
	sampler sampler_(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);

	if (myrank==0)
	{
		//output to screen
		cout<<"\n";
		cout<<"------------------------------------------------------\n";
		cout<<"\n";
		cout<<"num_cycles:     \t\t"<<num_cycles*nprocs<<"\n";
		cout<<"thermalization: \t\t"<<thermalization<<"\n";
		cout<<"delta_t:        \t\t"<<delta_t<<"\n";
		cout<<"num_part:       \t\t"<<num_part<<"\n";
		cout<<"nprocs:         \t\t"<<nprocs<<"\n";
		cout<<"omega:          \t\t"<<var_par[2]<<"\n";
		cout<<"------------------------------------------------------\n\n";
		cout<<"alpha values:\n";	
		for (int i=0;i<var_par_cyc[0];i++) 
		{ cout<<setprecision(10)<<setw(6)<<var_par[1]+(1.+(double)i)
			* var_par_inc[1]<<", "; }
		cout<<"\nbeta values:\n";	
		for (int i=0;i<var_par_cyc[0];i++) 
		{ cout<<setprecision(10)<<setw(6)<<var_par[0]+(1+(double)i) 
			* var_par_inc[0]<<", "; }
		cout<<"\n\n------------------------------------------------------\n\n";
	}	
	
	//storing the variational result
	double** variational_result_e = new double*[var_par_cyc[0]];
	for (int i=0;i<var_par_cyc[0];i++) 
		variational_result_e[i] = new double[var_par_cyc[1]];
	double** variational_result_e2 = new double*[var_par_cyc[0]];
	for (int i=0;i<var_par_cyc[0];i++) 
		variational_result_e2[i] = new double[var_par_cyc[1]];

	//temporary storage of energies and sigmas
	double* var_res_temp = new double[2];
	
	//loop over variational parameters
	int loop_v0, loop_v1;
	for (loop_v0=0; loop_v0<var_par_cyc[0]; loop_v0++)
	{
		var_par[0]+=var_par_inc[0];
		for (loop_v1=0; loop_v1<var_par_cyc[1]; loop_v1++)
		{
			var_par[1]+=var_par_inc[1];
			
			//start sampling
			sampler_.sample(num_cycles,thermalization,var_par,delta_t,var_res_temp);
			//collecting energy and energy squared from all processes
			MPI_Reduce(&var_res_temp[0], &variational_result_e[loop_v0][loop_v1],
					1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&var_res_temp[1], &variational_result_e2[loop_v0][loop_v1],
					1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			variational_result_e[loop_v0][loop_v1] /= nprocs;
			variational_result_e2[loop_v0][loop_v1] /= nprocs;
		}
		var_par[1]-=var_par_inc[1]*loop_v1;
	}
	var_par[0]-=var_par_inc[0]*loop_v0;
	
	//Write results to screen..
	if (myrank==0)
		writeToScreen(var_par, var_par_inc, var_par_cyc, 
				variational_result_e, variational_result_e2);
	
	//write to file
#if WRITEOFA
	if (myrank==0)
		writeToFile(variational_result_e,variational_result_e2, var_par,var_par_inc,
				var_par_cyc, num_cycles, thermalization, delta_t, num_part, nprocs);
#endif
	cout<<"\a";

}/*//endvimfold*/

void writeToScreen(double* var_par, double* var_par_inc, int* var_par_cyc, 
		double** variational_result_e, double** variational_result_e2)
{/*//startvimfold*/
	cout<<"\nenergy\naHOR,bVERT";	
	for (int i=0;i<var_par_cyc[1];i++) 
	{ 
		cout<<setprecision(10)<<setw(18)<<var_par[1]+(1.+(double)i)*var_par_inc[1]<<" "; 
	}
	cout<<"\n";	
	for (int i=0;i<var_par_cyc[0];i++) 
	{
		cout<<setprecision(10)<<setw(18)<<var_par[0]+(1+(double)i)*var_par_inc[0]<<" ";
		for (int j=0;j<var_par_cyc[1];j++) 
		{ 
			cout<<setprecision(10)<<setw(18)<<variational_result_e[i][j]; 
		}
	cout<<"\n";
	}
	cout<<"\nvariance\naHOR,bVERT";	
	for (int i=0;i<var_par_cyc[1];i++) 
	{ 
		cout<<setprecision(10)<<setw(18)<<var_par[1]+(1.+(double)i)*var_par_inc[1]<<" "; 
	}
	cout<<"\n";	
	for (int i=0;i<var_par_cyc[0];i++) 
	{
		cout<<setprecision(10)<<setw(18)<<var_par[0]+(1+(double)i)*var_par_inc[0]<<" ";
		for (int j=0;j<var_par_cyc[1];j++) 
		{ 
			cout<<setprecision(10)<<setw(18)<<variational_result_e2[i][j]; 
		}
	cout<<"\n";
	}
}/*//endvimfold*/

void writeToFile(double** variational_result_e, double** variational_result_e2,
		double* var_par,double* var_par_inc,int* var_par_cyc, int num_cycles, /*//startvimfold*/
		int thermalization, double delta_t, int num_part, int nprocs)
{
	ofstream ofile;
	ofile.open((OFPATHA));
		ofile << setiosflags(ios::showpoint | ios::uppercase);
		ofile<<"num_cycles*nproc: "<<num_cycles*nprocs<<"\n";
		ofile<<"thermalization:   "<<thermalization<<"\n";
		ofile<<"delta_t:          "<<delta_t<<"\n";
		ofile<<"num_part:         "<<num_part<<"\n";
		ofile<<"omega             "<<var_par[2]<<"\n"; 
		ofile<<"alpha_variations  "<<var_par_cyc[1]<<"\n"; 
		ofile<<"beta_variations   "<<var_par_cyc[0]<<"\n"; 
		ofile<<"------------------------------------------------------\n";
	//cout<<writing to file ((std::string)OFPATHB)((std::string)OFNAMEB);
		ofile<<"alpha values:\n";
		for (int i=0;i<var_par_cyc[0];i++)
        { 
			ofile<<setprecision(10)<<setw(16)<<var_par[1]+(1.+(double)i)*var_par_inc[1]; 
		}
        ofile<<"\nbeta values:\n";
		for (int i=0;i<var_par_cyc[0];i++)
		{ 
			ofile<<setprecision(10)<<setw(16)<<var_par[0]+(1+(double)i)*var_par_inc[0]; 
		}
		ofile<<"\n------------------------------------------------------\n";
		ofile<<"\nenergy (alpa horisontal axis, beta vertical axis)\n";	
		for (int i=0;i<var_par_cyc[0];i++) 
		{
			for (int j=0;j<var_par_cyc[1];j++) 
			{ 
				ofile<<setprecision(16)<<setw(28)<<variational_result_e[i][j]; 
			}
			ofile<<"\n";
		}
		ofile<<"\nvariance (alpa horisontal axis, beta vertical axis)";	
		ofile<<"\n";	
		for (int i=0;i<var_par_cyc[0];i++) 
		{
			for (int j=0;j<var_par_cyc[1];j++) 
			{ 
				ofile<<setprecision(16)<<setw(28)<<variational_result_e2[i][j]; 
			}
		ofile<<"\n";
		}
	ofile.close();
}/*//endvimfold*/
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold


