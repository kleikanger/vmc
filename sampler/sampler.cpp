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

double sampler::getEnergyGrad(int i)
{/*//startvimfold*/
	return energy_gradient[i];
}/*//endvimfold*/

void sampler::sample(int num_cycles, int thermalization, double* var_par, double delta_t, double* result)
{/*//startvimfold*/
#if WRITEOFB/*//startvimfold*/
	//Only for rank 0 proc
	//if (rank=0) ?
	//if (num_cycles>1.1e7)
    //ofile.open(OFPATHB);
	double* all_energies = new double[num_cycles + 1];
	//cout<<writing to file ((std::string)OFPATHB)((std::string)OFNAMEB);
#endif
#if WRITEOFC
	//if (num_cycles>1.1e7)
	//{
	//	
	//	cerr<<"Warning: sampler::sample(). Max size for num_cycles (=" 
	//		<<num_cycles<<")is set to 1e6 when WRITEOFC == true. Terminating";
	//	exit(1);
	//}
	ofstream ofilec;
	ofilec.open((OFPATHC));
#endif
#if CONJGRAD
	//dim: num of varpar-1, 2
	double e_grad_temp[2];
	double **e_grad_cum = (double**)matrix(2,2,sizeof(double));
	e_grad_cum[0][0] 	= e_grad_cum[1][0] = e_grad_cum[0][1]
						= e_grad_cum[1][1] = 0.0; //[num_var_par-1][2]
#endif/*//endvimfold*/
	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;

	int i,j,k,active_part;
	int accepted=0;
//sga test
#if 1
	
	ofstream ofilecga;
	ofilecga.open("cgadata/test.txt");
	
	int n_sga=1, m_sga=1;
	double len_energy_gradient_cum=0;
	double f_sga, e_local_sga=0;
	double var_par_cum[2];//num_of_var_par
	double var_par_cum2[2];//num_of_var_par
	
	double var_par_cum_t[2];//num_of_var_par
	double var_par_cum2_t[2];//num_of_var_par

	double sga_alph_c=0, sga_beta_c=0;
	double sga_alph_c_old=0, sga_beta_c_old=0;
	double sga_alph_c_oold=0, sga_beta_c_oold=0;

	double energy_gradient_old[2];
	energy_gradient_old[0]=energy_gradient_old[1]=1.;

	double m_v_sga[2];//num_of_var_par
	m_v_sga[0]=m_v_sga[1]=1;
	double m_sga_inc=1;

	var_par_cum[0]=var_par_cum[1]=0.0;
	var_par_cum2[0]=var_par_cum2[1]=0.0;
	var_par_cum_t[0]=var_par_cum[1]=0.0;
	var_par_cum2_t[0]=var_par_cum2[1]=0.0;
#endif

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
		e_local_temp = quantum_dot->calcLocalEnergy(var_par);
		e_local += e_local_temp;
		e_local_squared += e_local_temp*e_local_temp;	

#if CONJGRAD
		/*
		quantum_dot->getVarParGrad(e_grad_temp); //TODO NEW NAME
		e_grad_cum[0][0]+=e_grad_temp[0];
		e_grad_cum[0][1]+=e_grad_temp[1];
		e_grad_cum[1][0]+=e_grad_temp[0]*e_local_temp;
		e_grad_cum[1][1]+=e_grad_temp[1]*e_local_temp;
		*/
	/*	
		double de_l_da, de_l_db;
		//cout<<e_grad_temp[1]<<" "<<e_grad_temp[0]<<" \n";
		minimizeVarPar(de_l_da, de_l_db, e_grad_temp);
		//cout<<e_grad_temp[1]<<" "<<e_grad_temp[0]<<" \n";
		//cout<<de_l_da<<" "<<de_l_db<<" \n\n";
		
		de_l_db=de_l_db-e_grad_temp[0]*e_local_temp; //Wrong sign??
		de_l_da=de_l_da-e_grad_temp[1]*e_local_temp;
	
		de_l_da*=-1; //XXX XXX XXX XXX Correct sign. where is calculations wrong??
		//de_l_db*=-1;

		e_grad_cum[0][0]+=de_l_db;
		e_grad_cum[0][1]+=de_l_da;
		e_grad_cum[1][0]+=de_l_db*e_local_temp;
		e_grad_cum[1][1]+=de_l_da*e_local_temp;
	*/	
#endif
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
//SGA/*//startvimfold*/
#if 0 
		int sga_upd= 3000;//1500?
		int sga_therm=0;
		if (loop_c%100==0)
		//if ((loop_c%(sga_upd+sga_therm)>sga_therm)||(loop_c%(sga_upd+sga_therm)==0))
		{
			quantum_dot->getVarParGrad(e_grad_temp); //TODO NEW NAME
			e_grad_cum[0][0]+=e_grad_temp[0];
			e_grad_cum[0][1]+=e_grad_temp[1];
			e_grad_cum[1][0]+=e_grad_temp[0]*e_local_temp;
			e_grad_cum[1][1]+=e_grad_temp[1]*e_local_temp;
			e_local_sga+=e_local_temp;

		} 
		if (loop_c%(sga_upd+sga_therm)==0)
		{
			sga_upd=30;

			//energy minimization
			result[0] = (e_local_sga/(double)sga_upd);

			energy_gradient[0] = 2*(e_grad_cum[1][0]-e_grad_cum[0][0]*result[0])/(double)sga_upd;//!!
			energy_gradient[1] = 2*(e_grad_cum[1][1]-e_grad_cum[0][1]*result[0])/(double)sga_upd;//!!;
#if 1
			double len_energy_gradient = sqrt(cblas_ddot(2,energy_gradient,1,energy_gradient,1));	
			if (n_sga>30&&n_sga<70)
			{
				len_energy_gradient_cum+=len_energy_gradient;
			}
			else if (n_sga==70)
			{
				m_v_sga[0] = 
				m_v_sga[1] = pow(400.*(len_energy_gradient_cum/40.),.8);
				m_sga_inc  = m_v_sga[0]/10.; //1/20 gives faster convergence, but a more unstable algo
			}
			
			//setting max length to move to .05
			if ( (len_energy_gradient*pow((double)m_v_sga[0],-.8)>.01) //.5??
			|| (len_energy_gradient*pow((double)m_v_sga[1],-.8)>.01) )
			{
				energy_gradient[0]*=sqrt(0.01)/len_energy_gradient;
				energy_gradient[1]*=sqrt(0.01)/len_energy_gradient;
			}
			
			//can be proven to converge mathematically
			if (energy_gradient[0]/energy_gradient_old[0]<0) m_v_sga[0]+=m_sga_inc;
			if (energy_gradient[1]/energy_gradient_old[1]<0) m_v_sga[1]+=m_sga_inc;
			//m_v_sga[0]+=m_sga_inc;
			//m_v_sga[1]+=m_sga_inc;
			
			var_par[0]-=energy_gradient_old[0]*pow((double)m_v_sga[0],-.8);
			var_par[1]-=energy_gradient_old[1]*pow((double)m_v_sga[1],-.8);

			//prevent unphysical negative values of alpha,beta
			if (var_par[0]<0.01) var_par[0]=0.01;
			if (var_par[1]<0.01) var_par[1]=0.01;
#endif
			energy_gradient_old[0]=energy_gradient[0];
			energy_gradient_old[1]=energy_gradient[1];

			n_sga++;
			quantum_dot->wSetVarPar(var_par);
			//set variables to 0		
			e_grad_cum[0][0]=e_grad_cum[0][1]=e_grad_cum[1][0]=e_grad_cum[1][1]=0.0;
			result[0]=result[1]=0;
			e_local_sga=0.;

			int sga_out_upd=100; //Start collecting data after 100 cycles (therm..)
			if (n_sga>sga_out_upd) //make better test
			{
				//collect data
				var_par_cum[0]+=var_par[0];
				var_par_cum[1]+=var_par[1];
				var_par_cum2[0]+=var_par[0]*var_par[0];
				var_par_cum2[1]+=var_par[1]*var_par[1];
				
				//calculate error+std.dev
				if (n_sga%2==0)
				{
					cout<<"\r"<<"                                                                                                ";
					cout<<"\r"
						<<setprecision(8)//<<setw(10)
						<<"b_C: "<<var_par_cum[0]/(double)(n_sga-sga_out_upd)
						<<" a_C: "<<var_par_cum[1]/(double)(n_sga-sga_out_upd)
						<<"b: "<<var_par[0]
						<<" a: "<<var_par[1]
				    	<<" b_V: "<<sqrt(-pow(var_par_cum[0]/(double)(n_sga-sga_out_upd),2)+var_par_cum2[0]/(double)(n_sga-sga_out_upd))
						<<" a_V: "<<sqrt(-pow(var_par_cum[1]/(double)(n_sga-sga_out_upd),2)+var_par_cum2[1]/(double)(n_sga-sga_out_upd))
						<<" c: "<<(n_sga*sga_upd)
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
							<< n_sga*sga_upd << " " 							//loops
							<< var_par[0] << " " 								//alpha
							<< var_par_cum[0]/(double)(n_sga-sga_out_upd) <<" " //alpha mean
							<< var_par[1] << " " 								//beta
							<< var_par_cum[1]/(double)(n_sga-sga_out_upd) <<" "; //beta mean
					}
					ofilecga<<"\n";
				}
			}
		}/*//endvimfold*/



#endif
	} //************************** END OF MC sampling **************************

	ofilecga.close();


	//KEEP THE BELOW LINES WHILE DEVELOPING
	//cout<<setprecision(4)<<"alpha: "<<var_par[1];
	//cout<<setprecision(4)<<" beta: "<<var_par[0]<<"\t";
	//cout<<setprecision(10)<< "Local energy: "<<(e_local/(double)num_cycles);	
	//cout<<setprecision(5)<<" variance^2: "<<(e_local_squared-e_local*e_local/num_cycles)/num_cycles;
	//cout<<setprecision(5)<<" Acc.rate: "<<accepted/(double)(num_cycles*num_part)<<"\n";
	//return values
	result[0] = (e_local/(double)num_cycles);
	result[1] = (e_local_squared-e_local*e_local/(double)num_cycles)/(double)num_cycles;
#if CONJGRAD
	//for CGM minimization. see: sampler::getEnergyGrad()
	
	//energy minimization
	//energy_gradient[0] = 2*(e_grad_cum[1][0]-e_grad_cum[0][0]*result[0])/(double)num_cycles;
	//energy_gradient[1] = 2*(e_grad_cum[1][1]-e_grad_cum[0][1]*result[0])/(double)num_cycles;
	
	//variance minimization
	energy_gradient[0] = 2.*(e_grad_cum[1][0]-e_grad_cum[0][0]*result[0])/(double)num_cycles;
	energy_gradient[1] = 2.*(e_grad_cum[1][1]-e_grad_cum[0][1]*result[0])/(double)num_cycles;
#endif
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

//Experimental
#if 1
void sampler::minimizeVarPar(double &de_l_da, double &de_l_db, double* e_grad_temp)
{
	//quantum_dot->getdELdvElem(de_l_db, de_l_da, dbJ_ove_J, daD_ove_D)
	quantum_dot->getdELdvElem(de_l_db, de_l_da, e_grad_temp[0], e_grad_temp[1]);
}
#endif

// For vim users: Defining vimfolds.
/// vim:fdm=marker:fmr=//startvimfold,//endvimfold
