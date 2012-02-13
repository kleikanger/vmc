#include <cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "../lib/lib.h"
#include <iomanip>

#include "vmcsolver.h"
//#include "../QDslater/slaterMatrix.h"
//#include "../ipdist/ipdist.h"

int main()
{/*//startvimfold*/

	int num_cycles=1e8;
	int thermalization=0.3*num_cycles;//num_cycles*.5;
	int num_part=2;
	int spin_up_cutoff=1;
	int dimension=2;
	
	double ideal_step=2.4;
	//double omega = 1.0;

	//Number of variational parameters
	int num_of_var_par=1;	
	//Variational param
	double  var_par[num_of_var_par];
	
	var_par[0]=0.398;//beta (jastrow)	
	
	//Increase of variational parameters
	double var_par_inc[num_of_var_par];
	var_par_inc[0]=0.05;
	//Number of variations in each var. par.
	int var_par_cyc[num_of_var_par];
	var_par_cyc[0]=10;

	//construct vmcsolver object
	vmcsolver solver(num_part, spin_up_cutoff, dimension, num_of_var_par);	
	
	//output to screen
	cout<<"\n";
	cout<<"num_cycles: \t\t"<<num_cycles<<"\n";
	cout<<"thermalization:\t\t"<<thermalization<<"\n";
	cout<<"ideal_step: \t\t"<<ideal_step<<"\n";
	cout<<"beta: \t\t"<<var_par[0]<<"\n";

	//loop over variational parameters
	//start sampling
	solver.sample(num_cycles, thermalization, var_par);//pass on array returning e_local,variance.
	solver.sample(num_cycles, thermalization, var_par);//pass on array returning e_local,variance.
	solver.sample(num_cycles, thermalization, var_par);//pass on array returning e_local,variance.
	solver.sample(num_cycles, thermalization, var_par);//pass on array returning e_local,variance.
	solver.sample(num_cycles, thermalization, var_par);//pass on array returning e_local,variance.

}/*//endvimfold*/

vmcsolver::vmcsolver(int num_part, int spin_up_cutoff,int dimension, int num_of_var_par)
{/*//startvimfold*/
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;

	//init pos.vecs.
	//double* r = new (nothrow) double[n];
	//if (r==0) { cout<<"...\n"; exit(1);}
	r_old = new double*[num_part];
	for (int i=0; i<num_part; i++){
		r_old[i] = new double[dimension];
		//r_new[i] = new double[dimension];
		for (int j=0;j<dimension;j++)  
			r_old[i][j]=0.0;
	}
	//handling of the slater matrix
	slater = new slaterMatrix(num_part,spin_up_cutoff,num_of_var_par,dimension);
	//handling of the interparticle distances
	ipd = new ipdist(num_part,dimension,spin_up_cutoff);
	
	//Gradients of jastrow and laplacian
	jas_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) jas_grad[i]=new double[dimension];
	wf_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) wf_grad[i]=new double[dimension];
}/*//endvimfold*/

void vmcsolver::sample(int num_cycles, int thermalization, double* var_par)
{/*//startvimfold*/
	int accepted=0;
	double e_kinetic, r_12, r_12_new, e_potential, wf_R, r_squared, temp;
	double e_local=0.0;
	double e_local_squared=0.0;
	int i,j,k;
	int active_part;

	double ideal_step=2.2;
	
  	idum= - (time(NULL));//*(myrank));
	
	//''random'' startposition 
	for (int i = 0; i < num_part; i++) { 
		for (int j=0; j < dimension; j++) {
			//r_new[i][j]=
			r_old[i][j] = ideal_step*( ran2(&idum) - 0.5 );
		}
	}
	slater->updateVariationalParameters(var_par);
	slater->initSlaterMatrix(r_old);	
	slater->findInverse();
	ipd->init(r_old);
	
	//REMOVE
	int accepted_vec[num_part];
	for (i=0;i<num_part;i++) accepted_vec[i]=0;
	//interparticle distances updated.
	double* ipd_upd = new double[num_part-1];
	//storing new positions
	double* r_new = new double[dimension];
	//for summation using cblas ddot
	double* ones = new double[num_part];
	for (int i=0;i<num_part;i++) { ones[i]=1.0; }
	

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
		//Moving one particle at a time. (converges slowly?)
		//possible to try to move all part before collecting energy. 
		//or move two part. at a time, one in each SD?
		//*******
		active_part=loop_c%(num_part);
		//get new r_i (r_new) and r_{ij} (ipd_upd) for active particle
		getNewPos(active_part, ideal_step, r_new, ipd_upd);
	
		//jastrow ratio: 2 particle impl.
		//********************************
		//Summing over ipd_upd (sum r_12)
		r_12_new = cblas_ddot(num_part-1,ipd_upd,1,ones,1);// sum
		//corresponding elements in old vector
		r_12 = ipd->sumPart(active_part);
		double jas_R = exp( r_12_new/(1.0 + var_par[0]*r_12_new) - r_12/(1.0 + var_par[0]*r_12 ) );
		//********************************

		//wave func ratio
		wf_R=slater->waveFunction(r_new,active_part);
		//Metropolis test: update slater+ipd+positions if accepted.
		if ( ran2(&idum) <= (jas_R * jas_R * wf_R * wf_R ) ) 
		{
			//update positions	
			cblas_dcopy(dimension,&r_new[0],1,&r_old[active_part][0],1);
			//update inverse and slatermatrix
			slater->update(r_new,active_part);
			//update interparticle distances
			ipd->update(ipd_upd,active_part);
			//update acceptancerate
			if (loop_c<thermalization) { continue; } 
			//restart for loop or add one to accepted
				//DELETE accepted_vec
				accepted_vec[active_part]++;
				accepted++;
		}
		//if thermalization finished, start collecting data.
		else if (loop_c<thermalization) { continue; }
		//add local energy to variables e_local and e_local_squared.
		calcLocalEnergy(&e_local, &e_local_squared, var_par);
	}
	
	//************************** END OF MC sampling **************************

	cout<<setprecision(10)<<"Local energy:\t"<<(e_local/(double)num_cycles);	
	cout<<setprecision(10)<<"\tvariance:\t"<<(e_local_squared-e_local*e_local/num_cycles)/num_cycles;
	cout<<setprecision(10)<<"\tAcc.rate:\t"<<accepted/(double)(num_cycles)<<"\n";
	//cout<<"beta :\t\t"<<var_par[0]<<"\n";

	delete [] r_new;
	delete [] ones;
	delete [] ipd_upd;
}/*//endvimfold*/

void vmcsolver::calcLocalEnergy(double* e_local, double* e_local_squared, double* var_par)
{	/*//startvimfold*/
	
	int i;	
		
	double e_kinetic=0;
	//Laplacian: \nabla^2 \Psi_{\uparrow} + \nable^2 \Psi_{\downarrow}
	e_kinetic+=slater->lapl(r_old);
	//e_kinetic+=(cblas_ddot(dimension,r_old[1],1,r_old[1],1)+cblas_ddot(dimension,r_old[0],1,r_old[0],1) - dimension*num_part);
	//energy_s-=0.5*slater.lapl(r_old);
	e_kinetic+=ipd->jasLapl(var_par[0],r_old);
	//double beta = var_par[0];
	//e_kinetic +=  2*( 1.0 - beta*r_12_new ) / r_12_new / ( pow((1.0+r_12_new*beta),3) );
	//energy_j-=0.5*ipd.jasLapl(var_par[0], r_old);
	ipd->jasGrad(jas_grad, var_par[0],r_old);
	//temp = 1. / r_12_new / (1.+beta*r_12_new) / (1.+beta*r_12_new);
	//for (i=0;i<num_part;i++)
	//for (j=0;j<dimension;j++)
	//	jas_grad[i][j]=temp*(r_old[i][j]-r_old[(i+1)%2][j]);
	
	//TESTED: 2 part
	slater->grad(wf_grad,r_old);
	//for (i=0;i<num_part;i++)
	//for (j=0;j<dimension;j++)
	//	wf_grad[i][j]=-r_old[i][j];
	
	for (i=0;i<num_part;i++)
	{
		e_kinetic+=2.0*cblas_ddot(dimension,jas_grad[i],1,wf_grad[i],1);
		e_kinetic+=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
	}
	e_kinetic *= -0.5;

	double e_potential=0.0;
	//e - V_ext electrostatic pot
	//OPT: only |\vec r|^2 for one particle needs to be updated.
	for (i = 0; i < num_part; i++) 
	{
		e_potential+=cblas_ddot(dimension,r_old[i],1,r_old[i],1);
	}
	e_potential*=0.5;
	// e-e electrostatic interaction
	e_potential += ipd->sumInvlen();
	
	//for (i = 0; i < num_part-1; i++) 
	//{ 
	//	for (j = i+1; j < num_part; j++) 
	//	{
	//		r_12 = 0;  
	//		for (k = 0; k < dimension; k++) 
	//		{ 
	//			r_12 += (r_old[i][k]-r_old[j][k])*(r_old[i][k]-r_old[j][k]);
	//		}
	//		e_potential += 1/sqrt(r_12);         
	//		energy_ee +=1/sqrt(r_12); 
	//	}
	//}
	*e_local += e_potential+e_kinetic;
	*e_local_squared += pow((e_potential+e_kinetic),2); 
}/*//endvimfold*/

void vmcsolver::getNewPos(int active_part, double ideal_step, double* r_new, double* ipd_upd)
{/*//startvimfold*/
	int i,k;
	double temp;
	for (k = 0; k < dimension; k++)
	{ 
		r_new[k]=r_old[active_part][k]+ideal_step*(ran2(&idum)-0.5);
	}
	//New ipd : calculating new lengths betw particles.
	//only iterating over changed lengths.
	for (i=0; i<active_part; i++)
	{
		temp=0;
		for (k=0; k<dimension; k++)
		{
			temp+=pow((r_new[k]-r_old[i][k]),2);
		}
		//if (temp<1e-5) temp=1e-5;
		ipd_upd[i]=sqrt(temp);
	}
	for (i=active_part+1; i<num_part; i++)
	{
		temp=0;
		for (k=0; k<dimension; k++)
		{
			temp+=pow((r_new[k]-r_old[i][k]),2);
		}
		//if (temp<1e-5) temp=1e-5;
		ipd_upd[i-1]=sqrt(temp);
	}
}/*//endvimfold*/

vmcsolver::~vmcsolver()
{/*//startvimfold*/
	for (int i=0;i<num_part;i++)	
	{
		delete [] jas_grad[i];
		delete [] wf_grad[i];	
	}
	delete [] jas_grad;
	delete [] wf_grad;
	
	for (int i=0;i<num_part;i++)	
	{
		delete [] r_old[i];
	}
	delete[] r_old;
	
	//make destructor for these classes
	slater->clear();
	ipd->clear();
	delete slater;
	delete ipd;
}//End function /*//endvimfold*/
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
