#include <mkl_cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>
//write to file
#include <fstream>
//Ref: Jurgen A. Doornik 2005.
#include "../ziggurat/zigrandom.h"
#include "../ziggurat/zignor.h"
#include "walker.h"
//#include "../QDslater/slaterMatrix.h"
//#include "../ipdist/ipdist.h"

using std::cout;
using std::cerr;
//using std::ofile;
using std::ofstream;
using std::setprecision;
using std::setiosflags;
using std::setw;
using std::ios;

#ifndef OMG
#define OMG 1.0
#endif

#define RAN_NORM_SET RanNormalSetSeedZig32
#define RAN_NORM DRanNormalZig32
#define	RAN_UNI_SET RanSetSeed_MWC8222
#define RAN_UNI DRan_MWC8222

walker::walker(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank)
{/*//startvimfold*/
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;

	//init pos.vecs.
	//double* r = new (nothrow) double[n];
	//if (r==0) { cout<<"...\n"; exit(1);}

	r_old = new double*[num_part];
	for (int i=0; i<num_part; i++)
		r_old[i] = new double[dimension];
	r_new = new double*[num_part];
	for (int i=0; i<num_part; i++)
		r_new[i] = new double[dimension];

	//handling of the slater matrix
	slater = new slaterMatrix(num_part,spin_up_cutoff,num_of_var_par,dimension);
	//handling of the interparticle distances
	ipd = new ipdist(num_part,dimension,spin_up_cutoff);
	
	//Gradients of jastrow and laplacian
	jas_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) jas_grad[i]=new double[dimension];
	sla_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) sla_grad[i]=new double[dimension];
	//old gradients
	jas_grad_bu=new double*[num_part];
	for (int i=0;i<num_part;i++) jas_grad_bu[i]=new double[dimension];
	sla_grad_bu=new double*[num_part];
	for (int i=0;i<num_part;i++) sla_grad_bu[i]=new double[dimension];
	//quantum force (for importance sampling)
	q_force_new = new double*[num_part];
	for (int i=0;i<num_part;i++) q_force_new[i]=new double[dimension];
	q_force_old = new double*[num_part];
	for (int i=0;i<num_part;i++) q_force_old[i]=new double[dimension];
  
	//initialize random number generator	
	idum = (int)abs(time(NULL)*(myrank+1));
	int cseed=1;//diff sequence for different seed
	//init ran number generator. only necc once, move to constructor?
	RAN_NORM_SET(&idum,cseed);
	RAN_UNI_SET(&idum,cseed);

}/*//endvimfold*/

walker::~walker()
{/*//startvimfold*/
	for (int i=0; i<num_part; i++)
	{
		delete [] r_old[i];
		delete [] r_new[i];
		delete [] jas_grad[i];
		delete [] jas_grad_bu[i];
		delete [] sla_grad[i];
		delete [] sla_grad_bu[i];
		delete [] q_force_new[i];
		delete [] q_force_old[i];
	}
	delete [] r_old;
	delete [] r_new;
	delete [] jas_grad;
	delete [] jas_grad_bu;
	delete [] sla_grad;
	delete [] sla_grad_bu;
	delete [] q_force_new;
	delete [] q_force_old;
	delete slater;
	delete ipd;
}/*//endvimfold*/

void walker::initWalker(double* var_par, double delta_t)
{/*//startvimfold*/
	this-> delta_t=delta_t;
	//delta_t*diffusion constant
	dt_x_D=delta_t*0.5;
	int i,j;

	double init_sigma = 3.0;	
	//''random'' startposition 
	for (i = 0; i < num_part; i++) { 
		for (j=0; j < dimension; j++) {
			r_old[i][j] = init_sigma*(RAN_UNI()-0.5);
			r_new[i][j] = r_old[i][j];
		}
	}

	//update alpha in orbitals
	slater->setVarPar(var_par[1]);
	ipd->setBeta(var_par[0]);
	//initialize slatermatrix and ipd-matrix
	slater->initSlaterMatrix(r_old);	
	slater->findInverse();
	ipd->init(r_old);
	
	//init gradients of slatermatrix, (both spin up and spin down determinant)
	slater->grad(sla_grad,r_old,0);
	slater->grad(sla_grad,r_old,spin_up_cutoff);
	//init gradient of jastrow
	for (i=0;i<num_part;i++)
		ipd->jasGrad(jas_grad, r_new, r_old, i);
	//init Q-force
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
	{
		q_force_old[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
	}
	//init backup of jas_grad and slagrad.
	for (i=0; i<num_part;i++) for (j=0; j<dimension;j++)
	{
		jas_grad_bu[i][j]=jas_grad[i][j];
	}
	for (i=0; i<num_part;i++) for (j=0; j<dimension;j++)
	{
		sla_grad_bu[i][j]=sla_grad[i][j];
	}
}/*//endvimfold*/

bool walker::tryRandomStep(int active_part) 
{/*//startvimfold*/
	double jas_l_R, wf_R, greens_f;
	double* ipd_upd = new double[num_part];
	int i,j;

	//update r_i (r_new) and r_{ij} (ipd_upd) for the active particle
	getNewPos(active_part, ipd_upd);

	//add log of old jastrow ratio  		
	jas_l_R = -ipd->logJasR(active_part);
	//update interparticle distances
	ipd->update(ipd_upd,active_part);
	//subtract log of new jastrow ratio  		
	jas_l_R += ipd->logJasR(active_part);
	jas_l_R *=2.;

	//upd gradient of jastrow
	ipd->jasGrad(jas_grad,r_new,r_old, active_part);

	//update inverse determinant and determinant
	slater->update(r_new[active_part],active_part);
	//wave func ratio
	wf_R=slater->waveFunction(active_part);
	//update gradient of slatermatrix
	slater->grad(sla_grad,r_new,active_part);

	//calculate new qforce
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
	{
		q_force_new[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
	}

	//calculate log of greens ratio
	greens_f=0.0;
	for (i=0; i<num_part; i++) 
	for (j=0; j<dimension; j++)
	{
		if (i!=active_part)
		{
			greens_f+= 0.25 * dt_x_D *
				(q_force_old[i][j] + q_force_new[i][j]) *
				(q_force_old[i][j] - q_force_new[i][j]);
		} 
		else
		{
			//i=active_part;
			greens_f+=0.5 * (q_force_old[i][j] + q_force_new[i][j]) *
				( 0.5 * dt_x_D * (q_force_old[i][j] - q_force_new[i][j]) +
				r_old[i][j] - r_new[i][j] );
		}	
	}

	delete [] ipd_upd;
	//metropolis test	
	return (RAN_UNI()<=(exp(greens_f+jas_l_R)*wf_R*wf_R ));
}/*//endvimfold*/

void walker::acceptStep(int active_part)
{/*//startvimfold*/
	//cblas_dcopy(dimension,&r_new[active_part][0],1,&r_old[active_part][0],1);
	//accept new position
	int i,j;
	for (i=0;i<dimension;i++)
	{	
		r_old[active_part][i]=r_new[active_part][i]; 
	}
	//accept move in slatermatrix and jastrow
	slater->accept(active_part);
	ipd->accept(active_part);
	//update quantumforce, swap pointers q_force_new<->q_force_new
	double** s_p_temp = q_force_old;
	q_force_old = q_force_new;
	q_force_new = s_p_temp;
	//backup gradients
	for (i=0;i<num_part;i++) for (j=0;j<dimension;j++)
	{ 
		jas_grad_bu[i][j] = jas_grad[i][j]; 
	}
	//only active matr necc to update
	if (active_part<spin_up_cutoff)
	{
		for (i=0;i<spin_up_cutoff;i++) for (j=0;j<dimension;j++)
		{ 
			sla_grad_bu[i][j] = sla_grad[i][j]; 
		}
	}
	else
	{
		for (i=spin_up_cutoff;i<num_part;i++) for (j=0;j<dimension;j++)
		{ 
			sla_grad_bu[i][j] = sla_grad[i][j]; 
		}
	}
}/*//endvimfold*/

void walker::rejectStep(int active_part)
{/*//startvimfold*/
	//cblas_dcopy(dimension,&r_old[active_part][0],1,&r_new[active_part][0],1);
	int i,j;
	for (i=0;i<dimension;i++)
	{	
		r_new[active_part][i]=r_old[active_part][i];
	}
	//reject updates in slatermatrix and jastrow
	slater->reject(active_part);
	ipd->reject(active_part);
	//(ONLY NECC FOR UP OR DOWN MATR)
	//only active matr necc to update
	//reset slater gradient only if active_part is the last particle
	//in the determinant (spin_up or spin_down)
	if ( !((active_part+1)%spin_up_cutoff ))	
	{
	if (active_part<spin_up_cutoff)
	{
		for (i=0;i<spin_up_cutoff;i++) for (j=0;j<dimension;j++)
		{ 
			sla_grad[i][j] = sla_grad_bu[i][j]; 
		}
	}
	else
	{
		for (i=spin_up_cutoff;i<num_part;i++) for (j=0;j<dimension;j++)
		{ 
			sla_grad[i][j] = sla_grad_bu[i][j]; 
		}
	}
	}
	//Must reset every cycle because of new jastgrad routine
	//if (active_part==(num_part-1)) 
	{
		for (i=0;i<num_part;i++) for (j=0;j<dimension;j++) 
		{ 
			jas_grad[i][j] = jas_grad_bu[i][j]; 
		}
	}
}/*//endvimfold*/

double const walker::calcLocalEnergy(double* var_par)
{	/*//startvimfold*/
	int i;	
	double e_kinetic=0;
	//Laplacian of slatermatrices
	e_kinetic-=slater->lapl(r_old);
	//Laplacian of jastrow 
	e_kinetic-=ipd->jasLapl(r_old);
	for (i=0;i<num_part;i++)
	{
		e_kinetic-=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
	}
	//cross term of the total laplacian
	e_kinetic *= 0.5;	
	for (i=0;i<num_part;i++)
	{
		e_kinetic-=cblas_ddot(dimension,jas_grad[i],1,sla_grad[i],1);
	}
	double e_potential=0.0;
	//e - V_ext electrostatic pot
	for (i = 0; i < num_part; i++) 
	{
		e_potential+=cblas_ddot(dimension,r_old[i],1,r_old[i],1);
	}
	e_potential*=0.5*OMG*OMG;
	// e-e electrostatic interaction
	e_potential += ipd->sumInvlen();

	return e_potential+e_kinetic;
}/*//endvimfold*/

void walker::getNewPos(int active_part, double* ipd_upd)
{/*//startvimfold*/
	
	int i,k;
	double temp;

	for (k = 0; k < dimension; k++)
	{ 
		r_new[active_part][k] += dt_x_D * q_force_old[active_part][k] 
					+ RAN_NORM() * sqrt(delta_t);
	}

	//New ipd : calculating new lengths betw particles.
	//only iterating over changed lengths.
	for (i=0; i<active_part; i++)
	{
		temp=0;
		for (k=0; k<dimension; k++)
		{
			temp+=pow((r_new[active_part][k]-r_old[i][k]),2);
		}
		//if (temp<1e-5) temp=1e-5;
		ipd_upd[i]=sqrt(temp);
	}
	for (i=active_part+1; i<num_part; i++)
	{
		temp=0;
		for (k=0; k<dimension; k++)
		{
			temp+=pow((r_new[active_part][k]-r_old[i][k]),2);
		}
		//if (temp<1e-5) temp=1e-5;
		ipd_upd[i-1]=sqrt(temp);
	}
}/*//endvimfold*/

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
