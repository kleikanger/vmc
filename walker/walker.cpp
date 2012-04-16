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
#include "../newmatrix/newmatrix.h"
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

//Defining random number generators RAN_NORM,RAN_NORM_SET,RAN_UNI,RAN_UNI_SET
//File automatically generated in python script and looks approx like this:
//
// #define RAN_NORM DRanNormalZig32
// #define RAN_NORM_SET RanNormalSetSeedZig32
// #define RAN_UNI DRan_MWC8222
// #define RAN_UNI_SET RanSetSeed_MWC8222
#include "../definitions/randomNumberGenerators.h"

walker::walker(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank)
{/*//startvimfold*/
	this->num_part=num_part;
	this->spin_up_cutoff=spin_up_cutoff;
	this->dimension=dimension;
	this->num_of_var_par=num_of_var_par;

	r_old = (double**)matrix(num_part,dimension,sizeof(double));
	r_new = (double**)matrix(num_part,dimension,sizeof(double));
	jas_grad = (double**)matrix(num_part,dimension,sizeof(double));
	sla_grad = (double**)matrix(num_part,dimension,sizeof(double));
	jas_grad_bu = (double**)matrix(num_part,dimension,sizeof(double));
	sla_grad_bu = (double**)matrix(num_part,dimension,sizeof(double));
	q_force_new = (double**)matrix(num_part,dimension,sizeof(double));
	q_force_old = (double**)matrix(num_part,dimension,sizeof(double));
	
	//handling of the slater matrix
	slater = new slaterMatrix(num_part,spin_up_cutoff,num_of_var_par,dimension);
	//handling of the interparticle distances
	ipd = new ipdist(num_part,dimension,spin_up_cutoff);
	//initialize random number generator	
	idum = (int)abs(time(NULL)*(myrank+1));
	int cseed=1;//diff sequence for different seed
	//init ran number generator. only necc once, move to constructor?
	RAN_NORM_SET(&idum,cseed);
	RAN_UNI_SET(&idum,cseed);

}/*//endvimfold*/

walker::~walker()
{/*//startvimfold*/
	free_matrix((void **) r_old);
	free_matrix((void **) r_new);
	free_matrix((void **) jas_grad);
	free_matrix((void **) jas_grad_bu);
	free_matrix((void **) sla_grad);
	free_matrix((void **) sla_grad_bu);
	free_matrix((void **) q_force_new);
	free_matrix((void **) q_force_old);
	delete slater;
	delete ipd;
}/*//endvimfold*/

void walker::initWalker(double* var_par, double delta_t)
{/*//startvimfold*/
	this-> delta_t=delta_t;
	//delta_t*diffusion constant
	dt_x_D=delta_t*0.5;
	sq_delta_t=sqrt(delta_t);
	omega=var_par[2];
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
	slater->setVarPar(var_par[1],var_par[2]);
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

void walker::initEmptyWalker(double* var_par, double delta_t)
{/*//startvimfold*/
	this-> delta_t=delta_t;
	//delta_t*diffusion constant
	dt_x_D=delta_t*0.5;
	sq_delta_t=sqrt(delta_t);
	omega=var_par[2];
	//update alpha in orbitals
	slater->setVarPar(var_par[1],var_par[2]);
	ipd->setBeta(var_par[0]);
}/*//endvimfold*/

bool walker::tryRandomStep(int active_part) 
{/*//startvimfold*/
	double jas_l_R, greens_f; 
	//wf_R is a class variable. See h-file for explanation.
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
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
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

double walker::calcLocalEnergy(double* var_par) const 
{	/*//startvimfold*/
	int i;	
	double e_kinetic=0;
	//Laplacian of slatermatrices
	e_kinetic-=slater->lapl(r_old);
	//Laplacian of jastrow 
	e_kinetic-=ipd->jasLapl(r_old);
	for (i=0;i<num_part;i++)
	{
		//part of the exp for the jastrow laplacian
		e_kinetic-=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
	}
	e_kinetic *= 0.5;	
	//cross term of the total laplacian
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
	e_potential*=0.5*omega*omega;
	// e-e electrostatic interaction
	e_potential += ipd->sumInvlen();

	return e_potential+e_kinetic;
}/*//endvimfold*/

void walker::getNewPos(int const &active_part, double* ipd_upd)
{/*//startvimfold*/
	
	int i,k;
	double temp;

	for (k = 0; k < dimension; k++)
	{ 
		r_new[active_part][k] += dt_x_D * q_force_old[active_part][k] 
					+ RAN_NORM() * sq_delta_t;
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

void walker::getRi(int i_w, double* x)
{/*//startvimfold*/
	for (int i=0;i<dimension;i++)
	{
		x[i]=r_old[i_w][i];
	}
}/*//endvimfold*/

bool walker::nodeCrossed()
{/*//startvimfold*/
	if (wf_R<=0) // < ??
		return true;
	else
		return false;
}/*//endvimfold*/

void walker::getVarParGrad(double* grad_var_par) const
{/*//startvimfold*/
	grad_var_par[0] = ipd->getdPdA(); //TODO Change name to getdPdBoveP
	//opdim: only for one slatermatrix the gradient needs to be updated!
	grad_var_par[1] = slater->getdPdAoveP(r_old); //Change name to getdPdAoveP
}/*//endvimfold*/

//EXPERIMENTAL
#if 1
void walker::getdELdvElem(double &de_l_db, double &de_l_da, double &dbJ_ove_J, double &daD_ove_D)
{

	double da_lapl_D_ove_D=0; 
	double db_lapl_J_ove_J=0;
	double jas_lapl;
	int i,j,k,l;
	double wf_R;
	
	double** da_grad_D_ove_D=(double**)matrix(num_part,dimension,sizeof(double));
	double** db_grad_J_ove_J=(double**)matrix(num_part,dimension,sizeof(double));
	//double** da_grad_D_temp = (double**)matrix(num_part,dimension,sizeof(double));

	//JASTROW PART
	dbJ_ove_J = ipd->getdPdA();
	//Jas lapl (already calculated when calculating local energy)
	jas_lapl=-ipd->jasLapl(r_old);
	for (i=0;i<num_part;i++)
	{
		jas_lapl=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
	}
	//jas_lapl*=-.5;
	//***
	//get derivative w.r.t. beta
	//sum derivative

	db_lapl_J_ove_J=ipd->getdbJasLapl(r_old);
	//db_lapl_J_ove_J+=dbJ_ove_J*jas_lapl;
	for (int i=0; i<num_part; i++) for(int j=0; j<dimension; j++)
		db_grad_J_ove_J[i][j] = 0.0;	
	ipd->getdbJasGrad(db_grad_J_ove_J,r_old); 

	//for (int i=0; i<num_part; i++) for(int j=0; j<dimension; j++)
	//	 db_grad_J_ove_J+= dbJ_ove_J*jas_grad[i][j];
	de_l_db = db_lapl_J_ove_J; //XXX check exp. in ipdist XXX *-2??
	for (i=0;i<num_part;i++)   //XXX check exp. in ipdist
		{
			de_l_db +=  2.*cblas_ddot(dimension, db_grad_J_ove_J[i], 1, sla_grad[i], 1);
		}
	//de_l_db += dbJ_ove_J*jas_lapl;
	de_l_db += dbJ_ove_J*slater->lapl(r_old);
	de_l_db *= -.5;

	//SLATER PART
	daD_ove_D = slater->getdPdAoveP(r_old);

	da_lapl_D_ove_D=0.;
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
		da_grad_D_ove_D[i][j]=0.0;
	for (i=0; i<num_part; i++)
	{
		//change the i'th row of the slatermatrix to dpsi/da
		//update inverse
		slater->updateInvA(r_old[i],i);
		//calculate gradients and laplacian..
		wf_R=slater->waveFunction(i);
		//wf_R_inv=1./wf_R_inv;
		da_lapl_D_ove_D+=slater->lapl(r_old)*wf_R; //*-.5?
		//ok?
		
		//slater->grad(da_grad_D_temp,r_old,0); //use normal grads?
		//slater->grad(da_grad_D_temp,r_old,spin_up_cutoff);
		slater->grad(sla_grad,r_old,i);
		
		//multiply all new elements by wf_R_inv
		for (j=0; j<num_part; j++) for (k=0; k<dimension; k++)
			da_grad_D_ove_D[j][k]+=sla_grad[j][k]*wf_R;

		//reset determinant, inverse and gradient
		slater->reject(i);
		//reset slater gradient only if active_part is the last particle
		//in the determinant (spin_up or spin_down)
		if ( !((i+1)%spin_up_cutoff ))	
		{
			if (i<spin_up_cutoff)
			{
				for (k=0;k<spin_up_cutoff;k++) for (l=0;l<dimension;l++)
				{ 
					sla_grad[k][l] = sla_grad_bu[k][l]; 
				}
			}
			else
			{
				for (k=spin_up_cutoff;k<num_part;k++) for (l=0;l<dimension;l++)
				{ 
					sla_grad[k][l] = sla_grad_bu[k][l]; 
				}
			}
		}
	}
	//calculate the derivative of the local energy
	de_l_da = da_lapl_D_ove_D; 
	for (i=0;i<num_part;i++)
		{
			de_l_da += 2. * cblas_ddot(dimension, da_grad_D_ove_D[i], 1, jas_grad[i], 1);
		}
	//de_l_da += daD_ove_D*slater->lapl(r_old);
	de_l_da += daD_ove_D*jas_lapl;
	de_l_da *= -0.5; 

	//free_matrix((void **) da_grad_D_temp);
	free_matrix((void **) da_grad_D_ove_D);
	free_matrix((void **) db_grad_J_ove_J);

}/*//endvimfold*/
#endif

void walker::wSetVarPar(double* var_par)
{
	int i,j;
	
	//update alpha in orbitals & jastrow
	slater->setVarPar(var_par[1],var_par[2]);
	ipd->setBeta(var_par[0]);
	
	
	//Update all matrices vith new variational parameters
	//upd gradient of jastrow (could be done by one single loop)
	for (i=0;i<num_part;i++)
		ipd->jasGrad(jas_grad,r_new,r_new, i);
	//update inverse determinant and determinant
	for (i=0;i<num_part;i++)
		slater->update(r_old[i],i);
	slater->accept(0);
	slater->accept(spin_up_cutoff);
	//update gradient of slatermatrix
	slater->grad(sla_grad,r_new,0);
	slater->grad(sla_grad,r_new,spin_up_cutoff);
	//calculate old/new??? qforce
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
	{
		q_force_old[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
	}
	//necc??
//	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
//	{
//		q_force_new[i][j] = q_force_old[i][j];
//	}
	for (i=0;i<num_part;i++) for (j=0;j<dimension;j++)
	{ 
		jas_grad_bu[i][j] = jas_grad[i][j]; 
	}
	for (i=0;i<num_part;i++) for (j=0;j<dimension;j++)
	//only active matr necc to update
	{
		sla_grad_bu[i][j] = sla_grad[i][j]; 
	}
}

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
