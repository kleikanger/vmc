#include <cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <mpi.h>

//Ref: Jurgen A. Doornik 2005.
#include "../ziggurat/zigrandom.h"
#include "../ziggurat/zignor.h"
#include "vmcsolver.h"
//#include "../QDslater/slaterMatrix.h"
//#include "../ipdist/ipdist.h"

using std::cout;
using std::setprecision;

//see folder ziggurat. more generators+automatic benchmark
#define RAN_NORM_SET RanNormalSetSeedZig32
#define RAN_NORM DRanNormalZig32
#define	RAN_UNI_SET RanSetSeed_MWC8222
#define RAN_UNI DRan_MWC8222

#ifndef OMG
#define OMG 1.
#endif
#define OMG2 OMG*OMG

int main(int argc, char** argv)
{/*//startvimfold*/
	
	// ************ MPI INIT **************
	int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	// ************************************
	
	//Number of variational parameters
	int num_of_var_par=2;	
	//Variational param
	double  var_par[num_of_var_par];
	//Increase of variational parameters
	double var_par_inc[num_of_var_par];
	//Number of variations in each var. par.
	int var_par_cyc[num_of_var_par];
	
	int num_cycles=3.e6;
	int thermalization=0.3*num_cycles;//num_cycles*.5;
	int num_part=6;
	int spin_up_cutoff=3;
	int dimension=2;

	double delta_t = 0.05;

	//double ideal_step=.5;
	//double omega = 1.0;

	//beta (jastrow)	
	var_par[0]=0.55;
	var_par_inc[0]=0.00;
	var_par_cyc[0]=1;
	//alpha (orbitals)
	var_par[1]=.928;//0.92;
	var_par_inc[1]=0.0;
	var_par_cyc[1]=1;
	
	//construct vmcsolver object
	vmcsolver solver(num_part, spin_up_cutoff, dimension, num_of_var_par);	

	if (myrank==0)
	{
		//output to screen
		cout<<"\n\n";
		cout<<"------------------------------------------------------\n";
		cout<<"\n";
		cout<<"num_cycles:    \t\t"<<num_cycles<<"\n";
		cout<<"thermalization:\t\t"<<thermalization<<"\n";
		cout<<"delta_t:       \t\t"<<delta_t<<"\n";
		cout<<"num_part:      \t\t"<<num_part<<"\n\n";
		cout<<"------------------------------------------------------\n\n";
	}	

	//loop over variational parameters
	//start sampling
	for (int loop_v0=0; loop_v0<var_par_cyc[0]; loop_v0++)
	{
		var_par[0]+=var_par_inc[0];
		for (int loop_v1=0; loop_v1<var_par_cyc[1]; loop_v1++)
		{
			var_par[1]+=var_par_inc[1];
			solver.sample(num_cycles, thermalization, var_par, delta_t, myrank);
		}
	}

    MPI_Finalize();
    return 0;

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
	for (int i=0; i<num_part; i++)
		r_old[i] = new double[dimension];

	//handling of the slater matrix
	slater = new slaterMatrix(num_part,spin_up_cutoff,num_of_var_par,dimension);
	//handling of the interparticle distances
	ipd = new ipdist(num_part,dimension,spin_up_cutoff);
	
	//Gradients of jastrow and laplacian
	jas_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) jas_grad[i]=new double[dimension];
	sla_grad=new double*[num_part];
	for (int i=0;i<num_part;i++) sla_grad[i]=new double[dimension];
	//quantum force (for importance sampling)
	q_force_new = new double*[num_part];
	for (int i=0;i<num_part;i++) q_force_new[i]=new double[dimension];
	q_force_old = new double*[num_part];
	for (int i=0;i<num_part;i++) q_force_old[i]=new double[dimension];
}/*//endvimfold*/

void vmcsolver::sample(int num_cycles, int thermalization, double* var_par, double delta_t, int myrank)
{/*//startvimfold*/
	double e_kinetic, jaslR_new, jaslR_old, e_potential, wf_R, r_squared, temp;
	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;
	double greens_f;
	int i,j,k,active_part;
	int accepted=0;
	//diffusion const D * delta_t = 0.5 * delta_t
	double dt_x_D = 0.5 * delta_t;
	
	//updatevector interparticle distance.
	double* ipd_upd = new double[num_part-1];
	//new positions
	double** r_new = new double*[num_part];
	for (i=0;i<num_part;i++) { r_new[i]=new double[dimension]; }

  	idum = (int)time(NULL)*(myrank+1);
	//idum2 = (int)time(NULL);
	double cseed=5;//diff sequence for different seed
	//init ran number generator. only necc once, move to constructor?
	RAN_NORM_SET(&idum,cseed);
	RAN_UNI_SET(&idum,cseed);
	
	double init_sigma = 3.0;	
	//''random'' startposition 
	for (int i = 0; i < num_part; i++) { 
		for (int j=0; j < dimension; j++) {
			r_old[i][j] = init_sigma*(RAN_UNI()-0.5);
			r_new[i][j] = r_old[i][j];
		}
	}
	//update alpha in orbitals
	slater->setVarPar(var_par[1]);
	//initialize slatermatrix and ipd-matrix
	slater->updateVariationalParameters(var_par);
	slater->initSlaterMatrix(r_old);	
	slater->findInverse();
	ipd->init(r_old);
	
	
	//init q_force_old
	slater->grad(sla_grad,r_old);
	for (i=0; i<num_part; i++)
	{
		ipd->jasGrad(jas_grad, var_par[0],r_old,i);
	}
	for (i=0; i<num_part; i++)
	{
		for (j=0; j<dimension; j++)
			q_force_old[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
	}

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
	//	for (int loop_p=0;loop_p<num_part; loop_p++) //move all particles
	//	{
			//Moving one particle at a time. (converges slowly?)
			//possible to try to move all part. before collecting energy. 
			//or move two part. at a time, one in each SD?
			active_part=loop_c%(num_part);
			//active_part=loop_p;
			//get new r_i (r_new) and r_{ij} (ipd_upd) for active particle
			getNewPos(active_part, r_new, ipd_upd, dt_x_D, delta_t);
			//wave func ratio
			wf_R=slater->waveFunction(r_new[active_part],active_part);
			//log of old jastrow ratio  		
			jaslR_old = ipd->logJasR(active_part,var_par[0]);
			//update interparticle distances
			ipd->update(ipd_upd,active_part);
			//log of new jastrow ratio  		
			jaslR_new = ipd->logJasR(active_part,var_par[0]);
			//update inverse and slatermatrix (neccesary to find gradient??)
			//if not, wait, use and update() instead of accept().
			slater->update(r_new[active_part],active_part);
			//update gradients OPT: not necc to recalc for all part..
			slater->grad(sla_grad,r_new);//only necc to update one part (XXX necc to update inderse before taking gradients??)
			ipd->jasGrad(jas_grad,var_par[0],r_new,active_part);
			//*** void calcQF
			//calculate new qforce
			for (i=0; i<num_part; i++)
			{
				for (j=0; j<dimension; j++)
				{
					q_force_new[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
				}
			} 
			//double calcGreensf
			//calculate greensfunc
			greens_f=0.0;
			for (i=0; i<num_part; i++)
			for (j=0; j<dimension; j++)
			{
				//temp1=(q_force_old[i][j]+q_force_new[i][j]);
				//temp2=(q_force_old[i][j]-q_force_new[i][j]);
				if (i!=active_part)
				{
					greens_f+= 0.25 * dt_x_D *
						(q_force_old[i][j] + q_force_new[i][j]) *
						(q_force_old[i][j] - q_force_new[i][j]);
				} 
				else
				{
					//i=active_part
					greens_f+=0.5 * (q_force_old[i][j] + q_force_new[i][j]) *
					   ( 0.5 * dt_x_D * (q_force_old[i][j] - q_force_new[i][j]) +
						 r_old[i][j] - r_new[i][j] );
				}	
			}
			greens_f=exp(greens_f);
			//calculate jastrow ratio
			double jas_R = exp(jaslR_new-jaslR_old);
			//metropolis-hastings test
			if ( RAN_UNI() <= (greens_f * jas_R * jas_R * wf_R * wf_R ) ) 
			{
				//update positions (XXX loop faster)	
				cblas_dcopy(dimension,&r_new[active_part][0],1,&r_old[active_part][0],1);
				//accept updates
				slater->accept(active_part);
				ipd->accept(active_part);
				//update quantumforce, swap pointers q_force_new<->q_force_new
				double** s_p_temp;
				s_p_temp = q_force_old;
				q_force_old = q_force_new;
				q_force_new = s_p_temp;
				//restart for-loop in not thermalized
				if (loop_c<thermalization) { continue; } 
					//update acceptancerate
					accepted++;
					//update local energy
					e_local_temp=calcLocalEnergy(var_par);
			}
			else  
			{
			//reject updates
			slater->reject(active_part);
			ipd->reject(active_part);
			//update positions (XXX loop faster)	
			cblas_dcopy(dimension,&r_old[active_part][0],1,&r_new[active_part][0],1);
			//restart for-loop if not thermalized
			if (loop_c<thermalization) { continue; } 
		}
	//}//All particles moved

		//ADD LOOP (loop=0;loop<num_part;...)
		//Collect energy here (but acceptancerate as before)

		//if thermalization finished, start collecting data.
		//add local energy to variables e_local and e_local_squared.
		//e_local_temp = calcLocalEnergy(var_par);
		e_local += e_local_temp;
		e_local_squared += e_local_temp*e_local_temp; 	
	}
	e_local_temp=calcLocalEnergy(var_par);
	
	//************************** END OF MC sampling **************************

	cout<<"alpha: "<<var_par[1];
	cout<<" beta: "<<var_par[0]<<" ";
	cout<<setprecision(10)<<"Local energy:  "<<(e_local/(double)num_cycles);	
	cout<<setprecision(5)<<"\tvariance:  "<<(e_local_squared-e_local*e_local/num_cycles)/num_cycles;
	cout<<setprecision(5)<<"\tAcc.rate:  "<<accepted/(double)(num_cycles)<<"\n";
	//cout<<"beta :\t\t"<<var_par[0]<<"\n";

	for (i=0;i<num_part;i++) { delete [] r_new[i]; }
	delete [] r_new;
	//delete [] ones;
	delete [] ipd_upd;
}/*//endvimfold*/

double vmcsolver::calcLocalEnergy(double* var_par)
{	/*//startvimfold*/
	
	int i;	
		
	double e_kinetic=0;
	//Laplacian of slatermatrices: \nabla^2 \Psi_{\uparrow} + \nable^2 \Psi_{\downarrow}
	e_kinetic+=slater->lapl(r_old);
	//Laplacian of : 
	e_kinetic+=ipd->jasLapl(var_par[0],r_old);
	
	for (i=0;i<num_part;i++)
	{
		e_kinetic+=2.0*cblas_ddot(dimension,jas_grad[i],1,sla_grad[i],1);
		e_kinetic+=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
	}
	e_kinetic *= -0.5;

	double e_potential=0.0;
	//e - V_ext electrostatic pot
	//OPT: only |\vec r|^2 for one particle needs to be updated.
	for (i = 0; i < num_part; i++) 
	{
		e_potential+=OMG2*cblas_ddot(dimension,r_old[i],1,r_old[i],1);
	}
	e_potential*=0.5;
	// e-e electrostatic interaction
	e_potential += ipd->sumInvlen();
	//update local energy and local energy squared
	return e_potential+e_kinetic;
	//*e_local += e_potential+e_kinetic;
	//*e_local_squared += pow((e_potential+e_kinetic),2); 
}/*//endvimfold*/

void vmcsolver::getNewPos(int active_part, double** r_new, 
		double* ipd_upd, double dt_x_D, double delta_t)
{/*//startvimfold*/
	
	int i,k;
	double temp;
	for (k = 0; k < dimension; k++)
	{ 
		//r_new[k]=r_old[active_part][k]+ideal_step*(ran2(&idum)-.5);//gauDev(&idum,1.,0.);//ideal_step*(ran2(&idum)-0.5);
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

vmcsolver::~vmcsolver()
{/*//startvimfold*/
	for (int i=0;i<num_part;i++)	
	{
		delete [] jas_grad[i];
		delete [] sla_grad[i];	
		delete [] q_force_new[i];
		delete [] q_force_old[i];
	}
	delete [] jas_grad;
	delete [] sla_grad;
	delete [] q_force_old;
	delete [] q_force_new;
	
	for (int i=0;i<num_part;i++)	
	{
		delete [] r_old[i];
	}
	delete[] r_old;
	
	//make destructor for these classes
	ipd->clear();
	delete slater;
	delete ipd;
}//End function /*//endvimfold*/
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
