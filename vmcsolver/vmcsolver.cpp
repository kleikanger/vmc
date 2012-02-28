#include <cblas.h>
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
#include "vmcsolver.h"
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

//using namespace std;

//see folder ziggurat. more generators+automatic benchmark
#define RAN_NORM_SET RanNormalSetSeedZig32
#define RAN_NORM DRanNormalZig32
#define	RAN_UNI_SET RanSetSeed_MWC8222
#define RAN_UNI DRan_MWC8222

#ifndef OMG
#define OMG 1.0
#endif

//name and path of ofile
#ifndef OFNAMEB
#define OFNAMEB "test1.dat"
#endif
#ifndef OFPATHB
#define OFPATHB "/home/karleik/masterProgging/vmc/datafiles/zerotermalization.dat"
#endif
#ifndef OFNAMEB
#define OFNAMEB "test2.dat"
#endif

vmcsolver::vmcsolver(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par)
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

void vmcsolver::sample(int num_cycles, int thermalization, double* var_par, double delta_t, int myrank, double* result)
{/*//startvimfold*/

#if WRITEOFB
	//Only for rank 0 proc
	//if (rank=0) ?
	if (num_cycles>1.1e7)
	{
		//
		cerr<<"Warning: vmcsolver::sample(). Max size for num_cycles (=" 
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
	double e_kinetic, jaslR_new, jaslR_old, e_potential, wf_R, r_squared, temp;
	double e_local=0.0;
	double e_local_squared=0.0;
	double e_local_temp=0.0;
	double greens_f;
	int i,j,k,active_part;
	int accepted=0;
	//diffusion const D * delta_t = 0.5 * delta_t
	double dt_x_D = 0.5 * delta_t;
	//backup of jas_grad and sla_grad
	double jas_grad_bu[num_part][dimension];
	double sla_grad_bu[num_part][dimension];
	//updatevector interparticle distance.
	double* ipd_upd = new double[num_part-1];
	//temp for pointer swapping
	double** s_p_temp;
	//new positions
	double** r_new = new double*[num_part];
	for (i=0;i<num_part;i++) { r_new[i]=new double[dimension]; }

  	idum = (int)abs(time(NULL)*(myrank+1));
//	cout<<" "<<idum<< " "<<myrank<<" "<<time(NULL)<<"\n";
	//idum2 = (int)time(NULL);
	int cseed=1;//diff sequence for different seed
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
	slater->initSlaterMatrix(r_old);	
	slater->findInverse();
	ipd->init(r_old);
	
	//init q_force_old
	slater->grad(sla_grad,r_old);
	for (i=0; i<num_part; i++)
	{
		ipd->jasGrad(jas_grad, var_par[0],r_old,i);
	}
	for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
	{
		q_force_old[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
	}
	//init backup of jas_grad and slagrad.
    //not really necc since particles are thermalizing	
	for (i=0; i<num_part;i++) for (j=0; j<dimension;j++)
	{
		jas_grad_bu[i][j]=jas_grad[i][j];
	}
	for (i=0; i<num_part;i++) for (j=0; j<dimension;j++)
	{
		sla_grad_bu[i][j]=sla_grad[i][j];
	}
	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
		for (int loop_p=0;loop_p<num_part; loop_p++) //move all particles
		{
			//new position ++/*//startvimfold*/
			//possible to move two part. at a time, one in each SD?
			//active_part=loop_c%(num_part);
			active_part=loop_p;
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
			
			//upd gradient of jastrow
			//BEFORE UPD OF JASGRAD?
			ipd->jasGrad(jas_grad,var_par[0],r_new,active_part);
			
			//update inverse
			slater->update(r_new[active_part],active_part);
			//update gradients XXX OPT: ONLY ONE MATRIX AT A TIME. 2x. speedup
			slater->grad(sla_grad,r_new);//XXX necc to update inverse before taking gradients??)
			//*** void calcQF/*//endvimfold*/
			//calculate new qforce
			for (i=0; i<num_part; i++) for (j=0; j<dimension; j++)
			{
				//changing 2. will change the acceptance rate
				q_force_new[i][j] = 2.*(jas_grad[i][j]+sla_grad[i][j]);
			}
			//double calcGreensf
			//calculate greensfunc
			//only necc to calculate for 1 part??
			greens_f=0.0;/*//startvimfold*/
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
			//greens_f=exp(greens_f);/*//endvimfold*/


			//metropolis-hastings test
			if ( RAN_UNI() <= ( exp( greens_f+2.*(jaslR_new-jaslR_old) ) * wf_R * wf_R ) ) //XXX <! 
			{
				//accept updates/*//startvimfold*/
				//cblas_dcopy(dimension,&r_new[active_part][0],1,&r_old[active_part][0],1);//XXX loop faster
				for (i=0;i<dimension;i++)
				{	
					r_old[active_part][i]=r_new[active_part][i]; 
				}
				slater->accept(active_part);
				ipd->accept(active_part);
				//update quantumforce, swap pointers q_force_new<->q_force_new
				s_p_temp = q_force_old;
				q_force_old = q_force_new;
				q_force_new = s_p_temp;
				//backup gradients
				for (i=0;i<num_part;i++) for (j=0;j<dimension;j++)
				{ 
					jas_grad_bu[i][j] = jas_grad[i][j]; 
				}
				for (i=0;i<num_part;i++) for (j=0;j<dimension;j++)//only active matrix
				{ 
					sla_grad_bu[i][j] = sla_grad[i][j]; 
				}/*//endvimfold*/
				//update acceptancerate
				if (loop_c>thermalization) { accepted++; } 
			}
			else  
			{
				//reject updates/*//startvimfold*/
				slater->reject(active_part);
				ipd->reject(active_part);
				//revejt new position
				//cblas_dcopy(dimension,&r_old[active_part][0],1,&r_new[active_part][0],1);/*//endvimfold*/
				for (i=0;i<dimension;i++)
				{	
					r_new[active_part][i]=r_old[active_part][i];
				}
				//reset gradients if the energy are calculated after this cycle
				if (loop_p==(num_part-1)) 
				{
					for (i=0;i<num_part;i++) for (j=0;j<dimension;j++) 
					{ 
						jas_grad[i][j] = jas_grad_bu[i][j]; 
					}
					for (i=0;i<num_part;i++) for (j=0;j<dimension;j++) //XXX only active matrix necc
					{ 
						sla_grad[i][j] = sla_grad_bu[i][j]; 
					}
				}
			}
		}//All particles moved
		//restart for-loop if not thermalized
		if (loop_c<thermalization) { continue; }
		//if thermalization finished, start collecting data.
		e_local_temp = calcLocalEnergy(var_par);
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

	for (i=0;i<num_part;i++) { delete [] r_new[i]; }
	delete [] r_new;
	//delete [] ones;
	delete [] ipd_upd;
}/*//endvimfold*/

double vmcsolver::calcLocalEnergy(double* var_par)
{	/*//startvimfold*/
	int i;	
	
	double e_kinetic=0;
	//Laplacian of slatermatrices
	e_kinetic-=slater->lapl(r_old);
	//Laplacian of jastrow 
#if 1
	e_kinetic-=ipd->jasLapl(var_par[0],r_old); //x2??
#else/*//startvimfold*/
//mortens kode
double rij,a;
double beta = var_par[0];
int N=6;
int N2=3;
for ( i = 0 ; i < N-1; i ++) {
	for ( int j = i +1; j < N ; j ++) {
		if ( ( j < N2 && i < N2 ) || ( j >= N2 && i >= N2))
			a = 0.33333333;
		else
			a = 1.0;
		rij=0.0;
		for (int k = 0; k<dimension; k++)
			rij+=pow(r_old[i][k]-r_old[j][k],2);
		rij=sqrt(rij);
		e_kinetic -= a*(1.-rij*beta) / ( rij *pow(1.+ beta * rij , 3 ) );
	}
}/*//endvimfold*/
#endif
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
	//OPT: only |\vec r|^2 for one particle needs to be updated.
	for (i = 0; i < num_part; i++) 
	{
		e_potential+=cblas_ddot(dimension,r_old[i],1,r_old[i],1);
	}
	e_potential*=0.5*OMG*OMG;
	// e-e electrostatic interaction
	e_potential += ipd->sumInvlen();
	
	return e_potential+e_kinetic;
}/*//endvimfold*/

void vmcsolver::getNewPos(int active_part, double** r_new, double* ipd_upd, 
		double dt_x_D, double delta_t)
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
