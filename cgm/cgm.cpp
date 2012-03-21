#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include "cgm.h"

using namespace std;

//Algorithms from numerical recipes. Written by M. Hjorth-Jensen, slightly modified 
void lnsrch(int n, Vector &xold, double fold, Vector &g, Vector &p, Vector &x,
		 double *f, double stpmax, int *check, cgm *pCgm);
void dfpmin(Vector &p, int n, double gtol, int *iter, double *fret, cgm* pCgm);

cgm::cgm(int num_cyclesARG, int thermalizationARG, double delta_tARG, 
		double omegaARG, int num_partARG, int spin_up_cutoffARG,
		int dimensionARG, int num_of_var_parARG, int myrank_ARG, 
		int nprocs_ARG)
{/*//startvimfold*/
	num_cycles = num_cyclesARG;
	thermalization = thermalizationARG;
	delta_t = delta_tARG;
	omega = omegaARG;
	num_part = num_partARG;
	spin_up_cutoff = spin_up_cutoffARG;
	num_of_var_par = num_of_var_parARG;
	dimension = dimensionARG;
	myrank = myrank_ARG;
	nprocs = nprocs_ARG;
	//todo gtol
	//TODO num_of_var_par=2;
	dE_array = new double[2];
}/*//endvimfold*/

cgm::~cgm()
{/*//startvimfold*/
	delete dE_array;
}/*//endvimfold*/

void cgm::optimizeVarPar(double* initial_var_par)
{/*//startvimfold*/
	int  iter;
   	double gtol, fret;
    double alpha;
	//reserve space in memory for vectors containing the variational
	//parameters
		
	if (myrank==0)
	{
		cout<<"\nRunning CGM - minimization\n"
			<<num_part<<" particles, "
			<<num_cycles<<" cycles/proc, "
			<<nprocs<<" procs\n";
	}

	int n=2; //TODO num_of_var_par-1
	Vector g(n), p(n);
   	gtol = 1.0e-5; //TODO Should be imput variable

	//now call dfmin and compute the minimum
    p(0) = initial_var_par[0];
	p(1) = initial_var_par[1];

    dfpmin(p, n, gtol, &iter, &fret, this);

	if (myrank==0)
	{
		cout <<setprecision(8)<< "Value of energy minimum = " << fret << endl;
		cout << "Number of iterations = " << iter << endl;
		cout << "Value of alpha at minimum = " << p(1) << endl;
		cout << "Value of beta at minimum = " << p(0) << endl;
	} 	
	//TODO return alpha, beta
}/*//endvimfold*/

double cgm::E_function(Vector &variational_parameters)
{/*//startvimfold*/
	double* var_par = new double[3];
	double* var_res_temp = new double[2];
	double energy=0.0;

	//should be class variable
	sampler sampler_(num_part, spin_up_cutoff, dimension, num_of_var_par, myrank);
	
	//TESTING (abs) must find smarter method here TODO TODO TODO TODO
	//
	// Maybe return some negative (+sign) gradient or 0? if variational_parameters()<0 
	//
	//start sampling in all processes
	var_par[0]=fabs(variational_parameters(0));
	var_par[1]=fabs(variational_parameters(1));
	var_par[2]=omega;
	sampler_.sample(num_cycles,thermalization,var_par,delta_t,var_res_temp);

	//collecting energy from all processes, broadcasting mean 
	MPI_Reduce( &var_res_temp[0], &energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	energy /= (double)nprocs;
	MPI_Bcast (&energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//TODO 2 = num_var_par-1
	for (int i=0; i<2; i++)
	{
		//collecting gradient from all processes, broadcasting mean 
		double temp=sampler_.getEnergyGrad(i);
		MPI_Reduce(&temp, &dE_array[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		dE_array[i]/=(double)nprocs;
		MPI_Bcast (&dE_array[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	
	delete [] var_par;
	delete [] var_res_temp;
	
	if (myrank==0)
	{	
		cout<<setprecision(8)
			<<"energy="<<energy<<", (dE/da,dE/db)=("
			<<dE_array[1]<<","<<dE_array[0]<<")"
			<<" alpha: "<<variational_parameters(1) 
			<<" beta: "<<variational_parameters(0) 
			<<"\n";
	}
	return energy;
}/*//endvimfold*/

//variational_parameters not necc here. evt: remove from call to function.
void cgm::dE_function(Vector &variational_parameters, Vector &gradient)
{/*//startvimfold*/
	for (int i=0; i<2; i++)
	{
		gradient(i) = dE_array[i];
	}
}/*//endvimfold*/

//REF: numerical recipes
//TODO remove vector class
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

void dfpmin(Vector &p, int n, double gtol, int *iter, double *fret, cgm* pCgm)
{/*//startvimfold*/
  int check,i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  Vector dg(n), g(n), hdg(n), pnew(n), xi(n);
  Matrix hessian(n,n);

  fp=pCgm->E_function(p);
  pCgm->dE_function(p,g);

  for (i = 0;i < n;i++) {
    for (j = 0; j< n;j++) hessian(i,j)=0.0;
    hessian(i,i)=1.0;
    xi(i) = -g(i);
    sum += p(i)*p(i);
  }
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,pCgm);
    fp = *fret;
    for (i = 0; i< n;i++) {
      xi(i)=pnew(i)-p(i);
      p(i)=pnew(i);
    }
    test=0.0;
    for (i = 0;i< n;i++) {
      temp=fabs(xi(i))/FMAX(fabs(p(i)),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i);
	pCgm->dE_function(p,g);
    test=0.0;
    den=FMAX(*fret,1.0);
    for (i=0;i<n;i++) {
      temp=fabs(g(i))*FMAX(fabs(p(i)),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i)-dg(i);
    for (i=0;i<n;i++) {
      hdg(i)=0.0;
      for (j=0;j<n;j++) hdg(i) += hessian(i,j)*dg(j);
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) {
      fac += dg(i)*xi(i);
      fae += dg(i)*hdg(i);
      sumdg += SQR(dg(i));
      sumxi += SQR(xi(i));
    }
    if (fac*fac > EPS*sumdg*sumxi) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=0;i<n;i++) dg(i)=fac*xi(i)-fad*hdg(i);
      for (i=0;i<n;i++) {
	for (j=0;j<n;j++) {
	  hessian(i,j) += fac*xi(i)*xi(j)
	    -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
	}
      }
    }
    for (i=0;i<n;i++) {
      xi(i)=0.0;
      for (j=0;j<n;j++) xi(i) -= hessian(i,j)*g(j);
    }
  }
  cout << "too many iterations in dfpmin" << endl;
}/*//endvimfold*/
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

#define ALF 1.0e-4
#define TOLX 1.0e-7
void lnsrch(int n, Vector &xold, double fold, Vector &g, Vector &p, Vector &x,
	    double *f, double stpmax, int *check, cgm *pCgm)
{/*//startvimfold*/
  int i;
  double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  *check=0;
  for (sum=0.0,i=0;i<n;i++) sum += p(i)*p(i);
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=0;i<n;i++) p(i) *= stpmax/sum;
  for (slope=0.0,i=0;i<n;i++)
    slope += g(i)*p(i);
  test=0.0;
  for (i=0;i<n;i++) {
    temp=fabs(p(i))/FMAX(fabs(xold(i)),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=0;i<n;i++) x(i)=xold(i)+alam*p(i);
	*f=pCgm->E_function(x);
    if (alam < alamin) {
      for (i=0;i<n;i++) x(i)=xold(i);
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) cout << "Roundoff problem in lnsrch." << endl;
	  else tmplam=(-b+sqrt(disc))/(3.0*a);
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}/*//endvimfold*/
#undef ALF
#undef TOLX

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
