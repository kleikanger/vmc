/*
Implementation of class orbitals
 	*/

#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mkl_cblas.h>

using std::cout;
using std::cerr;

//const. for calc. of derivatives.

#define H 0.001
#define ONE_OVER_H (double)1000
#define ONE_OVER_H2 (double)1000000

#ifndef ANALYTIC_D1
#define ANALYTIC_D1 true
#endif
#ifndef ANALYTIC_D2
#define ANALYTIC_D2 true
#endif
#ifndef OMG
#define OMG 1.
#endif
#ifndef SQOMG
#define SQOMG 1.
#endif

//Orbitals
inline double psi_10(double r_sqrd, double omg_alp)
{
	return exp(-0.5 * omg_alp * r_sqrd);
}
//hermite polynomials. (h0=1);
inline double h1(double x, double sq_omg_alp)
{
	return 2.*x*sq_omg_alp; 
}
inline double h2(double x, double omg_alp)
{
	return 4.*x*x*omg_alp - 2.; 
}
inline double h3(double x, double sqalpha)
{
	return 8.*pow(x*SQOMG*sqalpha,3) - 12.; 
}
inline double h4(double x, double alpha)
{
	return 16.*pow(x,4)*pow(OMG*alpha,2) - 
		48.*pow(x,2)*OMG*alpha + 12.;
}
//Constructors
orbital::orbital(){}
orbital::orbital(int el,int am, bool su, int di){
//startvimfold
	energy_level=el;
	angular_momentum=am;
	spin_up=su;	
	dim=di;

	if (energy_level<0) {
	cout<<"\nerror: energy_level must be a positive integer.\n";
	exit(1);
	} else if (energy_level*energy_level<angular_momentum*angular_momentum) {
	cout<<"\nerror: |angular_momentum| larger than energy_level.\n";
	exit(1);
	} else {
	cout<<"particle initialized: energy_level: "<<energy_level
				<<", angular momentum: "<<angular_momentum
				<<", spin_up: "<<spin_up<<"\n"; 
	}
}
//endvimfold
void orbital::setValues(int el,int am, bool su, int di){
//startvimfold
	dim=di;
	energy_level=el;
	angular_momentum=am;
	spin_up=su;
	//dim should be input param.
	dim=di;

	if (energy_level<0) {
	cout<<"\nerror: energy_level must be a positive integer.\n";
	exit(1);
	} else if (energy_level*energy_level<angular_momentum*angular_momentum) {
	cout<<"\nerror: |angular_momentum| larger than energy_level.\n";
	exit(1);
//	} else {
//	cout<<"particle initialized: energy_level: "<<energy_level
//				<<", angular momentum: "<<angular_momentum
//				<<", spin_up: "<<spin_up<<"\n"; 
	}
}
//endvimfold
void orbital::setAlpha(double alph){
	alpha = alph;/*//startvimfold*/
	omg_alp = OMG*alpha;
	sqalpha = sqrt(alph);
	sq_omg_alp = sqrt(alpha*OMG);
}/*//endvimfold*/
double orbital::valueWF(double* dR) const {
//startvimfold
//CALCULATE: value of the orbital in some point dR.

	double r_sqrd=0.0;
	r_sqrd=cblas_ddot(dim,dR,1,dR,1);

	//Angular momentum can also be included 3.
	//if (spin_up) 
	switch (energy_level) 
	{
		case 1: return psi_10(r_sqrd,omg_alp);
		case 2: return h1(dR[0],sq_omg_alp)*psi_10(r_sqrd,omg_alp);//n_x=1,n_y=0
		case 3: return h1(dR[1],sq_omg_alp)*psi_10(r_sqrd,omg_alp);//n_y=1,n_x=0
		case 4: return h2(dR[1],omg_alp)*psi_10(r_sqrd,omg_alp);//n_y=2,n_x=0
		case 5: return h1(dR[0],sq_omg_alp)*h1(dR[1],sq_omg_alp)
				*psi_10(r_sqrd,omg_alp);//n_y=1,n_x=1
		case 6: return h2(dR[0],omg_alp)*psi_10(r_sqrd,omg_alp);//n_y=0,n_x=2
		default: 
				cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
					<<", energy_level= " <<energy_level<<"\n";
				exit(1);
	}
}//end of orbitalWavefunctions::orbitalWavefunctions()
//endvimfold
double orbital::D1(double* dR, const int &axis) const {
//startvimfold
//CALCULATE: gradient along one axis. axis E {0,1,2,..}
#if !ANALYTIC_D1
	//does not work when multiplying with psi( spin_..._matr ) in slaterMatrix
	double f_min;
	double f_plus;
	double dR_temp[dim];
	double result=0.0;

	cblas_dcopy(dim,dR,1,dR_temp,1);

	dR_temp[axis]+=H;
	f_plus=valueWF(dR_temp);
	dR_temp[axis]-=2*H;
	f_min=valueWF(dR_temp);
	result=(f_plus-f_min)*ONE_OVER_H*0.5;

	return result;
#else
	//Analytic expressions for the gradients
	//double abs_r = cblas_dnrm2(dim,dR,1);

//	if (spin_up) {
	// returns gradient/phi. Multiply with psi for the correct result
	switch (energy_level) 
	{
		case 1: return -dR[axis]*alpha*OMG;
		case 2:	//n_x=1 n_y=0
				if (axis==0) return 
				( 2.*sq_omg_alp/h1(dR[0],sq_omg_alp)-dR[0]*omg_alp );
				else return 
				( -dR[1]*omg_alp );	
		case 3: //n_x=0 n_y=1
				if (axis==0) return 
				( -dR[0]*omg_alp );
				else return
				( 2.*sq_omg_alp/h1(dR[1],sq_omg_alp)-dR[1]*omg_alp );	
		case 4: //x0 y2
				if (axis==0) return 
				( -dR[0]*omg_alp );
				else return
				( 4.*sq_omg_alp*h1(dR[1],sq_omg_alp)/h2(dR[1],omg_alp)-dR[1]*omg_alp );	
		case 5: //x1 y1
				if (axis==0) return 
				( 2.*sq_omg_alp/h1(dR[0],sq_omg_alp)-dR[0]*omg_alp );	
				else return
				( 2.*sq_omg_alp/h1(dR[1],sq_omg_alp)-dR[1]*omg_alp );	
		case 6: //x2 y0
				if (axis==0) return 
				( 4.*sq_omg_alp*h1(dR[0],sq_omg_alp)/h2(dR[0],omg_alp)-dR[0]*omg_alp );	
				else return
				( -dR[1]*omg_alp );
		default: 
				cerr<<"\n error in orbital::D1(): energy_level out of bounds\n"
					<<", energy_level= " <<energy_level<<"\n";
				exit(1);
	}	
//	} else {
#endif
}
//endvimfold
double orbital::D2(double* dR) const {
//startvimfold
//calc: laplacian in some point dR.
//Does not work as long as we are multiplying with psi in slatermatrix

#if !ANALYTIC_D2
	double f_min;
	double f;
	double f_plus;
	double dR_temp[dim];
	double result=0.0;
	cblas_dcopy(dim,dR,1,dR_temp,1);

	f=valueWF(dR);
	for (int i=0; i<dim; i++)
	{
		dR_temp[i]+=H;
		f_plus=valueWF(dR_temp);
		dR_temp[i]-=2*H;
		f_min=valueWF(dR_temp);
		dR_temp[i]=dR[i];
		result+=(f_plus+f_min-2*f)*ONE_OVER_H2;
	}	
	return result;
#else
	//Analytic expressions for the laplacians
	// returns laplacian / psi. Multiply with psi for the correct result
	
	double r_sq = cblas_ddot(dim,dR,1,dR,1);
	
	//if (spin_up) {
	switch (energy_level) 
	{
		case 1: return omg_alp*(omg_alp*r_sq -2.);
		case 2:	return //n_x=1 n_y=0 
				omg_alp*( omg_alp*r_sq-2.
						-4.*sq_omg_alp*dR[0]/h1(dR[0],sq_omg_alp)) ;
		case 3: return //n_x=0 n_y=1 
				omg_alp*( omg_alp*r_sq-2.
						-4.*SQOMG*sq_omg_alp*dR[1]/h1(dR[1],sq_omg_alp));
		case 4: return //x0,y2
				omg_alp*(8./h2(dR[1],omg_alp)+omg_alp*r_sq
						-8.*sq_omg_alp*dR[1]*h1(dR[1],sq_omg_alp)/h2(dR[1],omg_alp)-2.);
		case 5: return //x1,y1
				omg_alp*( omg_alp*r_sq-2.
						-4.*sq_omg_alp*dR[1]/h1(dR[1],sq_omg_alp)
						-4.*sq_omg_alp*dR[0]/h1(dR[0],sq_omg_alp) );
		case 6: return //x2,y0
				omg_alp*(8./h2(dR[0],omg_alp)+omg_alp*r_sq
						-8.*sq_omg_alp*dR[0]*h1(dR[0],sq_omg_alp)/h2(dR[0],omg_alp)-2.);
		default: 
				cerr<<"\n error in orbital::D2(): energy_level out of bounds\n"
					<<", energy_level= " <<energy_level<<"\n";
				exit(1);
	}
	//} else {
#endif
}
//endvimfold
int orbital::angularMomentum() const {
//startvimfold
	return angular_momentum;
}
//endvimfold
int orbital::energyLevel() const {
//startvimfold
	return energy_level;
}
//endvimfold
bool orbital::spinUp() const {
//startvimfold
//returns true if spin up
	return spin_up;
}
//endvimfold
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold

