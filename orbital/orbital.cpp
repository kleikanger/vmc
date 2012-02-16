/*
Implementation of class orbitals
 	*/

#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cblas.h>

using std::cout;
using std::cerr;

//const. for calc. of derivatives.

#define H 0.001
#define ONE_OVER_H 1000
#define ONE_OVER_H2 1000000

#ifndef ANALYTIC_D1
#define ANALYTIC_D1 true
#endif
#ifndef ANALYTIC_D2
#define ANALYTIC_D2 true
#endif
#ifndef OMG
#define OMG 1.
#endif
//squareroot of omega
#define SQOMG 1. 

//Orbitals
inline double psi_10(double r_sqrd, double alpha)
{
	return exp(-0.5 * alpha * OMG * r_sqrd);
}
//hermite polynomials. (h0=1);
inline double h1(double x)
{
	return 2.*x*SQOMG; 
}
inline double h2(double x)
{
	return 4.*pow(x*SQOMG,2) - 2.; 
}
inline double h3(double x)
{
	return 8.*pow(x*SQOMG,3) - 12.; 
}
inline double h4(double x)
{
	return 16.*pow(x*SQOMG,4) - 48.*pow(x*SQOMG,2) + 12.; 
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
	} else {
	cout<<"particle initialized: energy_level: "<<energy_level
				<<", angular momentum: "<<angular_momentum
				<<", spin_up: "<<spin_up<<"\n"; 
	}
}
//endvimfold
void orbital::setAlpha(double alph){
	alpha = alph;/*//startvimfold*/
}/*//endvimfold*/
const double orbital::valueWF(double* dR){
//startvimfold
//CALCULATE: value of the orbital in some point dR.

	// wf = h_(n_x)(dR[0])*h_(n_y)(dR[1])*psi_10;

	double r_sqrd=0.0;
	r_sqrd=cblas_ddot(dim,dR,1,dR,1);

	//Angular momentum can also be included 3.
	if (spin_up) 
	{
		switch (energy_level) 
		{
			case 1: return psi_10(r_sqrd,alpha);
			default: 
					cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	} 
	else 
	{
		switch (energy_level) 
		{
			case 1:	return psi_10(r_sqrd,alpha);
			default:
					cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}//end of orbitalWavefunctions::orbitalWavefunctions()
//endvimfold
const double orbital::D1(double* dR, int axis){
//startvimfold
//CALCULATE: gradient along one axis. axis E {0,1,2,..}
#if !ANALYTIC_D1
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
	double abs_r = cblas_dnrm2(dim,dR,1);
	if (spin_up) {

		switch (energy_level) {
			case 1: return -dR[axis]*alpha*OMG*valueWF(dR);
			case 2:	return 0;	
			case 3: return 0;
			case 4: return 0;
			default: 
					cerr<<"\n error in orbital::D1(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	} else {
		switch (energy_level) {
			case 1:	return -dR[axis]*alpha*OMG*valueWF(dR);	
			case 2:	return 0;	
			case 3: return 0;
			case 4: return 0;
			default:
					cerr<<"\n error in orbital::D1(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
#endif
}
//endvimfold
const double orbital::D2(double* dR){
//startvimfold
//calc: laplacian in some point dR.

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
	
	double r = cblas_dnrm2(dim,dR,1);
	
	if (spin_up) {
	
		switch (energy_level) {
			case 1: return alpha*OMG*(alpha*OMG*r*r -2.)*valueWF(dR);
			case 2:	return 0;	
			case 3: return 0;
			case 4: return 0;
			default: 
					cerr<<"\n error in orbital::D2(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	} else {
		switch (energy_level) {
			case 1:	return alpha*OMG*(alpha*OMG*r*r - 2.)*valueWF(dR);
			case 2:	return 0;	
			case 3: return 0;
			case 4: return 0;
			default:
					cerr<<"\n error in orbital::D2(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
#endif
}
//endvimfold
const int orbital::angularMomentum(){
//startvimfold
	return angular_momentum;
}
//endvimfold
const int orbital::energyLevel(){
//startvimfold
	return energy_level;
}
//endvimfold
const bool orbital::spinUp(){
//startvimfold
//returns true if spin up
	return spin_up;
}
//endvimfold
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold

