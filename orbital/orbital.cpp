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

#define ANALYTIC_D2 true
#define ANALYTIC_D1 true
#define alpha .987

//Orbitals:
inline double psi_10(double r_sqrd)
{
	return exp( - 0.5 * alpha * r_sqrd );
}
inline double psi_20(double r)
{
	return 1.; 
}
inline double psi_30(double r)
{
	return 1.; 
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
double orbital::valueWF(double* dR){
//startvimfold
//CALCULATE: value of the orbital in some point dR.

	double r_sqrd=0.0;
	r_sqrd=cblas_ddot(dim,dR,1,dR,1);

	//Angular momentum can also be included 3.
	if (spin_up) 
	{
		switch (energy_level) 
		{
			case 1: return psi_10(r_sqrd);
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
			case 1:	return psi_10(r_sqrd);
			default:
					cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}//end of orbitalWavefunctions::orbitalWavefunctions()
//endvimfold
double orbital::D1(double* dR, int axis){
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
			case 1: return -dR[axis]*alpha*valueWF(dR);
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
			case 1:	return -dR[axis]*alpha*valueWF(dR);	
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
double orbital::D2(double* dR){
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
			case 1: return alpha*(alpha*r*r -2.)*valueWF(dR);
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
			case 1:	return alpha*(alpha*r*r - 2.)*valueWF(dR);
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

