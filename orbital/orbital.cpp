/*
Implementation of class orbitals
Number of dimensions should be included.
 	*/

#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cblas.h>

using std::cout;
using std::cerr;

#define h 0.001
#define h2 100000

//double psi_10(double);
//double psi_20(double);
//double psi_30(double);

inline double psi_10(double r)
{
	return exp(-4.0 * r );
}
inline double psi_20(double r)
{
	return (2.0-4.0 * r) * exp(-2.0 * r ); 	
}
inline double psi_30(double r)
{
	return 1; 
}
//Constructors
orbital::orbital(){}
orbital::orbital(int el,int am, bool su){
//startvimfold
	energy_level=el;
	angular_momentum=am;
	spin_up=su;	

	if (energy_level<0) {
	cout<<"\nerror: energy_level must be a positive integer.\n";
	exit(1);
	} else if (energy_level<angular_momentum) {
	cout<<"\nerror: angular_momentum smaller then energy_level.\n";
	exit(1);
	} else {
	cout<<"particle initialized: energy_level: "<<energy_level
				<<", angular momentum: "<<angular_momentum
				<<", spin_up: "<<spin_up<<"\n"; 
	}
}
//endvimfold
void orbital::setValues(int el,int am, bool su){
//startvimfold
	energy_level=el;
	angular_momentum=am;
	spin_up=su;	

	if (energy_level<0) {
	cout<<"\nerror: energy_level must be a positive integer.\n";
	exit(1);
	} else if (energy_level<angular_momentum) {
	cout<<"\nerror: angular_momentum smaller then energy_level.\n";
	exit(1);
	} else {
	cout<<"particle initialized: energy_level: "<<energy_level
				<<", angular momentum: "<<angular_momentum
				<<", spin_up: "<<spin_up<<"\n"; 
	}
}
//endvimfold
//using const for better opimalization
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
//returns true if spin up
const bool orbital::spinUp(){
//startvimfold
	return spin_up;
}
//endvimfold
double orbital::valueWF(double* dR){
//startvimfold
	
	//Constants found in HF simulation. More decimals?
	double a=0.9955496248;
	double b=0.09423876105;
	//finding norm	
	double abs_r=0;
	abs_r=cblas_dnrm2(3,dR,1);
	abs_r=sqrt(abs_r);

	//Angular momentum can also be included.
	if (spin_up) {

		switch (energy_level) {
			case 1: return a*psi_10(abs_r)-b*psi_20(abs_r);
			case 2:	return -b*psi_10(abs_r)-a*psi_20(abs_r);	
			//XXX:Testing for larger systems
			case 3: return (2-4.0*abs_r)*exp(-4.0*abs_r);

			default: 
					cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	
	} else {

		switch (energy_level) {
			case 1:	return a*psi_10(abs_r)-b*psi_20(abs_r);
			case 2:	return -b*psi_10(abs_r)-a*psi_20(abs_r);	
			//XXX:Testing for larger systems
			case 3: return psi_10(abs_r)+psi_20(abs_r)+(2.0*abs_r*abs_r+2-4.0*abs_r)*exp(-4.0*abs_r);
			
			default:
					cerr<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}//end of orbitalWavefunctions::orbitalWavefunctions()
//endvimfold
double orbital::wFDeriv1(double* dR){
//startvimfold
	//Angular momentum can also be included.
	if (spin_up) {

		switch (energy_level) {
			case 1: return 1;
			case 2:	return 1;	
			case 3: return 1;
			case 4: return 1;
			default: 
					cerr<<"\n error in orbital::wFDeriv1(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	} else {
		switch (energy_level) {
			case 1:	return 1;	
			case 2:	return 1;	
			case 3: return 1;
			case 4: return 1;
			default:
					cerr<<"\n error in orbital::wFDeriv1(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}
//endvimfold
double orbital::wFDeriv2(double * dR){
//startvimfold
	//Angular momentum can also be included.
	if (spin_up) {

		switch (energy_level) {
			case 1: return 1;
			case 2:	return 1;	
			case 3: return 1;
			case 4: return 1;
			default: 
					cerr<<"\n error in orbital::wFderiv2(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	} else {
		switch (energy_level) {
			case 1:	return 1;	
			case 2:	return 1;	
			case 3: return 1;
			case 4: return 1;
			default:
					cerr<<"\n error in orbital::wFDeriv2(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}
//endvimfold
//Implementation example
//startvimfold
/*
int main(){

	orbital* B[4];
	orbital* C;
	
	C = new orbital[4];

	//orbital* pointer;
	for (int i = 0; i<4; i++){
		//pointer = new orbital(0,i,false);
		//B[i] = pointer; 
		B[i] = new orbital(i,2,false);	
		C[i].setValues(i,2,false);
	}
	
	double* a = new double[3];
	a[0]=a[1]=a[2]=0.2;
	
	for (int i = 0; i<2; i++){
		cout<<"\n"<<B[i]->angularMomentum()<<"\n";
		cout<<"\n"<<B[i]->valueWF(a)<<"\n";
		cout<<"\n"<<C[i].valueWF(a)<<"\n";
		
	}

     for (int i = 0; i < 4; i++){
         delete B[i];
         B[i] = NULL;//neccesary?
     }
	 delete a;
}
*/
//endvimfold
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold

