/*
Implementation of class orbitals
Number of dimensions should be included.
 	*/

#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cout;

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
//Functions that returns different object properties
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

const double orbital::valueWF(double* dR){
//startvimfold
	return orbitalWavefunctions(dR);
}
//endvimfold

 /***********************************************
  * 			 wavefunctions. 				*
  * 											*
  ***********************************************/ 
double orbital::orbitalWavefunctions(double* dR){
//startvimfold
	
	//Constants found in HF simulation. More decimals?
	double a=0.9955496248;
	double b=0.09423876105;
	
	//calculate distance to core
	double dAbsR=0;
	for (int i=0;i<3;i++){
		dAbsR+=dR[i]*dR[i];
	}
	dAbsR=sqrt(dAbsR);
	
	//n=1,l=0
	double psi_10 = 16.0*exp(-4.0*dAbsR);	//not normalized
	//n=2,l=0
	double psi_20 = 2.82842712474619*(2.0-4.0*dAbsR)*exp(-2.0*dAbsR); 	

	//Angular momentum can also be included.
	if (spin_up) {

		switch (energy_level) {
			case 1: return a*psi_10-b*psi_20;
			case 2: return -b*psi_10-a*psi_20;	
			
			default: 
					cout<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	
	} else {

		switch (energy_level) {
			case 0: return -a*psi_10+b*psi_20;	
			case 1: return b*psi_10+a*psi_20;	
			
			default:
					cout<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}//end of orbitalWavefunctions::orbitalWavefunctions()
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
