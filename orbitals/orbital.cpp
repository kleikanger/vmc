/*
Implementation of class orbitals
   */

#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using std::cout;

orbital::orbital(int el,int am, bool su){

	energy_level=el;
	angular_momentum=am;
	spin_up=su;	

	if (energy_level<0) {
	cout<<"\nerror: energy_level must be a positive integer.\n";
	exit(1);
	} else {
	cout<<"particle initialized: energy_level: "<<energy_level
				<<", angular momentum: "<<angular_momentum
				<<", spin_up: "<<spin_up<<"\n"; 
	}
}

//using const for better opimalization
//Functions that returns different object properties
const int orbital::angularMomentum(){
	return angular_momentum;
}

const int orbital::energyLevel(){
	return energy_level;
}

//returns true if spin up
const bool orbital::spinUp(){
	return spin_up;
}

const double orbital::valueWF(double** r){
	return orbitalWavefunctions(r);
}

double orbital::orbitalWavefunctions(double** r){


	//Angular momentum can also be included.
	if (spin_up) {

		switch (energy_level) {
			case 1: return 0;
			case 2: return 1;
			case 3: return 3;
			
			
			default: 
					cout<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	
	} else {

		switch (energy_level) {
			case 1: return 5;
			case 2: return 6;
			case 3: return 7;
			
			
			default:
					cout<<"\n error in orbital::orbitalWavefunctions(): energy_level out of bounds\n"
						<<", energy_level= " <<energy_level<<"\n";
					exit(1);
		}
	}	
}

int main(){

//implementation example;

	orbital* B[4];

	//orbital* pointer;
	for (int i = 0; i<4; i++){
		//pointer = new orbital(0,i,false);
		//B[i] = pointer; 
		B[i] = new orbital(0,i,false);
	}
	
	for (int i = 0; i<4; i++){
		cout<<"\n"<<B[i]->angularMomentum()<<"\n";
	}

     for (int i = 0; i < 4; i++){
         delete B[i];
         B[i] = NULL;
     }

}
