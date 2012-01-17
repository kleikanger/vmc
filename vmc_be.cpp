
/*
	Implementation of class slaterMatrix

   */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "lib/lib.h"
#include "./slaterMatrix/slaterMatrix.h"

using std::cout;

//Will include: MC sampling, statistics, write to outfile, paralellization.
//Calculates the local energy
double  localEnergy(double**,double, double,slaterMatrix*,int, int);

#define H 0.001
#define H2 1000000

int main(){
//startvimfold

	//variational parameter. Start value.
	double alpha=1.0;
	//Number of variations
	double number_of_alpha_variations=10;
	//increasing alpha with
	double alpha_increase=0.05;
	//number of particles
	int iNumPart=4;	
	//number of particles for separate spins.
	int iCutoff=2;
	//length of random walker step
	double ideal_step=0.6;
	int iNumber_of_iterations=500000;
	int iThermalization=(int)iNumber_of_iterations*.3;

	//allocating positions of the particles assuming 3 dimensions.
	//First <iCutoff> indexes: pos. of spin up particles. 
	//Next  <iCutoff> indexes: pos. of spin down particles.
	double** partPos = new double*[iNumPart];
	double** newPartPos = new double*[iNumPart];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[3];
		newPartPos[i] = new double[3];
		//neccesary?
		for (int j=0;j<3;j++){
			partPos[i][j]=0.0;
			newPartPos[i][j]=0.0;
		}
	}
	
	//Ecumulative: will be needed when parallellizing code.
	double cumulative_energy;
	double local_energy;
	//initialize slaterMatrix object
	slaterMatrix slater_Matrix(iNumPart,iCutoff);
	//
	double wfnew, wfold;

	long idum;
  	idum= - (1 +  time(NULL));//time(NULL)*(myrank));
	
	for (int iTemp=0; iTemp < number_of_alpha_variations; iTemp++){
	
		cumulative_energy=0;
		local_energy=0;
	
		//''random'' startposition 
    	for (int i = 0; i < iNumPart; i++) { 
    		for (int j=0; j < 3; j++) {
 				partPos[i][j] = ideal_step*(ran2(&idum) -0.5);
    	  	}
    	}

		//Variational parameter alpha is used to variate the jastrow factor
		wfold=slater_Matrix.waveFunction(partPos,alpha);	
	
		//mc sampling loop
		int iAccepted_jumps=0;
		int iIteration_count=0;

		cout<<"\nnumber of iterations = "<<iNumber_of_iterations;

		//cout<<"\n\niterations (.), accepted jumps (+):\n";

		while (iIteration_count < iNumber_of_iterations){
			iIteration_count++;

			int i,j;

			//Suggesting random step.
    		for (i = 0; i < iNumPart; i++) { 
      			for ( j=0; j < 3; j++) {
 					newPartPos[i][j] = partPos[i][j] + ideal_step*(ran2(&idum) -0.5);
      			}
    		}

			//Updating slatermatrix and finding wf (determinant*jastrow factor) in new position.
			slater_Matrix.updateSlaterMatrix(newPartPos,alpha);
			wfnew=slater_Matrix.waveFunction(newPartPos,alpha);
		
			//Metropolis algorithm
			if (wfnew*wfnew/wfold/wfold >= ran2(&idum)){
				for (i = 0; i < iNumPart; i++) { 
      				for ( j=0; j < 3; j++) {
 						partPos[i][j] = newPartPos[i][j];
      				}
    			}
				//Updating parameters and energy.
				wfold=wfnew;
				//Updating energy if system is thermalized.
				if (iIteration_count>iThermalization){
					local_energy+=localEnergy(partPos,alpha,wfnew,&slater_Matrix,3,iNumPart);
					iAccepted_jumps+=1;
				}
			}
		}//End while,

	cout << "\n";
	cout << "Energy = " << local_energy/(double)iAccepted_jumps<<"\n";
	cout << "Acceptance rate = " << iAccepted_jumps/(double)(iNumber_of_iterations-iThermalization)<<"\n";

	alpha+=alpha_increase;
	}//end for loop over variational parameters


}//End main()
//endvimfold

// Function to calculate the local energy 
double  localEnergy(double **r,	double beta , double wfold, 
					slaterMatrix* slater_Matrix, int dimension, int number_particles){
/*//startvimfold*/
int charge = 4;

  int i, j , k;
  double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12, r_single_particle;
  double **r_plus, **r_minus;

  // allocate matrices which contain the position of the particles  
  // the function matrix is defined in the progam library 
  r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
  r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
  for (i = 0; i < number_particles; i++) { 
    for ( j=0; j < dimension; j++) {
      r_plus[i][j] = r_minus[i][j] = r[i][j];
    }
  }
  // compute the kinetic energy  
  e_kinetic = 0;
  for (i = 0; i < number_particles; i++) {
    for (j = 0; j < dimension; j++) { 
      r_plus[i][j] = r[i][j]+H;
      r_minus[i][j] = r[i][j]-H;
	  slater_Matrix->updateSlaterMatrix(r_minus,beta);
	  wfminus = slater_Matrix->waveFunction(r_minus,beta);
	  slater_Matrix->updateSlaterMatrix(r_plus,beta); 
      wfplus  = slater_Matrix->waveFunction(r_plus,beta); 
      e_kinetic -= (wfminus+wfplus-2*wfold);
      r_plus[i][j] = r[i][j];
      r_minus[i][j] = r[i][j];
    }
  }
  // include electron mass and hbar squared and divide by wave function 
  e_kinetic = 0.5*H2*e_kinetic/wfold;
  // compute the potential energy 
  e_potential = 0;
  // contribution from electron-proton potential  
  for (i = 0; i < number_particles; i++) { 
    r_single_particle = 0;
    for (j = 0; j < dimension; j++) { 
      r_single_particle += r[i][j]*r[i][j];
    }
    e_potential -= charge/sqrt(r_single_particle);
  }
  // contribution from electron-electron potential  
  for (i = 0; i < number_particles-1; i++) { 
    for (j = i+1; j < number_particles; j++) {
      r_12 = 0;  
      for (k = 0; k < dimension; k++) { 
	r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
      }
      e_potential += 1/sqrt(r_12);          
    }
  }
  free_matrix((void **) r_plus); // free memory
  free_matrix((void **) r_minus);
  e_local = e_potential+e_kinetic;
  return e_local;
}//End function local energy  
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
