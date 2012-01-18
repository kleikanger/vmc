
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
double  localEnergy(double**,double, slaterMatrix*,int, int);

#define H 0.001
#define H2 100000

int main(){
//startvimfold
	
	//Number of variational parameters
	int iNumber_of_variational_parameters = 1;
	//Variational parameters
	double variational_parameters[iNumber_of_variational_parameters];
	//Number of variations
	int number_of_variations[iNumber_of_variational_parameters];
	//Increasing the variational parameter with
	double variational_parameter_increase[iNumber_of_variational_parameters];
	
	//variational parameter. Start value.
	variational_parameters[0]=0.08; //alpha_1 minima -14.5 at 0.13
		//variational_parameter[1]=...;
	//Number of variations
	number_of_variations[0]=10;
	//increasing alpha with
	variational_parameter_increase[0]=0.01;
	
	//number of particles
	int iNumPart=4;	
	//number of particles each spin direction.
	int iCutoff=2;
	//length of random walker step
	double ideal_step=0.6;
	//Number of iterations in mc sampling
	int iNumber_of_iterations=500000;
	//Length of thermalization phase
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

	long idum;
  	idum= - (1 +  time(NULL));//time(NULL)*(myrank));
	
	//XXX: The variables below must be changed to an array containing the 
	//energies from all the runs. Cumulative energy isn't needed before MC
	//sampling is introduced.
	double cumulative_energy;
	double local_energy;
	//initialize slaterMatrix object
	slaterMatrix slater_Matrix(iNumPart,iCutoff,iNumber_of_variational_parameters);
	//variables to store calculated energies
	double wfnew, wfold;
			

	//Loop over variational parameters	
	for (int iTemp=0; iTemp < number_of_variations[0]; iTemp++){
	
		//''random'' startposition 
    	for (int i = 0; i < iNumPart; i++) { 
    		for (int j=0; j < 3; j++) {
 				partPos[i][j] = ideal_step*(ran2(&idum) -0.5);
      		}
    	}
	
		//Initializing startposition.
		slater_Matrix.updateVariationalParameters(variational_parameters);
		slater_Matrix.updateSlaterMatrix(partPos);
		wfold=slater_Matrix.waveFunction()*slater_Matrix.jastrow(partPos);
	
		//XXX: The variables below must be changed to an array containing the 
		//energies from all the runs. Cumulative energy isn't needed before MC
		//sampling is introduced.
		cumulative_energy=0;
		local_energy=0;

		//mc sampling loop
		int iAccepted_jumps=0;
		int iIteration_count=0;

		cout<<"\nnumber of iterations = "<<iNumber_of_iterations;

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
			slater_Matrix.updateSlaterMatrix(newPartPos);
			wfnew=slater_Matrix.waveFunction()*slater_Matrix.jastrow(newPartPos);
		
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
					local_energy+=localEnergy(partPos,wfnew,&slater_Matrix,3,iNumPart);
					iAccepted_jumps+=1;
				}
			}
		}//End while, (MC sampling)

	cout << "\n";
	cout << "alpha_1="<<variational_parameters[0]<<"\n";
	cout << "Energy = " << local_energy/(double)iAccepted_jumps<<"\n";
	cout << "Acceptance rate = " << iAccepted_jumps/(double)(iNumber_of_iterations-iThermalization)<<"\n";

	//Changing variational parameter.
	variational_parameters[0]+=variational_parameter_increase[0];
	}//end for loop over variational parameters


}//End main()
//endvimfold

//XXX for larger systems (larger slatermatrixes then 3x3, the derivatives 
//should be taken by using the cofactor elements. Then we dont have to to 
//calculate the 'determinant' more then once (losely speaking.))
//O(n^4)+O(n)*n instead of O(n^4)*n!
//double derivative()
//double doublederivative()

// Function to calculate the local energy 
double  localEnergy(double **r, double wfold, 
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
	  
	  slater_Matrix->updateSlaterMatrix(r_minus);
	  wfminus = (slater_Matrix->waveFunction())*(slater_Matrix->jastrow(r_minus));
	  
	  slater_Matrix->updateSlaterMatrix(r_plus); 
      wfplus  = (slater_Matrix->waveFunction())*(slater_Matrix->jastrow(r_plus)); 
      
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
