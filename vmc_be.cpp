
/*
	Implementation of class slaterMatrix
	Find explicit expression for jastrow gradient (evt. see slides from lectures).
   */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "lib/lib.h"
#include "./slaterMatrix/slaterMatrix.h"

using std::cout;

//Will include: MC sampling, statistics, write to outfile, paralellization.
//Calculates the local energy
double localEnergy(double**,double, slaterMatrix*,int, int);
//Functions to calculate first and second derivatives.
//XXX: Not tested! 
double firstDerivative(double **,int,int,slaterMatrix*);
//XXX: not tested for matrixes larger then (2x2)(cofactor method)
double grad2(double **,int,int,double,slaterMatrix*);
//Modified gradient. Calculate q-force.
double qForce(double**,double**,double,int,int,slaterMatrix*);

#define H 0.001
#define H1 1000
#define H2 1000000

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
	// generates gaussian distributed random numbers
	double gaussian_random(long*);
	
	//variational parameter. Start value.
	variational_parameters[0]=0.08; //alpha_1 minima -14.48 (exact 14.6...) at 0.13
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
	int iNumber_of_iterations=100000;
	//Length of thermalization phase
	int iThermalization=(int)iNumber_of_iterations*.3;
	
	//allocating positions of the particles assuming 3 dimensions.
	//First <iCutoff> indexes: pos. of spin up particles. 
	//Next  <iCutoff> indexes: pos. of spin down particles.
	double** partPos = new double*[iNumPart];
	double** newPartPos = new double*[iNumPart];
	double** q_force = new double*[iNumPart];
	double** new_q_force = new double*[iNumPart];

	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[3];
		newPartPos[i] = new double[3];
		q_force[i] = new double[3];
		new_q_force[i] = new double[3];
		//neccesary?
		for (int j=0;j<3;j++){
			partPos[i][j]=0.0;
			newPartPos[i][j]=0.0;
			q_force[i][j]=0.0;
			new_q_force[i][j]=0.0;
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

			//METROPOLIS ALGO
			int i,j;
			//Suggesting random step.
    		for (i = 0; i < iNumPart; i++) { 
      			for ( j=0; j < 3; j++) {
 					newPartPos[i][j] = partPos[i][j] + ideal_step*(ran2(&idum) -0.5);
      			}
    		}
#if 0 //importance sampling

			//Importance sampling
			//Changing q_force.			
			//XXX XXX XXX Not tested!

			qForce(partPos,q_force,wfold,iNumPart,3,&slater_Matrix);
			int i,j;
			double deltat=0.1;
			//Suggesting random step.
    		for (i = 0; i < iNumPart; i++) { 
      			for ( j=0; j < 3; j++) {
 					newPartPos[i][j] = partPos[i][j] + gaussian_random(&idum) + q_force[i][j]*deltat;
      			}
    		}
			//cout<<q_force[3][2]<<"\n";
#endif
			//Updating slatermatrix and finding wf (determinant*jastrow factor) in new position.
			slater_Matrix.updateSlaterMatrix(newPartPos);
			wfnew=slater_Matrix.waveFunction()*slater_Matrix.jastrow(newPartPos);
#if 0	
			//Important that this happens after update of slater matrix.
			qForce(newPartPos,new_q_force,wfnew,iNumPart,3,&slater_Matrix);

			// Ratio between Green’sfunctions
			double D=0.5;//diffusion constanti
			double greensratio = 0;
			for ( int ii= 0;ii < iNumPart ; ii++){
				for (int j=0;j < 3;j++){
					greensratio+=.5*(new_q_force[ii][j]+q_force[ii][j])
						*(.5*D* deltat * (q_force[ii][j]-new_q_force[ii][j])+partPos[ii][j]-newPartPos[ii][j]);
					}
			}
			greensratio=exp(greensratio);
			
			//cout<<new_q_force[3][2]<<"\n";
			//cout<<greensratio<<"\n";
#endif	
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

//For larger slatermatrixes then 3x3, the derivatives 
//are taken by using the cofactor elements. Then we dont have to to 
//calculate the determinant more then once (losely speaking.)
//O(n^4)+O(n)*n instead of O(n^4)*n!

//First derivative
//implement new method O(N^2) method in lecnotes.
//f'(x)=(f(x+h)-f(x-h))/2h + O(h^2) h^2=10^-6
double firstDerivative(double** r, int number_particles, int dimension, slaterMatrix* slater_Matrix){
//startvimfold
	int i,j;
	double wfplus,wfminus;
	double first_derivative=0;
	double **r_plus, **r_minus;

	r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
	r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
	for (i = 0; i < number_particles; i++) { 
		for ( j=0; j < dimension; j++) {
			r_plus[i][j] = r_minus[i][j] = r[i][j];
	   } 
 	} 
	
	//Using cofactor method?
	if (!slater_Matrix->useCofact()){
	 	for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
			  	r_plus[i][j] = r[i][j]+H;
			  	r_minus[i][j] = r[i][j]-H;
			  	slater_Matrix->updateSlaterMatrix(r_minus);
			  	wfminus = (slater_Matrix->waveFunction());
				wfminus*=(slater_Matrix->jastrow(r_minus));
			  	slater_Matrix->updateSlaterMatrix(r_plus); 
		  		wfplus  = (slater_Matrix->waveFunction());
				wfplus *= (slater_Matrix->jastrow(r_plus)); 
			  	first_derivative += (wfplus-wfminus);
			  	r_plus[i][j] = r[i][j];
			  	r_minus[i][j] = r[i][j];
			}
		}
	 } else {
		for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
				r_plus[i][j] = r[i][j]+H;
				r_minus[i][j] = r[i][j]-H;
				wfminus = slater_Matrix->waveFunction(r_minus[i],i);
				wfminus *= slater_Matrix->jastrow(r_minus);
				wfplus = slater_Matrix->waveFunction(r_plus[i],i);
				wfplus *= slater_Matrix->jastrow(r_plus); 
				first_derivative += (wfplus-wfminus);
				r_plus[i][j] = r[i][j];
				r_minus[i][j] = r[i][j];
			}
		}
	}
  	free_matrix((void **) r_plus); // free memory
  	free_matrix((void **) r_minus);
	return first_derivative*H1;
}//end function grad
//endvimfold
//modified gradient
double qForce(double** r, double** q_force, double wfold, 
		int number_particles, int dimension, slaterMatrix* slater_Matrix){
//startvimfold
	int i,j;
	double wfplus,wfminus;
	double first_derivative=0;
	double **r_plus, **r_minus;

	r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
	r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
	for (i = 0; i < number_particles; i++) { 
		for ( j=0; j < dimension; j++) {
			r_plus[i][j] = r_minus[i][j] = r[i][j];
	   } 
 	} 
	
	//Using cofactor method?
	if (!slater_Matrix->useCofact()){
	 	for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
			  	r_plus[i][j] = r[i][j]+H;
			  	r_minus[i][j] = r[i][j]-H;
			  	slater_Matrix->updateSlaterMatrix(r_minus);
			  	wfminus = (slater_Matrix->waveFunction());
				wfminus*=(slater_Matrix->jastrow(r_minus));
			  	slater_Matrix->updateSlaterMatrix(r_plus); 
		  		wfplus  = (slater_Matrix->waveFunction());
				wfplus *= (slater_Matrix->jastrow(r_plus)); 
			  	q_force[i][j] = (wfplus-wfminus)*H1*2/wfold;
			  	r_plus[i][j] = r[i][j];
			  	r_minus[i][j] = r[i][j];
			}
		}
	 } else {
		for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
				r_plus[i][j] = r[i][j]+H;
				r_minus[i][j] = r[i][j]-H;
				wfminus = slater_Matrix->waveFunction(r_minus[i],i);
				wfminus *= slater_Matrix->jastrow(r_minus);
				wfplus = slater_Matrix->waveFunction(r_plus[i],i);
				wfplus *= slater_Matrix->jastrow(r_plus); 
				q_force[i][j] = (wfplus-wfminus)*H1*2/wfold;
				r_plus[i][j] = r[i][j];
				r_minus[i][j] = r[i][j];
			}
		}
	}
  	free_matrix((void **) r_plus); // free memory
  	free_matrix((void **) r_minus);
}//end function grad
//endvimfold

//second derivative
//Implement new method O(n^2) method in lecturenotes.
//f''=((f(x+h)+f(x-h)-2*f(x))/h^2
double grad2(double ** r, int number_particles, int dimension, double wfold, slaterMatrix* slater_Matrix){
//startvimfold
	int i,j;
	double wfplus,wfminus;
	double second_derivative=0;
	double **r_plus, **r_minus;

	r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
	r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
	for (i = 0; i < number_particles; i++) { 
		for ( j=0; j < dimension; j++) {
			r_plus[i][j] = r_minus[i][j] = r[i][j];
	   } 
 	} 
	
	//Using cofactor method?
	if (!slater_Matrix->useCofact()){
	 	for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
			  	r_plus[i][j] = r[i][j]+H;
			  	r_minus[i][j] = r[i][j]-H;
			  	slater_Matrix->updateSlaterMatrix(r_minus);
			  	wfminus = (slater_Matrix->waveFunction());
				wfminus*=(slater_Matrix->jastrow(r_minus));
			  	slater_Matrix->updateSlaterMatrix(r_plus); 
		  		wfplus  = (slater_Matrix->waveFunction());
				wfplus *= (slater_Matrix->jastrow(r_plus)); 
			  	second_derivative += (wfminus+wfplus-2*wfold);
			  	r_plus[i][j] = r[i][j];
			  	r_minus[i][j] = r[i][j];
			}
		}
	 } else {
		for (i = 0; i < number_particles; i++) {
			for (j = 0; j < dimension; j++) { 
				r_plus[i][j] = r[i][j]+H;
				r_minus[i][j] = r[i][j]-H;
				wfminus = slater_Matrix->waveFunction(r_minus[i],i);
				wfminus *= slater_Matrix->jastrow(r_minus);
				wfplus = slater_Matrix->waveFunction(r_plus[i],i);
				wfplus *= slater_Matrix->jastrow(r_plus); 
				second_derivative += (wfminus+wfplus-2*wfold);
				r_plus[i][j] = r[i][j];
				r_minus[i][j] = r[i][j];
			}
		}
	}
  	free_matrix((void **) r_plus); // free memory
  	free_matrix((void **) r_minus);
	return second_derivative*H2;
}//end function grad2();
//endvimfold

// Function to calculate the local energy 
double  localEnergy(double **r, double wfold, slaterMatrix* slater_Matrix, int dimension, int number_particles){
/*//startvimfold*/
int charge = 4; //XXX XXX XXX XXX Include in start of main.

	int i, j , k;
	double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12, r_single_particle;
	// compute the kinetic energy  
	// include electron mass and hbar squared and divide by wave function 
	e_kinetic=grad2(r,number_particles,dimension,wfold,slater_Matrix);
	e_kinetic = -0.5*e_kinetic/wfold;
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
			//XXX XXX XXX XXX ADD POTENTIAL CLASS
     		e_potential += 1/sqrt(r_12);         
		}
  	}
  	e_local = e_potential+e_kinetic;
  	return e_local;
}//End function local energy  
//endvimfold

#if 0
//calculates the quantum force for the importance sampling
double calculate_qforce(double**, double**, double, double, int, int);

//Finds Quantum-force (del Psi) / Psi
double calculate_qforce(double** r, double** q_force, double beta, double wf, int dimension, int number_particles){
//startvimfold
	double** r_plus;
   	double** r_minus;
  	r_plus = (double **) matrix(number_particles, dimension, sizeof(double));
	r_minus = (double **) matrix(number_particles, dimension, sizeof(double));
 	
	double wf_plus, wf_minus;
	double alpha=2.75;

	//int number_particles=3; //SHOULD ONYL BE NECC. TO DECLARE AT ONE PLACE
	//int dimension=3;
	
	//Moving only one particle at a time:
	for (int i=0; i<number_particles; i++){
		for (int j=0; j<number_particles; j++){
			r_minus[i][j] = r[i][j];
			r_plus[i][j] = r[i][j];
		}
	}

	for (int i=0; i<number_particles; i++){
		for (int j=0; j<number_particles; j++){
			r_plus[i][j]=r[i][j]+h;
			r_minus[i][j]=r[i][j]-h;
			wf_plus = wave_function(r_plus, beta, dimension, number_particles);
			wf_minus =  wave_function(r_minus, beta, dimension, number_particles);
			q_force[i][j]= 2*(wf_plus-wf_minus)/(h*wf);
			r_minus[i][j] = r[i][j];
			r_plus[i][j] = r[i][j];
		} 
	}

	2*grad/wf;

}//end function calculate_qforce();
//endvimfold
#endif

//returns gaussian distributed random numbers XXX what const?
double gaussian_random(long* idum){
//startvimfold	
	//double pi_sqrt = 1.7724538509055160;//sqrt(3.1415926535897932);
	double random=ran2(idum)-0.5;
	
	// (general gaussian  exp(-x^2/(2*pi*s^2)) / (sqrt(pi)* s) ))
	if (random<0){
		return exp(-random*random)/1.8;
	} else {
		return -exp(-random*random)/1.8;
	}

}//end: function gaussian_random()
//endvimfold


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
