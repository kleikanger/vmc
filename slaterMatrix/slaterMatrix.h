/*
Variational Monte Carlo class. Containing data only about the states (updated slater matrixes).
Doing operations on the slater matrixes. eg: finding determinants, calculating energies, ect.
 */

//should be written in such a manner that only the orbitals needs to be changed to change to an other 
//system..
	
#include "../orbital/orbital.h"

class slaterMatrix {


	private:

	//Class variables
		//Slater matrixes. Two separate matrixes containing spin up and spin 
		//down particles. (iNumPart x iNumPart) matrix.
		double** spinUpMatrix;
		double** spinDownMatrix;
	public:	
		//cofactormatrixes. cofactors used to speed up calc. of determinants bigger than 2x2.
		//When one row is changed. Ex: when calculating first and second derivatives.
		double** cofactorMatrix_up;
		double** cofactorMatrix_down;
	private:
		//initializes orbital* array. Each object one orbital with unique wf, quantumnumbers ect.
		orbital* orbital_;

		//Number of particles
		int iNumPart;

		//Number of particles each spin orientation
		int iCutoff;

		//if true; cofactormatrix is calculated and used to find determinants.
		bool bUse_cofactors;

		//keeping track of the determinants of spinUpMatrix and spinDownMatrix.
		//ensuring that they are not calculated more than once.
		//bool bDet_up;
		//bool bDet_down;
		double det_up;
		double det_down;

		//Vector containing all the variational parameters
		double* variational_parameters;
		int iNumber_of_variational_parameters;
	
	//class methods
		//Calculates determinant of the slater matrixes	
		double determinant(double**,int);
		
		//Calculate the cofactors
		void updateCofactors(double**, double**, int);
		//XXX
		

	public:

	//Class variables

	//class methods
		
		//constructor
		//Input: (number of particles, number of particles in each state (same number assumed),
		//iNumber_of_variational_parameters ).
		//Allocating determinant matrixes for spin ud and down states, and the cofactormatrix.
		//Initializing orbital objects. The list of orbitals (with quantumnumbers) can be changed 
		//here, wave-functions must be changed in class orbital.
		slaterMatrix(int, int, int);		
		
		//Update all variational parameters. Input: array of length iNumber_of_variational_parameters.
		void updateVariationalParameters(double*);
		
		//Use O(N^2) routine to update the inverse matrix when only one
		//particle r_i is moved. ((Updates det_up, det_down.))
		void updateInverse(double* dR, int iUp);
		
		//Initializing slater matrix. Updating elements of slatermatrix and elements of comatrix.
		//If the variational parameters are to be updated, call updateVariationalParameters() before
		//calling this function. This function has to be called before the waveFunction and
		//the jastrowfactor can be found.
		//ONLY NEEDED ONCE if we use algo for updating the inverse.
		void updateSlaterMatrix(double**);
		
		//returns: bool bUse_cofactors. If true (if iCutoff>3), cofactor method is used to
		//calculate the determinants
		const bool useCofact();
		
		//To be removed(?)
		//input: positions of all particles + variational parameter.
		//Brute force calculation of the determinant.
		double waveFunction();
		//Overloading the func.
		//Input: (position vector of moved particle, number of moved particle).
		//Using the cofactors to determine the determinant.
		//calculate determinant without updating cofactorMatrix if only one r_i is changed.
		//NB: Not neccesary to update slatermatrix when calculating derivatives.
		// Changing func !
		//double slaterMatrix::waveFunction(double * dR, int iCofac_column, int input_integer){
		double waveFunction(double*, int);
		
		//Calculates jastrowfactor. Explicit expression for jastrow in function.
		double jastrow(double**);
		
		//Delete matrixes and objects
		//void calculate_qforce(double** r, double** q_force, double beta){
};
