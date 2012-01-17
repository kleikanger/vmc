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
		
		//cofactormatrixes. cofactors used to speed up calc. of determinants bigger than 2x2.
		//When one row is changed. Ex: when calculating first and second derivatives.
		double** cofactorMatrix_up;
		double** cofactorMatrix_down;

		//initializes orbital* array. Each object one orbital with unique wf, quantumnumbers ect.
		orbital* orbital_;

		//Number of particles
		int iNumPart;

		//Number of particles each spin orientation
		int iCutoff;

		//if true; cofactormatrix is calculated and used to find determinants.
		bool bUse_cofactors;

		//keeping track of the determinants of spinUpMatrix and spinDownMatrix.
		//ensuring that they are not calculated more than one time.
		bool bDet_up;
		bool bDet_down;
		double det_up;
		double det_down;

	//class methods
		

	public:
		//XXX MOVE TO PRIVATE !!!
		//Calculates determinant of the slater matrixes	
		double determinant(double**,int);

		//Calculates jastrowfactor
		//move to private later
		double jastrow(double**, double);
		
		//Calculate the cofactors
		void updateCofactors(double**, double**, int);
		//XXX

	//Class variables

	//class methods
		
		//Initializing slater matrix
		void updateSlaterMatrix(double**, double);

		//constructor
		slaterMatrix(int, int);		
		
		//input: positions of all particles + variational parameter.
		//-Calculates wf in point dR (determinants) and multiplies with the jastrow factor with a variational parameter alpha.
		double waveFunction(double**, double);
		//Overloading the func.
		//Using the cofactors to determine the determinant
		//calculate determinant without updating cofactorMatrix if only one r_i is changed.
		double waveFunction(double*, double, int);
		
		//Delete matrixes and objects
		//void calculate_qforce(double** r, double** q_force, double beta){
};
