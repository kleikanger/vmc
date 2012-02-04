/*

Monte Carlo class. Containing data only about the updated slater matrices 
and the inverse of the slater matrices ant thd the jastrowfactor. Doing 
operations on the slater matrices. eg: Updating them, finding ratios between 
determinants, calculating i'th term in the laplacian or the gradient along 
some axis, ect.
 
 */

//should be written in such a manner that only the orbitals and their derivatives 
//in the orbital class needs to be changed to change to an other system..
	
#include "../orbital/orbital.h"

class slaterMatrix {
	private:
	//Class variables

		//Inverse matrixes (transpose for compputational reasons).
		double** inv_up_matr;
		double** inv_down_matr;
		//Slater matrixes. Two separate matrixes containing spin up and spin 
		//down particles. (iNumPart x iNumPart) matrix.
		double** spin_up_matr;
		double** spin_down_matr;
		//Number of particles of each spin orientation (of spin orientation up).
		int iCutoff;
		//initializes orbital* array. Each object one orbital with unique wf, 
		//quantumnumbers ect.
		orbital* orbital_;
		//Number of particles
		int iNumPart;
		//Number of dimensions
		int dim;

		// if true; cofactormatrix is calculated and used to find determinants.
		//bool bUse_cofactors;
		// keeping track of the determinants of inv_up_matr and inv_down_matr.
		// ensuring that they are not calculated more than once.
		//bool bDet_up;
		//bool bDet_down;
		
		//Remove? only used in waveFunction()
		double det_up;
		double det_down;
		
		////XXX Problem: Not initialized before iCutoff+1 first collecting cycles. XXX
		//double* grad_up;
		//double* grad_down;

		//Vector containing all the variational parameters
		double* variational_parameters;
		//length of this vector (remove?)
		int iNumber_of_variational_parameters;
	
	public:
	//Class variables

	//class methods
//slett();
		void print();
		/*
		Clear all malloc'ed variables.
		   */
		void clear();
		/*
		constructor
		Input: (number of particles, number of particles in each state (same number 
		assumed), Number of variational parameters, dimensions ).

		Allocating determinant matrixes for spin ud and down states, and the inverse 
		matrixes. Initializing orbital objects and class variables. The list of 
		orbitals (with quantumnumbers) can be changed here, wave-functions must be 
		changed in class orbital.
		 	*/
		slaterMatrix(int, int, int, int);
		/*
		Update all variational parameters. 
		Input: array of length iNumber_of_variational_parameters.
			*/
		void updateVariationalParameters(double*);
		/*
		Find inverse of both spin_..._matr. Store in both inv_..._matr.
		using lapack functions to invert inv_up_matr and inv_down_matr
		Storing the transpose of the inverse since many of the routines are 
		running through the last index
		   */
		void findInverse();
		/*
		Use O(N^2) routine to update the inverse matrix when only one
		particle i_upd is moved to coords r_new. 
		
		Updates the slater matrixes.
		*/
		void updateInverse(double* r_new, int i_upd);
		
		/*
		Initializing or updating slater matrix when all particles are moved to new coord's r.  
			*/
		void initSlaterMatrix(double** r);
		/*
		Updating slater matrix when particle i_upd are moved to coord r_new.  
		Find i
	  	
	 	MOVE TO PRIVATE LATER (when test not needed anymore).
			*/
		void updateSlaterMatrix(double* r_new, int i_upd);
		/*
		Input: (position vector of moved particle, number of moved particle).
		
		Using the inverse matrixes to determine the determinant. r_new is the position 
		of the new particle, i_upd the number of this particle.
		
		double waveFunction(double * dR, int iCofac_column,	int input_integer){
	 	input_integer=1
	 	input_integer=2
	 	input_integer=3
		 	*/
		double waveFunction(double*, int);
		/*	
		Calculates jastrowfactor. Explicit expression for jastrow in function.
		It is natural to have this function in this class since it contains the
		variational parameters.
			*/
		double jastrow(double**);
		/*
		Calculate one term of the gradient along axis 
		'axis' w.r.t particle number 'i_upd' over DET_old.
		Full gradient/DET_old = Sum over all i.
		EASY TO IMPLEMENT cblas_ddot() and temporary vector.
		 	*/		   
		double grad(double** dR, int axis);//, int i_upd);
		/*
		Calculate one term of the laplacian w.r.t. the 'i_upd''th 
		particle/DET_old (L_i f({r})). Full laplacian/DET/old = Sum_i L_i f({r}).
		EASY TO IMPLEMENT cblas_ddot() and temporary vector.
		   */
		double lapl(double** dR);

};
