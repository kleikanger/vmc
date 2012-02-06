/*

Slater matrix class. Containing data only about the updated slater matrices 
and the inverse of the slater matrices.

Methods: updating matrices, finding ratios between determinants, 
calculate laplacian or the gradient along 
some axis, ...
 
 */
	
#include "../orbital/orbital.h"

class slaterMatrix {
	private:
	//Class variables

		//Slater matrixes. Two separate matrixes containing spin up and spin 
		double** spin_up_matr;
		double** spin_down_matr;
		//Inverse matrixes
		double** inv_up_matr;
		double** inv_down_matr;
		//Number of particles of each spin orientation.
		int iCutoff;
		//Number of particles
		int iNumPart;
		//Number of dimensions
		int dim;
		//Remove? only used in waveFunction()
		double det_up;
		double det_down;
		
		//initializes orbital* array. Each object: orbital with unique wf, 
		//quantumnumbers ect.
		orbital* orbital_;
	
		//Vector containing all the variational parameters
		double* variational_parameters;
		//length of this vector (remove?)
		int iNumber_of_variational_parameters;

		//Experimental: store all variables that does not need to be updated 
		//each mc-cycle. ex: laplacian_up when down pareicle moved.
		// Problem: Not initialized before iCutoff+1 first collecting cycles.
		//double* grad_up;
		//double* grad_down;
	
	public:
	//Class variables

	//class methods
		/*
		   constructor
		   Input - (number of particles, number of particles in spin up/spin down state 
		   (same number assumed), Number of variational parameters, dimensions).
		   
		   Allocating determinant matrixes for spin ud and down states, and the inverse 
		   matrixes. Initializing orbital objects and class variables. The list of 
		   orbitals (with quantumnumbers) can be changed here, wave-functions must be 
		   changed in class orbital.
		 	*/
		slaterMatrix(int, int, int, int);
		/*
		   Input - position vector.
		   Initializing or updating slater matrix when more then one particle have new coords.	
			*/
		void initSlaterMatrix(double** r);
		/*
		   Find inverse (inv_..._matr) of both slater matrices (spin_..._matr).
		   using lapack functions to fint the inverse via lu-decomp
		   and backsubsitution.
		   
		   Storing the transpose of the inverse since most of the algos are 
		   running through the last index.
		   */
		void findInverse();
		/*
		   1 - Use O(N^2) routine to update the inverse matrices when
		   particle i_upd is moved to coords r_new. 
		   2 - Update the slater matrices.
		*/
		void update(double* r_new, int i_upd);
		/*
		   Input - (position vector of moved particle, index of moved particle).
		   
		   Using the inverse matrixes to calculate the determinant. r_new is the position 
		   of the new particle, i_upd the number of the particle.
		 	*/
		double waveFunction(double*, int);
		/*
		   Calculate ( \nabla_{axis,i_upd} \Psi ) / \Psi.
		   TEST!
		 	*/		   
		double grad(double* dR, int axis, int i_upd);
		/*
		   Calculate the laplacian of the slatermatrixes.
		   */
		double lapl(double** dR);
		/*
		   Clear all malloc'ed variables.
		   */
		void clear();
		/*
		   Print slatermatrices and the inverse.
		   */
		void print();
		/*
		   remove?
		   
		   Update all variational parameters. 
		   Input- array of length iNumber_of_variational_parameters.
			*/
		void updateVariationalParameters(double*);
		/*
		   Remove.

		   Updating slater matrix when particle i_upd are moved to coord r_new.  
			*/
		void updateSlaterMatrix(double* r_new, int i_upd);
};
