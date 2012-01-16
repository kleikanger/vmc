/*
Variational Monte Carlo class. Containing data only about the states (updated slater matrixes).
Doing operations on the slater matrixes. eg: finding determinants, calculating energies, ect.
 */


//This class shold contain the slatermatrixes and perform all operations on the slatermatrixes.
//eg: find energy, find quantum force, find determinant, ...
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
	
		//initializes orbital* array. Each object one orbital with unique wf, quantumnumbers ect.
		orbital* orbital_;

		//Number of particles
		int iNumPart;

		//Number of particles each spin orientation
		int iCutoff;

	//class methods
		
		//Initializing slater matrix
		void updateSlaterMatrix(double**, double);
		
		//Calculates determinant of the slater matrixes	
		double determinant(double**,int);

		//Calculates jastrowfactor
		//move to private later
		double jastrow(double**, double);
		
		
		

	public:

	//Class variables

	//class methods

		//constructor
		slaterMatrix(int, int);		
		
		//input: positions of all particles + variational parameter.
		//-Updates slater matrixes.
		//-Calculates wf in point dR (determinants) and multiplies with the jastrow factor with a variational parameter alpha.
		double waveFunction(double**, double);
		
		//Delete matrixes and objects
		//void calculate_qforce(double** r, double** q_force, double beta){
};

