/*
Variational Monte Carlo class. Containing data only about the states (updated slater matrixes).
Doing operations on the slater matrixes. eg: finding determinants, calculating energies, ect.
 */

#include "./orbital/orbital.h"

class slaterMatrix {

	private:

	//Class variables
		//Slater matrixes. Two separate matrixes containing spin up and spin 
		//down particles. (iNumPart x iNumPart) matrix.
		double** spinUpMatrix;
		double** spinDownMatrix;
		////Containing the positions to all patricles. (iNumPart x 3) matrix.
		//double** ppParticlePositions;
		
		orbital* orbital_;

		//Number of particles
		int iNumPart;
		//Number of particles each spin orientation
		int iCutoff;

	//class methods
		double orbitals(int, double*, double);

	public:

	//Class variables

	//class methods

		//constructor
		slaterMatrix(int, int);		
		
		//Initializing slater matrix
		void updateSlaterMatrix(double**, double);
	
		//Move to private later	
		double determinant(double**,int);

		//Calculates jastrowfactor
		double jastrow(double**, double);

		//Calculates wf in poinr dR with jastrow factor with variational parameter alpha.
		//slatermatrixes already calculated. (update())
		double waveFunction(double**, double);

};
#if 0 
   /*This class should contain all information about the different orbitals
	*eg: orbital wf_1s(quantumnumbers,...,).
	*methods: eg. wf_1s.value(r,alpha);
	*/
class orbital{

}
#endif
