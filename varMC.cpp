/*

	Variational Monte Carlo implementation.

   */
#include <cstdlib>
#include <cmath>
#include "varMC.h"
#include <iostream>
#include "lib/lib.h"

using std::cout;

//constructor
slaterMatrix::slaterMatrix(int a,int b) {
	/*//startvimfold*/

	//More elegant way to initialize class variables?
	iNumPart=a;
	iCutoff=b;
	
	//Allocating new matrices
	spinUpMatrix = new double*[iCutoff];
	spinDownMatrix = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){
		spinDownMatrix[i] = new double[iCutoff];
		spinUpMatrix[i] = new double[iCutoff];
		//neccesary?
		for (int j=0;j<iCutoff;j++){
			spinUpMatrix[i][j]=0.0;
			spinDownMatrix[i][j]=0.0;
		}
	}

	/*
	
	orbital* orbital_[iNumPart];
	
	orbital_[0] = new orbital(1,0,true);//n=1,l=0,spin up.
	orbital_[1] = new orbital(2,0,true);//n=2,l=0,spin up.
	orbital_[2] = new orbital(1,0.false);//n=1,l=0,spin down.
	orbital_[3] = new orbital(2,0,false);//n=2,l=0,spin down.
	
	orbital_[i]->valueWF(a)
		
	   */

}//end constructor varMC()
/*//endvimfold*/

//Sets up the slater determinant matrixes: spinUpMatrix** and spinUpMatrix**.
//Calls the varMC::orbitals() function find the matrix elements.
void slaterMatrix::updateSlaterMatrix(double** pParticlePositions, double dAlpha){
	/*//startvimfold*/
	
	int i,j;
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinUpMatrix[i][j] = orbitals(i,pParticlePositions[j],dAlpha);
		}
	}
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinDownMatrix[i][j] = orbitals(i,pParticlePositions[j+iCutoff],dAlpha);
		}
	}


}//End function varMC::DeterminantMatrix()
/*//endvimfold*/

//Calculates determinant of (iNxiN) matrix by recursion.
//Iterative function probably faster
//O(n^4)(??) For larger matrixes algo with LU decomp. much faster.
//Better if not in class varMC?
double slaterMatrix::determinant(double** ppM, int iN){
/*//startvimfold*/
	//if iN == 2, calculating smallest submatrix and breaking recursion.
	if (iN == 2){
		return  ppM[0][0]*ppM[1][1]-ppM[0][1]*ppM[1][0];
	}

	double dSum = 0.0;
	
	//declare and delete within forloop?
	double** ppSubMatrix = new double*[iN-1];

	//recursive loop
	for (int i=0;i<iN;i++){
		//Subtracting submatrix
		for (int j=0;j<i;j++){
			//1 holds adress to 2. element. New matrix (iN-1xiN-1).
			//last variable is iterated over later, so must take determinant along first coloumn.
			ppSubMatrix[j] = &ppM[j][1];
		}
		for (int j=i+1;j<iN;j++){
			ppSubMatrix[j-1] = &ppM[j][1];
		}

		if (i%2 == 0) {
			dSum += ppM[i][0] * determinant(ppSubMatrix,iN-1);
		} else {
			dSum -= ppM[i][0] * determinant(ppSubMatrix,iN-1);
		}
	}

	delete [] ppSubMatrix;
	return dSum;
}//end method varMC::determinant
/*//endvimfold*/

//Calculates wf including jastrow factor with variational parameter alpha.
//slatermatrixes already calculated. (updateSlaterMatrix())
double slaterMatrix::waveFunction(double ** dR, double alpha){
//startvimfold
	return determinant(spinDownMatrix,iCutoff)*determinant(spinUpMatrix,iCutoff)*jastrow(dR,alpha);
}
//endvimfold

//Returns value in position r for orbital with quantum number n,l.
//And variational parameter alpha. Alpha not used here.
double slaterMatrix::orbitals(int iN,  double* dR, double alpha){
/*//startvimfold*/
	
	//Constants found in HF simulation. More decimals?
	double a=0.9955496248;
	double b=0.09423876105;
	
	//calculate distance to core
	double dAbsR=0;
	for (int i=0;i<3;i++){
		dAbsR+=dR[i]*dR[i];
	}
	dAbsR=sqrt(dAbsR);
	
	//n=1,l=0
	double psi_10 = 16.0*exp(-4.0*dAbsR);	//not normalized
	//n=2,l=0
	double psi_20 = 2.82842712474619*(1-dAbsR)*exp(-2.0*dAbsR); 	



	//Best to insert explicit functions. faster!
	//XXX break; statement always neccesary?
	switch (iN){
	case 0:
			return a*psi_10-b*psi_20;
	case 1: 
   			return -b*psi_10-a*psi_20;			
	case 2:
			return -a*psi_10+b*psi_20;	
	case 3:
			return b*psi_10+a*psi_20;	
	default:
			cout<<"error: Non-matching quantum-number in function varMC::orbitals().";
			exit(1); 
	}
}//End varMC::orbitals()
/*//endvimfold*/

//Calculates the jastrow-factor
double slaterMatrix::jastrow(double** r, double alpha){
//startvimfold
	double length;
	int a;

    double argument = 0;
    for (int i = 0; i < iNumPart - 1; i++) {
        for (int j = i + 1; j < iNumPart; j++) {
   
            // Beregner avstand mellom partikkel i og j:
            length=0.0;
				for (a=0; a<3; a++){
					length+= (r[i][a]-r[j][a])*(r[i][a]-r[j][a]);
				}
			length=sqrt(length);

            // Dersom partikkel i og j har samme spinn:
            if ((i < iNumPart/2 && j < iNumPart/2) || (i >= iNumPart/2 && j >= iNumPart/2)) {

                argument += 0.25*length/(1 + alpha*length);
            }

            // Dersom partikkel i og j har ulike spinn:
            else {
                argument += 0.5*length/(1 + alpha*length);
            }
        }
    }

    return exp(argument);
}//End function slaterMatrix::jastrow()
//endvimfold

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

	//allocating positions of the particles assuming 3 dimensions.
	//First <iCutoff> indexes: pos. of spin up particles. 
	//Next  <iCutoff> indexes: pos. of spin down particles.
	double** partPos = new double*[iNumPart];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[3];
		//neccesary?
		for (int j=0;j<3;j++){
			partPos[i][j]=0.0;
		}
	}
	double** newPartPos = new double*[iNumPart];
	for (int i=0; i<iNumPart; i++){
		newPartPos[i] = new double[3];
		//neccesary?
		for (int j=0;j<3;j++){
			newPartPos[i][j]=0.0;
		}
	}
	
	//Ecumulative: will be needed when parallellizing code.
	double Ecumulative;
	double Elocal;
	
	//initialize slater matrix objects
	slaterMatrix newMatrix(iNumPart,iCutoff);
	slaterMatrix oldMatrix(iNumPart,iCutoff);
	
	double wfnew, wfold;

	long idum;
  	idum= - (1 +  time(NULL));//time(NULL)*(myrank));
	

	//''random'' startposition 
    for (int i = 0; i < iNumPart; i++) { 
    	for (int j=0; j < 3; j++) {
 			partPos[i][j] = ideal_step*(ran2(&idum) -0.5);
      	}
    }

	oldMatrix.updateSlaterMatrix(partPos,1.0);
	wfold=oldMatrix.waveFunction(partPos,alpha);	
	
	//mc sampling
	int iNumber_of_iterations=10;
	int iAccepted_jumps=0;

	int iIteration_count=0;	
	while (iIteration_count < iNumber_of_iterations){
		iIteration_count++;

		for (int i=0; i<iNumPart; i++){
			for (int j=0; j<3; j++){
				cout<<partPos[i][j]<<"\t";
			}
		cout<<"\n";
		}

	}//End while, 


#if 0
	//testing determinant function
	slaterMatrix A(4,2);
	int iN = 6;
	double** testmatrix = new double*[iN];
		for (int i=0; i<iN; i++){
			testmatrix[i] = new double[iN];
			for (int j=0; j<iN; j++){
				testmatrix[i][j] = 0;
				//((double)i+1.0)+((double)j+1.0);
			}
		testmatrix[i][i]=1;
		}
		//testmatrix[iN-2][iN-1]=45;
		//testmatrix[iN-1][iN-2]=47;


		for (int i=0; i<iN; i++){
			for (int j=0; j<iN; j++){
				std::cout<<testmatrix[i][j]<<"\t";
			}

		std::cout<<"\n";
		
		}
		double a = A.determinant(testmatrix,iN);
		std::cout<<"\n"<<" "<<a<<"\n";	
#endif
}//End main()
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
