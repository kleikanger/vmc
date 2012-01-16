/*
	Implementation of class slaterMatrix
	
	Contains:

	NB: A jastrow factor is included when calculating 
		the value of the wavefunction Psi(r). 

    NBB: XXX Only good for hamiltonians that isn't dependent of spin since D->Dup*Ddown.
   */
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"

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

	/****************************************
	 * 										*
	 *		Initializing orbital objects.	*
	 * 										*
	 ****************************************/	 
		
	orbital_ = new orbital[iNumPart];
	
	orbital_[0].setValues(1,0,true);//n=1,l=0,spin up.
	orbital_[1].setValues(2,0,true);//n=2,l=0,spin up.
	orbital_[2].setValues(1,0,false);//n=1,l=0,spin down.
	orbital_[3].setValues(2,0,false);//n=2,l=0,spin down.
		

	/****************************************
	 * 										*
	 * 										*
	 * 										*
	 ****************************************/

}//end constructor varMC()
/*//endvimfold*/

//Sets up the slater determinant matrixes: spinUpMatrix** and spinUpMatrix**.
//Calls the varMC::orbitals() function find the matrix elements.
void slaterMatrix::updateSlaterMatrix(double** partPos, double dAlpha){
	/*//startvimfold*/
	
	int i,j;
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinUpMatrix[i][j] = orbital_[i].valueWF(partPos[j]); //orbitals(i,pParticlePositions[j],dAlpha);
		}
	}
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinDownMatrix[i][j] = orbital_[i].valueWF(partPos[j+iCutoff]); //(i,pParticlePositions[j+iCutoff],dAlpha);
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
	updateSlaterMatrix(dR,1.0);
	return determinant(spinDownMatrix,iCutoff)*determinant(spinUpMatrix,iCutoff)*jastrow(dR,alpha);
}
//endvimfold

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
            if (orbital_[i].spinUp() == orbital_[j].spinUp()) {

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

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
