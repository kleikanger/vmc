
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"


///*TEST OF PROGRAM
//startvimfold

#include <iomanip>
int main(){
	
	long idum;
  	idum= - ( .2 );//time(NULL)*(myrank));

	int iNumPart=6;
	int iCutoff=3;
	
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

	double wfnew, wfold;
	double ideal_step=2.0;

	//for (int iTemp=0; iTemp < number_of_alpha_variations; iTemp++){
	//cumulative_energy=0;
	//local_energy=0;
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < 3; j++) {
			partPos[i][j] += ideal_step*(ran2(&idum) -0.5);
		}
	}
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < 3; j++) {
			partPos[i][j] += ideal_step*(ran2(&idum) -0.5);
		}
	}
	double  ar[1];
	ar[0]=0.13;	
	
	slaterMatrix a(iNumPart,iCutoff,1);
	slaterMatrix b(iNumPart,iCutoff,1);
	
	a.initSlaterMatrix(partPos);	
	b.initSlaterMatrix(partPos);	
	
	a.updateVariationalParameters(ar);
	b.updateVariationalParameters(ar);
	
	a.findInverse();
	b.findInverse();

#if 1	
	for (int g=0; g<10000000; g++){   
		
		double* test_qa = new double[3];
		int changed_part=0;
		for (int j = 0; j < 3; j++) { 
	 		partPos[changed_part][j] += ideal_step*(ran2(&idum) -0.5);
		}
		//a:probably as fast when working with few particles? TEST:when progr. working.
		//Uncomment lines below if one need to compare with unoptimalized method.
		//a.initSlaterMatrix(partPos);	
		//a.findInverse();
		
		//b:faster when systems are larger.
		b.updateInverse(partPos[changed_part],changed_part);
		//included in update inverse
		//b.updateSlaterMatrix(partPos[changed_part],changed_part);	
	}
#endif
	a.print();
	b.print();
	
	a.clear();
	b.clear();

}//End function 
//endvimfold
