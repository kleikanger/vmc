/*
	Implementation of class slaterMatrix
	
	Contains:

	NB: A jastrow factor is included when calculating 
	the value of the wavefunction Psi(r). 

	NBB: XXX Only good for hamiltonians that isn't dependent of spin since D->Dup*Ddown. New variable:
	bool bSeparate_up_down_matrixes.

	NBBB: Assumes that nr. of spin up particles equals nr of spin down particles. 
	If number of spin up particles does not equal number of spin down program needs to be 
	modified. New variables: iCutoff_up, iCutoff_down.
 	  
 	NBBBB: Variational parameters should be initialized in update, since slaterMatrix is updated between 
	each run. should be class variables, accessible in the entide class.
 
 	NBBBBB:Implement LU-fact. routine to calculate determinants for larger system.
 
 	NBBBBBB!:Implement lapack/cblas or lapack++/blas++ or blitz. 
 
 */
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"

using std::cout;

	
//constructor
slaterMatrix::slaterMatrix(int iNp,int iCo,int iNovp) {
	/*//startvimfold*/

	//More elegant way to initialize class variables?
	iNumPart=iNp;
	iCutoff=iCo;
	variational_parameters = new double[iNovp];
	iNumber_of_variational_parameters=iNovp;
	
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
	
	// *****************************************
#if 0		
	//States 0-(iCutoff-1) enters spinUpMatrix. 
	//States iCutoff-(iNumPart-1) enters spinDownMatrix 
	orbital_[0].setValues(1,0,true);//n=1,l=0,spin up.
	orbital_[1].setValues(2,0,true);//n=2,l=0,spin up.
	orbital_[2].setValues(1,0,false);//n=1,l=0,spin down.
	orbital_[3].setValues(2,0,false);//n=2,l=0,spin down.

	//testing for larger system
	//XXX Remove later
#endif
#if 1
	orbital_[0].setValues(1,0,true);//n=1,l=0,spin up.
	orbital_[1].setValues(2,0,true);//n=2,l=0,spin up.
	orbital_[2].setValues(3,0,true);//n=1,l=0,spin down.
	orbital_[3].setValues(4,0,true);//n=2,l=0,spin down.
	if (iNumPart>4)
	orbital_[4].setValues(1,0,false);//n=1,l=0,spin up.
	if (iNumPart>5)
	orbital_[5].setValues(2,0,false);//n=2,l=0,spin up.
	if (iNumPart>6)
	orbital_[6].setValues(3,0,false);//n=1,l=0,spin down.
	if (iNumPart>7)
	orbital_[7].setValues(4,0,false);//n=2,l=0,spin down.
#endif
	// ****************************************

	//Use cofactormatrix only if determinantmatrix are bigger than (2x2)
	if (iCutoff<3){
		bUse_cofactors=false;
		cout<<"\nSlatermatrix: iCutoff<3. Cofactormatrixes not allocated.\n";
	} else {
		bUse_cofactors=true;
		cout<<"\nSlatermatrix: iCutoff>2. Allocating cofactormatrixes.\n";
	}

	//if matrix>2x2. We use cofactors to calculate the determinants 
	//change name to inverse_up and inverse_down.
	if (bUse_cofactors) {
		//declare and allocate matrixes
		cofactorMatrix_up = new double*[iCutoff];
		cofactorMatrix_down = new double*[iCutoff];
		for (int i=0; i<iCutoff; i++){
			cofactorMatrix_up[i] = new double[iCutoff];
			cofactorMatrix_down[i] = new double[iCutoff];
			//neccesary?
			for (int j=0;j<iCutoff;j++){
				cofactorMatrix_up[i][j]=0.0;
				cofactorMatrix_down[i][j]=0.0;
			}
		}
	} //else {
	//	delete  cofactorMatrix_down;
	//	delete  cofactorMatrix_up;
	//}

}//end constructor varMC()
/*//endvimfold*/

const bool slaterMatrix::useCofact(){
//startvimfold
	return bUse_cofactors;
}
//endvimfold

//Sets up the slater determinant matrixes: spinUpMatrix** and spinUpMatrix**.
//Calls the varMC::orbitals() function find the matrix elements.
void slaterMatrix::updateSlaterMatrix(double** partPos){
	/*//startvimfold*/
	
	int i,j;
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinUpMatrix[i][j] = orbital_[i].valueWF(partPos[j]); //orbitals(i,pParticlePositions[j],dAlpha);
		}
	}
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spinDownMatrix[i][j] = orbital_[i+iCutoff].valueWF(partPos[j+iCutoff]); //(i,pParticlePositions[j+iCutoff],dAlpha);
		}
	}
	//Update cofactor matrixes (maybe not necc. to do this each update, 
	//only if MC jump accepted and we need to take derivatives??)
	if (bUse_cofactors){
	 	updateCofactors(cofactorMatrix_up,spinUpMatrix,iCutoff);
	 	updateCofactors(cofactorMatrix_down,spinDownMatrix,iCutoff);
	}
}//End function varMC::DeterminantMatrix()
/*//endvimfold*/

//Calculate the elements of the cofactormatrix.
//Examplecall: updateCofactors(cofactorMatrix_up,spinUpMatrix,iCutoff)
//XXX XXX XXX XXX
//CHANGE TO UPDATE INVERSE MATRIX maybe better to uptate inverse matrix
//NEW METHOD: CALCULATEINVERSEMATRIX: No need for cofactors. better to invert the matrix directly. 
void slaterMatrix::updateCofactors(double** cofactor_matrix, double** determinant_matrix, int iDim_determinant){
//startvimfold
	//Keeping track of the determinants det_up,det_down (spin up and spin down) 
	//so that we dont have to calculate any determinants more than once. 
	//Using variables in waveFunction(double*,double,int).
	
	int i,j,k,l;
	double sign;
	//allocate minor matrix
	double** minorMatrix = new double*[(iDim_determinant-1)];
	for (i=0;i<(iDim_determinant-1);i++){
		minorMatrix[i] = new double[iDim_determinant-1];
		for (j=0;j<(iDim_determinant-1);j++){
			minorMatrix[i][j]=0.0;
		}
	}
	for (i=0;i<iDim_determinant;i++){
		for (j=0;j<iDim_determinant;j++){
		//Find minor matrix for entry i,j
			for (k=0;k<i;k++){
				for (l=0; l<j; l++){
					minorMatrix[k][l]=determinant_matrix[k][l];
				}
				for (l=j+1;l<iDim_determinant; l++){
					minorMatrix[k][l-1]=determinant_matrix[k][l];
				}
			}
			for (k=i+1;k<iDim_determinant;k++){
				for (l=0; l<j; l++){
					minorMatrix[k-1][l]=determinant_matrix[k][l];
				}
				for (l=j+1;l<iDim_determinant; l++){
					minorMatrix[k-1][l-1]=determinant_matrix[k][l];
				}
			}
		//Finding the cofactors
		if ((i+j)%2==0) sign=1.0; 
		else sign=-1.0;
		cofactor_matrix[i][j]=sign*determinant(minorMatrix,(iDim_determinant-1));	
		}
	}
		//Need'nt calculate both determinants each time. eg: when calculating derivatives.
		//XXX: Smart to always take determinant along first row?
		det_up=0.0;
		for (i=0;i<iCutoff;i++){
			det_up+=spinUpMatrix[0][i]*cofactorMatrix_up[0][i];
		}
		det_down=0.0;
		for (i=0;i<iCutoff;i++){
			det_down+=spinDownMatrix[0][i]*cofactorMatrix_down[0][i];
		}

	//delete minor matrix XXX XXX
}//End function slatermatrix::updateCofactors()
//endvimfold

//Should maybe be done directly in orbital class. Not inportant,
//not frequently used. UPDATE VAR PARAM IN ORBITAL CLASS AS WELL from this func. 
void slaterMatrix::updateVariationalParameters(double* vp){
//startvimfold
	for (int i=0;i<iNumber_of_variational_parameters;i++){
		variational_parameters[i]=vp[i];
	}
}
//endvimfold

//Calculates determinant of (iNxiN) matrix by recursion.
//Iterative function probably faster
//O(n^4)(??) For larger matrixes algo with LU decomp. much faster.
//XXX XXX XXX XXX IMPLEMENT LU-algo. Evt use lapack or other library. 
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
//Brute force calculation of the determinants.
//slatermatrixes already calculated. (updateSlaterMatrix())
double slaterMatrix::waveFunction(){
//startvimfold
	//cout<<"\nJastrow not included\n";
	//updateSlaterMatrix(r);
	if (!bUse_cofactors){
		return determinant(spinDownMatrix,iCutoff)*determinant(spinUpMatrix,iCutoff);//*jastrow(r);
	} else {
		return det_down*det_up;
	}
}
//endvimfold

//double slaterMatrix::waveFunction(double * dR, int iCofac_column, int input_integer){
double slaterMatrix::waveFunction(double * dR, int iCofac_column){
	//startvimfold
	if (!bUse_cofactors){
		cout<< "\nerror in waveFunction(). Since iCutoff=2 the determinant must be calculated directly."
			<<"\nUse function call waveFunction(double**,double).";
		exit(1);
	}
	int i;
	//double funct_[i][iNumPart];
	//switch (input_integer):
	//case 1: 	for (i=0;i<iNumPart;i++)
	//				funct_[i]=orbital_[i].valueWF(); 
	//			break;//now the ratio wf_new/wf_old is returned
	//case_2: 	for (i=0;i<iNumPart;i++)
	//				funct_[i]=orbital_[i].wFDeriv1(); //now the gradient is returned 
	//			break;
	//case 3: 	for (i=0;i<iNumPart;i++)
	// 				funct_[i]=orbital_[i].wFDeriv2(); //now the laplacian (Del^2 wf_new)/wf_old is returned
	// 			break;
	//default: 
	//			cout<<"\nerror slaterMatrix::wavefunction(). Value of input_integer"; 
	// 			cout<<"does not fit any case in switch statement.\n" 
	// 			exit(1);
	//
	// ********************************************************************************
	// 	1: change orbital_ to funct_ below 2: make sure det_up and det_down are updated in upd.func.
	// ********************************************************************************
	double temp=0.0;
	if (iCofac_column<iCutoff){//evt if(orbital_[iCofac_column].spinUp();)
		for (i=0;i<iCutoff;i++){
			temp+=orbital_[i].valueWF(dR)*cofactorMatrix_up[i][iCofac_column];
		}
		return temp*det_down;
	} else if ((iCofac_column>=iCutoff)&&(iCofac_column<iNumPart)){ //XXX REMOVE TEST LATER 
		for (i=iCutoff;i<iNumPart;i++){
			temp+=orbital_[i].valueWF(dR)*cofactorMatrix_down[i-iCutoff][iCofac_column-iCutoff];
		}
		return temp*det_up;
	} else {
		cout<< "error in function waveFunction: iCofac_column out of bounds"; 
		exit(1);
	}
}//End function waveFunction();
//endvimfold

//Update C(cofactor_matrix->inverse). O(N**2) algo.
//*********************************************************
// 				XXX XXX	TEST XXX XXX
//*********************************************************
//Test by doing updates here, and calculating the matrix directly 
//in two different objects.
void slaterMatrix::updateInverse(double* dR, int iUp){
	
	//New determinant
	double RR = waveFunction(dR,iUp);
	int k,j,l;
	double temp_factor=0;

	for (k=0;k<iCutoff;k++){
		for (j=0;j<iCutoff;j++){
			
			if (iUp<iCutoff){//evt if(orbital_[iCofac_column].spinUp();)
				for (l=0;l<iCutoff;l++){
					temp_factor+=orbital_[l].valueWF(dR)*cofactorMatrix_up[l][j];
				}
				temp_factor*=det_down;
			} else if ((iUp>=iCutoff)&&(iUp<iNumPart)){ //XXX REMOVE TEST LATER
				for (l=iCutoff;l<iNumPart;l++){
					temp_factor+=orbital_[l].valueWF(dR)*cofactorMatrix_down[l-iCutoff][j-iCutoff];
				}
				temp_factor*= det_up;
			} else {
				 cout<< "error in function updateInverse: iUp out of bounds"; 
				 exit(1);
			}
			//Careful here!
			if (iUp!=j&&iUp<iCutoff){
				temp_factor*=cofactorMatrix_up[k][iUp]/RR;
				temp_factor+=cofactorMatrix_up[k][j];
			} else if (iUp!=j&&iUp>=iCutoff){
				temp_factor*=cofactorMatrix_down[k][iUp-iCutoff]/RR;
				temp_factor+=cofactorMatrix_down[k][j];
			} else if (iUp=j&&iUp<iCutoff){
				temp_factor*=cofactorMatrix_up[k][iUp]/RR;
			} else if (iUp=j&&iUp>=iCutoff){
				temp_factor*=cofactorMatrix_down[k][iUp-iCutoff]/RR;
			}
			if (k<iCutoff){
				cofactorMatrix_up[k][j]=temp_factor;
			} else {
				cofactorMatrix_down[k][j]=temp_factor;
	 		}
		temp_factor=0;
		}
	}
}

//Calculates the jastrow-factor
double slaterMatrix::jastrow(double** r){
//startvimfold
	double length;
	int a;
	
	//variational parameter.
	double alpha = variational_parameters[0];

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

///*TEST OF PROGRAM
//startvimfold
#include <iomanip>
int main(){

	int iNumPart=8;
	int iCutoff=4;
	
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

	slaterMatrix slater_Matrix(iNumPart,iCutoff,1);
	slaterMatrix slater_Matrix2(iNumPart,iCutoff,1);
	
	double wfnew, wfold;
	double ideal_step=2.0;

	long idum;
  	idum= - ( .6 );//time(NULL)*(myrank));
	
	//for (int iTemp=0; iTemp < number_of_alpha_variations; iTemp++){
	
		//cumulative_energy=0;
		//local_energy=0;
	
		//''random'' startposition 
	for (int g=0; g<1; g++){   

		for (int i = 0; i < iNumPart; i++) { 
   			for (int j=0; j < 3; j++) {
 				partPos[i][j] += ideal_step*(ran2(&idum) -0.5);
   		  	}
   		}
		
		double  ar[1];
	   	ar[0]=0.13;	
		
		slater_Matrix.updateVariationalParameters(ar);
		slater_Matrix.updateSlaterMatrix(partPos);	
		
		slater_Matrix2.updateVariationalParameters(ar);
		slater_Matrix2.updateSlaterMatrix(partPos);	
		
		int changed_part=5;
		for (int j = 0; j < iNumPart; j++) { 
 				partPos[changed_part][j] += ideal_step*(ran2(&idum) -0.5);
   		}
		
		slater_Matrix.updateVariationalParameters(ar);
		slater_Matrix.updateSlaterMatrix(partPos);	
		
		slater_Matrix2.updateInverse(partPos[changed_part],changed_part);
		
		cout<<"\n\n";
		for (int i = 0; i < iCutoff; i++) { 
			for (int j = 0; j < iCutoff; j++ ){
				cout<<slater_Matrix.cofactorMatrix_down[i][j]<<"\t";
			}
			cout<<"\n";
		}

		cout<<"\n\n";
		for (int i = 0; i < iCutoff; i++) { 
			for (int j = 0; j < iCutoff; j++ ){
				cout<<slater_Matrix2.cofactorMatrix_down[i][j]<<"\t";
			}
			cout<<"\n";
		}

#if 0
		for (int i=0; i<iNumPart;i++){
			cout<<setw(8)<<setprecision(16)<<"WF direct  :"<<slater_Matrix.waveFunction()<<"\n";
			cout<<setw(8)<<setprecision(16)<<"WF cofactor:"<<slater_Matrix.waveFunction(partPos[i],i)<<"\n";
		}
		
		for (int i=0; i<iNumPart;i++){
			cout<<setw(8)<<setprecision(16)<<"WF direct  :"<<slater_Matrix2.waveFunction()<<"\n";
			cout<<setw(8)<<setprecision(16)<<"WF cofactor:"<<slater_Matrix2.waveFunction(partPos[i],i)<<"\n";
		}
#endif

	}

}//End function 
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold