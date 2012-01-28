/*
	Implementation of class slaterMatrix
	
	Contains:

	NB: XXX Only good for hamiltonians that isn't dependent of spin since D->Dup*Ddown. New variable:
	bool bSeparate_up_down_matrixes.

	NBB: Assumes that nr. of spin up particles equals nr of spin down particles. 
	If number of spin up particles does not equal number of spin down program needs to be 
	modified. New variables: iCutoff_up, iCutoff_down.
 
 */
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"
#include <cblas.h>

using std::cout;

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int*, int*, double*, int*, int*, int*);
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int*, double*, int*, int*, double*, int*, int*);
};
	
//constructor
slaterMatrix::slaterMatrix(int iNp,int iCo,int iNovp) {
	/*//startvimfold*/

	//More elegant way to initialize class variables?
	iNumPart=iNp;
	iCutoff=iCo;
	variational_parameters = new double[iNovp];
	iNumber_of_variational_parameters=iNovp;
	
	//Allocating new matrices
	inv_up_matr = new double*[iCutoff];
	inv_down_matr = new double*[iCutoff];
	spin_up_matr = new double*[iCutoff];
	spin_down_matr = new double*[iCutoff];

	for (int i=0; i<iCutoff; i++){
		inv_down_matr[i] = new double[iCutoff];
		inv_up_matr[i] = new double[iCutoff];
		spin_up_matr[i] = new double[iCutoff];
		spin_down_matr[i] = new double[iCutoff];
	//	for (int j=0;j<iCutoff;j++){
	//		inv_up_matr[i][j]=0.0;
	//		inv_down_matr[i][j]=0.0;
	//		spin_up_matr[i][j]=0.0;
	//		spin_down_matr[i][j]=0.0;
	//	}
	}
	
	//Initializing orbital objects
	orbital_ = new orbital[iNumPart];
	//States 0-(iCutoff-1) enters inv_up_matr. 
	//States iCutoff-(iNumPart-1) enters inv_down_matr 
	for (int i=0; i<iCutoff; i++)
		orbital_[i].setValues(i+1,0,true);//n=i,l=0,spin up.
	for (int i=iCutoff; i<iNumPart; i++)
		orbital_[i].setValues(i-iCutoff+1,0,false);//n=i,l=0,spin down.

}//end constructor varMC()
/*//endvimfold*/
void slaterMatrix::initSlaterMatrix(double** partPos){
	/*//startvimfold*/
	int i,j;
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spin_up_matr[i][j] = orbital_[j].valueWF(partPos[i]); // Transpose
		}
	}
	for (i=0; i<iCutoff; i++){
		for (j=0; j<iCutoff; j++){
			spin_down_matr[i][j] = orbital_[j+iCutoff].valueWF(partPos[i+iCutoff]); //(i,pParticlePositions[j+iCutoff],dAlpha);
		}
	}
}//End function varMC::DeterminantMatrix()
/*//endvimfold*/
//NOT NEEDED ANYMORE.
void slaterMatrix::updateSlaterMatrix(double* partPos, int i_upd){
	/*//startvimfold*/
	int i;
	if (i_upd<iCutoff){
		for (i=0; i<iCutoff; i++){
			spin_up_matr[i_upd][i] = orbital_[i].valueWF(partPos); // Transpose
		}
	} else {
		for (i=0; i<iCutoff; i++){
			spin_down_matr[i_upd-iCutoff][i] = orbital_[i+iCutoff].valueWF(partPos); // Transpose
		}
	}	
}//End function varMC::DeterminantMatrix()
/*//endvimfold*/
void slaterMatrix::findInverse(){
//startvimfold
	int n = iCutoff;
	int *ipiv = new int[n];//n+1?
	int lwork = n*n;
	double *work = new double[lwork];
	int info;
	double* B = new double[n*n];
	
	for (int h=0; h<2; h++){
		//Lapack using coloumn major!
		//Matrixes must be transposed before they are sent to inverse
		int k=0;
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				if (h==0) B[k+j] = spin_up_matr[j][i];
				else B[k+j] = spin_down_matr[j][i];
			}
		k+=n;
		}
			
		dgetrf_(&n,&n,B,&n,ipiv,&info);
		//Even if we know that that the determinant is !=0, maybe we should check if matrix is singular or near to...
		if (info==0) {
			dgetri_(&n,B,&n,ipiv,work,&lwork,&info);
			}
		if (info!=0) {
			cout<<"error slaterMatrix::findInverse in lapack dgetri_ or degetrf_ info="<<info<<"\n";
			exit(1);
			}
		k=0;
		//inverse
		for (int i=0; i<n; i++){
			for (int j=0; j<n; j++){
				if (h==0) inv_up_matr[i][j]=B[k+j];
				else inv_down_matr[i][j]=B[k+j];
			}
		k+=n;
		}
	}    
	delete B;
    delete ipiv;
    delete work;

}//End function slatermatrix::updateCofactors()
//endvimfold
void slaterMatrix::updateInverse(double* d_R, int i_upd){
//startvimfold	
	
	int k,j;
	double d_upd[iCutoff];
	double R, temp;
    
	if (i_upd<iCutoff){
	//SPINUP 
		//find D^new/D^old
		for (k=0;k<iCutoff;k++){
				d_upd[k]=orbital_[k].valueWF(d_R);
		}
		//Updating ingerse Blas routine : d_up*inv_up_matr[j][:], 
		//remember inverse M are transposed
		R = cblas_ddot(iCutoff,d_upd,1,inv_up_matr[i_upd],1);

		for (j=0;j<iCutoff;j++){
			if (i_upd!=j){
				temp = cblas_ddot(iCutoff,d_upd,1,inv_up_matr[j],1);
				cblas_daxpy( iCutoff,-temp/R,inv_up_matr[i_upd],1,inv_up_matr[j], 1);
			}
		}
		//inv_up_matr[i_upd][:]<-inv_up_matr[i_upd][:]/R
		cblas_dscal(iCutoff,1/R,inv_up_matr[i_upd],1);
		//Copy d_upd to spin_down_matr, updating spinupmatr.
		cblas_dcopy(iCutoff, d_upd, 1, spin_up_matr[i_upd], 1);

	} else {
	//Repeat all for SPINDOWN if i_upd>iCutoff.
		for (k=0;k<iCutoff;k++){
				d_upd[k]=orbital_[k+iCutoff].valueWF(d_R);
		}
		R = cblas_ddot(iCutoff,d_upd,1,inv_down_matr[i_upd-iCutoff],1);

		for (j=0;j<iCutoff;j++){
			if ((i_upd-iCutoff)!=j){
				temp = cblas_ddot(iCutoff,d_upd,1,inv_down_matr[j],1);
				cblas_daxpy(iCutoff,-temp/R,inv_down_matr[i_upd-iCutoff],1,inv_down_matr[j],1);
			}
		}
		cblas_dscal(iCutoff,1/R,inv_down_matr[i_upd-iCutoff],1);
		cblas_dcopy(iCutoff, d_upd, 1, spin_down_matr[i_upd-iCutoff], 1);
	}	
}
//endvimfold
double slaterMatrix::waveFunction(double* dR, int iCofac_column){
	//startvimfold
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
	// 	1: change orbital_ to funct_ below 2: make sure det_up and det_down are updated .
	// ********************************************************************************
	det_up=0.0;
	det_down=0.0;
	//EASY TO USE cblas ddot()
	//if (iCofac_column<iCutoff){//evt if(orbital_[iCofac_column].spinUp();)
		for (i=0;i<iCutoff;i++){
			det_up+=orbital_[i].valueWF(dR)*inv_up_matr[iCofac_column][i];
		}
	//EASY TO USE cblas_ddot()
	//} else if ((iCofac_column>=iCutoff)&&(iCofac_column<iNumPart)){ //XXX REMOVE TEST LATER 
		for (i=iCutoff;i<iNumPart;i++){
			det_down+=orbital_[i].valueWF(dR)*inv_down_matr[iCofac_column-iCutoff][i-iCutoff];
		}
		return det_up*det_down;
		//switch (input_integer==1):
		//case 1:
		//RATIO DET_new/DET_old
		//case 2: 
		//1. derived
		//case 3:
		//2. derived
		//	return det_up + det_down + jastrow();
	
		
	//} else {
	//	cout<< "error in function waveFunction: iCofac_column out of bounds"; 
	//	exit(1);
	//}
}//End function waveFunction();
//endvimfold
void slaterMatrix::updateVariationalParameters(double* vp){
//startvimfold
	for (int i=0;i<iNumber_of_variational_parameters;i++){
		variational_parameters[i]=vp[i];
	}
}
//endvimfold
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
void slaterMatrix::clear(){
//startvimfold
	delete variational_parameters;
	for (int i=0; i<iCutoff; i++){
		delete inv_down_matr[i];
		delete inv_up_matr[i];
		delete spin_up_matr[i];
		delete spin_down_matr[i];
	}
	delete inv_down_matr;
	delete inv_up_matr;
	delete spin_up_matr;
	delete spin_down_matr;
	delete orbital_; //DELETE ORBITAL[i]
}
//endvimfold
void slaterMatrix::print(){
//startvimfold	
	cout<<"\ninvUp\n";
	for (int i = 0; i < iCutoff; i++) { 
		for (int j = 0; j < iCutoff; j++ ){
			cout<<setprecision(16)<<inv_up_matr[i][j]<<"\t";
		}
		cout<<"\n";
	}

	cout<<"\nspinUp:\n";
	for (int i = 0; i < iCutoff; i++) { 
		for (int j = 0; j < iCutoff; j++ ){
			cout<<setprecision(16)<<spin_up_matr[i][j]<<"\t";
		}
		cout<<"\n";
	}
	cout<<"\ninvDown\n";
	for (int i = 0; i < iCutoff; i++) { 
		for (int j = 0; j < iCutoff; j++ ){
			cout<<setprecision(16)<<inv_down_matr[i][j]<<"\t";
		}
		cout<<"\n";
	}

	cout<<"\nspinDown:\n";
	for (int i = 0; i < iCutoff; i++) { 
		for (int j = 0; j < iCutoff; j++ ){
			cout<<setprecision(16)<<spin_down_matr[i][j]<<"\t";
		}
		cout<<"\n";
	}
}
//endvimfold


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
