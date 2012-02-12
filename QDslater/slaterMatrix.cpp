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

//LAPACK functions.
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int*, int*, double*, int*, int*, int*);
    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int*, double*, int*, int*, double*, int*, int*);
};
	
//constructor
slaterMatrix::slaterMatrix(int iNp,int iCo,int inovp, int di){
	/*//startvimfold*/

	dim=di;
	iNumPart=iNp;
	iCutoff=iCo;
	variational_parameters = new double[inovp];
	iNumber_of_variational_parameters=inovp;
	
	//Allocating new matrices and arrays
	//grad_up = new double[dim];
	//grad_down = new double[dim];

	inv_up_matr = new double*[iCutoff];
	inv_down_matr = new double*[iCutoff];
	spin_up_matr = new double*[iCutoff];
	spin_down_matr = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){
		inv_down_matr[i] = new double[iCutoff];
		inv_up_matr[i] = new double[iCutoff];
		spin_up_matr[i] = new double[iCutoff];
		spin_down_matr[i] = new double[iCutoff];
	}
	
	//Initializing orbital objects
	orbital_ = new orbital[iNumPart];
	//States 0-(iCutoff-1) enters inv_up_matr. 
	//States iCutoff-(iNumPart-1) enters inv_down_matr 
	for (int i=0; i<iCutoff; i++)
		orbital_[i].setValues(i+1,0,true,dim);//n=i,l=0,spin up.
	for (int i=iCutoff; i<iNumPart; i++)
		orbital_[i].setValues(i-iCutoff+1,0,false,dim);//n=i,l=0,spin down.

}//end constructor varMC()
/*//endvimfold*/
void slaterMatrix::initSlaterMatrix(double** partPos){
	/*//startvimfold*/
	int i,j;
	for (i=0; i<iCutoff; i++)
	{
		for (j=0; j<iCutoff; j++)
		{
			spin_up_matr[i][j] = orbital_[j].valueWF(partPos[i]); // Transpose
		}
	}
	for (i=0; i<iCutoff; i++)
	{
		for (j=0; j<iCutoff; j++)
		{
			spin_down_matr[i][j] = orbital_[j+iCutoff].valueWF(partPos[i+iCutoff]); 
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
	
	for (int h=0; h<2; h++)
	{
		//Lapack using coloumn major
		//Matrixes must be transposed before they are sent to inverse
		//Better to initiate transpose matrix at once !!
		int k=0;
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				if (h==0) B[k+j] = spin_up_matr[j][i];
				else B[k+j] = spin_down_matr[j][i];
			}
		k+=n;
		}
			
		dgetrf_(&n,&n,B,&n,ipiv,&info);
		//maybe check if matrix is singular...
		if (info==0) 
		{
			dgetri_(&n,B,&n,ipiv,work,&lwork,&info);
		}
		if (info!=0) 
		{
			cout<<"error slaterMatrix::findInverse in lapack dgetri_ or degetrf_ info="<<info<<"\n";
			exit(1);
		}
		k=0;
		//inverse
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
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
void slaterMatrix::update(double* d_R, int i_upd){
//startvimfold	
	
	int k,j;
	double d_upd[iCutoff];
	double R, temp;
    
	if (i_upd<iCutoff)
	{
	//SPINUP 
		//find D^new/D^old
		for (k=0;k<iCutoff;k++)
		{
				d_upd[k]=orbital_[k].valueWF(d_R);
		}
		//Updating ingerse Blas routine : d_up*inv_up_matr[j][:], 
		//remember inverse M are transposed
		R = cblas_ddot(iCutoff,d_upd,1,inv_up_matr[i_upd],1);

		for (j=0;j<iCutoff;j++)
		{
			if (i_upd!=j)
			{
				temp = cblas_ddot(iCutoff,d_upd,1,inv_up_matr[j],1);
				cblas_daxpy( iCutoff,-temp/R,inv_up_matr[i_upd],1,inv_up_matr[j], 1);
			}
		}
		//inv_up_matr[i_upd][:]<-inv_up_matr[i_upd][:]/R
		cblas_dscal(iCutoff,1/R,inv_up_matr[i_upd],1);
		//Copy d_upd to spin_down_matr, updating spin_up_matr.
		cblas_dcopy(iCutoff, d_upd, 1, spin_up_matr[i_upd], 1);
	} 
	else 
	{
	//repeat all for SPINDOWN if i_upd>iCutoff.
		for (k=0;k<iCutoff;k++)
		{
				d_upd[k]=orbital_[k+iCutoff].valueWF(d_R);
		}
		R = cblas_ddot(iCutoff,d_upd,1,inv_down_matr[i_upd-iCutoff],1);

		for (j=0;j<iCutoff;j++)
		{
			if ((i_upd-iCutoff)!=j)
			{
				temp = cblas_ddot(iCutoff,d_upd,1,inv_down_matr[j],1);
				cblas_daxpy(iCutoff,-temp/R,inv_down_matr[i_upd-iCutoff],1,inv_down_matr[j],1);
			}
		}
		cblas_dscal(iCutoff,1/R,inv_down_matr[i_upd-iCutoff],1);
		cblas_dcopy(iCutoff, d_upd, 1, spin_down_matr[i_upd-iCutoff], 1);
	}	
}
//endvimfold
const double slaterMatrix::waveFunction(double* dR, int i_upd){
	//startvimfold
	//returns ratio between new and old determinant when one particle moved.
	int i;
	double new_vec[iCutoff];
	
	//RATIO DET_new/DET_old
	if (i_upd<iCutoff)
	{
		for (i=0;i<iCutoff;i++)
		{
			new_vec[i]=orbital_[i].valueWF(dR);
		}
		//Det down not changed: ratio =1;
		return cblas_ddot(iCutoff,inv_up_matr[i_upd],1,new_vec,1);
	} 
	else 
	{
		for (i=iCutoff;i<iNumPart;i++)
		{
			new_vec[i-iCutoff]=orbital_[i].valueWF(dR);
		}
		//Det up not changed: ratio =1;
		return cblas_ddot(iCutoff,inv_down_matr[i_upd-iCutoff],1,new_vec,1);
	}
}
//endvimfold
//test! 
const void slaterMatrix::grad(double** ret_vec, double** dR)//, int axis, int i_upd)
{//startvimfold
//OPTIMALIZATION: only one some grads needs to be updated.
//Grad_{i,axis}. Full gradient/DET_old = Sum_axis Sum_i Grad_{i,axis} \vec e_{i,axis}.
//OPT two loops instead of one

//TEST

	int i, axis, i_upd;
	double d_upd[iCutoff];

	for (i_upd=0;i_upd<iNumPart;i_upd++)
	{
		if (i_upd<iCutoff)
		{
			for (axis=0;axis<dim;axis++)
			{
				for (i=0;i<iCutoff;i++)
				{
					//temp1+=orbital_[i].D1(dR,axis)*inv_up_matr[i_upd][i];
					d_upd[i]=orbital_[i].D1(dR[i_upd],axis);
				}
			//grad_up[axis]=cblas_ddot(iCutoff,d_upd,1,inv_up_matr[i_upd],1);
			ret_vec[i_upd][axis] = cblas_ddot(iCutoff,d_upd,1,inv_up_matr[i_upd],1);
			}
		}	
		else 
		{
			for (axis=0;axis<dim;axis++)
			{
				for (i=iCutoff;i<iNumPart;i++)
				{
					//temp2+=orbital_[i].D1(dR,axis)*inv_down_matr[i_upd-iCutoff][i-iCutoff];
					d_upd[i-iCutoff]=orbital_[i].D1(dR[i_upd],axis);
				}
			//grad_down[axis]=cblas_ddot(iCutoff,d_upd,1,inv_down_matr[i_upd-iCutoff],1);
			ret_vec[i_upd][axis] = cblas_ddot(iCutoff,d_upd,1,inv_down_matr[i_upd-iCutoff],1);
			}
		}
	}
}
//endvimfold
const double slaterMatrix::lapl(double** dR){
//startvimfold
	//OPTIMALIZATION: only one of the matrices needs to be calculated.
	double d_upd[iCutoff];
	double temp=0.;
	int i,j;
	for (j=0;j<iCutoff;j++)
	{		
		for (i=0;i<iCutoff;i++)
		{
			d_upd[i]=orbital_[i].D2(dR[j]);
		}
		temp += cblas_ddot(iCutoff,d_upd,1,inv_up_matr[j],1);
	}	
	for (j=0;j<iCutoff;j++)
	{
		for (i=0;i<iCutoff;i++)
		{
			d_upd[i]=orbital_[i+iCutoff].D2(dR[j+iCutoff]);
		}
		temp += cblas_ddot(iCutoff,d_upd,1,inv_down_matr[j],1);
	}
	return temp;
}//End function 
//endvimfold
void slaterMatrix::clear(){
//startvimfold
	//delete variational_parameters;
	for (int i=0; i<iCutoff; i++){
		delete inv_down_matr[i];
		delete inv_up_matr[i];
		delete spin_up_matr[i];
		delete spin_down_matr[i];
	}
	//delete grad_up;
	//delete grad_down;
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
//NOT USED. CAN BE REMOVED.
void slaterMatrix::updateSlaterMatrix(double* partPos, int i_upd){
	/*//startvimfold*/
	int i;
	if (i_upd<iCutoff)
	{
		for (i=0; i<iCutoff; i++)
		{
			spin_up_matr[i_upd][i] = orbital_[i].valueWF(partPos); // Transpose
		}
	} 
	else 
	{
		for (i=0; i<iCutoff; i++)
		{
			spin_down_matr[i_upd-iCutoff][i] = orbital_[i+iCutoff].valueWF(partPos); // Transpose
		}
	}	
}//End function varMC::DeterminantMatrix()
/*//endvimfold*/
//MOVauE TO ORBITAL CLASS
void slaterMatrix::updateVariationalParameters(double* vp){
//startvimfold
	for (int i=0;i<iNumber_of_variational_parameters;i++){
		variational_parameters[i]=vp[i];
	}
	//for..
	//orbital_[i].updateVarPar(...)
}
//endvimfold


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
