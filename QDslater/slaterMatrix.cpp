/*
	Implementation of class slaterMatrix
	
	Contains:

	NB: XXX Only good for hamiltonians that isn't dependent of spin since D->Dup*Ddown. New variable:
	bool bSeparate_up_down_matrixes.

	NBB: Assumes that nr. of spin up particles equals nr of spin down particles. 
	If number of spin up particles does not equal number of spin down program needs to be 
	modified. New variables: iCutoff_up, iCutoff_down.
 
	NBBB: later initiate orbital objects in solver class, and pass pointers to this
	class as argument when initializing the class. Then only one set of orbital obj.
	needs to be initialized and stored. Aldernative: heritage??
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
	spin_up_matr = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){spin_up_matr[i] = new double[iCutoff];}
	spin_down_matr = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){spin_down_matr[i] = new double[iCutoff];}
	inv_up_matr = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){inv_up_matr[i] = new double[iCutoff];}
	inv_down_matr = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){inv_down_matr[i] = new double[iCutoff];}

	spin_up_backup = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){spin_up_backup[i] = new double[iCutoff];}
	spin_down_backup = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){spin_down_backup[i] = new double[iCutoff];}
	inv_up_backup = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){inv_up_backup[i] = new double[iCutoff];}
	inv_down_backup = new double*[iCutoff];
	for (int i=0; i<iCutoff; i++){inv_down_backup[i] = new double[iCutoff];}
		
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
slaterMatrix::~slaterMatrix(){
//startvimfold
	//delete variational_parameters;
	for (int i=0; i<iCutoff; i++){
		delete [] inv_down_matr[i];
		delete [] inv_up_matr[i];
		delete [] spin_up_matr[i];
		delete [] spin_down_matr[i];
		delete [] inv_down_backup[i];
		delete [] inv_up_backup[i];
		delete [] spin_up_backup[i];
		delete [] spin_down_backup[i];
	}
	//delete grad_up;
	//delete grad_down;
	delete [] inv_down_matr;
	delete [] inv_up_matr;
	delete [] spin_up_matr;
	delete [] spin_down_matr;
	delete [] inv_down_backup;
	delete [] inv_up_backup;
	delete [] spin_up_backup;
	delete [] spin_down_backup;
	delete [] orbital_; //DELETE ORBITAL[i]
}
//endvimfold
void slaterMatrix::setVarPar(double alpha){
	for (int i=0; i<iNumPart; i++)/*//startvimfold*/
	{
		orbital_[i].setAlpha(alpha); 
	}
}/*//endvimfold*/
void slaterMatrix::initSlaterMatrix(double** partPos){
	/*//startvimfold*/
	//XXX OPT: init as 1D array.
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
	//copying initialized matrices to backupmatr.
	for (i=0; i<iCutoff; i++)
	{
		for (j=0; j<iCutoff; j++)
		{
			spin_up_backup[i][j] = spin_up_matr[i][j];
		}
	}
	for (i=0; i<iCutoff; i++)
	{
		for (j=0; j<iCutoff; j++)
		{
			spin_down_backup[i][j] = spin_down_matr[i][j]; 
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
	delete [] B;
    delete [] ipiv;
    delete [] work;

	//init backupmatr
	accept(0);
	accept(iCutoff);

}//End function slatermatrix::updateCofactors()
//endvimfold
void slaterMatrix::update(double* d_R, int i_upd){
//startvimfold	
//XXX OPT: dcopy & daxpy much slower than 2D loop, dscal same, ddot 2x faster	
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
void slaterMatrix::reject(int i_upd)
{/*//startvimfold*/
	//copying entire matrix : inv_up_matr<-inv_up_backup
	int i,j;
	//copying the last updated row : spin_up_matr <- spin_up_backup
	if (i_upd<iCutoff)
	{
		for (i=0;i<iCutoff;i++)
		{
			for (j=0;j<iCutoff;j++)
			{
				inv_up_matr[i][j]=inv_up_backup[i][j];
			}
		}
	//if (i_upd<iCutoff)
		for (i=0;i<iCutoff;i++)
		{
			spin_up_matr[i_upd][i]=spin_up_backup[i_upd][i];
		}
	}
	else
	{
		for (i=0;i<iCutoff;i++)
		{
			for (j=0;j<iCutoff;j++)
			{
				inv_down_matr[i][j]=inv_down_backup[i][j];
			}
		}
		i_upd-=iCutoff;
	//if (i_upd>iCutoff)
		for (i=0;i<iCutoff;i++)
		{
			spin_down_matr[i_upd][i]=spin_down_backup[i_upd][i];
		}
	}
}/*//endvimfold*/
void slaterMatrix::accept(int i_upd)
{/*//startvimfold*/
	int i,j;
	if (i_upd<iCutoff)
	{
	//copying entire matrix : inv_up_backup <- inv_up_matr
		for (i=0;i<iCutoff;i++)
		{
			for (j=0;j<iCutoff;j++)
			{
				inv_up_backup[i][j]=inv_up_matr[i][j];
			}
		}
	//copying the last updated row : spin_up_backup <- spin_up_matr
		for (i=0;i<iCutoff;i++) 
		{
			spin_up_backup[i_upd][i]=spin_up_matr[i_upd][i];
		}
	}
	else
	{	
		for (i=0;i<iCutoff;i++)
		{
			for (j=0;j<iCutoff;j++)
			{
				inv_down_backup[i][j]=inv_down_matr[i][j];
			}
		}
		i_upd-=iCutoff;
		for (i=0;i<iCutoff;i++)
		{
			spin_down_backup[i_upd][i]=spin_down_matr[i_upd][i];
		}
	}
}/*//endvimfold*/
double const slaterMatrix::waveFunction(int i_upd){
	//startvimfold
	//returns ratio between new and old determinant when one particle moved.
	int i;
	//double new_vec[iCutoff];
	
	//double sum=0;
	if (i_upd<iCutoff)
	{
		//for (i=0;i<iCutoff;i++)
		//{
			//sum+=orbital_[i].valueWF(dR)*inv_up_matr[i_upd][i];
			//sum+=spin_up_matr[i_upd][i]*inv_up_backup[i_upd][i];
		//}
		return cblas_ddot(iCutoff,spin_up_matr[i_upd],1,inv_up_backup[i_upd],1);
		//Det down not changed: ratio =1;
		//return sum;
	} 
	else 
	{
		//for (i=iCutoff;i<iNumPart;i++)
		//{
			//sum+=orbital_[i].valueWF(dR)*inv_down_matr[i_upd-iCutoff][i-iCutoff];
			//sum+=spin_down_matr[i_upd-iCutoff][i-iCutoff]*inv_down_backup[i_upd-iCutoff][i-iCutoff];
		//}
		return cblas_ddot(iCutoff,spin_down_matr[i_upd-iCutoff],1,inv_down_backup[i_upd-iCutoff],1);
		//Det up not changed: ratio =1;
		//return sum;
	}

}
//endvimfold
void const slaterMatrix::grad(double** ret_vec, double** dR, int active_part)//, int axis, int i_upd)
{//startvimfold
//OPTIMALIZATION: only one some grads needs to be updated.
//Grad_{iaxis}. Full gradient/DET_old = Sum_axis Sum_i Grad_{i,axis} \vec e_{i,axis}.
//OPT two loops instead of one

//TEST

	int i, axis, i_upd;
	double d_upd[iCutoff];

	//Only update 'active' matrix
	if (active_part<iCutoff)
	{
		for (i_upd=0;i_upd<iCutoff;i_upd++)
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
	}	
	else
	{	
		for (i_upd=iCutoff;i_upd<iNumPart;i_upd++)
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
double const slaterMatrix::lapl(double** dR){
//startvimfold
	//OPTIMALIZATION: only one of the matrices needs to be calculated.
	//probably faster method than using double loops
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


// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
