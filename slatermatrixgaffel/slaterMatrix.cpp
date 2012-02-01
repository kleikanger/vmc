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

	//More elegant way to initialize class variablesi?
	dim=3;
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
			spin_down_matr[i][j] = orbital_[j+iCutoff].valueWF(partPos[i+iCutoff]); //(i,pParticlePositions[j+iCutoff],dAlpha);
		}
	}

			//inv_up_matr[0][0]=1.0;
			//inv_down_matr[0][0]=1.0;
			//spin_up_matr[0][1]=2.0;
			//spin_down_matr[0][1]=2.0;
			//inv_up_matr[1][0]=3.0;
			//inv_down_matr[1][0]=3.0;
			//spin_up_matr[1][1]=4.0;
			//spin_down_matr[1][1]=4.0;

}//End function varMC::DeterminantMatrix()
/*//endvimfold*/
//ONLY FOR TESTING.
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
		//Lapack using coloumn major!
		//Matrixes must be transposed before they are sent to inverse
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
		//Even if we know that that the determinant is !=0, maybe we should check if matrix is singular or near to...
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
void slaterMatrix::updateInverse(double* d_R, int i_upd){
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
		//Copy d_upd to spin_down_matr, updating spinupmatr.
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
//returns ratio between new and old determinant when one particle moved.
//test
double slaterMatrix::waveFunction(double* dR, int i_upd){
	//startvimfold
	int i;
	double det=0;
	double new_vec[iCutoff];
	
	//RATIO DET_new/DET_old
	if (i_upd<iCutoff)//evt if(orbital_[iCofac_column].spinUp();)
	{
		for (i=0;i<iNumPart;i++)
		{
			new_vec[i]=orbital_[i].valueWF(dR);
		}
		det=cblas_ddot(iCutoff,inv_up_matr[i_upd],1,new_vec,1);
		//Det down not changed: ratio =1;
		return det;
		//det*=cblas_ddot(iCutoff,inv_down_matr[i_upd],1,spin_down_matr[i_upd],1);
	} 
	else 
	{
		for (i=iCutoff;i<iNumPart;i++)
		{
			new_vec[i-iCutoff]=orbital_[i].valueWF(dR);
		}
		det=cblas_ddot(iCutoff,inv_down_matr[i_upd-iCutoff],1,new_vec,1);
		//Det up not changed: ratio =1;
		return det;
		//det_up=cblas_ddot(iCutoff,inv_up_matr[i_upd-iCutoff],1,spin_up_matr[i_upd-iCutoff],1);
	}
}
//endvimfold
//test
double slaterMatrix::grad(double** dR, int axis){
//startvimfold
//Gradient along some axis w.r.t one particle / DET_old.
//Grad_i. Full gradient/DET_old = Sum_i Grad_i.
//EASY TO IMPLEMENT cblas_ddot() and temporary vector.
	double temp1=0.;
	//double temp2=0.;
	int i,j;
	double d_upd[iCutoff];

	//SINCE ONLY ONE PARTICLE IS UPDATED, IT IS ENOUGH TO CALCULATE ONE OF THE GRADS
	for (j=0;j<iCutoff;j++)
	{
	//if (i_upd<iCutoff){
		for (i=0;i<iCutoff;i++)
		{
			//temp1+=orbital_[i].D1(dR[j],axis)*inv_up_matr[j][i];
			d_upd[i]=orbital_[i].D1(dR[j],axis);
		}
		temp1=cblas_ddot(iCutoff,d_upd,1,inv_up_matr[j],1);
	}
	for (j=iCutoff; j<iNumPart; j++)
	{
		for (i=iCutoff;i<iNumPart;i++)
		{
			//temp2+=orbital_[i].D1(dR[j],axis)*inv_down_matr[j-iCutoff][i-iCutoff];
			d_upd[i-iCutoff]=orbital_[i].D1(dR[j],axis);
		}
		temp1*=cblas_ddot(iCutoff,d_upd,1,inv_down_matr[j-iCutoff],1);
	//XXX change to + when including jastrow	
	}
	return temp1;//*temp2;
//	}
}
//endvimfold
double slaterMatrix::lapl(double** dR){
//startvimfold

//SUM OF THE LAPLACIANS FOR up and down natr
//Laplacian w.r.t. i'th particle/DET_old (L_i f({r})).
//Full laplacian/DET/old = Sum_i L_i f({r}).
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
    for (int i = 0; i < iNumPart - 1; i++) 
	{
        for (int j = i + 1; j < iNumPart; j++) 
		{
            // Beregner avstand mellom partikkel i og j:
            length=0.0;
				for (a=0; a<3; a++)
				{
					length+= (r[i][a]-r[j][a])*(r[i][a]-r[j][a]);
				}
			length=sqrt(length);
            // Dersom partikkel i og j har samme spinn:
            if (orbital_[i].spinUp() == orbital_[j].spinUp()) 
			{
                argument += 0.25*length/(1 + alpha*length);
            }
            // Dersom partikkel i og j har ulike spinn:
            else 
			{
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
	//for (int i=0; i<iNumPart; i++)
	//	delete orbital_[i];
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
