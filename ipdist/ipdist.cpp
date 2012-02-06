#include "ipdist.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cblas.h>

using std::cout;

ipdist::ipdist(int n, int di, int iC)
{/*//startvimfold*/
	//storing class variables
	n_min_one=n-1;
	dim=di;
	iCutoff=iC;
	//allocating arrays of length's between particles
	ip_len = new double*[n_min_one];
	for (int i=0; i<n_min_one; i++) ip_len[i] = new double[i+1];
	ip_invlen = new double*[n_min_one];
	for (int i=0; i<n_min_one; i++) ip_invlen[i] = new double[i+1];
}/*//endvimfold*/
void ipdist::init(double** r)
{/*//startvimfold*/
	double temp;
	int i,j,k;
	//Initializing lower triag. array of r_ij
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{
			temp=0;
			for (k=0; k<dim; k++)
			{
				temp+=(r[i+1][k]-r[j][k])*(r[i+1][k]-r[j][k]);
			}
			temp=sqrt(temp);
			ip_len[i][j]=temp;
		}
	}
	//Initializing lower triag. array of 1/r_ij
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{	
			ip_invlen[i][j]=1/ip_len[i][j];
		}
	}
}/*//endvimfold*/
void ipdist::update(double* r, int i_upd)
{/*//startvimfold*/
	int i, i_upd_mo;
	//i_upd minus one
	i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be updated
	for (i=0;i<i_upd;i++)
	{
		ip_len[i_upd_mo][i]=r[i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		ip_len[i][i_upd]=r[i+1];
	}
	//n-1 elements in ip_invlen needs to be updated
	for (i=0;i<i_upd;i++)
	{
		ip_invlen[i_upd_mo][i]=1/r[i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		ip_invlen[i][i_upd]=1/r[i+1];
	}
}/*//endvimfold*/
const double ipdist::sumPart(int i_upd)
{/*//startvimfold*/
	double sum=0.0;
	int i, i_upd_mo;
	//i_upd minus one
	i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be summed
	for (i=0;i<i_upd;i++)
	{
		sum+=ip_len[i_upd_mo][i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		sum+=ip_len[i][i_upd];
	}
	return sum;
}/*//endvimfold*/
const double ipdist::sumInvlen()
{/*//startvimfold*/
	double sum=0;
	int i,j;
	//summing all inverse elements 1/r_ij
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{	
			sum+=ip_invlen[i][j];
		}
	}
	return sum;
}/*//endvimfold*/
const void ipdist::jasGrad(double** ret_vec, double beta, double** r)
{//startvimfold
	double r_kj;
	int j,k,J;
   	//summation variable
	double sum = 0.0;
	//temporary array
	double temp_arr[dim];
	//jastrow par
	double a=0.3333333333333333;
	double temp;

//XXX XXX XXX XXX
//NB: CALC GRAD FOR SUB DIAG ELEMENTS.
//THE TOTAL GRAD SHOULD BE MULT. BY 2????

	//set all enements in ret_vec = 0
	//for (int i=0;i<((n_min_one+1)*dim);i++) { cout<<(*(ret_vec)+i)<<"\n"; }
	for (k=0;k<n_min_one+1;k++)	
		for (j=0;j<dim;j++)
			ret_vec[k][j]=0;	
	
	//Sum over all particles.	
	for (k=0;k<=n_min_one;k++)
	{
		for (j=0;j<k;j++)
		{       
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[k-1][j]; 
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],1,temp_arr,1);
			//temp_arr< r[k]-r[j]
			cblas_daxpy(dim,-1.0,r[j],1,temp_arr,1);
			//Particles with parallel spins:
			if (j<iCutoff==k<iCutoff) 
			{
				temp = a /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//cblas_dscal(dim, temp, temp_arr, 1);
			//cblas_daxpy(dim, 1, temp_arr, 1, ret_vec[k], 1);
			cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1);
		}
		for (j=k;j<n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j][k];
			J=j+1; 
			
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],1,temp_arr,1);
			//temp_arr< r(k)-r[j]
			cblas_daxpy(dim,-1.0,r[J],1,temp_arr,1);
				
			//Particles with parallel spins:
			if (k<iCutoff==J<iCutoff) 
			{
				temp = a /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//cblas_dscal(dim, temp, temp_arr, 1);
			//cblas_daxpy(dim, 1, temp_arr, 1, ret_vec[k], 1);
			cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1);
		}
	}
}//endvimfold
const double ipdist::jasLapl(double beta, double** r)
{//startvimfold
	//summation variables
	int i,j,k,ii,jj;
	//temporary variables
	double temp;
	double r_kj, r_ki;
	double ove_ki, ove_kj;
   	//summation variable
	double sum_1 = 0;
	double sum_2 = 0;
	//temporary array
	double temp_arr_i[dim];
	double temp_arr_j[dim];
	//jastrow parameter
	double a=0.3333333333333333;

//XXX XXX XXX
//CALC LAPL FOR SUB DIAG ELEMENTS.
//TOTAL LAPL MULT. BY 2?

	for (k=0;k<=n_min_one;k++)
	{
		for (j=0;j<k;j++)
		{   
			//pick the correct matrix element	
			r_kj=ip_len[k-1][j];
			
			for (i=0;i<k;i++)
			{	
				r_ki=ip_len[k-1][i];

				temp=0.;
				//temp_arr<-r[]
				cblas_dcopy(dim,r[k],1,temp_arr_i,1);
				cblas_dcopy(dim,r[k],1,temp_arr_j,1);
				//temp_arr<r[k]-temp
				cblas_daxpy(dim,-1.0,r[i],1,temp_arr_i,1);
				cblas_daxpy(dim,-1.0,r[j],1,temp_arr_j,1);
				//(r_k-r_j).(r_k-r_i)
				temp=cblas_ddot(dim,temp_arr_i,1,temp_arr_j,1);
				temp/=r_ki*r_kj; //*ipd_invlen
				//temporary vars
				ove_ki=1./(1+beta*r_ki);
				ove_kj=1./(1+beta*r_kj);
				temp*=ove_ki*ove_ki*ove_kj*ove_kj;
				sum_1+=temp;
			}
			//Particles with parallel spins:
			if (i<iCutoff==j<iCutoff) 
			{
				sum_2 +=   a*ove_kj * ove_kj / r_kj ;
				sum_2 -=   a*beta * ove_kj * ove_kj * ove_kj;
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum_2 +=   ove_kj * ove_kj / r_kj ;
				sum_2 -=   beta * ove_kj * ove_kj * ove_kj;
			}
		}
		for (j=k;j<n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j][k];
			jj=j+1; 
	
			for (i=k;i<n_min_one;i++)
			{	
				//pick the correct element from the ip_len matrix
				r_ki=ip_len[i][k];
				ii=i+1; 
				
				temp=0;
				//temp_arr<-r[]
				cblas_dcopy(dim,r[k],1,temp_arr_i,1);
				cblas_dcopy(dim,r[k],1,temp_arr_j,1);
				//temp_arr<r[k]-temp
				cblas_daxpy(dim,-1.0,r[ii],1,temp_arr_i,1);
				cblas_daxpy(dim,-1.0,r[jj],1,temp_arr_j,1);
				//(r_k-r_j).(r_k-r_i)
				temp=cblas_ddot(dim,temp_arr_i,1,temp_arr_j,1);
				temp/=r_ki*r_kj;
				//temporary vars
				ove_ki=1./(1+beta*r_ki);
				ove_kj=1./(1+beta*r_kj);
				temp*=ove_ki*ove_ki*ove_kj*ove_kj;
				sum_1+=temp;
			}
			//Particles with parallel spins:
			if (i<iCutoff==j<iCutoff) 
			{
				sum_2 +=   a * ove_kj * ove_kj / r_kj ;
				sum_2 -=   a * beta * ove_kj * ove_kj * ove_kj;
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum_2 +=   ove_kj * ove_kj / r_kj ;
				sum_2 -=   beta * ove_kj * ove_kj * ove_kj;
			}
		}
	}

	sum_1*=a*a;
	sum_2*=2.;

    return sum_1+sum_2;
}//End function slaterMatrix::jastrow()
//endvimfold
void ipdist::clear()
{/*//startvimfold*/
	for (int i=0; i<n_min_one; i++) delete ip_len[i];
	delete ip_len;
	for (int i=0; i<n_min_one; i++) delete ip_invlen[i];
	delete ip_invlen;
}/*//endvimfold*/
void ipdist::print()
{/*//startvimfold*/
	int i,j;
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{	
			//cout<<"\t\t"<<ip_invlen[i][j];
			cout<<"\t\t"<<ip_len[i][j];
		}
		cout<<"\n";
	}
	cout<<"\n";
}/*//endvimfold*/
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold

