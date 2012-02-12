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
		for (j=0; j<=i; j++)
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
		for (j=0; j<=i; j++)
		{	
			ip_invlen[i][j]=1.0/ip_len[i][j];
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
		ip_len[i][i_upd]=r[i];
	}
	//n-1 elements in ip_invlen needs to be updated
	for (i=0;i<i_upd;i++)
	{
		ip_invlen[i_upd_mo][i]=1.0/r[i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		ip_invlen[i][i_upd]=1.0/r[i];
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
	double sum=0.0;
	int i,j;
	//summing all inverse elements 1/r_ij
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<=i; j++)
		{	
			sum+=ip_invlen[i][j];
		}
	}
	return sum;
}/*//endvimfold*/
const void ipdist::jasGrad(double** ret_vec, double beta, double** r)
{//startvimfold
	int j,k,J,kk;
	//temporary variable
	double r_kj=0;
	//temporary array
	double temp_arr[dim];
	//jastrow par
	double a=0.3333333333333333;
	double temp;

	//set all enements in ret_vec = 0
	//for (int i=0;i<((n_min_one+1)*dim);i++) { cout<<(*(ret_vec)+i)<<"\n"; }
	for (k=0;k<=n_min_one;k++)	
		for (j=0;j<dim;j++)
			ret_vec[k][j]=0.0;	
	
	//XXX OPT XXX ONLY ONE PART NEEDS TO BE UPDATED
	//INPUT i_upd, set k=i_upd

	//Sum over all particles.	
	for (k=0;k<=n_min_one;k++)
	{
		for (j=0;j<k;j++)
		{       
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[k-1][j]; //k>=1 since 0<=j<k 
			
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],1,temp_arr,1);
			//temp_arr< r[k]-r[j]
			cblas_daxpy(dim,-1.0,r[j],1,temp_arr,1);
			//Particles with parallel spins:
			if ((j<iCutoff)==(k<iCutoff)) 
			{
				temp =  a /r_kj /(1.+beta*r_kj) /(1.+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /r_kj /(1.+beta*r_kj) /(1.+beta*r_kj);
			}
			cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1);
		}
		for (j=k+1;j<=n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j-1][k];

			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],1,temp_arr,1);
			//temp_arr< r(k)-r[j]
			cblas_daxpy(dim,-1.0,r[j],1,temp_arr,1);
				
			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				temp =  a /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1.0 /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1);
		}
	}
}//endvimfold
const double ipdist::jasLapl(double beta, double** r)
{//startvimfold
	//summation indexes
	int j,k;
	//temporary variable
	double r_kj;
   	//summation variable
	double sum = 0.0;
	//jastrow parameter
	double a=0.3333333333333333;

	for (k=0;k<=n_min_one;k++)
	{
		for (j=0;j<k;j++)
		{   
			//pick the correct matrix element	
			r_kj=ip_len[k-1][j];
			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				sum += a * ( 1.0 - beta*r_kj ) / r_kj / ( pow((1.0+r_kj*beta),3) );
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum +=  ( 1.0 - beta*r_kj ) / r_kj / ( pow((1.0+r_kj*beta),3) );
			}
		}
		for (j=k+1;j<=n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j-1][k];
			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				sum +=  a*( 1.0 - beta*r_kj ) / r_kj / ( pow((1.0+r_kj*beta),3) );
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum +=  ( 1.0 - beta*r_kj ) / r_kj / ( pow((1.0+r_kj*beta),3) );
			}
		}
	}
  return sum;
	
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

