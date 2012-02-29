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
	//and corresponding backupmatr	
	ip_len_backup = new double*[n_min_one];
	for (int i=0; i<n_min_one; i++) ip_len_backup[i] = new double[i+1];
}/*//endvimfold*/
ipdist::~ipdist()
{/*//startvimfold*/
	for (int i=0; i<n_min_one; i++) delete ip_len[i];
	delete ip_len;
	for (int i=0; i<n_min_one; i++) delete ip_len_backup[i];
	delete ip_len_backup;
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
	//Initializing backup matrices
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<=i; j++)
		{
			ip_len_backup[i][j]=ip_len[i][j];
		}
	}
}/*//endvimfold*/
void ipdist::setBeta(double ibeta)
{/*//startvimfold*/
	beta=ibeta;
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
}/*//endvimfold*/
void ipdist::accept(int i_upd)
{/*//startvimfold*/
	int i, i_upd_mo;
	//i_upd minus one
	i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be updated
	for (i=0;i<i_upd;i++)
	{
		ip_len_backup[i_upd_mo][i]=ip_len[i_upd_mo][i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		ip_len_backup[i][i_upd]=ip_len[i][i_upd];
	}
}/*//endvimfold*/
void ipdist::reject(int i_upd)
{/*//startvimfold*/
	int i, i_upd_mo;
	//i_upd minus one
	i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be updated
	for (i=0;i<i_upd;i++)
	{
		ip_len[i_upd_mo][i]=ip_len_backup[i_upd_mo][i];
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		ip_len[i][i_upd]=ip_len_backup[i][i_upd];
	}
}/*//endvimfold*/
double const ipdist::sumInvlen()
{/*//startvimfold*/
	double sum=0.0;
	int i,j;
	//summing all inverse elements 1/r_ij
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<=i; j++)
		{	
			//sum+=ip_invlen[i][j];
			sum+=1./ip_len[i][j];
		}
	}
	return sum;
}/*//endvimfold*/
void const ipdist::jasGrad(double** ret_vec, double** r)
{//startvimfold
	int j,k,i;
	//temporary variable
	double r_kj=0;
	//jastrow par
	double a=0.3333333333333333;
	double temp, sum;
	//set all enements in ret_vec = 0
	//XXX OPT: take away this summation. 
	//sum as in mortens code
	for (k=0;k<=n_min_one;k++)	
		for (j=0;j<dim;j++)
		{
			ret_vec[k][j]=0.0;
		}

	//Sum over all particles in j & k. j!=k	
	for (k=0;k<=n_min_one;k++)
	{
		for (j=0;j<k;j++)
		{       
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[k-1][j]; 

			//Particles with parallel spins:
			if ((j<iCutoff)==(k<iCutoff)) 
			{
				temp =  a /(r_kj * pow(1.+beta*r_kj,2));
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1.0 /(r_kj * pow(1.+beta*r_kj,2));
			}
			for (i=0;i<dim;i++)
			{
				ret_vec[k][i]+=temp*(r[k][i]-r[j][i]);
			}
				
		}
		for (j=k+1;j<=n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j-1][k];

			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				temp =  a / (r_kj * pow(1+beta*r_kj,2));
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1.0 / (r_kj * pow(1+beta*r_kj,2));
			}
			for (i=0;i<dim;i++)
			{
				ret_vec[k][i]+=temp*(r[k][i]-r[j][i]);
			}
		}
	}
	//EVT set up temp matrix and use level 3 blas func.
}//endvimfold
double const ipdist::jasLapl(double** r)
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
#if 1
		for (j=0;j<k;j++)
		{   
			//pick the correct matrix element	
			r_kj=ip_len[k-1][j];
			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				sum += a * ( 1.0-beta*r_kj ) / ( r_kj * pow((1.0+r_kj*beta),3) );
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum +=  ( 1.0-beta*r_kj ) / ( r_kj * pow((1.0+r_kj*beta),3) );
			}
		}
#endif
		for (j=k+1;j<=n_min_one;j++)
		{
			//pick the correct element from the ip_len matrix
			r_kj=ip_len[j-1][k];
			//Particles with parallel spins:
			if ((k<iCutoff)==(j<iCutoff)) 
			{
				sum +=  a*( 1.0-beta*r_kj ) / ( r_kj * pow((1.0+r_kj*beta),3) );
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum +=  ( 1.0-beta*r_kj ) / ( r_kj * pow((1.0+r_kj*beta),3) );
			}
		}
	}
  	return sum;
#if 0
	//will give same results if we do not sum the particle twice
	return 2.*sum;
#endif
	
}//End function slaterMatrix::jastrow()
//endvimfold
double const ipdist::logJasR(int i_upd)
{/*//startvimfold*/
	double r_12, a;
	int i;
	double sum=0.0;
	double con=0.3333333333333333;
	bool i_upd_spin_up=(i_upd<iCutoff);
	int i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be summed
	for (i=0;i<i_upd;i++)
	{
		(i_upd_spin_up==i<iCutoff) ? a=con : a=1.0; 
		//if (i_upd_spin_up==(i<iCutoff))  a=con; else a=1.0; 
		r_12=ip_len[i_upd_mo][i];
		sum+=a*r_12/(1.0 + beta*r_12);
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		(i_upd_spin_up==((i+1)<iCutoff)) ? a=con : a=1.0;
		//if (i_upd_spin_up==((i+1)<iCutoff))  a=con; else a=1.0;
		r_12=ip_len[i][i_upd];
		sum+=a*r_12/(1.0 + beta*r_12);
	}
	return sum;
}/*//endvimfold*/
void ipdist::print()
{/*//startvimfold*/
	int i,j;
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{	
			cout<<"\t\t"<<ip_len[i][j];
		}
		cout<<"\n";
	}
	cout<<"\n";
}/*//endvimfold*/
// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold

