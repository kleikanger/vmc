#include "ipdist.h"
#include "../newmatrix/newmatrix.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <mkl_cblas.h>

using std::cout;

ipdist::ipdist(int n, int di, int iC)
{/*//startvimfold*/
	//storing class variables
	n_min_one=n-1;
	dim=di;
	iCutoff=iC;
	//allocating arrays of length's between particles and corresponding backupmatr	
	//triangular matrix with matrixdimension (n-1)
	ip_len = (double**)tria_matrix(n_min_one,sizeof(double));
	ip_len_backup = (double**)tria_matrix(n_min_one,sizeof(double));
}/*//endvimfold*/
ipdist::~ipdist()
{/*//startvimfold*/
	free_matrix((void **) ip_len);
	free_matrix((void **) ip_len_backup);
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
void ipdist::setBeta(double betaARG)
{/*//startvimfold*/
	beta=betaARG;
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
double ipdist::sumInvlen() const 
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
void ipdist::jasGrad(double** ret_vec, double** r, double** r_old, const int & i_upd) const 
{//startvimfold
	int j,i; //k
	//temporary variable
	double r_kj=0;
	//jastrow par
	double a=0.3333333333333333;
	double temp;
	
	//Only changing the elements that has nev values changed
	int i_upd_mo;
	//i_upd minus one
	i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be updated
	
	double r_kj_old, temp_old;
	
	
	for (j=0;j<dim;j++)
	{
		ret_vec[i_upd][j]=0.0;
	}

	bool spin_up_particle = (i_upd<iCutoff);
	for (i=0;i<i_upd;i++)
	{
		r_kj=ip_len[i_upd_mo][i];
		r_kj_old=ip_len_backup[i_upd_mo][i];

		if ((i<iCutoff)==(spin_up_particle)) 
		{
			temp =  a /(r_kj * pow(1.+beta*r_kj,2));
			temp_old = a /(r_kj_old * pow(1.+beta*r_kj_old,2));
		}
		//particles with antiparallel spins
		else	
		{
			//a=1.0
			temp = 1.0 /(r_kj * pow(1.+beta*r_kj,2));
			temp_old = 1.0 /(r_kj_old * pow(1.+beta*r_kj_old,2));
		}
/*
		a = ((i<iCutoff)==(spin_up_particle)) ? .3333333333333333 : 1;
		temp =  a /(r_kj * pow(1.+beta*r_kj,2));
		temp_old = a /(r_kj_old * pow(1.+beta*r_kj_old,2));
*/
		for (j=0;j<dim;j++)
		{
			ret_vec[i][j]-=temp_old*(r_old[i][j]-r_old[i_upd][j]);
			ret_vec[i][j]+=temp*(r[i][j]-r[i_upd][j]);
		}

		for (j=0;j<dim;j++)
			ret_vec[i_upd][j]+=temp*(r[i_upd][j]-r[i][j]);
		
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		r_kj=ip_len[i][i_upd];
		r_kj_old=ip_len_backup[i][i_upd];
		
		if (((i+1)<iCutoff)==(spin_up_particle)) 
		{
			temp =  a /(r_kj * pow(1.+beta*r_kj,2));
			temp_old = a /(r_kj_old * pow(1.+beta*r_kj_old,2));

		}
		//particles with antiparallel spins
		else	
		{
			//a=1.0
			temp = 1.0 /(r_kj * pow(1.+beta*r_kj,2));
			temp_old = 1.0 /(r_kj_old * pow(1.+beta*r_kj_old,2));
		}
		/*
		a = (((i+1)<iCutoff)==(spin_up_particle)) ? .3333333333333333 : 1;
		temp =  a /(r_kj * pow(1.+beta*r_kj,2));
		temp_old = a /(r_kj_old * pow(1.+beta*r_kj_old,2));
		*/
		for (j=0;j<dim;j++)
		{
			ret_vec[i+1][j]-=temp_old*(r_old[i+1][j]-r_old[i_upd][j]);
			ret_vec[i+1][j]+=temp*(r[i+1][j]-r[i_upd][j]);
		}

		for (j=0;j<dim;j++)
			ret_vec[i_upd][j]+=temp*(r[i_upd][j]-r[i+1][j]);
	}
	//EVT set up temp matrix and use level 3 blas func.
}//endvimfold
double ipdist::jasLapl(double** r) const 
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
#if 0
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
  //	return sum;
#if 1
	//will give same results if we do not sum the particle twice
	return 2.*sum;
#endif
	
}//End function slaterMatrix::jastrow()
//endvimfold
double ipdist::logJasR(const int &i_upd) const 
{/*//startvimfold*/
	double r_12;
	int i;
	double sum=0.0;
	double a=0.3333333333333333;
	bool i_upd_spin_up=(i_upd<iCutoff);
	int i_upd_mo=i_upd-1;
	//n-1 elements in ip_len needs to be summed
	for (i=0;i<i_upd;i++)
	{
		r_12=ip_len[i_upd_mo][i];
		if (i_upd_spin_up==i<iCutoff) //? a=con : a=1.0; 
		//if (i_upd_spin_up==(i<iCutoff))  a=con; else a=1.0; 
		{
			sum+=a*r_12/(1.0 + beta*r_12);
		}
		else
		{
			sum+=r_12/(1.0 + beta*r_12);
		}
	}
	for (i=i_upd;i<n_min_one;i++)
	{
		r_12=ip_len[i][i_upd];
		if (i_upd_spin_up==((i+1)<iCutoff)) //? a=con : a=1.0;
		//if (i_upd_spin_up==((i+1)<iCutoff))  a=con; else a=1.0;
		{
			sum+=a*r_12/(1.0 + beta*r_12);
		}
		else
		{
			sum+=r_12/(1.0 + beta*r_12);
		}

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

