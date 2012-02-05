#include "ipdist.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

using std::cout;

ipdist::ipdist(int n, int di)
{/*//startvimfold*/
	//storing class variables
	n_min_one=n-1;
	dim=di;
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
	//n-1 elements in ip_len needs to be updated
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
