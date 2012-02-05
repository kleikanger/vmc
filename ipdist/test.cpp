#include "ipdist.h"
#include "cmath"

#include <iostream>
using std::cout;

main()
{

	//TEST CLASS ipdist
	//calculating distances directly and via update
	int n=11;
	int dim=1;
	int i_upd=0;

	//Initializing a matrix
	int i,k;
	double ** r = new double*[n];
	for (int i=0;i<n;i++)
	{
		r[i] = new double[dim];
	   for (int j=0;j<dim;j++)
	   {
		   r[i][j] = (i+3)*(j+1);
	   }	   
	}

	//init ipdist objects
	ipdist a(n,dim);
	ipdist b(n,dim);

	a.init(r);
	b.init(r);

	//test:
	cout<<"\nprint init matrix\n";
	a.print();
	cout<<"\ncharge r, particle"<<i_upd<<"\n";

	//updating r, position i_upd
	for (int i=0;i<dim;i++)
	{
	   r[i_upd][i] = (i_upd+1)/(i+1);
	}	   

	//calculating new lengths betw particles
	double upd_vec[n-1];
	for (i=0; i<i_upd; i++)
	{
		double temp=0;
		for (k=0; k<dim; k++)
		{
			temp+=(r[i_upd][k]-r[i][k])*(r[i_upd][k]-r[i][k]);
		}
		upd_vec[i]=sqrt(temp);
	}
	for (i=i_upd+1; i<n; i++)
	{
		double temp=0;
		for (k=0; k<dim; k++)
		{
			temp+=(r[i_upd][k]-r[i][k])*(r[i_upd][k]-r[i][k]);
		}
		upd_vec[i]=sqrt(temp);
	}

	a.update(upd_vec,i_upd);
	b.init(r);

	cout<<"\nprint updated matrix\n";
	a.print();
	cout<<"\nprint init matrix\n";
	b.print();

	cout<<"sumInvlen: "<<a.sumInvlen()<<"\n";
	for (int s=0; s<n; s++)
	cout<<"sumPart("<<s<<"):"<<a.sumPart(s)<<"\n";

}
