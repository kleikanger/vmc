/*

Test Laplacian and grad.

   */




#include "orbital.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cblas.h>

using std::cout;
using std::cerr;
using std::setprecision;

int main(){

	orbital a(4,0,true);
	double* r = new double[3];

	for (int i=0; i<3; i++)
	{ 
		r[i]=i+1;
	}
	
	cout<<"grad x "<<setprecision(16)<<a.D1(r,0)<<"\n";
	cout<<"grad y "<<setprecision(16)<<a.D1(r,1)<<"\n";
	cout<<"grad z "<<setprecision(16)<<a.D1(r,2)<<"\n";
	cout<<"lapl. "<<setprecision(16)<<a.D2(r)<<"\n";

	delete r;
}
