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
	
	cout<<"x "<<setprecision(16)<<a.D1(r,0)<<"\n";
	cout<<"y "<<setprecision(16)<<a.D1(r,1)<<"\n";
	cout<<"z "<<setprecision(16)<<a.D1(r,2)<<"\n";

	cout<<"x "<<setprecision(16)<<a.D2(r)<<"\n";

	delete r;
}
