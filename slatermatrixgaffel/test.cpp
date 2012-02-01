#include <cblas.h>
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"


///*TEST OF PROGRAM
//startvimfold

#include <iomanip>
int main(){

	
	long idum;
  	idum= - (time(NULL));//*(myrank));

	int iNumPart=4;
	int iCutoff=2;
	
	double** partPos = new double*[iNumPart];
	double** newPartPos = new double*[iNumPart];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[3];
		newPartPos[i] = new double[3];
		//neccesary?
		//for (int j=0;j<3;j++){
		//	partPos[i][j]=0.0;
		//	newPartPos[i][j]=0.0;
		//}
	}

	double ideal_step=0.9;

	//for (int iTemp=0; iTemp < number_of_alpha_variations; iTemp++){
	//cumulative_energy=0;
	//local_energy=0;
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < 3; j++) {
			partPos[i][j] += 4.*(ran2(&idum) -0.5);
		}
	}
	double  ar[1];
	ar[0]=0.13;	
	
	slaterMatrix a(iNumPart,iCutoff,1);
	slaterMatrix b(iNumPart,iCutoff,1);
	
	a.updateVariationalParameters(ar);
	b.updateVariationalParameters(ar);
	
	a.initSlaterMatrix(partPos);	
	b.initSlaterMatrix(partPos);	
	
	a.findInverse();
	b.findInverse();
//test of update of inverse
#if 0	

	for (int g=0; g<10000000; g++){   
		
		int active_part=0;
		for (int j = 0; j < 3; j++) { 
	 		partPos[active_part][j] += ideal_step*(ran2(&idum) -0.5);
		}
		//a:probably as fast when working with few particles? TEST:when progr. working.
		//Uncomment lines below if one need to compare with unoptimalized method.
		a.initSlaterMatrix(partPos);	
		a.findInverse();
		
		//b:faster when systems are larger.
		b.updateInverse(partPos[active_part],active_part);
		//included in update inverse
		//b.updateSlaterMatrix(partPos[active_part],active_part);	
	}
	a.print();
	b.print();
	
#endif

#if 0
		//Test of wavefunction. Both determinants should be one
		cout<<"wavefunction test:\n"<<setprecision(16)<<"correct ans: 1 : ans: "<<a.waveFunction(partPos[0],0)<<"\n";	
		cout<<setprecision(16)<<"correct ans: 1 : ans: "<<a.waveFunction(partPos[iCutoff],iCutoff)<<"\n";	
		cout<<"\n";	
#endif	
			
	double e_kinetic, r_single_particle, r_12, e_potential, wf_R;
	int dimension=3;
	int i,j,k;
	int charge = 4; //Include in start of main.
	int active_part;
	bool test;
	ideal_step=0.80;
	
	double e_local = 0.0;
	
	int num_cycles=30000000;
	int accepted=0;
	int thermalization=num_cycles*.3;
	
	double accepted_vec[iNumPart];
	for (i=0;i<iNumPart;i++) accepted_vec[i]=0;
	
	cout<<"\n";
	cout<<"num_cycles: \t\t"<<num_cycles<<"\n";
	cout<<"thermalization:\t\t"<<thermalization<<"\n";
	cout<<"ideal_step: \t\t"<<ideal_step<<"\n";

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int re=0;re<num_cycles;re++)
	{
		//Moving one particle at a time. (converges slowly)
		//possible to move two part. at a time?
		active_part=re%(iNumPart);
		for (j = 0; j < dimension; j++) 
		{ 
			newPartPos[active_part][j] += ideal_step*(ran2(&idum) -0.5);
		}		
		//Metropolis test
		wf_R=b.waveFunction(newPartPos[active_part],active_part);
		test = (ran2(&idum) < wf_R*wf_R);
		if (!test)
		{
			//start over with the same coords	
			for (j = 0; j < dimension; j++) 
			{ 
				newPartPos[active_part][j] = partPos[active_part][j];
			}
			continue;
		}
		//continue with updated position
		for (j = 0; j < dimension; j++) 
		{ 
			partPos[active_part][j] = newPartPos[active_part][j];
		}//cblas_dcopy

		b.updateInverse(partPos[active_part],active_part);
		
		//****** if thermalization finished, start collecting data. ***********
		if (re<thermalization){ continue; }
		
		accepted_vec[active_part]++;
		accepted++;
		//Find and collect energy
		//startvimfold
		e_kinetic=0.0;
		//Returns the sum of the laplacians of det_up and det_down
		e_kinetic+=b.lapl(partPos);
		//incl. jastrow later. returns \nabla|D_{up}|\cdot\nabla|D_{down}|
		for (k = 0; k < iCutoff; k++)
		{
			e_kinetic+=2*b.grad(partPos,k);
		}
		e_kinetic *= -0.5;
		
		e_potential=0;
		// contribution from electron-proton potential  
		for (i = 0; i < iNumPart; i++) 
		{ 
			r_single_particle = 0;
			r_single_particle+=cblas_dnrm2(dimension,partPos[i],1);
			e_potential -= charge/r_single_particle;
		}
	
		// contribution from electron-electron potential  
		for (i = 0; i < iNumPart-1; i++) 
		{ 
			for (j = i+1; j < iNumPart; j++) 
			{
				r_12 = 0;  
				for (k = 0; k < dimension; k++) 
				{ 
					r_12 += (partPos[i][k]-partPos[j][k])*(partPos[i][k]-partPos[j][k]);
				}
				e_potential += 1/sqrt(r_12);         
			//cblas_daxpy...(partPos[i],partPos[j])
			//cblas_dnrm2...(partPos[j])
			}
		}
	//endvimfold
	e_local+=e_potential+e_kinetic;
	}	

	cout<<"Local energy:\t\t"<<(e_local/(double)accepted)<<"\n";	
	cout<<"Acc.rate:\t\t"<<accepted/(double)(num_cycles-thermalization)<<"\n";
	for (i=0;i<iNumPart;i++)
		cout<<"accepted_vec["<<i<<"]:\t"<<accepted_vec[i]<<"\n";

	a.initSlaterMatrix(partPos);	
	a.findInverse();
		
	//a.print();
	//b.print();

	for (i=0;i<iNumPart;i++)
	{
		delete partPos[i];
		delete newPartPos[i];
	}
	delete partPos;
	delete newPartPos;
	
	a.clear();
	b.clear();

}//End function 
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
