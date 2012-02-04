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

	int iNumPart=2;
	int iCutoff=1;
	int dimension=3;
	
	double** partPos = new double*[iNumPart];
	double** newPartPos = new double*[iNumPart];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[3];
		newPartPos[i] = new double[3];
	}

	double ideal_step=3.1;

	//for (int iTemp=0; iTemp < number_of_alpha_variations; iTemp++){
	//cumulative_energy=0;
	//local_energy=0;
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < 3; j++) {
			newPartPos[i][j]=partPos[i][j] += 2.*( ran2(&idum) - 0.5 );
		}
	}

	double  ar[1];
	ar[0]=0.13;	
	
	slaterMatrix a(iNumPart,iCutoff,1, dimension);
	slaterMatrix b(iNumPart,iCutoff,1, dimension);
	
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
	int i,j,k;
	double charge = 2; //Include in start of main.
	int active_part;
	bool test;
	ideal_step= 1.30;
	
	double e_local=0.;
	
	int num_cycles=70000000;
	int accepted=0;
	int thermalization=num_cycles*.5;
	
	int accepted_vec[iNumPart];
	for (i=0;i<iNumPart;i++) accepted_vec[i]=0;
	
	cout<<"\n";
	cout<<"num_cycles: \t\t"<<num_cycles<<"\n";
	cout<<"thermalization:\t\t"<<thermalization<<"\n";
	cout<<"ideal_step: \t\t"<<ideal_step<<"\n";

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int re=0;re<num_cycles;re++)
	{
		//Moving one particle at a time. (converges slowly?)
		//possible to move two part. at a time, one in each SD?
		active_part=re%(iNumPart);
		for (j = 0; j < dimension; j++) 
		{ 
			newPartPos[active_part][j] += ideal_step*(ran2(&idum) -0.5);
		}		
		//Metropolis test
		wf_R=b.waveFunction(newPartPos[active_part],active_part);
		
		//double rra=cblas_dnrm2(3,partPos[0],1);
		//double rrb=cblas_dnrm2(3,partPos[1],1);
		//double rrc=cblas_dnrm2(3,newPartPos[0],1);
		//double rrd=cblas_dnrm2(3,newPartPos[1],1);
		//wf_R=exp(1.69*(rra+rrb-rrc-rrd));
		
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
		//Returns the sum of the laplacians of det_up and det_down
		
			
		//double rr0=cblas_dnrm2(3,partPos[0],1);
		//double rr1=cblas_dnrm2(3,partPos[1],1);
		//e_kinetic=0.;
		//e_kinetic+=1.69*(1.69-2/rr0);//b.lapl(partPos);
		//e_kinetic+=1.69*(1.69-2/rr1);//
		
		e_kinetic=b.lapl(partPos);
		//incl. jastrow later. 
		//returns \nabla|D_{up}|\cdot\nabla|D_{down}|
		for (k = 0; k < dimension; k++)
		{
			e_kinetic+=2*b.grad(partPos,k);//,active_part);
		}
		//double reQ=cblas_ddot(3,partPos[0],1,partPos[1],1);
		//double reX=cblas_dnrm2(3,partPos[0],1);
		//double reZ=cblas_dnrm2(3,partPos[1],1);
		//e_kinetic+=1.69*reQ*exp(-1.69*(reX+reZ));

		e_kinetic *= -0.5;
		
		e_potential=0;
		// contribution from electron-proton potential  
		for (i = 0; i < iNumPart; i++) 
		{ 
			r_single_particle=cblas_dnrm2(dimension,partPos[i],1);
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
			}
		}
		
		//cblas_daxpy...(partPos[i],partPos[j])
		//cblas_dnrm2...(partPos[j])
	
	//endvimfold
	e_local+=e_potential+e_kinetic;
//	double tr1=cblas_dnrm2(3,partPos[1],1);
//	double tr0=cblas_dnrm2(3,partPos[0],1);

//	r_12 = (partPos[0][1]-partPos[1][1])*(partPos[0][1]-partPos[1][1]);
//	r_12 += (partPos[0][2]-partPos[1][2])*(partPos[0][2]-partPos[1][2]);
//	r_12 += (partPos[0][0]-partPos[1][0])*(partPos[0][0]-partPos[1][0]);
	

//	cout<<e_potential+e_kinetic<<", ";
//	cout<<(1.69-2.)*(1./tr1+1./tr0)+1/r_12-1.69*1.69<<"\n";


	}	

	//************************** END OF MC sampling **************************

	cout<<"Local energy:\t\t"<<(e_local/(double)accepted)<<"\n";	
	cout<<"Acc.rate:\t\t"<<accepted/(double)(num_cycles-thermalization)<<"\n";
	for (i=0;i<iNumPart;i++)
		cout<<"accepted_vec["<<i<<"]:\t"<<accepted_vec[i]<<"\n";

	a.initSlaterMatrix(partPos);	
	a.findInverse();
		
	a.print();
	b.print();

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
