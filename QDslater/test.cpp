#include <cblas.h>
#include <cstdlib>
#include <cmath>
#include "slaterMatrix.h"
#include <iostream>
#include "../lib/lib.h"
#include <iomanip>


///*TEST OF PROGRAM
//startvimfold


int main(){

	int num_cycles=7000000;
	int thermalization=num_cycles*.3;
	int accepted=0;
	int iNumPart=2;
	int iCutoff=1;
	int dimension=2;
	
	double ideal_step= 2.0;
	//double omega = 1.0;

	//Number of variational parameters
	int number_of_var_par;	
	//Variational param
	double  var_par[1];
	var_par[0]=0.398;//beta (jastrow)	
	//Increase of variational parameters
	double var_par_inc[1];
	var_par_inc[0]=0.1;
	//Number of variations in each var. par.
	int var_par_cyc[1];
	var_par_cyc[0]=10;
	
	//init pos.vecs.
	double** partPos = new double*[iNumPart];
	double* newPartPos = new double[dimension];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[dimension];
		//newPartPos[i] = new double[dimension];
	}
	
	long idum;
  	idum= - (time(NULL));//*(myrank));
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < dimension; j++) {
			//newPartPos[i][j]=
			partPos[i][j] += 2.*( ran2(&idum) - 0.5 );
		}
	}

	//slaterMatrix a(iNumPart,iCutoff,1, dimension);
	slaterMatrix b(iNumPart,iCutoff,number_of_var_par, dimension);
	//a.updateVariationalParameters(ar);
	b.updateVariationalParameters(var_par);
	//a.initSlaterMatrix(partPos);	
	b.initSlaterMatrix(partPos);	
	//a.findInverse();
	b.findInverse();

	double e_kinetic, r_12, e_potential, wf_R, r_squared;
	double e_local=0.;
	int i,j,k;
	int active_part;
	bool test;
	int accepted_vec[iNumPart];
	for (i=0;i<iNumPart;i++) accepted_vec[i]=0;
	
	cout<<"\n";
	cout<<"num_cycles: \t\t"<<num_cycles<<"\n";
	cout<<"thermalization:\t\t"<<thermalization<<"\n";
	cout<<"ideal_step: \t\t"<<ideal_step<<"\n";

	//******** START Monte Carlo SAMPLING LOOP ***********
	for (int loop_c=0;loop_c<num_cycles;loop_c++)
	{
		//Moving one particle at a time. (converges slowly?)
		//possible to move two part. at a time, one in each SD?
		active_part=loop_c%(iNumPart);
		for (j = 0; j < dimension; j++) 
		{ 
			newPartPos[j] = partPos[active_part][j] + ideal_step*(ran2(&idum) -0.5);
		}		
		//Metropolis test: update position if accepted, try new .
		wf_R=b.waveFunction(newPartPos,active_part);
		if (ran2(&idum) < wf_R*wf_R) 
		{
			cblas_dcopy(dimension,newPartPos,1,partPos[active_part],1);
		}
		else continue;

		//update inverse and slatermatrix
		b.update(partPos[active_part],active_part);
		
		//****** if thermalization finished, start collecting data. ***********
		if (loop_c<thermalization){ continue; }
		accepted_vec[active_part]++;
		accepted++;
		
		//collect local energy
//startvimfold

		//Laplacian: \nabla^2 \Psi_{\uparrow} + \nable^2 \Psi_{\downarrow}
		e_kinetic=b.lapl(partPos);
		//incl. jastrow here.
		//returns \nabla \Psi_{\uparrow} + \nabla Psi_{\downarrow}
		/*for (k = 0; k < dimension; k++)
		{
			e_kinetic-=2*b.grad(partPos,k);//,active_part);
		}
		*/
		e_kinetic *= -0.5;
		
		//ANALYTIC expr. for \nabla^2 \Psi \over \Psi	
		/*double rR = cblas_ddot(dimension,partPos[0],1,partPos[0],1)
					+ cblas_ddot(dimension,partPos[1],1,partPos[1],1);
		e_kinetic=  - 0.5 * ( rR - 4. );
		*/

		e_potential=0;
		//e - V_ext electrostatic pot
		//OPT: only |\vec r|^2 for one particle needs to be updated.
		for (i = 0; i < iNumPart; i++) 
		{ 
			r_squared=cblas_ddot(dimension,partPos[i],1,partPos[i],1);
			e_potential += r_squared;
		}
		e_potential*=0.5;
	
		// e-e electrostatic interaction
		//OPT: only N-1 of the interpart. dist. r_{12} needs to be updated.
		//Ex: n*(n-1)/2 unique elements. 
		//update n-1 elements. ex ip part. k moved then update row k and col k exept diagonal elements. sum 1/r_{ij}.
		//O(n-1) instead of O(n*(n-1)/2) algo.
		/*
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
		*/
		//cblas_daxpy...(partPos[i],partPos[j])
		//cblas_dnrm2...(partPos[j])
	
	//endvimfold
	e_local+=e_potential+e_kinetic;
	}	

	//************************** END OF MC sampling **************************

	cout<<setprecision(16)<<"Local energy:\t\t"<<(e_local/(double)accepted)<<"\n";	
	cout<<"Acc.rate:\t\t"<<accepted/(double)(num_cycles-thermalization)<<"\n";
	for (i=0;i<iNumPart;i++)
		cout<<"accepted_vec["<<i<<"]:\t"<<accepted_vec[i]<<"\n";

	//a.initSlaterMatrix(partPos);	
	//a.findInverse();
	//a.print();
	//b.print();

	for (i=0;i<iNumPart;i++)
	{
		delete partPos[i];
		//delete newPartPos[i];
	}
	delete partPos;
	delete newPartPos;
	
	//a.clear();
	b.clear();

}//End function 
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
