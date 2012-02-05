#include <cblas.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "../lib/lib.h"
#include <iomanip>

#include "slaterMatrix.h"
#include "../ipdist/ipdist.h"

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
	int num_of_var_par=1;	
	//Variational param
	double  var_par[num_of_var_par];
	var_par[0]=0.398;//beta (jastrow)	
	//Increase of variational parameters
	double var_par_inc[num_of_var_par];
	var_par_inc[0]=0.1;
	//Number of variations in each var. par.
	int var_par_cyc[num_of_var_par];
	var_par_cyc[0]=10;
	
	//init pos.vecs.
	double** partPos = new double*[iNumPart];
	double* newPartPos = new double[dimension];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[dimension];
		//newPartPos[i] = new double[dimension];
	}

	//interparticle distances updated.
	//1 particle upd, n-1 new parameters (r_ij)
	double* ipd_upd = new double[iNumPart-1];
	double* ones = new double[iNumPart-1];
	for (int i=0;i<iNumPart-1;i++) ones[i]=1.0;
	
	long idum;
  	idum= - (time(NULL));//*(myrank));
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < dimension; j++) {
			//newPartPos[i][j]=
			partPos[i][j] += ideal_step*( ran2(&idum) - 0.5 );
		}
	}

	//slaterMatrix a(iNumPart,iCutoff,1, dimension);
	slaterMatrix b(iNumPart,iCutoff,num_of_var_par,dimension);
	//a.updateVariationalParameters(ar);
	b.updateVariationalParameters(var_par);
	//a.initSlaterMatrix(partPos);	
	b.initSlaterMatrix(partPos);	
	//a.findInverse();
	b.findInverse();

	//ipl: keeping track of len betw part r_ij, and 1/r_ij
	ipdist ipd(iNumPart,dimension);
	ipd.init(partPos);

	double e_kinetic, r_12, r_12_new, e_potential, wf_R, r_squared;
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

		//New ipd : calculating new lengths betw particles
		for (i=0; i<active_part; i++)/*//startvimfold*/
		{
			r_12=0;
			for (k=0; k<dimension; k++)
			{
				r_12+=(newPartPos[k]-partPos[i][k])*(newPartPos[k]-partPos[i][k]);
			}
			ipd_upd[i]=sqrt(r_12);
		}
		for (i=active_part+1; i<iNumPart; i++)
		{
			r_12=0;
			for (k=0; k<dimension; k++)
			{
				r_12+=(newPartPos[k]-partPos[i][k])*(newPartPos[k]-partPos[i][k]);
			}
			ipd_upd[i]=sqrt(r_12);
		}/*//endvimfold*/
	   	
		//jastrow ratio: brute force impl.
		//********************************
		double beta = var_par[0];
		
		//Summing the objects in ipd_upd (sum r_12)
		r_12_new = cblas_ddot(iNumPart-1,ipd_upd,1,ones,1);// sum
		//corresponding elements in old vector
		r_12 = ipd.sumPart(active_part);	
		
		double jas_R=exp(0.5*( r_12_new/(1 + beta*r_12_new) -  0.5*( r_12/(1 + beta*r_12) )));
		//********************************

		//Metropolis test: update position if accepted, try new .
		wf_R=b.waveFunction(newPartPos,active_part);
	
		if (ran2(&idum) < jas_R*jas_R*wf_R*wf_R) 
		{
			cblas_dcopy(dimension,newPartPos,1,partPos[active_part],1);
		}
		else continue;

		//update inverse and slatermatrix
		b.update(partPos[active_part],active_part);
		//update interparticle distances
		ipd.update(ipd_upd,active_part);

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
		
		e_potential += ipd.sumInvlen();
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

	//cout<<"a\n";


	//OBS: pointer to VAR_PAR (malloced) 

	for (i=0;i<iNumPart;i++)
	{
		delete partPos[i];
		//delete newPartPos[i];
	}
	delete partPos;
	delete newPartPos;
	
	//a.clear();
	b.clear();
	ipd.clear();

}//End function 
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
