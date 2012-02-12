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

	int num_cycles=1e7;
	int thermalization=0.2*num_cycles;//num_cycles*.5;
	int iNumPart=2;
	int iCutoff=1;
	int dimension=2;
	
	int accepted=0;
	double ideal_step=2;
	//double omega = 1.0;

	//Number of variational parameters
	int num_of_var_par=1;	
	//Variational param
	double  var_par[num_of_var_par];
	var_par[0]=0.398;//beta (jastrow)	
	//Increase of variational parameters
	double var_par_inc[num_of_var_par];
	var_par_inc[0]=0.05;
	//Number of variations in each var. par.
	int var_par_cyc[num_of_var_par];
	var_par_cyc[0]=10;
	
	//init pos.vecs.
	double** partPos = new double*[iNumPart];
	double* newPartPos = new double[dimension];
	for (int i=0; i<iNumPart; i++){
		partPos[i] = new double[dimension];
		//newPartPos[i] = new double[dimension];
		for (int j=0;j<dimension;j++)  
			partPos[i][j]=0.0;
	}

	//interparticle distances updated.
	//1 particle upd, n-1 new parameters (r_ij)
	double* ipd_upd = new double[iNumPart-1];
	double* ones = new double[iNumPart-1];
	for (int i=0;i<iNumPart-1;i++) { ones[i]=1.0; }
	
	long idum;
  	idum= - (time(NULL));//*(myrank));
	//''random'' startposition 
	for (int i = 0; i < iNumPart; i++) { 
		for (int j=0; j < dimension; j++) {
			//newPartPos[i][j]=
			partPos[i][j] = ideal_step*( ran2(&idum) - 0.5 );
		}
	}

	//slaterMatrix a(iNumPart,iCutoff,num_of_var_par,dimension);
	slaterMatrix b(iNumPart,iCutoff,num_of_var_par,dimension);
	//a.updateVariationalParameters(var_par);
	b.updateVariationalParameters(var_par);
	//a.initSlaterMatrix(partPos);	
	b.initSlaterMatrix(partPos);	
	//a.findInverse();
	b.findInverse();

	//ipl: keeping track of len betw part r_ij, and 1/r_ij
	ipdist ipd(iNumPart,dimension,iCutoff);
	ipd.init(partPos);

	double e_kinetic, r_12, r_12_new, e_potential, wf_R, r_squared, temp;
	double e_local=0.0;
	double e_local_squared=0.0;
	int i,j,k;
	int active_part;
	int accepted_vec[iNumPart];
	for (i=0;i<iNumPart;i++) accepted_vec[i]=0;

	//Gradients of jastrow and laplacian
	double** jas_grad=new double*[iNumPart];
	for (i=0;i<iNumPart;i++) jas_grad[i]=new double[dimension];
	double** wf_grad=new double*[iNumPart];
	for (i=0;i<iNumPart;i++) wf_grad[i]=new double[dimension];
	
//**********************
// DELETE 
double energy_ee=0.0;
double energy_ev=0.0;
double energy_j=0.0;
double energy_s=0.0;
//int cct;
//**********************

	cout<<"\n";
	cout<<"num_cycles: \t\t"<<num_cycles<<"\n";
	cout<<"thermalization:\t\t"<<thermalization<<"\n";
	cout<<"ideal_step: \t\t"<<ideal_step<<"\n";

	//******** START Monte Carlo SAMPLING LOOP ***********
	for(int inc = 0; inc<1; inc++){
		var_par[0] += 0.0;// var_par_inc[0];
		e_local=0.0;
		e_local_squared=0.0;
		accepted=0;
	for (int loop_c=0;loop_c<num_cycles+thermalization;loop_c++)
	{
		//Moving one particle at a time. (converges slowly?)
		//possible to try to move all part before collecting energy. 
		//or move two part. at a time, one in each SD?
		//*******
		//function : get_new_pos(active_part, newPartPos, ipd_upd..)
		active_part=loop_c%(iNumPart);
		for (k = 0; k < dimension; k++) 
		{ 
			newPartPos[k]=partPos[active_part][k]+ideal_step*(ran2(&idum)-0.5);
		}
		//New ipd : calculating new lengths betw particles.
		//only iterating over changed lengths.
	
		for (i=0; i<active_part; i++)/*//startvimfold*/
		{
			temp=0;
			for (k=0; k<dimension; k++)
			{
				temp+=pow((newPartPos[k]-partPos[i][k]),2);
			}
			//if (temp<1e-5) temp=1e-5;
			ipd_upd[i]=sqrt(temp);
		}
		for (i=active_part+1; i<iNumPart; i++)
		{
			temp=0;
			for (k=0; k<dimension; k++)
			{
				temp+=pow((newPartPos[k]-partPos[i][k]),2);
			}
			//if (temp<1e-5) temp=1e-5;
			ipd_upd[i-1]=sqrt(temp);
		}/*//endvimfold*/
	
		//Summing over ipd_upd (sum r_12)
		r_12_new= cblas_ddot(iNumPart-1,ipd_upd,1,ones,1);// sum
		//corresponding elements in old vector
		r_12 = ipd.sumPart(active_part);
		
		//jastrow ratio: 2 particle impl.
		//********************************
		double jas_R = exp( r_12_new/(1.0 + var_par[0]*r_12_new) - r_12/(1.0 + var_par[0]*r_12 ) );
		//********************************

		//wave func ratio
		wf_R=b.waveFunction(newPartPos,active_part);
		//Metropolis test: update slater+ipd+positions if accepted.
		if ( ran2(&idum) <= (jas_R * jas_R * wf_R * wf_R ) ) 
		{ //function update
			//update positions	
			cblas_dcopy(dimension,&newPartPos[0],1,&partPos[active_part][0],1);
			//update inverse and slatermatrix
			b.update(newPartPos,active_part);
			//update interparticle distances
			ipd.update(ipd_upd,active_part);
			//update acceptancerate
			if (loop_c>thermalization) 
			{
				accepted_vec[active_part]++;
				accepted++;
			}
		}
		//****** if thermalization finished, start collecting data. ***********
		if (loop_c<thermalization) { continue; }
		//collect local energy function -> get energy
		e_kinetic=0;
		//Laplacian: \nabla^2 \Psi_{\uparrow} + \nable^2 \Psi_{\downarrow}
		e_kinetic+=b.lapl(partPos);
		//e_kinetic+=(cblas_ddot(dimension,partPos[1],1,partPos[1],1)+cblas_ddot(dimension,partPos[0],1,partPos[0],1) - dimension*iNumPart);
		//energy_s-=0.5*b.lapl(partPos);
		e_kinetic+=ipd.jasLapl(var_par[0],partPos);
		//double beta = var_par[0];
		//e_kinetic +=  2*( 1.0 - beta*r_12_new ) / r_12_new / ( pow((1.0+r_12_new*beta),3) );
		//energy_j-=0.5*ipd.jasLapl(var_par[0], partPos);
		ipd.jasGrad(jas_grad, var_par[0],partPos);
		//temp = 1. / r_12_new / (1.+beta*r_12_new) / (1.+beta*r_12_new);
		//for (i=0;i<iNumPart;i++)
		//for (j=0;j<dimension;j++)
		//	jas_grad[i][j]=temp*(partPos[i][j]-partPos[(i+1)%2][j]);
		
		//TESTED: 2 part
		b.grad(wf_grad,partPos);
		//for (i=0;i<iNumPart;i++)
		//for (j=0;j<dimension;j++)
		//	wf_grad[i][j]=-partPos[i][j];
		
		for (i=0;i<iNumPart;i++)
		{
			e_kinetic+=2.0*cblas_ddot(dimension,jas_grad[i],1,wf_grad[i],1);
			e_kinetic+=cblas_ddot(dimension,jas_grad[i],1,jas_grad[i],1);
		}
		e_kinetic *= -0.5;

		e_potential=0.0;
		//e - V_ext electrostatic pot
		//OPT: only |\vec r|^2 for one particle needs to be updated.
		for (i = 0; i < iNumPart; i++) 
		{
			e_potential+=cblas_ddot(dimension,partPos[i],1,partPos[i],1);
			energy_ev += e_potential*.5;
		}
		e_potential*=0.5;
		// e-e electrostatic interaction
		e_potential += ipd.sumInvlen();
		energy_ee += ipd.sumInvlen();//ipd.sumInvlen();
		
		//for (i = 0; i < iNumPart-1; i++) 
		//{ 
		//	for (j = i+1; j < iNumPart; j++) 
		//	{
		//		r_12 = 0;  
		//		for (k = 0; k < dimension; k++) 
		//		{ 
		//			r_12 += (partPos[i][k]-partPos[j][k])*(partPos[i][k]-partPos[j][k]);
		//		}
		//		e_potential += 1/sqrt(r_12);         
		//		energy_ee +=1/sqrt(r_12); 
		//	}
		//}
		//cblas_daxpy...(partPos[i],partPos[j])
		//cblas_dnrm2...(partPos[j])
	e_local += e_potential+e_kinetic;
	e_local_squared += pow((e_potential+e_kinetic),2); 
	}
	//endvimfold
	cout<<"variance:\t\t"<<(e_local_squared-e_local*e_local/num_cycles)/num_cycles<<"\n";
	cout<<setprecision(16)<<"Local energy:\t\t"<<(e_local/(double)num_cycles)<<"\n";	
	cout<<"ev "<<energy_ev/(double)num_cycles<<"\n";
	cout<<"ee "<<energy_ee/(double)num_cycles<<"\n";
	cout<<"Acc.rate:\t\t"<<accepted/(double)(num_cycles)<<"\n";
	cout<<"beta :\t\t"<<var_par[0]<<"\n";
	energy_ee=0.0;
	energy_ev=0.0;
	}	

	//************************** END OF MC sampling **************************

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
		delete jas_grad[i];
		delete wf_grad[i];	
	}
	delete partPos;
	delete newPartPos;
	delete jas_grad;
	delete wf_grad;
	
	delete ipd_upd;
	delete ones;
	
	//a.clear();
	b.clear();
	ipd.clear();

}//End function 
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold,//endvimfold
