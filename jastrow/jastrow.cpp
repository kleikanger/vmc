

double jastrow::jastRat(double r_12_new, double r_12, double* variational_parameters){
	double beta = variational_parameters[0];
	return jas_R=exp(0.5*( r_12_new/(1 + beta*r_12_new) -  0.5*( r_12/(1 + beta*r_12) )));
}

double jastrow::jast(double** ip_len, double* variational_parameters){
//startvimfold
	double r_12;
	int i,j,k;
   	//summation variable
	double value = 0;
	//variational parameter.
	double beta = variational_parameters[0];
  	
   	n_min_one=iNumPart-1;	
	
	for (i=0; i<n_min_one; i++)
	{
		for (j=0; j<i+1; j++)
		{	
            r_12=ip_len[i][j];
            //Particles with parallel spins:
            if (i<iCutoff==j<iCutoff) 
			{
                value += 0.3333333333333333 * r_12/(1 + beta*r_12);
            }
            //particles with antiparallel spins
            else 
			{
                value += r_12/(1 + beta*r_12);
            }
        }
    }
    return exp(value*0.5);
}//End function slaterMatrix::jastrow()
//endvimfold

double jastrow::lapl(double** ip_len, double* variational_parameters, double** r){
//startvimfold
	double r_12;
	int i,j,k;
   	//summation variable
	double sum_1 = 0;
	double sum_2 = 0;
	//temporary array
	double* temp_arr_i[dim];
	double* temp_arr_j[dim];
	//variational parameter.
	double beta = variational_parameters[0];
  	//temp vars
	double ove_ki;
	double ove_ki;

	double a=0.3333333333333333;
	double aa=a*a;
//SUM OVER K
//NB: THIS CALC LAPL FOR SUB DIAG ELEMENTS.
//THE TOTAL LAPL SHOULD BE MULT. WITH 2!

   	n_min_one=iNumPart-1;	
	
	for (k=0; k<n_min_one; k++) //sum over k here?
	{
		for (j=0; j<k+1; j++)
		{
			//CAN BE DONE MORE EFFICIENTLY: 2 LOOPS, LESS IFTESTS.
			if (j==k) continue;
			for (i=0; i<k+1; i++)
			{	
				if (i==k) continue;
				temp=0;
				//temp_arr<-r[]
				cblas_dcopy(dim,r[i],temp_arr_i,1);
				cblas_dcopy(dim,r[j],temp_arr_j,1);
				//temp_arr<r[k]-temp
				cblas_daxpy(dim,-1.0,r[k],1,temp_arr_i,1);
				cblas_daxpy(dim,-1.0,r[k],1,temp_arr_j,1);
				
				//(r_k-r_j).(r_k-r_i)
				temp=cblas_ddot(dim,temp_arr_i,1,temp_arr_j,1);
				
				temp/=ipd_len[k][i]*ipd_len[k][j]; //*ipd_invlen
			
				//temporary vars
				ove_ki=(1+beta*ipd_len[k][i]);
				ove_kj=(1+beta*ipd_len[k][j]);

				temp*=ove_ki*ove_ki*ove_kj*ove_kj;
				sum_1+=temp;
			}
		
			//Particles with parallel spins:
			if (i<iCutoff==j<iCutoff) 
			{
				sum_2 +=   a*ove_kj * ove_kj / ip_len[k][j] ;
				sum_2 -=   a*beta * ove_kj * ove_kj * ove_kj;
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				sum_2 +=   ove_kj * ove_kj / ip_len[k][j] ;
				sum_2 -=   beta * ove_kj * ove_kj * ove_kj;
			}
		}
	}

	sum_1*=a*a;
	sum_2*=2;

    return sum_1+sum_2;
}//End function slaterMatrix::jastrow()

//Changing elements in ret_vec to gradient
double jastrow::grad(double** ret_vec, double** ip_len, double* variational_parameters, double** r){
//startvimfold
	double r_12;
	int j,k,J;
   	//summation variable
	double sum = 0;
	//temporary array
	double* temp_arr[dim];
	//variational parameter.
	double beta = variational_parameters[0];

	double a=0.3333333333333333;

//NB: THIS CALC LAPL FOR SUB DIAG ELEMENTS.
//THE TOTAL LAPL SHOULD BE MULT. WITH 2????

	

   	n_min_one=iNumPart-1;	
	double sum=0.0;
	//i_upd minus one

	
	for (k=0;k<n;k++)
	{
		for (j=0;j<k;j++)
		{       
			//k=n-1: J=0,1,2,...,n-2, K=n-1
			r_kj=ip_len[k-1][j]; 
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],temp_arr,1);
			//temp_arr< r(k)-r[j]
			cblas_daxpy(dim,-1.0,r[j],1,temp_arr,1);
			//Particles with parallel spins:
			if (J<iCutoff==K<iCutoff) 
			{
				temp = a /r_kj /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /ipd_len*r_kj) /(1+beta*r_kj) /(1+beta*r_kj);
			}
			cblas_dscal(dim, temp, temp_arr, 1);
			cblas_daxpy(dim, 1, temp_arr, 1, ret_vec[k], 1)
			//same as?? cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1)
		}
		for (j=k;j<n_min_one;j++)
		{
			//k=0: J=1,2,3,4,5,6,..,n-1
			r_kj=ip_len[j][k];
			J=j+1; 
			
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],temp_arr,1);
			//temp_arr< r(k)-r[j]
			cblas_daxpy(dim,-1.0,r[J],1,temp_arr,1);
				
			//Particles with parallel spins:
			if (K<iCutoff==J<iCutoff) 
			{
				temp = a /r_kj) /(1+beta*r_kj) /(1+beta*r_kj);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /r_kj) /(1+beta*r_kj) /(1+beta*r_kj);
			}
			cblas_dscal(dim, temp, temp_arr, 1);
			cblas_daxpy(dim, 1, temp_arr, 1, ret_vec[k], 1)
			//same as?? cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1)
		}
	}
/*
	for (k=0; k<n_min_one; k++) //sum over k here?
	{
		for (j=0; j<k+1; j++)
		{
			//CAN BE DONE MORE EFFICIENTLY: 2 LOOPS, LESS IFTESTS.
			if (j==k) continue;

			temp=0;
			
			//temp_arr<-r[k]
			cblas_dcopy(dim,r[k],temp_arr,1);
			//temp_arr< r(k)-r[j]
			cblas_daxpy(dim,-1.0,r[j],1,temp_arr,1);
				
			//Particles with parallel spins:
			if (i<iCutoff==j<iCutoff) 
			{
				temp = a /ipd_len[k][j]) /(1+beta*ipd_len[k][j]) /(1+beta*ipd_len[k][j]);
			}
			//particles with antiparallel spins
			else	
			{
				//a=1.0
				temp = 1. /ipd_len[k][j]) /(1+beta*ipd_len[k][j]) /(1+beta*ipd_len[k][j]);
			}
			cblas_dscal(dim, temp, temp_arr, 1);
			cblas_daxpy(dim, 1, temp_arr, 1, ret_vec[k], 1)
			//same as?? cblas_daxpy(dim, temp, temp_arr, 1, ret_vec[k], 1)
		}
	}
*/

}//End function slaterMatrix::jastrow()
//endvimfold
