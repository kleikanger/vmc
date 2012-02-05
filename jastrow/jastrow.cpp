

double jastrow::jastRati(double r_12_new, double r_12, double* variational_parameters){
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
