/*

   -Generalize to more dimensions.
   -Include Jastrow factor.
   -Make classes?
		
   */

//Define USE_MPI true if using MPI.
#define USE_MPI false

#if USE_MPI
	#include "mpi.h"
#endif

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib/lib.h"
#include <time.h>

using namespace  std;
// output file as glbal variable
ofstream ofile;  
// the step length and its squared inverse for the second derivative 
#define h 0.001
#define h2 1000000

// declaraton of functions 

// The Mc sampling for the variational Monte Carlo 
void  mc_sampling(int, int, int, int, int, int, int, int, double *, double *);

// The variational wave function 
double  wave_function(double **, double, int, int);

// The local energy 
double  local_energy(double **, double, double,  int, int);

// prints to screen the results of the calculations  
void  output(int, int, int, double *, double *);

// generates gaussian distributed random numbers
double gaussian_random(long *idum);

//calculates the quantum force for the importance sampling
double calculate_qforce(double**, double**, double, double, int, int);

//Calculates the determinant




// Begin of main program   

//int main()
int main(int argc, char* argv[]){ 
//startvimfold
    char *outfilename;
   	int max_variations, thermalization, charge;
  	int dimension, number_particles; 
   	double *cumulative_e, *cumulative_e2;
   	double *local_cum_e,* local_cum_e2;
   	int numprocs, my_rank, number_cycles;  
   	double local_sum;
   	
#if (!USE_MPI)
	my_rank=0;	
	numprocs=1;
#endif

	//Initital values
   	number_particles = 3;
   	charge = 3;
   	dimension = 3;
    
   	max_variations = 20;
 
  	cumulative_e = new double[max_variations];
   	cumulative_e2 = new double[max_variations];
   	local_cum_e = new double[max_variations];
   	local_cum_e2 = new double[max_variations];

   	//   MPI initializations
#if USE_MPI
   	MPI_Init (&argc, &argv);
   	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
   	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
   	MPI_Status status;

   	// Read from screen a possible new vaue of n
   	// Read in output file, abort if there are too few command-line arguments
   	if(my_rank == 0){
#endif
   		if( argc <= 1 ){
    		cout << "Bad Usage: " << argv[0] <<  " read also output file on same line" << endl; 
        	exit(1); 
	
		} else {
	
    		outfilename=argv[1];
      		ofile.open(outfilename);
       		cout << "# MC steps= ";
      		cin >> number_cycles;
			//cout << "steplength= ";
        	//cin >> ;
		}
#if USE_MPI 
	}
  	//broadcast
  	MPI_Bcast (&number_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
   	
	thermalization = 0.10 * number_cycles;

  	//  Do the mc sampling  
	mc_sampling(my_rank, numprocs, dimension, number_particles, charge, 
              max_variations, thermalization, number_cycles, cumulative_e, cumulative_e2);
  	
#if USE_MPI
  	// Calculations are done. Time to send and recive. First cumulative_e
  	if (my_rank==0) {
    	for (int i=1; i<numprocs; i++) {
       		MPI_Recv(local_cum_e, max_variations,MPI_DOUBLE,MPI_ANY_SOURCE,500,MPI_COMM_WORLD,&status);
#endif 	
			for(int n=0;n < max_variations; n++){
   				cumulative_e[n] += local_cum_e[n];
  			}
#if USE_MPI			
    	}
   	} else {
    	MPI_Send(cumulative_e,max_variations,MPI_DOUBLE,0,500,MPI_COMM_WORLD);
  	}

 	//Sending and reciving cumulative_e2
  	if (my_rank==0) {
    	for (int i=1; i<numprocs; i++) {
       	MPI_Recv(local_cum_e2,max_variations,MPI_DOUBLE,MPI_ANY_SOURCE,501,MPI_COMM_WORLD,&status);
#endif
		for(int n=0; n < max_variations; n++){
   			cumulative_e2[n] += local_cum_e2[n];
  		}
#if USE_MPI
	}
#endif
    output(max_variations, number_cycles, charge, cumulative_e, cumulative_e2);
#if USE_MPI
  } else {
      MPI_Send(cumulative_e2,max_variations,MPI_DOUBLE,0,501,MPI_COMM_WORLD);
  //    cout << my_rank << " : done!" << endl;
  }
#endif

  delete [] cumulative_e; delete [] cumulative_e2; 
  delete [] local_cum_e; delete [] local_cum_e2;
  ofile.close();  // close output file
  
#if USE_MPI  
  	MPI_Finalize ();
#endif
  return 0;
}
//endvimfold

// Monte Carlo sampling with the Metropolis algorithm  
void mc_sampling(int my_rank, 
	int numprocs, int dimension, int number_particles, int charge, 
                 int max_variations, int thermalization, int number_cycles, 
    			double *cumulative_e, double *cumulative_e2){
//startvimfold
  	int cycles, variate, accept, dim, i, j;
  	long idum;
  	double wfnew, wfold, alpha, beta, energy, energy2, delta_e;
  	double **r_old, **r_new;
  	double ideal_step;

  	beta = 1.25;
  	alpha = 2.75;

  	idum=-((my_rank +1) * time(NULL) + 1); // my_rank+1<-my_rank (MC sampl orig. prog.)
  	// allocate matrices which contain the position of the particles  
  	r_old = (double **) matrix( number_particles, dimension, sizeof(double));
  	r_new = (double **) matrix( number_particles, dimension, sizeof(double));
  		for (i = 0; i < number_particles; i++) { 
    		for ( j=0; j < dimension; j++) {
      			r_old[i][j] = r_new[i][j] = 0;
    		}
  		}
  	// loop over variational parameters  
  	for (variate=0; variate < max_variations; variate++){

  	  	// initialisations of variational parameters and energies 
  	  	beta  += 0.025; 
  	  	//alpha += 0.025;                                                           ////////////  
  	  	energy = energy2 = 0; accept =0; delta_e=0;

    	//Find the ideal step legth for given alpha (- ln(0.25) / (2 alpha))
    	ideal_step =0.6; //65(0.25) / (2* alpha));// 0.34657/ alpha ; //XXX Right place to declare this var?
   
    	//  initial trial position, note calling with alpha 
    	//  and in three dimensions 

    	for (i = 0; i < number_particles; i++) { 
      		for ( j=0; j < dimension; j++) {
 				r_old[i][j] = ideal_step*(ran1(&idum) -0.5);
      		}
    	}

   		wfold = wave_function(r_old, beta, dimension, number_particles);
    	
		// loop over monte carlo cycles 
    	for (cycles = 1; cycles <= number_cycles+thermalization; cycles++){ 
      		// new position 
      		for (i = 0; i < number_particles; i++) { 
 				for ( j=0; j < dimension; j++) {
   					//r_new[i][j] = r_old[i][j]+ideal_step*(ran1(&idum) -0.5);
   					r_new[i][j] = r_old[i][j]+gaussian_random(&idum); //XXX + qforce
 				}
      		}
			//XXX
      
			wfnew = wave_function(r_new, beta, dimension, number_particles); 

   			// Metropolis test 
    		if(ran1(&idum) <= wfnew*wfnew/wfold/wfold ) { 
 				for (i = 0; i < number_particles; i++) { 
   					for ( j=0; j < dimension; j++) {
     					r_old[i][j]=r_new[i][j];
   					}
 				}
			wfold = wfnew;
 			accept = accept+1;
  			}

 			// compute local energy  
 
   	 		if ( cycles > thermalization ) {
 				delta_e = local_energy(r_old, beta, wfold, dimension, number_particles);

				// update energies  
        		energy += delta_e;
        		energy2 += delta_e*delta_e;
      		}
    	}   // end of loop over MC trials   

	cout << "variational parameter= " << beta 
    << " accepted steps= " << accept << endl;
    // Calculations are done. Time to send and recive. First cumulative_e
    
	if (my_rank==0) {
    	// update the energy average and its squared 
    	cumulative_e[variate] = energy / number_cycles   /numprocs;
    	cumulative_e2[variate] = energy2/ number_cycles /numprocs;
	}
  }    // end of loop over variational  steps 
  
  free_matrix((void **) r_old); // free memory
  free_matrix((void **) r_new); // free memory

}// end mc_sampling function  
//endvimfold

// Function to compute the (unsquared) wave function, simplest form 
double  wave_function(double **r,// double alpha,
		double beta,int dimension, int number_particles){
/*//startvimfold*/
int i, j, k;                                
double  wf_1=0, wf_2=0 ,wf_3=0, toS=0, toP=0, toQ=0;
double  r_p1=0, r_p2=0, r_p3=0, r_ento=0, r_totre=0, r_entre=0;

double alpha=2.75;

//Find r_p1,r_p2, length of vectors
	for (j = 0; j < dimension; j++) { 
		r_p1  += r[0][j]*r[0][j];
		r_p2  += r[1][j]*r[1][j];
		r_p3  += r[2][j]*r[2][j];
    		}
	r_p1 = sqrt(r_p1);
	r_p2 = sqrt(r_p2);
	r_p3 = sqrt(r_p3);

	//find r_ento, distance between particles 1,2
	r_ento=  pow(r[0][0]-r[1][0],2)
   		    +pow(r[0][1]-r[1][1],2)
        	+pow(r[0][2]-r[1][2],2);
		r_ento=sqrt(r_ento);
	//find r_totre, distance between particles 2,3
	r_totre= pow(r[1][0]-r[2][0],2)
        	+pow(r[1][1]-r[2][1],2)
        	+pow(r[1][2]-r[2][2],2);
	r_totre=sqrt(r_totre);
	//find r_entre, distance between particles 2,3
	r_entre= pow(r[0][0]-r[2][0],2)
        	+pow(r[0][1]-r[2][1],2)
        	+pow(r[0][2]-r[2][2],2);
	r_entre=sqrt(r_entre);

	//find 2 S  
	toS=2*(1.+beta*r_ento);
	//finr 2 P
	toP=2*(1.+beta*r_totre);
	//find 2 Q
	toQ=2*(1.+beta*r_entre);

	//wawefunction  1
	wf_1 = (2.-alpha*r_p1)*exp(-alpha*(r_p1/2. + r_p2));
	//wawefunction 2
	wf_2 = (2.-alpha*r_p2)*exp(-alpha*(r_p2/2. + r_p1));
	//wawefunction 3
	wf_3 = exp(alpha*(-r_p3+r_ento/toS+r_totre/toP+r_entre/toQ));

	//XXX This is the place where we need to include an slater determinant
	//return full wawefunction (antisymmetric spatial part)
	return (wf_1-wf_2)*wf_3;
}
/*//endvimfold*/

// Function to calculate the local energy 
double  local_energy(double **r,
		double beta , double wfold, int dimension, 
                        int number_particles){
/*//startvimfold*/
int charge =3;

  int i, j , k;
  double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12, 
    r_single_particle;
  double **r_plus, **r_minus;
 	
 // double alpha = 2.75; //XX to be removed

  // allocate matrices which contain the position of the particles  
  // the function matrix is defined in the progam library 
  r_plus = (double **) matrix( number_particles, dimension, sizeof(double));
  r_minus = (double **) matrix( number_particles, dimension, sizeof(double));
  for (i = 0; i < number_particles; i++) { 
    for ( j=0; j < dimension; j++) {
      r_plus[i][j] = r_minus[i][j] = r[i][j];
    }
  }
  // compute the kinetic energy  
  e_kinetic = 0;
  for (i = 0; i < number_particles; i++) {
    for (j = 0; j < dimension; j++) { 
      r_plus[i][j] = r[i][j]+h;
      r_minus[i][j] = r[i][j]-h;
      wfminus = wave_function(r_minus,beta,dimension, number_particles); 
      wfplus  = wave_function(r_plus,beta,dimension, number_particles); 
      e_kinetic -= (wfminus+wfplus-2*wfold);
      r_plus[i][j] = r[i][j];
      r_minus[i][j] = r[i][j];
    }
  }
  // include electron mass and hbar squared and divide by wave function 
  e_kinetic = 0.5*h2*e_kinetic/wfold;
  // compute the potential energy 
  e_potential = 0;
  // contribution from electron-proton potential  
  for (i = 0; i < number_particles; i++) { 
    r_single_particle = 0;
    for (j = 0; j < dimension; j++) { 
      r_single_particle += r[i][j]*r[i][j];
    }
    e_potential -= charge/sqrt(r_single_particle);
  }
  // contribution from electron-electron potential  
  for (i = 0; i < number_particles-1; i++) { 
    for (j = i+1; j < number_particles; j++) {
      r_12 = 0;  
      for (k = 0; k < dimension; k++) { 
	r_12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
      }
      e_potential += 1/sqrt(r_12);          
    }
  }
  free_matrix((void **) r_plus); // free memory
  free_matrix((void **) r_minus);
  e_local = e_potential+e_kinetic;
  return e_local;
}    
/*//endvimfold*/

//write data to outfile
void output(int max_variations, int number_cycles, int charge, 
            double *cumulative_e, double *cumulative_e2){
//startvimfold
  int i;
  double alpha, variance, error;
  alpha = 1.25;
  for( i=0; i < max_variations; i++){
    alpha += 0.025;  
    variance = cumulative_e2[i]-cumulative_e[i]*cumulative_e[i];
    error=sqrt(variance/number_cycles);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << alpha;
    ofile << setw(15) << setprecision(8) << cumulative_e[i];
    ofile << setw(15) << setprecision(8) << variance;
    ofile << setw(15) << setprecision(8) << error << endl;

//    fprintf(output_file, "%12.5E %12.5E %12.5E %12.5E \n",
//alpha,cumulative_e[i],variance, error);
  }
//  fclose (output_file);
}  // end of function output         
//endvimfold

//returns gaussian distributed random numbers XXX what const?
double gaussian_random(long *idum){
//startvimfold	
	//double pi_sqrt = 1.7724538509055160;//sqrt(3.1415926535897932);
	//any nr ok
	double random=ran3(idum)-0.5;
	
	// (general gaussian  exp(-x^2/(2*pi*s^2)) / (sqrt(pi)* s) ))
	if (random<0){
		return exp(-random*random)/1.8;
	} else {
		return -exp(-random*random)/1.8;
	}

}//end: function gaussian_random()
//endvimfold

//Finds Quantum-force (del Psi) / Psi
double calculate_qforce(double** r, double** q_force, double beta, double wf, int dimension, int number_particles){
//startvimfold
	double** r_plus;
   	double** r_minus;
  	r_plus = (double **) matrix(number_particles, dimension, sizeof(double));
	r_minus = (double **) matrix(number_particles, dimension, sizeof(double));
 	
	double wf_plus, wf_minus;
	double alpha=2.75;

	//int number_particles=3; //SHOULD ONYL BE NECC. TO DECLARE AT ONE PLACE
	//int dimension=3;
	
	//Moving only one particle at a time:
	for (int i=0; i<number_particles; i++){
		for (int j=0; j<number_particles; j++){
			r_minus[i][j] = r[i][j];
			r_plus[i][j] = r[i][j];
		}
	}

	for (int i=0; i<number_particles; i++){
		for (int j=0; j<number_particles; j++){
			r_plus[i][j]=r[i][j]+h;
			r_minus[i][j]=r[i][j]-h;
			wf_plus = wave_function(r_plus, beta, dimension, number_particles);
			wf_minus =  wave_function(r_minus, beta, dimension, number_particles);
			q_force[i][j]= 2*(wf_plus-wf_minus)/(h*wf);
			r_minus[i][j] = r[i][j];
			r_plus[i][j] = r[i][j];
		} 
	}
}//end function calculate_qforce();
//endvimfold

// For vim users: Defining vimfolds.
// vim:fdm=marker:fmr=//startvimfold//endvimfold
