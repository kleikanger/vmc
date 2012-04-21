
#ifndef WALKER_H
	#include "../walker/walker.h"
#endif
#ifndef SGA_H
#define SGA_H

#include "mpi.h"

class sga {

	private:
	int num_part;
	int spin_up_cutoff;
	int dimension;
	int num_of_var_par;
	int myrank;
	int nprocs;
	walker** quantum_dot;
	MPI_Status status;
	popControl* popCtr;

//	double *energy_gradient;
	
	public:	
	sga(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank, int nprocs, MPI_Status status);
	~sga();
	/*
	   main sampling loop. 
	   */
	void SGAMin(
		int num_cycles, 
		int thermalization, 
		double* var_par, 
		double delta_t, 
		int num_c_dmc_main_loop,
		int	num_c_dmc_inner_loop,
		int initial_number_of_walkers
		);
	void initWalkers(int initial_number_of_walkers, int thermalization, int corr_length, double *var_par);
};
#endif
