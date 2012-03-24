
#ifndef WALKER_H
	#include "../walker/walker.h"
#endif
#ifndef DMCSAMPLER_H
#define DMCSAMPLER_H

#include "mpi.h"

class dmcsampler {

	private:
	int num_part;
	int spin_up_cutoff;
	int dimension;
	int num_of_var_par;
	int myrank;
	int nprocs;
	walker** quantum_dot;
	popControl* popCtr;
	MPI_Status status;

	double *energy_gradient;
	
	public:	
	dmcsampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank, int nprocs, MPI_Status status);
	~dmcsampler();
	/*
	   main sampling loop. 
	   */
	void sampleDMC(
		int num_cycles, 
		int thermalization, 
		double* var_par, 
		double delta_t, 
		double e_trial,
		int num_c_dmc_main_loop,
		int	num_c_dmc_inner_loop,
		int	num_c_dmc_equilibriation_loop,
		int initial_number_of_walkers
		);
};
#endif
