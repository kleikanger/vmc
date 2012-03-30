
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

//	double *energy_gradient;
	
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
		void sortWalkers(int &num_alive, int killsd, int num_resurrected, bool *occupancy, double *e_local_old);	
		void redistributeWalkers(int myrank, int nprocs, int &num_alive, bool *occupancy, double *e_local_old, double renorm_threshh);
		void initializeSys(int initial_number_of_walkers, int thermalization, int corr_length, double *var_par, double *e_local_old);
};
#endif
