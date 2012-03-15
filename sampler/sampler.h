#ifndef WALKER_H
	#include "../walker/walker.h"
#endif
#ifndef SAMPLER_H
#define SAMPLER_H


class sampler {

	private:
	int num_part;
	int spin_up_cutoff;
	int dimension;
	int num_of_var_par;
	int myrank;
	walker* quantum_dot;
	
	double *energy_gradient;
	
	public:	
	sampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank);
	~sampler();
	/*
		Returning energy_gradient[i]. For use in the cgm - optimization class.
	   */
	double getEnergyGrad(int i);
	/*
	   main sampling loop. metropolis-hastings test and collecting of energy samples.
	   */
	void sample(int num_cycles, int thermalization, double* var_par, double delta_t, double* result);

};

#endif
