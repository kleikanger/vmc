#include "../walker/walker.h"
class sampler {

	private:
	int num_part;
	int spin_up_cutoff;
	int dimension;
	int num_of_var_par;
	int myrank;
	walker* quantum_dot;

	public:	
	sampler(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank);
	~sampler();

	/*
	   main sampling loop. metropolis-hastings test and collecting of energy samples.
	   */
	void sample(int num_cycles, int thermalization, double* var_par, double delta_t, double* result);

};
