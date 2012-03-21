
#ifndef SAMPLER_H
	#include "../sampler/sampler.h"
#endif
#include "vectormatrixclass.h"

#ifndef CGM_H
#define CGM_H
class cgm {

	private:
		int num_cycles;
		int thermalization;
		double delta_t;
		double omega;
		int num_of_var_par;
		int myrank;
		int nprocs;
		int num_part;
		int spin_up_cutoff;
		int dimension;
		int argc_c;
		char **argv_c;
		double *dE_array;

	sampler* sampler_;

	public:
		cgm(int num_cyclesARG, 
			int thermalizationARG, 
			double delta_tARG, 
			double omegaARG, 
			int num_partARG, 
			int spin_up_cutoffARG,
			int dimensionARG, 
			int num_of_var_parARG, 
			int myrank_ARG, 
			int nprocsARG);
		~cgm();
		void optimizeVarPar(double* initial_var_par);
		double E_function(Vector &variational_parameters);
		void dE_function(Vector &variational_parameters, Vector &gradient);
};
#endif
