#include "../QDslater/slaterMatrix.h"
#include "../ipdist/ipdist.h"

class vmcsolver{

	private:
			
		int num_part;
		int spin_up_cutoff;
		int dimension;
		int num_of_var_par;
		
		long idum; 

		double** r_old;
		double** jas_grad;
		double** wf_grad;
		
		//handling of the slater matrix
		slaterMatrix* slater;
		
		//handling of the interparticle distances
		ipdist* ipd;
				
	public:

		vmcsolver(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par);
		~vmcsolver();

		void sample(int num_cycles, int thermalization, double* var_par);
		//void calcEnergy()
		
		void getNewPos(int active_part, double ideal_step, double* r_new, double* ipd_upd);
		
		void calcLocalEnergy(double* e_local, double* e_local_squared,double* var_par);
};
