#include "../QDslater/slaterMatrix.h"
#include "../ipdist/ipdist.h"

class vmcsolver{

	private:
			
		int num_part;
		int spin_up_cutoff;
		int dimension;
		int num_of_var_par;
	
		//seed for the random number generator	
		long idum; 
		
		//particle positions 
		double** r_old;
		//gradient of the jastrow function
		double** jas_grad;
		//gradient of the slater wavefunction
		double** sla_grad;
		
		//handling of the slater matrix
		slaterMatrix* slater;
		
		//handling of the interparticle distances
		ipdist* ipd;
				
	public:

		//constructor, destructor
		vmcsolver(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par);
		~vmcsolver();
		
		//Main function, mc sampling
		void sample(int num_cycles, int thermalization, double* var_par);
		//void calcEnergy()
		
		//random update of the position of active_part.
		void getNewPos(int active_part, double ideal_step, double* r_new, double* ipd_upd);
		
		//calculate the local energy in the present position.
		void calcLocalEnergy(double* e_local, double* e_local_squared,double* var_par);
};
