#include "../QDslater/slaterMatrix.h"
#include "../ipdist/ipdist.h"

class vmcsolver{

	private:
			
		int num_part;
		int spin_up_cutoff;
		int dimension;
		int num_of_var_par;
	
		//seed for the random number generator	
		int idum;
	   	//int idum2;	

		//particle positions 
		double** r_old;
		//gradient of the jastrow function
		double** jas_grad;
		//gradient of the slater wavefunction
		double** sla_grad;
		//new and old quantum force
		double** q_force_new;
		double** q_force_old;

		//handling of the slater matrix
		slaterMatrix* slater;
		
		//handling of the interparticle distances
		ipdist* ipd;
				
	public:

		//constructor, destructor
		vmcsolver(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par);
		~vmcsolver();
		
		//Main function, mc sampling
		void sample(int num_cycles, int thermalization, double* var_par, double delta_t);
		//void calcEnergy()
		
		//random update of the position of active_part.
		void getNewPos(int active_part, double** r_new, 
				double* ipd_upd, double dt_x_D, double diff_const);
		
		//calculate the local energy in the present position.
		double calcLocalEnergy(double* var_par);
};
