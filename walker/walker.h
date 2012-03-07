
#include "../QDslater/slaterMatrix.h"
#include "../ipdist/ipdist.h"

class walker {

	private:
			
		int num_part;
		int spin_up_cutoff;
		int dimension;
		int num_of_var_par;

		//delta_t not necc. DELETE
		double delta_t;
		double dt_x_D;
		double sq_delta_t;

		//seed for the random number generator	
		int idum;
	   	//int idum2;

		//DELETE XXX REMOVE
		double* var_par;

		//new and old particle positions 
		double** r_new;
		double** r_old;
		//gradient of the jastrow function
		double** jas_grad;
		double** jas_grad_bu;
		//gradient of the slater wavefunction
		double** sla_grad;
		double** sla_grad_bu;
		//new and old quantum force
		double** q_force_new;
		double** q_force_old;
		//handling of the slater matrix
		slaterMatrix* slater;
		//handling of the interparticle distances
		ipdist* ipd;
				
	public:

		//constructor, destructor
		walker(int num_part, int spin_up_cutoff, int dimension, int num_of_var_par, int myrank);
		~walker();
		/*
		   Initialization of walker, random startposition etc. Initialising
		   and slaterMatrix classes with variational parameters and initial
		   configuration.
		   */
		void initWalker(double* var_par, double delta_t);
		/*
		   Moving one particle and returning result of metropolis-hastings test.
		   */
		bool tryRandomStep(int active_part);
		/*
		   Updating all relevant quantities in this class and in slaterMatrix and ipdist.
		   */
		void acceptStep(int active_part);
		/*
		   Reset all relevant quantities in this class and in slaterMatrix and ipdist.
		   */
		void rejectStep(int active_part);
		/*
		   Return local energy in the current position
		 */
		double calcLocalEnergy(double* var_par) const; 
		/*
		   Find a new random position for the given particle. Update position
		   and calculate all new interparticle distances
		   */
		void getNewPos(const int &active_part, double* ipd_upd);
		/*
		   Copy position of particle i_w to x.
		   */
		void getRi(int i_w, double* x);
};  
