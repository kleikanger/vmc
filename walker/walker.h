#ifndef POPCONTROL_H
	#include "../popControl/popControl.h"
#endif
#ifndef SLATERMATRIX_H
	#include "../QDslater/slaterMatrix.h"
#endif
#ifndef IPDIST_H
	#include "../ipdist/ipdist.h"
#endif

#ifndef WALKER_H
#define WALKER_H
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
		double omega;

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
				
		//Must be able to return the sign of this value for the fixed node approw in dmc
		double wf_R;

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
		   Initialize walker with no start values. For DMC initialization.
		   */
		void initEmptyWalker(double* var_par, double delta_t);
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
		   Overloaded function. The second variant returns 
		   double* energies: add pot, osc and kin E to vector
		 */
		double calcLocalEnergy() const; 
		double calcLocalEnergy(double* energies) const; 
		/*
		   Find a new random position for the given particle. Update position
		   and calculate all new interparticle distances
		   */
		void getNewPos(const int &active_part, double* ipd_upd);
		/*
		   Copy position of particle i_w to x.
		   */
		void getRi(int i_w, double* x);
		/*
		   Returning ... for calculating the minima via the conjugate gradient method
		   */
		void getVarParGrad(double* grad_var_par) const;
		/*
		   Returning the value of the sign of wf_R to check whether a walker has crossed a node.
		   if true : node crossed
		   if false: node not crossed
		   */
		bool nodeCrossed();
		/*
		   Resets the variational parameters
		   */
		void wSetVarPar(double* varPar);
		/*
		   Returns the length of the last move
		   */
		double getLenOfMoveSqrd();
	
	friend class popControl;
};
#endif
