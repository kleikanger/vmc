/*
	keeps track of all interparticle distanses. for n particles n*(n-1)/2 elements 
	are needed, and n-1 elements needs to be updated every time a patricle is moved.
	update O(n-1) instead of O(n(n-1)/2).

	class keeping two matrices one for the interparticle distances r_ij, and
   	one	for corresponding elements for the old matrix. 

	returns values for the gradient, the laplacian and the log of the jastrow function.
	Maybe it will pay off to init a new matrix g_ij = r_ij /(1+beta*rij) ??

	does not check if 1/r_ij is singular. that could in principle be a problem since 
	a jump to a configuration where r_i=r_j will always be accepted by the metropolis-algo
	since \Psi^{new}(r_i=r_j) = 0. iwhen calculating the new r_{i_upd,j} vec, one could
	check if any elements are 0. Chances are very small for this to happen \sim 10^-16 
	/ mc-cycle.

*/
#ifndef POPCONTROL_H
	#include "../popControl/popControl.h"
#endif
#ifndef IPDIST_H
#define IPDIST_H
class ipdist{

	private: 

		//number of particles -1
		int n_min_one;
		//dimension
		int dim;
		//spin down particles: only for calculating jastrow grad and lapl
		int iCutoff;
		//variational parameter
		double beta;
		
		/*
		   lower triangular with no diag elements n(n-1)/2 elem.

		   ( 	0 		0 		0 		0 		 0)
		   ( 	r_10 	0 		0 		0 		 0)
		   ( 	r_20 	r_21 	0 		0 		 0)
		   ( 	. 		. 		... 	0 		 0)
		   ( 	. 		.       ... 	0 		 0)
		   ( 	r_n0 	r_n1 	... 	r_n(n-1) 0)
			
		   r_{ij} = |\vec r_i - \vec r_j|
		 
		   element (i,j) corresponds to r_(i+1,j)   
		 
		 */

		// elements (r_ij)
		double** ip_len;
		// elements (r_ij). old iplen.
		double** ip_len_backup;

	public: 
		/*
		   declare and alloc ip_len and ip_invlen.
		   set n_min_one to n-1, dim to di and iCutoff to iC 
		 */
		ipdist(int n, int di, int iC);
		/*
		   clear all malloced vars
		 */
		~ipdist();
		/*
		   init all elements in ip_len and ip_invlen
		   NB! one should check if any of these elements are 0.
		 */
		void init(double** x);
		/*
		   Set variational parameter beta
		   */
		void setBeta(double beta);
		/*
		   updates n-1 elements of ip_len and ip_invlen
		   u_upd column and row.
		  
		   input =
		   		( 
		 		x_{i_upd,0}
				,...
				,x_{i_upd,i_upd-1}
				,x_{i_upd,i_upd+1}
				,...
				,x_{i_upd,n-1} 
				) 	

		 */
		void update(double* r, int i_upd);
		/*
		   Reset matrices to last value before update (backup'matr)
		   */
		void reject(int active_part);
		/*
		   Accept last update (update backup'matr)
		   */
		void accept(int i_upd);
		/*
		   Return sum of inverse elements. e-e electrostatic energy = ee/sum
		 */	   
		double sumInvlen() const;
		/*
		   returns R. exp(R_new-R_old) = the jastrow ratio when only one 
		   particle r_{i_upd}  is moved.
		   */
		double logJasR(const int &i_upd) const;
		/*
		   returns gradient of jastrow. i_upd part. in ret_vec will be changed to new 
		   gradient. beta is the variational parameter, and r is the positionvector.
		*/	   
		void jasGrad(double** ret_vec, double** r, double** r_old, int const &i_upd) const;
		/*
		   Calculate laplacian. beta is the variational parameter.
		   */
		double jasLapl(double** r) const;
		/*
			Returns the derivative of the jastrow w.r.t. beta.
			Analytical expression:
			\sum_{i<j} \frac{ a r_{12}^2 }{ ( 1 + \beta r_{ij} )^2 }
		   */
		double getdPdA();
		//Testing
		void print();
		
		//const double* iLen(int i_upd);
		//const double** allLen();
		//const double* iInvlen(int i_upd);
		//friend: jast?

		friend class popControl;
};
#endif
