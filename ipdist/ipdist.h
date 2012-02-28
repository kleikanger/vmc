/*
	keeps track of all interparticle distanses. for n particles n*(n-1)/2 elements 
	are needed, and n-1 elements needs to be updated every time a patricle is moved.
	update O(n-1) instead of O(n(n-1)/2).

	class keeping two matrices one for the interparticle distances r_ij, 
	and one for the inverse elements 1/r_ij. 

	does not check if 1/r_ij is singular. that could in principle be a problem since 
	a jump to a configuration where r_i=r_j will always be accepted by the metropolis-algo
	since \Psi^{new}(r_i=r_j) = 0. iwhen calculating the new r_{i_upd,j} vec, one could
	check if any elements are 0. Chances are very small for this to happen \sim 10^-16 
	/ mc-cycle.

*/
class ipdist{

	private: 

		//number of particles -1
		int n_min_one;
		//dimension
		int dim;
		//spin down particles: only for calculating jastrow grad and lapl
		int iCutoff;
		
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
		// XXX remove ip_invlen !!
		// inverse elements (1/r_ij)
		double** ip_invlen;
		// elements (r_ij)
		double** ip_len_backup;
		// inverse elements (1/r_ij)
		double** ip_invlen_backup;

	public: 
		/*
		   declare and alloc ip_len and ip_invlen.
		   set n_min_one to n-1, dim to di and iCutoff to iC 
		 */
		ipdist(int n, int di, int iC);
		/*
		   init all elements in ip_len and ip_invlen
		   NB! one should check if any of these elements are 0.
		 */
		void init(double** x);
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
		double const sumInvlen();
		/*
		   Sum all elements where one of the indices = i_upd
		   		( 
		 		r_{i_upd,0}
				,...
				,r_{i_upd,i_upd-1}
				,r_{i_upd,i_upd+1}
				,...
				,r_{i_upd,n-1} 
				) 	
		 */
		double const sumPart(int i_upd);
		/*
		   returns R. exp(R_new-R_old) = the jastrow ratio when only one 
		   particle r_{i_upd}  is moved.
		   */
		double const logJasR(int i_upd, double beta);
		/*
		   returns gradient of jastrow. i_upd part. in ret_vec will be changed to new 
		   gradient. beta is the variational parameter, and r is the positionvector.
		*/	   
		void const jasGrad(double** ret_vec, double beta, double** r, int i_upd);
		/*
		   Calculate laplacian. beta is the variational parameter.
		   */
		double const jasLapl(double beta, double** r);
		/*
		   clear all malloced vars
		 */
		void clear();
		//Testing
		void print();
		
		//const double* iLen(int i_upd);
		//const double** allLen();
		//const double* iInvlen(int i_upd);
		//friend: jast?
};
