/*
	keeps track of all interparticle distanses. for n particles n*(n-1)/2 elements 
	are needed, and n-1 elements needs to be updated every time a patricle is moved.
	update O(n-1) instead of O(n(n-1)/2).

	class keeping two matrices one for the interparticle distances r_ij, 
	and one for the inverse 1/r_ij. 

	does not check if 1/r_ij is singular. that could be a problem since 
	a jump to a configuration where r_i=r_j will be accepted by the metropolis-algo
	since \Psi^{new}(r_i=r_j) = 0. when finding the new r_{i_upd,j} vec, one should
	check if any elements are 0.

*/
class ipdist{

	private: 

		//n-1
		int n_min_one;
		//dimension
		int dim;

	//protected: //private within this class and friends
		
		//lower triangular with no diag elements n(n-1)/2 elem.
		//ip_dist = new double[n_min_one];
		//for (int i=0; i<n_min_one; i++)) ip_dist[i] = new double[i+1];
		/*

		   ( 	0 		0 		0 		0 		 0)
		   ( 	r_10 	0 		0 		0 		 0)
		   ( 	r_20 	r_21 	0 		0 		 0)
		   ( 	. 		. 		... 	0 		 0)
		   ( 	. 		.       ... 	0 		 0)
		   ( 	r_n0 	r_n1 	... 	r_n(n-1) 0)
			
		   r_{ij} = |\vec r_i - \vec r_j|
		 
		   OBS: element i,j corresponds to particle r_(i+1,j)   
		 
		 */

		// (r_ij)
		double** ip_len;
		// ip_len with inverse elements (1/r_ij)
		double** ip_invlen;

	public: 
		/*
		   declares and alloc. ip_len and ip_invlen.
		   set n_min_one to n-1, and dim to di
		 */
		ipdist(int n, int di);
		/*
		   init all elements in ip_len and ip_invlen
		   
		   OBS Not testing if elements are 0 !!
		 
		 */
		void init(double** x);
		/*
		   updates n-1 elements of ip_len and ip_invlen
		   u_upd column and row.
		  
		   input( 
		 		x_{i_upd,0}
				,...
				,x_{i_upd,i_upd-1}
				,x_{i_upd,i_upd+1}
				,...
				,x_{i_upd,n-1} 
				) 	

			NB! one should check if any of these elements are 0.
		 */
		void update(double* r, int i_upd);
		/*
		   Return sum of inverse elements. e-e electrostatic energy = ee/sum
		 */	   
		const double sumInvlen();
		/*
		   Sum all elements
		   		( 
		 		r_{i_upd,0}
				,...
				,r_{i_upd,i_upd-1}
				,r_{i_upd,i_upd+1}
				,...
				,r_{i_upd,n-1} 
				) 	
		 */
		const double sumPart(int i_upd);
		/*
		   clear all malloced vars
		 */
		void clear();
		//Testing
		void print();
		
		//returns pointers to ...
		//const double* iLen(int i_upd);
		//const double** allLen();
		//const double* iInvlen(int i_upd);

};
/*
class jastrow{
	private: 
		//
		int i_cutoff
	
	public:
		//
};
*/
