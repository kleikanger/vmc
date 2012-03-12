/*
   Contains the wave functions's of the orbitals.
   Each object an orbital with quantum nr. n,ml and spin+-1.
   methods that returns laplacian + grad + object properties.
*/

class orbital {
	
	private:
	//variables;
		int energy_level;
		int angular_momentum;
		int dim;
		bool spin_up;
		double alpha;
		double omega;
		double sq_omg_alp;
		double omg_alp;
		//int n_x
		//int n_y

	public:
		orbital();
		//(int energy_level, int angular_momentum, bool spin_up, int dim)
		orbital(int,int,bool,int);
		//setValues(int energy_level, int angular_momentum, bool spin_up, int dim)
		void setValues(int,int,bool,int);
		//Set variational parameter alpha
		void setOmgAlp(double,double);
		//input: position r of particle
		//returns value of orbital in r.
		double valueWF(double*) const;
		//input: position r of particle, axis
		//returns gradient in r along axis
		double D1(double*, const int&) const;
		//input: position r of particle
		//returns laplacian in r.
		double D2(double*) const;
		//Return wavefunction derived w.r.t. alpha
		double valuedPdA(double* dR);
		//setting sq_omg_alp and omg_alp
		inline void updOmgAlp(double alphaARG, double omegaARG);
		//Functions that returns different object properties
		int angularMomentum() const;
		int energyLevel() const;
		bool spinUp() const;
};
