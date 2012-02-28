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

		//int n_x
		//int n_y

	public:
		orbital();
		//(int energy_level, int angular_momentum, bool spin_up, int dim)
		orbital(int,int,bool,int);
		//setValues(int energy_level, int angular_momentum, bool spin_up, int dim)
		void setValues(int,int,bool,int);
		//Set variational parameter alpha
		void setAlpha(double);
		//input: position r of particle
		//returns value of orbital in r.
		double const valueWF(double*);
		//input: position r of particle, axis
		//returns gradient in r along axis
		double const D1(double*, int);
		//input: position r of particle
		//returns laplacian in r.
		double const D2(double*);
		//Functions that returns different object properties
		int const angularMomentum();
		int const energyLevel();
		bool const spinUp();
};

