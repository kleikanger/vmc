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

	public:
		orbital();
		//(int energy_level, int angular_momentum, bool spin_up, int dim)
		orbital(int,int,bool,int);
		//setValues(int energy_level, int angular_momentum, bool spin_up, int dim)
		void setValues(int,int,bool,int);
		//input: position r of particle
		//returns value of orbital in r.
		double valueWF(double*);
		//input: position r of particle, axis
		//returns gradient in r along axis
		double D1(double*, int);
		//input: position r of particle
		//returns laplacian in r.
		double D2(double*);
		//Functions that returns different object properties
		const int angularMomentum();
		const int energyLevel();
		const bool spinUp();
		
	//Make method that returns energies
};

