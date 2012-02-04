/*
Contains orbitals that is energy eigenfunctions.
Each object an orbital with quantum nr. n,ml and spin+-1.
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
		//Functions that returns different object properties
		const int angularMomentum();
		const int energyLevel();
		const bool spinUp();
		double valueWF(double*);
		//tested for f = x.*x and f = sin (x.'x)
		double D1(double*, int);
		double D2(double*);
		
	//Make method that returns energies
};

