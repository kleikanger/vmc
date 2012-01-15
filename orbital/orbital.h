/*
Contains orbitals that is energy eigenfunctions.
Each object an orbital with quantum nr. n,ml and spin+-1.
*/


class orbital {
	
	private:
		//variables;
		int energy_level;
		int angular_momentum;
		bool spin_up;

		//methods
		double orbitalWavefunctions(double*);

	public:
		//orbitals(int <energy_level>, int <angular_momentum>, bool <spin_up>
		orbital(int,int,bool);
		//using const for better opimalization
		//Functions that returns different object properties
		const double valueWF(double*);
		const int angularMomentum();
		const int energyLevel();
		const bool spinUp();
};

