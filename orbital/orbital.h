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

	public:
	//Constructors
		orbital();
		//(int <energy_level>, int <angular_momentum>, bool <spin_up>)
		orbital(int,int,bool);
	
		//methods	
		//setValues(int <energy_level>, int <angular_momentum>, bool <spin_up>
		void setValues(int,int,bool);
		//Functions that returns different object properties
		//using const for better opimalization
		const int angularMomentum();
		const int energyLevel();
		const bool spinUp();
		double valueWF(double*);
		//XXX NOT TESTED  XXX
		double wFDeriv1(double*);
		double wFDeriv2(double*);

		
	//Make method that returns the energy of the orbital
};

