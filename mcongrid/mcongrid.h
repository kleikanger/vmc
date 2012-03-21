class mcongrid {
	private:
	public:
		mcongrid();
		~mcongrid();
		void runMCongrid(
			int dimension,
			int num_cycles,
			int thermalization,
			int num_part,
			int spin_up_cutoff,
			double delta_t,
			double* var_par,
			double* var_par_inc,
			int* var_par_cyc,
			int num_of_var_par,
			int myrank, 
			int nprocs );
};
