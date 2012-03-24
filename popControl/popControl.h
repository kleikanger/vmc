
#ifndef POPCONTROL_H
#define POPCONTROL_H

#ifndef SLATERMATRIX_H
	#include "../QDslater/slaterMatrix.h"
#endif
#ifndef WALKER_H
	#include "../walker/walker.h"
#endif
#ifndef IPDIST_H
	#include "../ipdist/ipdist.h"
#endif
#include "mpi.h"

class popControl {
	private:
		int num_part;
		int n_min_one;
		int spin_up_cutoff;
		int dimension;
		MPI_Status status;
	public: 
		popControl(int num_part, int spin_up_cutoff, int dimension, MPI_Status status);
		~popControl();
		void cloneWalker(walker *parent, walker *child);
		void transmitWalker(walker *parent, int parent_rank, int child_rank, int myrank);

		//testing
		void print(walker* quantum_dot);
};
#endif
