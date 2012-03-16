
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

class popControl {
	private:
		int num_part;
		int n_min_one;
		int spin_up_cutoff;
		int dimension;
	public: 
		popControl(int num_part, int spin_up_cutoff, int dimension);
		~popControl();
		void cloneWalker(walker *parent, walker *child);

		//testing
		void print(walker* quantum_dot);
};
#endif
