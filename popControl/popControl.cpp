#include "popControl.h"

popControl::popControl(int num_part_ARG, int spin_up_cutoff_ARG, int dimension_ARG)
{
	num_part=num_part_ARG;
	n_min_one=num_part_ARG-1;
	spin_up_cutoff=spin_up_cutoff_ARG;
	dimension=dimension_ARG;
}
popControl::~popControl(){}
void popControl::cloneWalker(walker *parent, walker *child)
{

	//NB: Only send spin down if spin_up_cutoff<num_part

	int i,j;
	//copy ipdist matrix from parent to child
	for (i=0; i<n_min_one; i++)
		for (j=0; j<=i; j++)
			child->ipd->ip_len[i][j]=parent->ipd->ip_len[i][j];
	//copy ipdist backup matrix from parent to child
	for (i=0; i<n_min_one; i++)
		for (j=0; j<=i; j++)
			child->ipd->ip_len_backup[i][j]=parent->ipd->ip_len_backup[i][j];
	//copy slater matrixes from parent to child
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->spin_up_matr[i][j]=parent->slater->spin_up_matr[i][j];
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->spin_down_matr[i][j]=parent->slater->spin_down_matr[i][j];
	//copy slater backup matrixes from parent to child
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->spin_up_backup[i][j]=parent->slater->spin_up_backup[i][j];
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->spin_down_backup[i][j]=parent->slater->spin_down_backup[i][j];
	//copy inverse slater matrixes from parent to child
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->inv_up_matr[i][j]=parent->slater->inv_up_matr[i][j];
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->inv_down_matr[i][j]=parent->slater->inv_down_matr[i][j];
	//copy inverse slater backup matrixes from parent to child
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->inv_up_backup[i][j]=parent->slater->inv_up_backup[i][j];
	for (i=0; i<spin_up_cutoff; i++)
		for (j=0; j<spin_up_cutoff; j++)
			child->slater->inv_down_backup[i][j]=parent->slater->inv_down_backup[i][j];
	//copy position vectors from parent to child
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->r_new[i][j]=parent->r_new[i][j];
	//copy position vectors from parent to child
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->r_old[i][j]=parent->r_old[i][j];
	//copy gradient of the jastrow function from parent to child
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->jas_grad[i][j]=parent->jas_grad[i][j];
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->jas_grad_bu[i][j]=parent->jas_grad_bu[i][j];
	//copy quantum force from parent to child
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->q_force_new[i][j]=parent->q_force_new[i][j];
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->q_force_old[i][j]=parent->q_force_old[i][j];
	//copy gradient of the slater matrix from parent to child
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->sla_grad[i][j]=parent->sla_grad[i][j];
	for (i=0;i<num_part;i++)
		for (j=0;j<dimension;j++)
			child->sla_grad_bu[i][j]=parent->sla_grad_bu[i][j];
	
	//copy orbitals from parent to child would be necc for open shell model
	//copy alpha, beta from parent to child would be necc for open shell model
}
//void popControl::transmitWalker(walker *parent, int parent_rank, walker *child, child_rank)
//MPI_Send and MPI_Recv all parameters. 

//Testing
void popControl::print(walker* quantum_dot)
{
	quantum_dot->slater->print();
	quantum_dot->ipd->print();
}	
