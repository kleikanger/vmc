#include "popControl.h"
#include <cmath>

popControl::popControl(int num_part_ARG, int spin_up_cutoff_ARG, int dimension_ARG, MPI_Status status_ARG)
{
	num_part=num_part_ARG;
	n_min_one=num_part_ARG-1;
	spin_up_cutoff=spin_up_cutoff_ARG;
	dimension=dimension_ARG;
	status=status_ARG;
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
void popControl::transmitWalker(walker *w_trans, int parent_rank, int child_rank, int myrank)
{
//MPI_Send and MPI_Recv all parameters. 
	int i,j;
	
	//copy ipdist matrix from parent to child
	int n_ipd = num_part*n_min_one/2;

	if (myrank==parent_rank)
	MPI_Send(w_trans->ipd->ip_len[0], n_ipd, MPI_DOUBLE, child_rank, 500, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->ipd->ip_len[0], n_ipd, MPI_DOUBLE, parent_rank, 500, MPI_COMM_WORLD, &status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->ipd->ip_len_backup[0], n_ipd, MPI_DOUBLE, child_rank, 501, MPI_COMM_WORLD);
	if (myrank==child_rank)
	MPI_Recv(w_trans->ipd->ip_len_backup[0], n_ipd, MPI_DOUBLE, parent_rank, 501, MPI_COMM_WORLD ,&status);

	//copy slater matrixes from parent to child
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->spin_up_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 502, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->spin_up_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 502, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->spin_down_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 503, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->spin_down_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 503, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->spin_up_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 504, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->spin_up_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 504, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->spin_down_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 505, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->spin_down_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 505, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->inv_up_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 506, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->inv_up_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 506, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->inv_down_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 507, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->inv_down_matr[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 507, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->inv_up_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 508, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->inv_up_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 508, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->slater->inv_down_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, child_rank, 509, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->slater->inv_down_backup[0], pow(spin_up_cutoff,2), MPI_DOUBLE, parent_rank, 509, MPI_COMM_WORLD ,&status);
	
	//copy position vectors from parent to child
	if (myrank==parent_rank)
	MPI_Send(w_trans->r_new[0], num_part*dimension, MPI_DOUBLE, child_rank, 510, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->r_new[0], num_part*dimension, MPI_DOUBLE, parent_rank, 510, MPI_COMM_WORLD, &status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->r_old[0], num_part*dimension, MPI_DOUBLE, child_rank, 511, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->r_old[0], num_part*dimension, MPI_DOUBLE, parent_rank, 511, MPI_COMM_WORLD , &status);
	
	//copy gradient of the jastrow function from parent to child
	if (myrank==parent_rank)
	MPI_Send(w_trans->jas_grad[0], num_part*dimension, MPI_DOUBLE, child_rank, 512, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->jas_grad[0], num_part*dimension, MPI_DOUBLE, parent_rank, 512, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->jas_grad_bu[0], num_part*dimension, MPI_DOUBLE, child_rank, 513, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->jas_grad_bu[0], num_part*dimension, MPI_DOUBLE, parent_rank, 513, MPI_COMM_WORLD ,&status);
	
	//copy quantum force from parent to child
	if (myrank==parent_rank)
	MPI_Send(w_trans->sla_grad[0], num_part*dimension, MPI_DOUBLE, child_rank, 514, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->sla_grad[0], num_part*dimension, MPI_DOUBLE, parent_rank, 514, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->sla_grad_bu[0], num_part*dimension, MPI_DOUBLE, child_rank, 515, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->sla_grad_bu[0], num_part*dimension, MPI_DOUBLE, parent_rank, 515, MPI_COMM_WORLD ,&status);
	
	//copy gradient of the slater matrix from parent to child
	if (myrank==parent_rank)
	MPI_Send(w_trans->q_force_new[0], num_part*dimension, MPI_DOUBLE, child_rank, 516, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->q_force_new[0], num_part*dimension, MPI_DOUBLE, parent_rank, 516, MPI_COMM_WORLD ,&status);
	if (myrank==parent_rank)
	MPI_Send(w_trans->q_force_old[0], num_part*dimension, MPI_DOUBLE, child_rank, 517, MPI_COMM_WORLD);	
	if (myrank==child_rank)
	MPI_Recv(w_trans->q_force_old[0], num_part*dimension, MPI_DOUBLE, parent_rank, 517, MPI_COMM_WORLD ,&status);
}
//Testing
void popControl::print(walker* quantum_dot)
{
	quantum_dot->slater->print();
	quantum_dot->ipd->print();
}	
