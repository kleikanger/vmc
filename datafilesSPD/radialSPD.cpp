
/*

   Blocking analysis of vmc data.
   ./<blockingAnalysis>.out min_bs num_bs int_bs n_proc infilename 
   
	min_bs 		- Blocksize of smallest block 	
	num_bs 		- Number of blocks
	int_bs 		- Increment in blocksize for each block
	n_proc 		- Number of infiles
	ifilename 	- infilename : ex: for blocking_<i>.dat write blocking 
	
	Writes to outfile with the name
	ofilename = <infilename>.txt

   */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpi.h>
#include <cblas.h>

using std::cout;
using std::endl;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;

int main (int argc, char *argv[])
{
	//
	// XXX MUCH BETTER to use openMP with shared memory
	//
	// ************ MPI INIT **************
	int myrank, nprocs;
    MPI_Init(&argc, &argv);
   	MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	// ************************************
	
	if (argc != 6)
	{
		cout << "Bad Usage: " << argv[0] << " " << endl;
	}
	double *x, *y, *w, *r, *bins, *n;
	double binsize, max_elem;
	char *ifilename, *ofilename;
	struct stat result;
	int n_bins, i, j;
	ostringstream ost, ost_t;
	ifstream infile;
	ofstream outfile;

	//read variables from commandline
	n_bins = atof(argv[1]);
	ifilename = argv[2];
	ofilename = ifilename;//argv[6];

	//infile size
	ost_t << ifilename << "0.dat" ;
	if (stat((ost_t.str().c_str()),&result) == 0)
	{
		num_c=result.st_size/sizeof(double);
	}
	ost_t.clear();
	
	x = new double[num_c];
	y = new double[num_c];
	//XXX w[j] = 1 for all j if VMC, if DMC weights are stored during simulation
	w = new double[num_c];
	r = new double[num_c];
	bins = new double[num_c];
	n = new int[num_c];
	
	double **temp = new double[3];
	temp[0] = x;
	temp[1] = y;
	temp[2] = w;
	
	//loading datafiles to indata array
	indata = new double[num_c];
	for (int i=0; i<3; i++)
	//for (int i=0+myrank; i<n_proc; i+=2)
	{
		ost << ifilename << i << ".dat";
		cout<<"proc "<<myrank<<" reading data from: "<<ost.str().c_str()<<endl;;
		infile.open(ost.str().c_str(),ios::in|ios::binary);
		infile.read((char*)&(temp[i]),result.st_size);
		infile.close();
		ost.clear(); ost.str("");
	}
	delete [] temp;

	//calculate and store results (in outfile)
	//calculate r vector
	for (i=0;i<num_c;i++)
		r[i]=x[i]*x[i]+y[i]*y[i];
	for (i=0;i<num_c;i++)
		r[i]=sqrt(r);
	//find largest element in r
	i=cblas_idamax(num_c,r,1);
	max_elem=r[i];
	//bins: first from 0 to binsize, second from binsize to 2*binsize, ...
	binsize=(double)i*max_elem/n_bins;
	for (i=1;i<=bins;i++) //write this array to file for plotting
		bins[i]=(double)i*binsize-binsize/2.;//take midfoint in bin
	//find density
	for (i=1;i<=num_c;i++)
		n[i]=0;
	for (i=1;i<=num_c;i++)
	{
		j=(int)r[i]/binsize;
		n[j]+=w[j];
	}
#define pi 3.141592//......
	//find density ac function of |r|
	for (i=1;i<=num_c;i++)
		n[j]=n[j]/(bins[j]*2.*PI)
	//normalize
	for (i=1;i<=num_c;i++)
		n[j]=n_bins/(cblas_sum(n)*max_elem);

	//Write results to file
	if (myrank==0)
	{
		ost << ofilename  << ".txt";
		outfile.open(ost.str().c_str());
		cout<<"\nwriting blocking results to: "<<(ost.str().c_str())<<endl;
		for (int i = 0; i<num_bs; i++) 
		{
			cout<<num_bs<<" "<<arr_size_bs[i]<<" "<<arr_blc_res[i]<<" "<<i<<"\n";
			outfile<<setw(16)<<setprecision(16)
				<<arr_size_bs[i]<<"\t"<< arr_mean[i] 
				<< "\t"<<arr_blc_res[i]<< endl;
		}
		cout<<"\n";
		outfile.close();
	}
	MPI_Finalize();
}

