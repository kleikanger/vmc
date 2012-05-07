
/*

   Write MC results to file
   ./<blockingAnalysis>.out n_bins infilename 
   
	n_bins 		- resolution of plot 	
	ifilename 	- infilename : ex: for SPD0.dat write SPD
	
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

#define PI 3.1415926535897932

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
	double binsize, max_elem, scale;
	char *ifilename, *ofilename;
	struct stat result;
	int n_bins, i, j;
	ostringstream ost, ost_t;
	ifstream infile;
	ofstream outfile;

	//read variables from commandline
	n_bins = atoi(argv[1]);
	ifilename = argv[2];
	ofilename = ifilename;

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
	for (i=1;i<=bins;i++)
		bins[i]=(double)i*binsize-binsize/2.;//take midfoint in bin
	//find density
	for (i=0;i<n_bins;i++)
		n[i]=0;
	for (i=0;i<n_bins;i++)
	{
		j=(int)r[i]/binsize;
		n[j]+=w[j];
	}
	//find density ac function of |r|: divide by \pi(r_{i}^2-r_{i+1}^2)
	for (i=0;i<n_bins;i++)
		n[j]=n[j]/((double) PI*binsize*binsize*(2*i+1) );
	//normalize 
	scale=cblas_dasum(n_bins,r,1);
	scale=n_bins/(scale*max_elem)
	cblas_dscal(n_bins,scale,n,1)
	
	//MPI_Reduce(n)
	//MPI_Reduce(bins)

	//Write results to file
	if (myrank==0)
	{
		ost << ofilename  << ".txt";
		outfile.open(ost.str().c_str());
		cout<<"\nwriting SPD results to: "<<(ost.str().c_str())<<endl;
		for (int i = 0; i<num_bs; i++) 
		{
			outfile<<setw(16)<<setprecision(16)
				<<n[i]<<"\t"<< bins[i] << endl;
		}
		outfile.close();
	}
	MPI_Finalize();
}

