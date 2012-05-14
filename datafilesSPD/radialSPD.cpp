
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
//#include <mpi.h>
#include <cblas.h>

using namespace std;
/*
using std::cout;
using std::endl;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
*/
#define PI 3.1415926535897932

void readFromFile(int i, char *ifilename, double *in_arr, int pos, struct stat result);
//void writeToFile(char *ofilename, int n_bins, double *n, double *bins);

int main (int argc, char *argv[])
{
	//
	// XXX MUCH BETTER to use openMP with shared memory
	//
	// ************ MPI INIT **************
	//int myrank, nprocs;
    //MPI_Init(&argc, &argv);
   	//MPI_Status status;
    //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	// ************************************
	
	if (argc != 3)
	{
		cout << "Bad Usage: " << argv[0] << " " << endl;
	}
	double *x, *y, *w, *r, *bins, *n;
	double binsize, max_elem, scale;
	char *ifilename, *ofilename;
	struct stat result;
	int n_bins, i, j;
	long int num_c;
	ostringstream ost_t;

	//read variables from commandline
	ifilename = argv[1];
	n_bins = atoi(argv[2]);
	ofilename = ifilename;

	//open ofilestream	
	ofstream ofile;
	ostringstream ost;
	ost << ofilename << ".txt";
	cout<<"\nwriting SPD results to:"<<(ost.str().c_str())<<".\n";
	ofile.open((ost.str().c_str()));

	//infile size
	ost_t << ifilename << "0.dat" ;
	if (stat((ost_t.str().c_str()),&result) == 0)
	{
		num_c=result.st_size/sizeof(double);
	} 
	else
	{
		cout<<"\nerror reading from file\n";
	}
	ost_t.clear(); ost_t.str("");

	//XXX XXX XXX XXX
	num_c-=1;

	x = new double[num_c];
	y = new double[num_c];
	bins = new double[n_bins];
	n = new double[n_bins];
	
	//loading datafiles to x,y array
	readFromFile(0, ifilename, x, 0, result);
	readFromFile(1, ifilename, y, 0, result);

	//calculate and store results (in outfile)
	//calculate |r| vector (x<-sqrt(x^2+y^2))

#if 0
	for (i=0;i<num_c;i++)
		x[i]=x[i]*x[i];
	for (i=0;i<num_c;i++)
		x[i]+=y[i]*y[i];
	for (i=0;i<num_c;i++)
		x[i]=sqrt(x[i]);

	//*********************
	// FIND RADIAL DENSITY
	//*********************
	//rename wectors for a tidyer code, save memory by reusing the allocated memory
	r = x;
	//XXX w[j] = 1 for all j if VMC, if DMC weights are stored during simulation
	w = y;
	//x=NULL,y=NULL;
	//loading datafile to w array
	readFromFile(2, ifilename, w, 0 ,result);
	//find largest element in r
	i=cblas_idamax(num_c,r,1);
	max_elem=r[i];
	//bins: first from 0 to binsize, second from binsize to 2*binsize, ...
	binsize=max_elem/(double)n_bins;
	for (i=0;i<n_bins;i++)
		bins[i]=(double)(i+1)*binsize-binsize/2.;//take midpoint in bin
	//find density
	for (i=0;i<n_bins;i++)
		n[i]=0;
	for (i=0;i<num_c;i++)
	{
		j=(int)(r[i]/binsize);
		n[j]+=w[i];
	}
	//find density ac function of |r|: divide by \pi(r_{i}^2-r_{i+1}^2)
	for (i=0;i<n_bins;i++)
		n[i]=n[i]/( PI*binsize*binsize*((double)2.*i+1.) );
	//normalize 
	scale=cblas_dasum(n_bins,n,1);
	scale=n_bins/(scale*max_elem);
	cblas_dscal(n_bins,scale,n,1);
	//Write to file
	for (i = 0; i<n_bins; i++) 
	{
		ofile<<setw(16)<<setprecision(16)
			<<n[i]<<"\t"<< bins[i] << "\n";
	}
	ofile.close();
#endif

	//*************************
	// FIND DENSITY IN X,Y
	//*************************
#if 1
	w = new double[num_c];
	readFromFile(2, ifilename, w, 0 ,result);
	//find min,max elem
	/*
	i=cblas_idamax(num_c,x,1);
	double x_max=x[i];
	//cout<<x_max;
	//i=cblas_idamin(num_c,x,1);
	//double x_min=x[i];
	double x_min=x_max;
	i=cblas_idamax(num_c,y,1);
	double y_max=y[i];
	//cout<<y_max;
	//i=cblas_idamin(num_c,y,1);
	//double y_min=y[i];
	double y_min=y_max;
*/

	double x_max,y_max,x_min,y_min;
	x_max=y_max=x_min=y_min=57;
	//find bin intervals in x and y coords
	double x_binsize = (fabs(x_min)+fabs(x_max))/(double)n_bins;
	double y_binsize = (fabs(y_min)+fabs(y_max))/(double)n_bins;

	//initialize axis and grid
	double x_axis[n_bins];
	double y_axis[n_bins];
	double grid[n_bins][n_bins];
	for (i=0;i<n_bins;i++)
		x_axis[i]=(double)(i+1)*x_binsize-x_binsize/2.-fabs(x_min);//take midpoint in bin
	for (i=0;i<n_bins;i++)
		y_axis[i]=(double)(i+1)*y_binsize-y_binsize/2.-fabs(y_min);//take midpoint in bin
	for (i=0;i<n_bins;i++)
		for (j=0;j<n_bins;j++)
		{
			grid[i][j]=0.0;
		}
	//sort data
	for (int k=0;k<num_c;k++)
	{
		i=(int)((x[k]+fabs(x_min))/x_binsize);
		j=(int)((y[k]+fabs(y_min))/y_binsize);
		grid[i][j]+=w[k];
	}
	//normalize
	scale=0.0;
	for (i=0;i<n_bins;i++)
		for (j=0;j<n_bins;j++)
		{
			scale+=grid[i][j];
		}
	for (i=0;i<n_bins;i++)
		for (j=0;j<n_bins;j++)
		{
			grid[i][j]/=scale;
		}
	
	//Write to file
	for (i = 0; i<n_bins; i++) 
		ofile<<x_axis[i]<<"\t";
//		ofile<<"\n";
	for (i = 0; i<n_bins; i++) 
//		ofile<<y_axis[i]<<"\t";
//		ofile<<"\n";
	for (i = 0; i<n_bins; i++) 
	{
		for (j=0;j<n_bins;j++)
		{
			if (grid[i][j]<.001)
			ofile<<setw(16)<<setprecision(16)
				//<<x_axis[i]<<"\t"<<y_axis[j]<<"\t"
				<<grid[i][j]<<"\t";
		}
		ofile<<"\n";
	}
	ofile.close();

#endif
	//MPI_Reduce(n)
	//MPI_Reduce(bins)
	//Write results to file
	//if (myrank==0)
	//writeToFile(ofilename, n_bins, n, bins); 
	//MPI_Finalize();

	delete [] x;
	delete [] y;
	//delete [] w;
	//delete [] r;
	//delete [] bins;
	//delete [] n;
}
void readFromFile(int i, char *ifilename, double *in_arr, int pos, struct stat result) 
{
	ifstream infile;
	ostringstream ost;
	ost << ifilename << i << ".dat";
	cout<<"reading data from: "<<ost.str().c_str()<<endl;;
	infile.open(ost.str().c_str(),ios::in|ios::binary);
	infile.read((char*)&(in_arr[pos]),result.st_size); //TODO w+1??
	infile.close();
	//ost.clear(); ost.str("");
}
/*
void writeToFile(char *ofilename, int n_bins, double *n, double *bins) 
{
	for (int i = 0; i<n_bins; i++) 
	{
		cout<<setw(16)<<setprecision(16)
			<<n[i]<<"\t"<< bins[i] << endl;
	}
	cout<<"\n";
	ofstream ofile;
	ostringstream ost;
	ost << ofilename << ".txt";
	cout<<"\nwriting SPD results to:"<<(ost.str().c_str())<<".\n";
	ofile.open((ost.str().c_str()));
if (ofile.is_open())
{
	for (int i = 0; i<n_bins; i++) 
	{
		ofile<<setw(16)<<setprecision(16)
			<<n[i]<<"\t"<< bins[i] << "\n";
	}
	ofile.close();
}
else cout<<"could not open file\n";
}*/
