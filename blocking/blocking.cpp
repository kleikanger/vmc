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

using std::cout;
using std::endl;
using std::ostringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;

void blocking(double* indata , int num_cxnum_bs , int size_bs, double* res);
void meanvar(double *blockindata, int nblocks, double *res);
double mean(double* res,int size_bs);

int main (int argc, char *argv[])
{
	if (argc != 6)
	{
		cout << "Bad Usage: " << argv[0] << " " << endl;
	}
	double *indata, *res;
	int min_bs, num_bs, int_bs, num_c, n_proc, size_bs;
	double mean, sigma;
	char *ifilename, *ofilename;
	struct stat result;
	ostringstream ost, ost_t;
	ifstream infile;
	ofstream outfile;

	//read variables from commandline
	indata = new double[num_c];
	min_bs = atoi(argv[1]);	
	num_bs = atoi(argv[2]);
	int_bs = atoi(argv[3]);
	n_proc = atoi(argv[4]);
	ifilename = argv[5];
	ofilename = ifilename;//argv[6];

	//assuming that the infiles have the same size
	ost_t << ifilename << "0.dat" ;
	if (stat((ost_t.str().c_str()),&result) == 0)
	{
		num_c=result.st_size/sizeof(double);
	}
	ost_t.clear();
	//loading datafiles to indata array
	indata = new double[num_c*n_proc];
	for (int i=0; i<n_proc; i++)
	//for (int i=0+myrank; i<n_proc; i+=2)
	{
		ost << ifilename << i << ".dat";
		cout<<"reading data from: "<<ost.str().c_str()<<endl;;
		infile.open(ost.str().c_str(),ios::in|ios::binary);
		infile.read((char*)&(indata[i*num_c]),result.st_size);
		infile.close();
		ost.clear(); ost.str("");
	}
	
	
	//calculate and store results (in outfile)
	res = new double[2];
	ost << ofilename  << ".txt";
	outfile.open(ost.str().c_str());
	cout<<"writing blocking results to: "<<(ost.str().c_str())<<endl;
	for (int i = 0; i<num_bs ; i++) 
	//for (int i = 0+myrank; i<num_bs ; i++) //WILL WORK IF MYRANK IS EVEN
	{
		size_bs = min_bs+i*int_bs;
		blocking(indata, num_c * n_proc, size_bs, res);
		fflush(stdout);
		cout<<"\r"<<"                                        ";
		cout<<"\rprocessing block "<<i+1<<" of "<<num_bs;
		mean = res[0];
		sigma = res[1];
//MPI SAVE TO
		outfile<<setw(16)<<setprecision(16)
			<<size_bs<<"\t"<< mean << "\t"
			<< sqrt(sigma / ((num_c*n_proc) / ((double)size_bs)-1.0))
			<< endl;
	}
	outfile.close();
}	
//	MPI_Reduce(&in, &out,
//					1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


void blocking(double *indata, int size_indata, int size_bs, double *res)
{
	int nblocks = size_indata / size_bs;
	double* blockindata = new double[nblocks];
	for ( int i =0; i <nblocks ; i ++)
		blockindata[i] = mean(indata+i*size_bs,size_bs);
	meanvar(blockindata,nblocks,res);
	delete [] blockindata;
}

//calculate mean energy in one block.
double mean(double* indata,int size_bs)
{
	double mean_ = 0.0;
	for (int i=0; i<size_bs; i++)
		mean_ += indata[i];
	mean_ /= (double)size_bs;
	return mean_;
}

//calculate mean of variance and sigma for all blocks.  
void meanvar(double *blockindata, int nblocks, double *res)
{
	int i;
	double m=0, s=0;
	for (i=0; i<nblocks; i++)
		m += blockindata[i];
	m /= (double)nblocks;
	for (i=0; i<nblocks; i++)
		s += pow(blockindata[i],2);
	s /= ((double)nblocks);
	s = s-pow(m,2);
	res[0]=m, res[1]=s;
//	cout<<res[0]<<" \n";
}






