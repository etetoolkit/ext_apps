#include "phylo.h"

int main(int argc, char* argv[])	{
	
	int N;
	ifstream is(argv[1]);
	string filename = argv[2];
	double obs = atof(argv[3]);
	int ncat = 20;

	is >> N;
	double* val = new double[N];
	for (int i=0; i<N; i++)	{
		is >> val[i];
	}
	MakeHisto(val,N,filename,ncat);
	double mean = 0;
	double var = 0;
	double pp = 0;
	for (int i=0; i<N; i++)	{
		mean += val[i];
		var += val[i] * val[i];
		if (val[i] > obs) pp++;
	}
	mean /= N;
	var /= N;
	var -= mean*mean;
	pp /= N;
	cout << '\n';
	cout << "mean : " << mean << '\n';
	cout << "z    : " << (mean - obs) / sqrt(var) << '\n';
	cout << "pv   : " << pp << '\n';
	cout << '\n';
}

