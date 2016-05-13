#include "phylo.h"

int main(int argc, char* argv[])	{

	if (argc==1)	{
		cerr << "subtaxa datafile N outputfile\n";
		cerr << "makes a new datamatrix by selecting a random subsample of N taxa\n";
		exit(1);
	}

	// initialise random 
	Random::Random();

	string datafile = argv[1];
	MCParameters* mParam = new MCParameters;
	mParam->ReadDataFromFile(datafile);
	
	int N;
	int contnchar;
	ifstream is(argv[2]);
	is >> N >> contnchar;
	string names[N];
	double contdata[N][contnchar];
	for (int i=0; i<N; i++)	{
		is >> names[i];
		for (int k=0; k<contnchar; k++)	{
			is >> contdata[i][k];
		}
	}

	int* mask = new int[mParam->Ntaxa];
	for (int i=0; i<mParam->Ntaxa; i++)	{
		mask[i] = 0;
	}
	for (int i=0; i<mParam->Ntaxa; i++)	{	
		int k = 0;
		while ((k<N) && (mParam->SpeciesNames[i] != names[k]))	{
			k++;
		}
		if (k<N)	{
			cerr << mParam->SpeciesNames[i] << '\t' << names[k];
			for (int k=0; k<contnchar; k++)	{
				cerr << '\t' << contdata[i][k];
			}
			cerr << '\n';
			mask[i] = 1;
		}
	}
	string output = argv[3];
	mParam->WriteDataToFile(output, mask);	

}

