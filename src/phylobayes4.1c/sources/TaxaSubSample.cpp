#include "phylo.h"

int main(int argc, char* argv[])	{

	if (argc==1)	{
		cerr << "subtaxa datafile taxlist outputfile\n";
		exit(1);
	}

	// initialise random 
	Random::Random();

	string datafile = argv[1];
	MCParameters* mParam = new MCParameters;
	mParam->ReadDataFromFile(datafile);
	
	int N;
	ifstream is(argv[2]);
	is >> N;
	// cerr << "number of taxa in datafile : " << mParam->Ntaxa << '\n';
	// cerr << "number of taxa : " << N << '\n';
	string names[N];
	for (int i=0; i<N; i++)	{
		is >> names[i];
	}
 	if ((N<=0) || (N>mParam->Ntaxa)) {
		cerr << "error : bad number of taxa\n";
		exit(1);
	}
	int* mask = new int[mParam->Ntaxa];
	for (int i=0; i<mParam->Ntaxa; i++)	{
		mask[i] = 0;
	}
	for (int i=0; i<N; i++)	{	
		int k = 0;
		while ((k<mParam->Ntaxa) && (names[i] != mParam->SpeciesNames[k]))	{
			k++;
		}
		if (k==mParam->Ntaxa)	{
			cerr << "error : do not find " << names[i] << '\n';
			exit(1);
		}
		mask[k] = 1;
	}
	string output = argv[3];
	mParam->WriteDataToFile(output, mask);	
}

