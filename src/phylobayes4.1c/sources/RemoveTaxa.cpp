#include "phylo.h"

int main(int argc, char* argv[])	{

	if (argc==1)	{
		cerr << "rmtaxa datafile name1 name2 ... namep outputfile\n";
		cerr << "makes a new datamatrix by selecting a random subsample of N taxa\n";
		exit(1);
	}

	int P = argc- 3;
	string taxa[P];
	for (int j=0; j<P; j++)	{
		taxa[j] = argv[j+2];
	}

	// initialise random 
	Random::Random();

	string datafile = argv[1];
	MCParameters* mParam = new MCParameters;
	mParam->ReadDataFromFile(datafile);
	
	int* mask = new int[mParam->Ntaxa];
	for (int i=0; i<mParam->Ntaxa; i++)	{
		mask[i] = 1;
	}
	for (int i=0; i<P; i++)	{	
		int j=0; 
		while ((j<mParam->Ntaxa) && (mParam->SpeciesNames[j] != taxa[i]))	{
			j++;
		}
		if (j==mParam->Ntaxa)	{
			cerr << "error : did not find taxon : " << taxa[i] << '\n';
			exit(1);
		}
		mask[j] = 0;
	}
	string output = argv[argc-1];
	mParam->WriteDataToFile(output, mask);	

}

