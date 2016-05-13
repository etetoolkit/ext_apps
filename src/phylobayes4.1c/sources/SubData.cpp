#include "phylo.h"

const char* alphabet = AminoAcids;

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "makeconcat list output\n";
		exit(1);
	}
	string datafile = argv[1];
	int min = atoi(argv[2]);
	int max = atoi(argv[3]);
	string outname = argv[4];

	MCParameters mParam;
	mParam.ReadDataFromFile(datafile);
	int Nsite = mParam.Nsite;
	int* mask = new int[Nsite];
	for (int i=0; i<Nsite; i++)	{
		mask[i] = 0;
	}
	for (int i=min; i<max; i++)	{
		mask[i] = 1;
	}

	mParam.WriteDataToFileSiteMask(outname,mask);
}

	

