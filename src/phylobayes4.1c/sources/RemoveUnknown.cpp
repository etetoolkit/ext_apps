#include "phylo.h"

int main(int argc, char* argv[])	{

	string data = argv[1];
	string taxon = argv[2];
	string output = argv[3];

	MCParameters* mParam = new MCParameters();

	mParam->ReadDataFromFile(data);
	int taxindex = 0;
	while ((taxindex < mParam->Ntaxa) && (mParam->SpeciesNames[taxindex] != taxon)) taxindex++;
	if (taxindex == mParam->Ntaxa)	{
		cerr << "error: did not find taxon\n";
		exit(1);
	}

	int* Mask = new int[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		Mask[i] = (mParam->Data[taxindex][i] != unknown);
	}

	mParam->WriteDataToFileSiteMask(output,Mask);
}

