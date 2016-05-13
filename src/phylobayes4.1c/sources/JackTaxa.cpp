#include "phylo.h"

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	int min = atoi(argv[2]);
	int max = atoi(argv[3]);
	int step = atoi(argv[4]);
	int nrep = atoi(argv[5]);
	string outfile = argv[6];

	MCParameters* mParam = new MCParameters();
	mParam->ReadDataFromFile(datafile);
	int* mask = new int[mParam->Ntaxa];

	Random::Random();
	
	for (int n=min; n<=max; n+=step)	{

		int* index = new int[n];
		for (int rep=0; rep<nrep; rep++)	{
	
			Random::DrawFromUrn(index,n,mParam->Ntaxa);
			for (int j=0; j<mParam->Ntaxa; j++)	{
				mask[j] = 0;
			}
			for (int i=0; i<n; i++)	{
				mask[index[i]] = 1;
			}
			ostringstream s;
			s << outfile << n << "_" << rep << ".ali";
			string name = s.str();
			mParam->WriteDataToFile(name,mask);
		}
		delete[] index;
	}
}

