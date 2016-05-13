#include  "phylo.h"


int main(int argc, char* argv[])	{

	if (argc != 5)	{
		cerr << "jackboot datafile N Nrep outfile\n";
		exit(1);
	}

	// Random::Random();
	string name = argv[1];
	int N = atoi(argv[2]);
	int Nrep = atoi(argv[3]);
	string outfile = argv[4];
	
	MCParameters* mParam = new MCParameters();
	mParam->ReadDataFromFile(name);
	
	if (N > mParam->Nsite)	{
		cerr << "error : target size should be < Nsite\n";
		exit(1);
	}

	int* mask = new int[mParam->Nsite];
	int* index = new int[N];
	
	for (int rep=0; rep<Nrep; rep++)	{
		for (int i=0; i<mParam->Nsite; i++)	{
			mask[i] = 0;
		}
		Random::DrawFromUrn(index,N,mParam->Nsite);
		for (int i=0; i<N; i++)	{
			mask[index[i]] = 1;
		}
		ostringstream s;
		s << outfile << rep << ".ali";
		mParam->WriteDataToFileSiteMask(s.str(),mask);
	}

}

