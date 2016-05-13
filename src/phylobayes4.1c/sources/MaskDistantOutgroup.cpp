#include  "phylo.h"


int main(int argc, char* argv[])	{

	if (argc != 6)	{
		cerr << "maskog datafile taxfile outfile\n";
		exit(1);
	}

	// Random::Random();
	string name = argv[1];
	string outfile = argv[3];
	
	ifstream is(argv[2]);
	int Ntaxa;
	is >> Ntaxa;
	string name[Ntaxa];
	for (int i=0; i<Ntaxa; i++)	{
		is >> name[i];
	}

	MCParameters* mParam = new MCParameters();
	mParam->ReadDataFromFile(name);
	
	int* mask = new int[mParam->Nsite];
	for (int i=0; i<mParam->Nsite; i++)	{
		mask[i] = 0;
	}

	for (int i=0; i<mParam->Nsite; i++)	{
		int present = 0;
		for (int j=0; j<mParam->Ntaxa; j++)	{
			for (int k=0; k<Ntaxa; k++)	{
				if 

		Random::DrawFromUrn(index,N1+N2,mParam->Nsite);
		for (int i=0; i<N1; i++)	{
			mask[index[i]] = 1;
		}
		ostringstream s;
		s << outfile << rep << "train.ali";
		mParam->WriteDataToFileSiteMask(s.str(),mask);

		for (int i=0; i<mParam->Nsite; i++)	{
			mask[i] = 0;
		}
		for (int i=N1; i<N1+N2; i++)	{
			mask[index[i]] = 1;
		}
		ostringstream s2;
		s2 << outfile << rep << "test.ali";
		mParam->WriteDataToFileSiteMask(s2.str(),mask);
	}

}

