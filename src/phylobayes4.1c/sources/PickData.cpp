#include "phylo.h"

int main(int argc, char* argv[])	{

	Random::Random();

	ifstream is(argv[1]);
	int ntaxamin = atoi(argv[2]);
	int ntaxamax = atoi(argv[3]);
	int nsitemin = atoi(argv[4]);
	int nsitemax = atoi(argv[5]);
	int p = atoi(argv[6]);

	string dir = argv[7];

	int n = 0;
	is >> n;
	string list[n];
	int avail[n];
	for (int i=0; i<n; i++)	{
		avail[i] = 0;
	}	

	int nav = 0;
	for (int i=0; i<n; i++)	{
		is >> list[i];
		cerr << list[i] << '\n';
		MCParameters mParam;
		mParam.ReadDataFromFile(list[i]);
		cerr << mParam.Ntaxa << '\t' << mParam.Nsite << '\n';
		if ((mParam.Nsite >= nsitemin) && (mParam.Nsite <= nsitemax) && (mParam.Ntaxa >= ntaxamin) && (mParam.Ntaxa<=ntaxamax))	{
			avail[i] = 1;
			nav ++;
		}
	}
	
	if (nav < p)	{
		cerr << "error : cannot choose " << p << " among " << nav << "\n";
		exit(1);
	}
	for (int i=0; i<p; i++)	{
		int choose = (int) (Random::Uniform() * nav);
		int k = 0;
		while ((k<n) && (! avail[k])) k++;
		if (k == n)	{
			cerr << "error : overflow\n";
			exit(1);
		}
		int q = 0;
		while (q < choose)	{
			k++;
			while ((k<n) && (! avail[k])) k++;
			if (k == n)	{
				cerr << "error : overflow\n";
				exit(1);
			}
			q++;
		}
		avail[k] = 0;
		nav--;
		cout << list[k] << '\n';
		system(("cp " + list[k] + " " + dir).c_str());
	}
}

