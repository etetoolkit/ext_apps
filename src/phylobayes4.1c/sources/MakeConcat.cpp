#include "phylo.h"

const char* alphabet = AminoAcids;

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "makeconcat list output\n";
		exit(1);
	}
	string outname = argv[2];

	cerr << "reading list\n";
	ifstream is(argv[1]);
	int N;
	is >> N;
	cerr << N << " genes\n";
	cerr.flush();
	string name[N];
	for (int i=0; i<N; i++)	{
		is >> name[i];
		cerr << name[i] << '\n';
	}
	is.close();
	cerr << "ok\n\n";
	cerr.flush();

	int GeneSize[N];
	int GeneFirstSite[N];
	int Nsite = 0;	

	int Nmax = 1000;
	string SpeciesNames[Nmax];
	int Ntaxa = 0;

	for (int i=0; i<N; i++)	{
		MCParameters mParam;
		mParam.ReadDataFromFile(name[i]);
		GeneSize[i] = mParam.Nsite;
		GeneFirstSite[i] = Nsite;
		for (int j=0; j<mParam.Ntaxa; j++)	{
			int k = 0;
			int found = 0;
			while ((! found) && (k<Ntaxa))	{
				if (mParam.SpeciesNames[j] == SpeciesNames[k])	{
					found = 1;
				}
				else	{
					k++;
				}
			}
			if (! found)	{
				SpeciesNames[Ntaxa++] = mParam.SpeciesNames[j];
			}
		}
		Nsite += mParam.Nsite;
	}

	int** Data = new int*[Ntaxa];
	for (int j=0; j<Ntaxa; j++)	{
		Data[j] = new int[Nsite];
		for (int i=0; i<Nsite; i++)	{
			Data[j][i] = -1;
		}
	}

	int total = 0;
	for (int i=0; i<N; i++)	{
		MCParameters mParam;
		mParam.ReadDataFromFile(name[i]);
		
		for (int j=0; j<mParam.Ntaxa; j++)	{
			// find species
			int k=0;
			int found = 0;
			while ((! found) && (k<Ntaxa))	{
				if (mParam.SpeciesNames[j] == SpeciesNames[k])	{
					found = 1;
				}
				else	{
					k++;
				}
			}
			if (! found)	{
				cerr << "error : cannot identify species " <<mParam.SpeciesNames[j] << '\n';
				cerr << "gene : " << i << " of name " << name[i] << '\n';
				exit(1);
			}

			for (int l=0; l<mParam.Nsite; l++)	{
				Data[k][total + l] = mParam.Data[j][l];
			}
		}
		
		total += mParam.Nsite;
	}


	ofstream pos((outname + ".partition").c_str());
	pos << "Partition   " << N << '\n';
	for (int i=0; i<N; i++)	{
		pos <<  GeneSize[i] << '\t';
	}
	pos << '\n';
	pos.close();

	ofstream os((outname + ".phy").c_str());

	os << Ntaxa << '\t' << Nsite << '\n';

	unsigned int maxsize = 0;
	for (int i=0; i<Ntaxa; i++)	{
		if (maxsize < SpeciesNames[i].length())	{
			maxsize = SpeciesNames[i].length();
		}
	}

	maxsize += 5;

	for (int i=0; i<Ntaxa; i++)	{
		os << SpeciesNames[i];
		for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
			os << ' ';
		}

		for (int j=0; j<Nsite; j++)	{
			if (Data[i][j] == unknown)	{
				os << '-';
			}
			else	{
				os << alphabet[Data[i][j]];
			}
		}
		os << '\n';
	}
}

