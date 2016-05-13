#include "phylo.h"

int MakeConcat(int N, string* name, string outname, int ntaxa, string excludedtaxon, string includedtaxon);

const double cutoff = 0.70;

int main(int argc, char* argv[])	{

	if (argc == 1)	{
		cerr << "jack <genelist> <ngene> <ntaxa> <nrep> <basename>\n";
		cerr << "\t<genelist>: list of genes; preceded by their total number\n";
		cerr << "\t<ngene>   : number of genes to be included in each jacknife replicate\n";
		cerr << "\t<ntaxa>   : number of taxa to be included in each jacknife replicate\n";
		cerr << "\t<nrep>    : number of replicates\n";
		cerr << "\t<basename>: prefix for all output files\n";
		cerr << "\t<excluded>\n";
		exit(1);
	}
	ifstream is(argv[1]);
	int ngene = atoi(argv[2]);
	int ntaxa = atoi(argv[3]);
	int nrep = atoi(argv[4]);
	string basename = argv[5];
	string excludedtaxon = "";
	if (argc >= 7)	{
		excludedtaxon = argv[6];
	}
	string includedtaxon = "";
	if (argc == 8)	{
		includedtaxon = argv[7];
	}

	int N;
	is >> N;
	cerr.flush();
	string name[N];
	for (int i=0; i<N; i++)	{
		is >> name[i];
	}
	is.close();

	/*
	int Nex = 0;
	string* exname;
	if (argc == 7)	{
		ifstream isex(argv[1]);
		isex >> Nex;
		exname = new string[Nex];
		for (int i=0; i<Nex; i++)	{
			is >> exname[i];
			cerr << exname[i] << '\n';
		}
		is.close();
	}
	*/

	Random::Random();
	string* subname = new string[ngene];
	
	int rep = 0;
	int miss = 0;
	while (rep < nrep)	{
		int* index = new int[ngene];
		Random::DrawFromUrn(index,ngene,N);
		for (int i=0; i<ngene; i++)	{
			subname[i] = name[index[i]];
		}
		ostringstream s;
		s << basename << ngene << '_' << ntaxa << '_' << rep;
		if (MakeConcat(ngene,subname,s.str(),ntaxa,excludedtaxon, includedtaxon))	{
			rep++;
		}
		else	{
			miss++;
		}
		delete[] index;
	}
	cerr << "missed: " << miss << '\n';
	delete[] subname;
}


int MakeConcat(int N, string* name, string outname, int ntaxa, string excludedtaxon, string includedtaxon) {

	int GeneSize[N];
	int GeneFirstSite[N];
	int Nsite = 0;	

	int Nmax = 1000;
	string SpeciesNames[Nmax];
	int Ntaxa = 0;

	int ret = 1;
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

	if (ntaxa > Ntaxa)	{
		ret = 0;
	}

	ofstream pos((outname + ".ali.partition").c_str());
	pos << N << '\n';
	for (int i=0; i<N; i++)	{
		pos <<  GeneSize[i] << '\t';
	}
	pos << '\n';
	for (int i=0; i<N; i++)	{
		pos <<  name[i] << '\n';
	}
	pos.close();

	ofstream os((outname + ".ali").c_str());

	int* index = new int[ntaxa];
	int* mask = new int[Ntaxa];
	Random::DrawFromUrn(index,ntaxa,Ntaxa);
	for (int j=0; j<Ntaxa; j++)	{
		mask[j] = 0;
	}

	for (int i=0; i<ntaxa; i++)	{
		mask[index[i]] = 1;
		if (SpeciesNames[index[i]] == excludedtaxon)	{
			cerr << "missed : " << excludedtaxon << "\n";
			ret = 0;
		}
	}
	if (includedtaxon != "")	{
		int found = 0;
		for (int i=0; i<ntaxa; i++)	{
			if (SpeciesNames[index[i]] == includedtaxon)	{
				found = 1;
			}
		}
		if (! found) ret = 0;
	}

	delete[] index;

	int max = 0;
	int imax = 0;
	for (int j=0; j<Ntaxa; j++)	{
		int missing = 0;
		for (int i=0; i<Nsite; i++)	{
			if (Data[j][i] == unknown)	{
				missing++;
			}
		}
		if (max < missing)	{
			max = missing;
			imax = j;
		}
	}
	if (((double) max) / Nsite > cutoff)	{
		cerr << "too many missing : " << SpeciesNames[imax] << '\t' << ((double) max) / Nsite << "\n";
		ret = 0;
	}
	cerr << max << '\t' << Nsite << '\t' << ((double) max) / Nsite << '\t' << SpeciesNames[imax] << '\t' << ret << '\n';
	if (max == Nsite)	{
		cerr << "one taxon with purely missing; skip\n";
		ret = 0;
	}
	
	os << ntaxa << '\t' << Nsite << '\n';

	unsigned int maxsize = 0;
	
	for (int i=0; i<Ntaxa; i++)	{
		if (mask[i])	{
			if (maxsize < SpeciesNames[i].length())	{
				maxsize = SpeciesNames[i].length();
			}
		}
	}

	maxsize += 5;

	for (int i=0; i<Ntaxa; i++)	{
		if (mask[i])	{
			os << SpeciesNames[i];
			for (unsigned int j=0; j<maxsize-SpeciesNames[i].length(); j++)	{
				os << ' ';
			}

			for (int j=0; j<Nsite; j++)	{
				if (Data[i][j] == unknown)	{
					os << '-';
				}
				else	{
					os << AminoAcids[Data[i][j]];
				}
			}
			os << '\n';
		}
	}
	delete[] mask;
	for (int i=0; i<Ntaxa; i++)	{
		delete[] Data[i];
	}
	delete[] Data;
	return ret;
}

