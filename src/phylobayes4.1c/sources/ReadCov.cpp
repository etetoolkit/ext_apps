
#include "phylo.h"

double** Stationary;

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	string SampleName;
	string ChainName;

	int oes = 0;

	int extract = 0;
	int burnin = -1;
	int every = 1;
	int until = -1;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-oes")	{
				oes = 1;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				extract = 1;
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				ChainName = argv[i];
			}
			i++;
		}
		if ((SampleName == "") && (ChainName == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readcov  [-oes] [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		cerr << "\tcomputes the variance covariance matrix\n";
		cerr << "\t-oes: outgroup is excluded (as is done by estbranch)\n";
		cerr << "\n";
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}

	try	{

		Sample* sample = new Sample(ChainName,burnin,every,until);

		cout << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();

		if (oes)	{
			sample->CovWoOutgroup();
		}
		else	{
			sample->Cov();
		}

		// delete sample;
	}
	catch(...)	{
		exit(1);
	}
}

