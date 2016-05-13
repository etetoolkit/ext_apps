
#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	string SampleName;
	string ChainName;
	string OutGroupFile;

	string Path = "";

	int rates = 0;
	int modes = 0;
	int sitestat = 0;

	int extract = 0;
	int burnin = 0;
	int every = 1;
	int until = -1;
	int ncat = 10;

	int clock = 0;
	int ps = 0;

	int clustermodes = 0;

 	double mindist = 0.03;
 	int minsize = 10;

	double cutoff = 0.5;

	// read arguments

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-c")	{
				i++;
				if (i == argc) throw(0);
				if (! IsFloat(argv[i])) throw(0);
				cutoff = atof(argv[i]);
			}
			else if (s == "-div")	{
				clock = 1;
			}
			else if (s == "-p")	{
				i++;
				Path = argv[i];
			}
			else if (s == "-ps")	{
				ps = 1;
			}
			else if (s == "-cl")	{
				clustermodes = 1;
			}
			else if (s == "-sz")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				minsize = atoi(argv[i]);
			}
			else if (s == "-ds")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsFloat(s))	{
					throw(0);
				}
				mindist = atof(argv[i]);
			}
			else if (s == "-ncat")	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				ncat = atoi(argv[i]);
			}
			else if ( (s == "-m") || (s == "-modes") )	{
				modes = 1;
			}
			else if ( (s == "-r") || (s == "-rates") )	{
				rates = 1;
			}
			else if ( (s == "-ss") || (s == "-sitestat") )	{
				sitestat = 1;
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
		cerr << "readsample  [-c cutoff] [-m -ss -cl] [-ms <minsize>] [-md <mindist>] [-r -ncat <n>] [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		cerr << "\tdefaults : burnin = 0, every = 1, until the end\n";
		cerr << "\t-c cutoff : collapses all groups with posterior probability lower than cutoff\n"; 
		cerr << "\t-m : posterior distribution of the number of modes\n";
		cerr << "\t-ss: mean posterior site-specific stationaries\n";
		cerr << "\t-r : mean posterior site-specific rates (continuous gamma only)\n";
		cerr << '\n';
		cerr << "\t-ncat <n> : defines number of bins for rate histogram (default 20)\n";
		cerr << '\n';
		cerr << "\t-cl: mode clustering\n";
		cerr << "\t\t-ms: cluster min size (default : 1)\n";
		cerr << "\t\t-md: aggregating distance threshold (default : 0.01)\n";
		cerr << '\n';
		cerr << "\t-ps   : postscript output for tree (latex must be installed)\n";
		cerr << '\n';
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}

	try	{
		ifstream true_is((ChainName + ".chain").c_str());
		MCParameters* mParam = new MCParameters();
		ifstream Param_is((ChainName+".param").c_str());
		if (! Param_is)	{
			cerr << "?? non existing chain : " << ChainName << '\n';
			exit(1);
		}
		Param_is >> *mParam;
		mParam->Update();
		PhyloBayes* truepb = mParam->GetCurrentState();
		if (nrep == -1)	{
			nrep = mParam->HowManySaved;
		}
		cerr << nrep << '\n';

		int histo[ncat];
		for (int rep=0; rep<nrep; rep++)	{
			true_is >> *truepb;
			Tree tree(truepb);
			BipartitionList truebplist(&tree);
			if (verbose)	{
				cerr << "rep : " << rep << "\n";
			}
			ostringstream s;
			s << ChainName << rep << ".phy.treelist";
			string name = s.str();
			if (! ifstream(name.c_str()))	{
				name = ChainName;
				if (! ifstream(name.c_str()))	{
					cerr << "error: non existing chain\n";
					exit(1);
				}
			}
			BipartitionList bplist(name,burnin,every,until);
			if (!bplist.Ntree)	{
				cerr << "empty tree list\n";
				exit(1);
			}
			cout << bplist.Ntree << " trees were read\n";
			cout.flush();
			
			for (int i=0; i<truebplist.GetSize(); i++)	{
				Bipartition& truebp = truebpllist.GetBipartition(i);
				int found=0;
				int k= 0;
				while ((k<bplist.GetSize()) && (truebp != )) k++;
			}
			
		}

		Sample* sample = new Sample(ChainName,burnin,every,until,Path);

		MCParameters* mParam = sample->GetParameters();
		// initialising

		int Nsite = mParam->Nsite;
		int Ntaxa = mParam->Ntaxa;

		cout << '\n';
		cout << "Nsite : " << Nsite << '\n';
		cout << "Ntaxa : " << Ntaxa << '\n';
		cout << sample->GetSize() << " points to read\n";
		cout << '\n';
		cout.flush();


		if (mParam->NormalApprox)	{
			if (clustermodes)	{
				cerr << "no mode clustering under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (modes)	{
				cerr << "no mode analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (rates)	{
				cerr << "no rate analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
			if (sitestat)	{
				cerr << "no site-specific profile analysis under normal approx\n";
				cerr << '\n';
				exit(1);
			}
		}

		if (clock)	{
			sample->Dating(ps);
			sample->Reset();
		}
			
		else if (clustermodes)	{
			if (! mParam->SaveAll)	{
				cerr << "error : biochemical profiles were not saved. cannot cluster them\n";
				exit(1);
			}
			sample->ClusterModes(minsize, mindist, 1);
		}
		else {
			sample->Read(rates, modes, sitestat, cutoff, ncat, ps);
		}
	}
	catch(...)	{
		exit(1);
	}
}

