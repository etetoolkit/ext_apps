#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	int n = 0;
	int em = 1;
	int pp = 0;
	int ss = 0;
	int countvector = 0;
	int mccountvector = 0;
	int nrep = 10000;
	double cutoff = 0.01;
	int fixedlengths = 0;
	int fixedalpha = 0;
	int fixedweights = 0;
	double initalpha = 1;
	string initalphafile = "";
	int randominitalpha = 1;
	int randominitlength = 1;

	int every = 1;
	int until = -1;
	int deleteconstant = 0;

	int autorestart = 0;
	string datafile = "";
	string seplist = "";
	string initfile = "";
	string treefile = "";
	string name = "";
	string path = "";
	string directory = "";
	string partition = "";
	string calibration = "";
	string statfix = "";
	string outputfile = "";
	int saveall = 0;

	RRMode rr = poisson;
	RASMode ras = gam;
	int discrate = 4;
	int empfreq = 0; // 1 : global freq, 2: site specific freq
	int ncat = -2; // -2 : catfix 0 : cat, -1 max, otherwise, fixed number of modes
	int statcenter = 1;
	int fixtopo = 0;
	int force = 0;

	int forcecat = 0;
	int qmode = 0;

	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if (s == "-d")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -d <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				datafile = argv[i];
			}
			else if (s == "-em")	{
				em = 1;
			}
			else if (s == "-n")	{
				i++;
				n = atoi(argv[i]);
			}
			else if (s == "-pp")	{
				pp = 1;
			}
			else if (s == "-ss")	{
				ss = 1;
			}
			else if (s == "-cv")	{
				countvector = 1;
			}
			else if (s == "-mccv")	{
				mccountvector = 1;
			}
			else if (s == "-fl")	{
				fixedlengths = 1;
				randominitlength = 0;
			}
			else if (s == "-l")	{
				randominitlength = 0;
			}
			else if (s == "-fa")	{
				cerr << "fixed alpha\n";
				fixedalpha = 1;
				randominitalpha = 0;
				i++;
				initalphafile = argv[i];
			}
			else if (s == "-a")	{
				randominitalpha = 0;
				i++;
				initalphafile = argv[i];
			}
			else if (s == "-fw")	{
				fixedweights = 1;
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-o")	{
				i++;
				outputfile = argv[i];
			}
			else if (s == "-c")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -c <float>\n";
					cerr << '\n';
					exit(1);
				}
				cutoff = atof(argv[i]);
			}
			else if (s == "-dc")	{
				deleteconstant = 1;
			}
			else if (s == "-T")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -T <treefile>\n";
					cerr << '\n';
					exit(1);
				}
				treefile = argv[i];
				fixtopo = 1;
			}
			else if (s == "-s")	{
				saveall = 1;
			}
			else if (s == "-i")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -i <initfile>\n";
					cerr << '\n';
					exit(1);
				}
				initfile = argv[i];
			}

			else if (s == "-poisson")	{
				rr = poisson;
			}
			else if (s == "-jtt")	{
				rr = jtt;
			}
			else if (s == "-wag")	{
				rr = wag;
			}
			else if (s == "-mtrev")	{
				rr = mtrev;
			}
			else if (s == "-gtr")	{
				rr = gtr;
			}
			else if (s == "-empfreq")	{
				empfreq = 1;
			}
			else if (s == "-siteempfreq")	{
				empfreq = 2;
			}
			else if (s == "-uni")	{
				ras = uni;
				discrate = 0;
			}
			else if (s == "-gamma")	{
				ras = gam;
			}
			else if (s == "-invgamma")	{
				ras = invgam;
			}
			else if (s == "-cgam")	{
				discrate = 0;
			}
			else if (s == "-dgam")	{
				discrate =  4;
				i++;
				if (i == argc) 	{
					cerr << "error in command: -dgam <# disc cat>\n";
					cerr << '\n';
					exit(1);
				}
				string temp = argv[i];
				if (IsInt(temp))	{
					discrate = Int(temp);
				}
				else	{
					cerr << "error in command: -dgam <# disc cat>\n";
					cerr << '\n';
					exit(1);
				}
			}
			else if (s == "-ncat")	{
				i++;
				forcecat = 1;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -ncat <integer>\n";
					cerr << '\n';
					exit(1);
				}
				ncat = atoi(argv[i]);
			}
			else if (s == "-catfix")	{
				ncat = -2;
				i++;
				statfix = argv[i];
				forcecat = 1;
			}
			else if (s == "-statfree")	{
				statcenter = 1;
			}
			else if (s == "-statflat")	{
				statcenter = 0;
			}
			else if (s == "-part")	{
				i++;
				partition = argv[i];
			}
			else if (s == "-R")	{
				autorestart = 1;
			}
			else if (s == "-f")	{
				force = 1;
				autorestart = 0;
			}
			else if (s == "-p")	{
				i++;
				if (i == argc) {
					cerr << "error in command: -p <path>\n";
					cerr << '\n';
					exit(1);
				}
				path = argv[i];
			}
			else if ((s == "-qsub") || (s == "-qprep"))	{
				qmode = 2;
			} 
			else if (s == "-bg")	{
				qmode = 1;
			}

			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (! IsInt(s))	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				every = atoi(argv[i]);
				i++;
				if (i == argc)	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
				s = argv[i];
				if (IsInt(s))	{
					until = atoi(argv[i]);
				}
				else	{
					cerr << "error in command: -x <every> <until>\n";
					cerr << '\n';
					exit(1);
				}
			}
			else	{
				if (i != (argc -1))	{
					cerr << "error in command: unrecognized " << argv[i] << '\n';
					cerr << "assume chain name is " << argv[argc -1] << '\n';
					cerr << '\n';
					exit(1);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			cerr << "error in command: should specifiy a name for the chain\n";
			cerr << '\n';
			exit(1);
		}
	}
	catch(...)	{
		cerr << "\n";
		cerr << "em [options] <chainname>\n";
		cerr << "\tcreates a new chain, sampling from the posterior distribution, conditional on specified data\n";
		cerr << "\n";
		cerr << "\t-d <filename>       : file containing an alignment in phylip format; only amino acids\n";
		cerr << "\t-dc                 : deletes constant columns of the alignment\n";
		cerr << "\t-T <treefile>       : chain run under fixed, specified tree\n"; 
		cerr << "\t-catfix <profiles>  : specifying a fixed set of pre-defined categories\n";
		cerr << "\t                      file should be formatted as follows:\n";
		cerr << "\t                      <ncat>\n";
		cerr << "\t                      <weight1> <freq1_1> ... <freq1_20> \n";
		cerr << "\t                      <weight2> <freq2_1> ... <freq2_20> \n";
		cerr << "\t                      ...\n";
		cerr << "\t-uni                : uniform rates across sites\n";
		cerr << "\t-dgam <ncat>        : discrete gamma. ncat = number of categories (4 by default)\n";
		cerr << "\t-c <cutoff>         : stops EM when improvement in log likelihood is < cutoff\n"; 
		cerr << "\t-n <ncycle>         : stops EM after <ncycles> cycles\n";
		cerr << "\t-w                  : optimises the weights\n";
		cerr << "\t-cv                 : computes the count vectors\n";
		cerr << "\t-mccv               : computes the count vectors by monte carlo\n";
		cerr << "\t-nrep <nrep>        : specifies number of replicates for the monte carlo\n";
		cerr << "\t-p                  : makes a posterior predictive diversity test\n";
		cerr << "\t-f                  : force (overwrite already existing output file)\n";
		cerr << "\t-o <output>         : in case -mc is activated, appends count vectors to specified output file\n";
		cerr << "\n";
		cerr << "\t-fw                 : fix weights to those given in the input file\n";
		cerr << "\t-fa <file>          : fix alpha to the value specified in <file>\n";
		cerr << "\t-a <file>           : start with alpha equal to the value specified in <file>\n";
		cerr << "\t-fl                 : fix lengths to those specified in the input tree\n";
		cerr << "\t-l                  : start from lengths specified in the input tree\n";
		exit(1);
	}

	cerr << '\n';

	if (qmode > 1)	{
		ofstream os((name + ".launch").c_str());
		// os << "# init\n";
		// os << "#$ -v PATH\n";
		int i = 0;
		while (i<argc)	{
			string tmp = argv[i];
			if ((tmp == "-qsub") ||	(tmp == "-qprep"))	{
			}
			else if (tmp == "-path")	{
				i++;
			}
			else	{ 
				os << argv[i] << ' ';
			}
			i++;
		}
		os << '\n';
		os.close();
		string s = "qsub -cwd -o " + name + ".out -e " + name + ".err " + name + ".launch";
		if (qmode == 2)	{
			system(s.c_str());
		}
		exit(1);
	}

	if ((autorestart) || ((initfile == "") && (datafile == "")))	{		// assumes this is an already existing chain
		ifstream Param_is((path + name + ".param").c_str());
		if (! Param_is)	{
			if (! autorestart)	{
				cerr << "?? non existing chain : " <<  path + name << '\n';
				cerr << "to create a new chain, a datafile must be specified\n";
				cerr << '\n';
				exit(1);
			}
		}
		else	{
			Chain chain(name, 0, path);
			MCParameters* mParam = chain.mParam;
			cout << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\t';
			cerr << "how many saved : " << chain.mParam->HowManySaved << '\n';
			cout << "\nchain started\n\n";
			chain.Start();
			exit(1);
		}
	}

	if ((!force) && (ifstream((path + name + ".log").c_str())))	{
		cerr << "Chain " << name << " seems to already exist\n";
		cerr << "use \"-f\" option to override\n";
		exit(1);
	}
	MCParameters* mParam = new MCParameters();
	mParam->ActivateEM = 1;

	mParam->SaveEvery = every;
	mParam->StopAfter = until;
	mParam->Path = path;
	mParam->SaveAll = saveall;
	ofstream os((name + ".log").c_str());

	if (deleteconstant)	{
		mParam->DeleteConstant = Yes;
		os << "constant sites deleted\n";
	}
	if (datafile != "")	{
		mParam->ReadDataFromFile(datafile);
		os << "data file : " << datafile << '\n';
		os << "number of taxa : " << mParam->Ntaxa << '\n';
		os << "number of sites: " << mParam->Nsite << '\n';
	}

	if (fixtopo)	{
		if (treefile == "")	{
			cerr << "error: no topology file was specified\n";
			exit(1);
		}
		os << "fixed tree topology: " << treefile << '\n';
	}
	os << '\n';

	if (fixedweights)	{
		mParam->FixedWeights = 1;
	}
	if (fixedlengths)	{
		mParam->FixedLengths = 1;
	}
	if (fixedalpha)	{
		mParam->FixedAlpha = 1;
	}
	if (randominitlength)	{
		mParam->RandomInitLength = 1;
	}
	else	{
		mParam->RandomInitLength = 0;
	}

	
	if ((rr != poisson) && (! forcecat))	{
		ncat = 1;
	}
	if (ncat == -2)	{
		os << "fixed profiles taken from " << statfix << '\n';
		if (statfix != "")	{
			mParam->ReadStatFix(statfix);
		}
		else	{
			cerr << "error : should specify profile mixture\n";
			exit(1);
		}
	}
	else if (ncat == 0)	{
		os << "cat model: mixture of a free number of profiles\n";
		if (statcenter)	{
			os << "flexible hyperprior on profiles\n";
		}
		else	{
			os << "flat prior on profiles\n";
		}
	}
	else if (ncat == -1)	{
		os << "max model: site-specific profiles\n";
		if (statcenter)	{
			os << "flexible hyperprior on profiles\n";
		}
		else	{
			os << "flat prior on profiles\n";
		}
	}
	else if (ncat == 1)	{
		statcenter = 0;
		os << "site homogeneous model: ";
		if (rr == poisson)	{
			os << "poisson\n";
		}
		else if (rr == gtr)	{
			os << "gtr\n";
		}
		else if (rr == wag)	{
			os << "wag\n";
		}
		else if (rr == jtt)	{
			os << "jtt\n";
		}
		else if (rr == mtrev)	{
			os << "mtrev\n";
		}
		if (empfreq == 0)	{
			os << "free equilibrium frequencies (inferred from data)\n";
		}
		else if (empfreq == 1)	{
			os << "empirical equilibrium frequencies\n";
		}
		else	{
			os << "site-specific equilibrium frequencies\n";
		}
	}
	else	{
		os << "mixture of a fixed (" << ncat << ") number of profiles\n";
		if (statcenter)	{
			os << "flexible hyperprior on profiles\n";
		}
		else	{
			os << "flat prior on profiles\n";
		}
	}
		
	if (ras == uni)	{
		os << "uniform rates across sites\n";
	}
	else if (ras == gam)	{
		os << "gamma distributed rates: ";
		if (discrate)	{
			os << discrate << " discrete categories\n";
		}
		else	{
			os << "continuous distribution\n";
		}
	}
	os << '\n';
	if (initfile != "")	{
		mParam->InitFromFile(initfile);
	}
	else	{
		mParam->Init(rr,empfreq, ncat, statcenter, ras, discrate);
	}

	if (! mParam->DataOK)	{
		cerr << "error: should specify a datafile\n";
		exit(1);
	}

	if (treefile != "")	{
		ifstream is((path + treefile).c_str());
		mParam->GetCurrentState()->SetTree(is);
		if (fixtopo)	{
			mParam->FixTopo = 1;
		}
	}

	if (initfile == "")	{
		mParam->InitMove();
	}

	if (randominitalpha)	{
		mParam->GetCurrentState()->gamma = Random::sExpo();
	}
	else	{
		ifstream is(initalphafile.c_str());
		if (!is)	{
			cerr << "cannot find " << initalphafile << '\n';
			exit(1);
		}
		is >> initalpha;
		mParam->GetCurrentState()->gamma = initalpha;
	}

	mParam->Reset();

	int Nsite = mParam->Nsite;
	int Ntaxa = mParam->Ntaxa;
	cout << '\n';
	cout << "Nsite : " << Nsite << '\n';
	cout << "Ntaxa : " << Ntaxa << '\n';
	cout << "\nem started\n\n";
	cout.flush();

		
	PhyloBayes* pb = mParam->GetCurrentState();
	os << "initial tree\n";
	pb->Phylip(os,1,0,1,0);
	os << "initial alpha : " << pb->gamma << "\n";

	// HERE IS THE EM 
	double lnL = pb->EM(cutoff, n);

	// here you can change the topology
	// and re-call EM
 
	// pb->NNI_EM()

	ofstream emos((name + ".em").c_str());
	emos << "lnL: " << lnL << '\n';
	emos << '\n';
	pb->Phylip(emos,1,0,1,0);
	emos << '\n';
	emos << "alpha: " << pb->gamma << '\n';
	emos << '\n';
	emos << "weights:\n";
	for (int k=0; k<pb->Nmode; k++)	{
		emos << pb->ModeWeight[k] << '\n';
	}
	emos << '\n';
	emos << "log prob inf count : " << mParam->LogProbInfCount << "\n";
	emos.close();
	string echo = "cat " + name + ".em";
	system(echo.c_str());

	if (countvector)	{
		pb->ComputeCountVector();
		/*
		if (outputfile == "")	{
			ofstream count_os((name + ".counts").c_str());
			for (int i=0; i<mParam->Nsite; i++)	{
				count_os << 0 << '\t' << 0 << '\t' << Decimal(pb->TotalCountVector[i],5);
				for (int k=0; k<mParam->Nstate; k++)	{
					count_os << '\t' << Decimal(pb->CountVector[i][k],5);
				}
				count_os << '\n';
			}
		}
		else	{
			ofstream count_os(outputfile.c_str(),APPEND);
			for (int i=0; i<mParam->Nsite; i++)	{
				count_os << 0 << '\t' << 0 << '\t' << Decimal(pb->TotalCountVector[i],5);
				for (int k=0; k<mParam->Nstate; k++)	{
					count_os << '\t' << Decimal(pb->CountVector[i][k],5);
				}
				count_os << '\n';
			}
		}
		*/
		ofstream count_os((name + ".counts").c_str());
		pb->OutputEM(count_os);
	}	
	if (mccountvector)	{
		ofstream mcmccount_os((name + ".mcmccounts").c_str());
		pb->MCMCComputeCountVector(nrep);
		for (int i=0; i<mParam->Nsite; i++)	{
			mcmccount_os << 0 << '\t' << 0 << '\t' << pb->MCMCTotalCountVector[i];
			for (int k=0; k<mParam->Nstate; k++)	{
				mcmccount_os << '\t' << pb->MCMCCountVector[i][k];
			}
			mcmccount_os << '\n';
		}
	}

	if (pp)	{
		pb->Diversity(name,nrep);
	}
	if (ss)	{
		ofstream site_os((name + ".sslogL").c_str());
		pb->OutputEMSiteLogL(site_os);
	}
}


