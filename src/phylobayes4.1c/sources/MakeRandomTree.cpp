#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	int every = 1;
	int until = -1;
	int deleteconstant = 0;

	int clock = 0;
	int normalapprox = 0;
	int separate = 0; // 0 : concat 1 : separate 2 : separate with multipliers		
	int autorestart = 0;

	int covarion = 0;
	int external = 0;
	int conjugate = 1;
	int tempering = 0;
	double initbeta = 0;
	double betastep = 0;
	int burnin = 0;

	int nh = 0;
	int nnh = 10;
	int nhprior = 0;

	// clock models : KT, CIR, Strict WhiteNoise
	// Normal Approx : rigid for KT Strict and CIR. Flexible for WhiteNoise
	// Pruning: no rigid model. Only flexible, CIR or White noise
	// if genemode = 2: underlying global model is CIR

	string datafile = "";
	string seplist = "";
	string initfile = "";
	string treefile = "";
	string name = "";
	string path = "";
	string directory = "";
	string partition = "";
	string outgroup = "";
	string calibration = "";
	string statfix = "";
	
	int nrep = 10;

	int timeprior = 0; // 0: uniform, 1:birth death, 2:dirichlet

	int saveall = 1;
	// int saveall = 0;
	RecodingMode rec = NoRecoding;
	string recfile = "";

	RRMode rr = poisson;
	RASMode ras = gam;
	int discrate = 4;
	int empfreq = 0; // 1 : global freq, 2: site specific freq
	int ncat = 0; // -2 : catfix 0 : cat, -1 max, otherwise, fixed number of modes
	int statcenter = 1;
	int fixtopo = 0;
	int withpartial = 1;
	int force = 0;

	int forcecat = 0;
	int zipgtr = 0;
	int zipprior = 0;

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
			else if (s == "-fast")	{
				withpartial = 1;
			}
			else if (s == "-lm")	{
				withpartial = 0;
			}
			else if (s == "-cl")	{
				clock = 1;
			}
			else if (s == "-ln")	{
				clock = 2;
			}
			else if (s == "-cir")	{
				clock = 3;
			}
			else if (s == "-wn")	{
				clock = 4;
			}
			else if (s == "-ugam")	{
				clock = 5;
			}
			else if (s == "-bd")	{
				timeprior = 1;
			}
			else if (s == "-dir")	{
				timeprior = 2;
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			/*
			else if (s == "-ncat")	{
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -ncat <integer>\n";
					cerr << '\n';
					exit(1);
				}
				meanscale = atof(argv[i]);
			}
			*/
			else if (s == "-covarion")	{
				covarion = 1;
				external = 0;	
			}
			else if (s == "-covext")	{
				covarion = 1;
				external = 1;	
			}
			else if (s == "-zipgtr")	{
				i++;
				zipgtr = atoi(argv[i]);
			}
			else if (s == "-zipprior")	{
				i++;
				zipprior = atoi(argv[i]);
			}
			else if (s == "-mh")	{
				conjugate = 0;
			}
			else if (s == "-tmp")	{
				conjugate = 0;
				tempering = 1;
				i++;
				if (IsFloat(argv[i]))	{
					initbeta = atof(argv[i]);
					i++;
					if (IsFloat(argv[i]))	{
						betastep = atof(argv[i]);
						i++;
						if (IsInt(argv[i]))	{
							burnin = atoi(argv[i]);
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
					i--;
				}
			}
			else if (s == "-nh")	{
				nh = 1;
				i++;
				if (IsInt(argv[i]))	{
					nnh = atoi(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-bnh")	{
				nh = 2;
			}
			else if (s == "-fnh")	{
				nh = 3;
			}
			else if (s == "-anh")	{
				nh = 4;
			}
			else if (s == "-cnh")	{
				nh = 5;
			}
			
			else if (s == "-dnh")	{
				nh = 6;
			}
			else if (s == "-nhvar")	{
				nhprior = 1;
			}
			else if (s == "-nhcov")	{
				nhprior = 2;
			}
			else if (s == "-cov")	{
				i++;
				datafile = argv[i];
				normalapprox = 1;
				if (! clock)	{
					clock = 2;
				}
			}
			else if (s == "-sepcov")	{
				normalapprox = 1;
				if (!clock)	{
					clock = 2;
				}
				i++;
				seplist = argv[i];
			}
			else if (s == "-sep")	{
				separate = 1;
			}
			else if (s == "-mulsep")	{
				separate = 2;
			}
			else if (s == "-conc")	{
				separate = 0;
			}
			else if (s == "-cal")	{
				i++;
				calibration = argv[i];
				if (!clock)	{
					clock = 2;
				}
			}	
			else if (s == "-dc")	{
				deleteconstant = 1;
			}
			else if ((s == "-og") || (s == "-r"))	{
				i++;
				outgroup = argv[i];
			}
			else if ((s == "-rec") || (s == "-recode"))	{
				i++;
				s = argv[i];
				if (s == "hp")	{
					rec = HP;
				}
				else if (s == "dayhoff6")	{
					rec = Dayhoff6;
				}
				else if (s == "dayhoff4")	{
					rec = Dayhoff4;
				}
				else {
					rec = Custom;
					recfile = argv[i];
					cerr << "loading recoding table from " << recfile << '\n';	
				}
			}
			else if (s == "-t")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -t <treefile>\n";
					cerr << '\n';
					exit(1);
				}
				treefile = argv[i];
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
			else if (s == "-ratecat")	{
				ras = dprate;
			}
			else if (s == "-catfix")	{
				ncat = -2;
				i++;
				statfix = argv[i];
				forcecat = 1;
			}
			else if (s == "-cat")	{
				ncat = 0;
				forcecat = 1;
			}
			else if (s == "-max")	{
				ncat = -1;
				forcecat = 1;
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
		cerr << "randtree -d <data> -nrep <nrep> <output>\n";
		cerr << "makes <nrep> random instances of phylogenetic trees\n";
		cerr << "with taxa taken from dataset <data>\n";
		cerr << "list of tree output into <output>\n";
		cerr << "default nrep = 10\n";
		exit(1);
	}

	cerr << '\n';

	// make a MCParameters object, and prepare it according to the options
	MCParameters* mParam = new MCParameters();
	mParam->SaveEvery = every;
	mParam->StopAfter = until;
	mParam->Path = path;
	if (clock || fixtopo)	{
		saveall = 1;
	}
	mParam->SaveAll = saveall;
	mParam->Conjugate = conjugate;
	mParam->ZipGTR = zipgtr;
	mParam->ZipPrior = zipprior;
	if (zipprior == 2)	{
		mParam->ZipGTRDP = 0;
	}

	if (tempering)	{
		mParam->Tempering = 1;
		mParam->MSMode = Thermo;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->BetaStep = betastep;
	}
	mParam->ReadDataFromFile(datafile);
	cerr << "data file : " << datafile << '\n';
	cerr << "number of taxa : " << mParam->Ntaxa << '\n';
	cerr << "number of sites: " << mParam->Nsite << '\n';
	
	if ((rr != poisson) && (! forcecat))	{
		ncat = 1;
	}
	if (ncat == -2)	{
	}
	else if (ncat == 0)	{
		if (rr != poisson)	{
			statcenter = 0;
			mParam->AlphaPrior = Exponential;
		}
	}
	else if (ncat == 1)	{
		statcenter = 0;
	}

	mParam->Init(rr,empfreq, ncat, statcenter, ras, discrate);

	if (! mParam->DataOK)	{
		cerr << "error: should specify a datafile\n";
		exit(1);
	}

	mParam->InitMove();

	// create a chain ...
	Chain chain(mParam,name);

	for (int rep=0; rep<nrep; rep++)	{
		ostringstream s;
		s << name << rep << ".tre";
		ofstream os(s.str().c_str());
		mParam->GetCurrentState()->MakeRandomTree();
		mParam->GetCurrentState()->DrawLengths();
		mParam->GetCurrentState()->Phylip(os,1,0,1,0);
		os.close();
	}
}


