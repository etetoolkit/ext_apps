#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	GeneticCodeType gencode = Universal;

	int rrprior = 0;

	int prior = 0;
	string rrfile = "";

	int x_option = 0;
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
	int nchain = 1;
	int burnin = 1;
	double cutoff = 0.1;

	int tempering = 0;
	double initbeta = 0;
	double betastep = 0;

	int popeff = 0;

	int nh = 0;
	int nnh = 10;
	int nhprior = 0;

	int genebl = 0;
	string genepartition = "";
	int mbl = 0;
	int mutmode = 0;
	int selmode = 0;
	double gwf = 1;

	double meanscale = 0;

	double meanlength = 0.1;
	int fixmeanlength = 0;

	// clock models : KT, CIR, Strict WhiteNoise
	// Normal Approx : rigid for KT Strict and CIR. Flexible for WhiteNoise
	// Pruning: no rigid model. Only flexible, CIR or White noise
	// if genemode = 2: underlying global model is CIR

	string datafile = "";
	string contdatafile = "";
	string seplist = "";
	string initfile = "";
	string treefile = "";
	string name = "";
	string path = "";
	string directory = "";
	string outgroup = "";
	string calibration = "";
	string statfix = "";
	
	int timeprior = 0; // 0: uniform, 1:birth death, 2:dirichlet
	double chi = -1;
	double chi2 = -1;

	int nstart = 0;

	int saveall = 0;
	// int saveall = 0;
	RecodingMode rec = NoRecoding;
	string recfile = "";

	RRMode rr = poisson;
	int qmode = 0;
	RASMode ras = gam;
	int discrate = 4;
	int empfreq = 0; // 1 : global freq, 2: site specific freq
	int ncat = 0; // -2 : catfix 0 : cat, -1 max, otherwise, fixed number of modes
	int statcenter = 1;
	int fixtopo = 0;
	int withpartial = 1;
	int force = 0;
	int freeweight = 0;

	int forcecat = 0;
	int zipgtr = 0;
	int zipprior = 0;

	int qmod = 0;
	int lengthprior = 1; // gamma
	// 0: exponential	

	int nrep = 10;

	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if ((s == "-h") || (s == "--h") || (s == "help") || (s == "-help")) throw(0);
			else if (s == "-meanlength")	{
				fixmeanlength = 1;
				i++;
				meanlength = atof(argv[i]);
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-d")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -d <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				datafile = argv[i];
			}
			else if (s == "-cd")	{
				i++;
				if (i == argc) 	{
					cerr << "error in command: -cd <datafile>\n";
					cerr << '\n';
					exit(1);
				}
				contdatafile = argv[i];
			}
			else if (s == "-lm")	{
				withpartial = 0;
			}
			else if (s == "-nmode")	{
				nstart = 2;
			}
			else if (s == "-1mode")	{
				nstart = 1;
			}
			else if (s == "-universalcode")	{
				gencode = Universal;
			}
			else if (s == "-mtmamcode")	{
				gencode = MtMam;
			}
			else if (s == "-rrexp")	{
				rrprior = 0;
			}
			else if (s == "-rrgam")	{
				rrprior = 1;
			}
			else if (s == "-rp")	{
				i++;
				if ((i == argc) || (! IsFloat(argv[i])))	{
					cerr << "error in command: -ncat <integer>\n";
					cerr << '\n';
					exit(1);
				}
				meanscale = atof(argv[i]);
			}
			else if (s == "-prior")	{
				prior = 1;
				nstart = 1;
			}
			else if (s == "-nchain")	{
				i++;
				if ((i == argc) || (! IsInt(argv[i])))	{
					cerr << "error in command: -nchain <integer>\n";
					cerr << '\n';
					exit(1);
				}
				nchain = atoi(argv[i]);
				if (nchain<1)	{
					cerr << "error : nchain should be at least 1\n";
					exit(1);
				}
				i++;
				if (i==argc)	{
					cerr << "error in command\n";
					exit(1);
				}
				string s = argv[i];
				if (IsInt(argv[i]))	{
					burnin = atoi(argv[i]);
					i++;
					if (i==argc)	{
						cerr << "error in command\n";
						exit(1);
					}
					s = argv[i];
					if (IsFloat(s))	{
						cutoff = atof(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else if (s == "-lexp")	{
				lengthprior = 0;
			}
			else if (s == "-lgam")	{
				lengthprior = 1;
			}
			else if (s == "-gbl")	{
				genebl = 1;
				i++;
				genepartition = argv[i];
				forcecat = 1;
			}
			else if (s == "-mbl")	{
				mbl = 2;
			}
			else if (s == "-mmbl")	{
				mbl = 1;
			}
			else if (s == "-gw")	{
				selmode = 1;
				if (mutmode == 0)	{
					mutmode = 1;
				}
			}
			else if (s == "-hb")	{
				selmode = 0;
				if (mutmode == 0)	{
					mutmode = 1;
				}
			}
			else if (s == "-chky")	{
				mutmode = 2;
			}
			else if (s == "-cgtr")	{
				mutmode = 3;
			}
			else if (s == "-gwf")	{
				i++;
				gwf = atof(argv[i]);
			}
			else if (s == "-qmod")	{
				qmod = 1;
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
			else if (s == "-ar")	{
				clock = 6;
			}
			else if (s == "-bd")	{
				timeprior = 1;
				i++;
				if (IsFloat(argv[i]))	{
					chi = atof(argv[i]);
					i++;
					if (! IsFloat(argv[i]))	{
						cerr << "error : -bd <chi1> >chi2>\n";
						exit(1);
					}
					chi2 = atof(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-dir")	{
				timeprior = 2;
				i++;
				if (IsFloat(argv[i]))	{
					chi = atof(argv[i]);
				}
				else	{
					i--;
				}
			}
			else if (s == "-covarion")	{
				covarion = 1;
				external = 0;	
			}
			else if (s == "-covext")	{
				covarion = 1;
				external = 1;	
			}
			else if (s == "-fast")	{
				zipgtr = 4;
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
			else if (s == "-popeff")	{
				popeff = 1;
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
				//empfreq = 3;
			}
			else if (s == "-waghssp")	{
				rr = waghssp;
				//empfreq = 3;
			}
			else if (s == "-lg")	{
				rr = lg;
				//empfreq = 3;
			}
			else if (s == "-mtrev")	{
				rr = mtrev;
				//empfreq = 3;
			}
			else if (s == "-gtr")	{
				rr = gtr;
			}
			else if (s == "-rr")	{
				i++;
				rrfile = argv[i];
				rr = custom;
			}
			else if (s == "-freefreq")	{
				empfreq = 0;
			}
			else if (s == "-matfreq")	{
				empfreq = 3;
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
				discrate = 0;
				ras = dprate;
			}
			else if (s == "-catfix")	{
				ncat = -2;
				i++;
				statfix = argv[i];
				forcecat = 1;
				if ((statfix=="ul2")||(statfix == "ul3")||(statfix == "UL2")||(statfix == "UL3")){
					rr = gtr;
					qmod = 1;
				}
			}
			else if (s == "-fw")	{
				freeweight = 1;
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
				x_option = 1;
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
		cerr << "pb [options] <chainname>\n";
		cerr << "\tcreates a new chain, sampling from the posterior distribution, conditional on specified data\n";
		cerr << "\n";
		cerr << "data options:\n";
		cerr << "\t-d <filename>      : file containing an alignment in phylip or nexus format; dna, rna or amino acids\n";
		cerr << "\t-cov <varcovar>     : file containing a variance-covariance matrix (e.g. output of readcov or estbranches)\n";
		cerr << "\t-dc                 : deletes constant columns of the alignment\n";
		cerr << "\t-rec <recoding>     : data recoding. <recoding> is a pre-specified recoding (dayhoff6, dayhoff4, hp),\n";
		cerr << "\t                      or a file containing recoding table\n"; // , written as follows:\n";
		// cerr << "\t\tAA   <Letter>\n";
		// cerr << "\t\tAA   <Letter>\n";
		// cerr << "\t\t....\n";
		cerr << "tree:\n";
		cerr << "\t-t <treefile>       : starts from specified tree\n"; 
		cerr << "\t-T <treefile>       : chain run under fixed, specified tree\n"; 
		cerr << "\t-r <outgroup>       : re-root the tree (useful under clock models)\n";
		cerr << "substitution model:\n";
		cerr << "\tmixture configuration\n";
		cerr << "\t\t-cat         : free number of categories (Dirichlet Process, Lartillot and Philippe 2004)\n";
		cerr << "\t\t-ncat <ncat> : fixed number of categories\n";
		cerr << "\t\t-catfix <profiles>    : specifying a fixed set of pre-defined categories\n";
		// cerr << "\t\t\t<profiles> : name of file containing the profiles.\n";
		// cerr << "\t\t\tshould be written like:\n";
		// cerr << "\t\t\t<ncat>\n";
		// cerr << "\t\t\t<freq1.1> <freq1.2> ... <freq1.20>\n";
		// cerr << "\t\t\t<freq2.1> <freq2.2> ... <freq2.20>\n";
		// cerr << "\t\t\t...\n";
		// cerr << "\t\t-max         : one profile per site\n";
		cerr << "\tchoice of relative rates of substitution\n";
		cerr << "\t\t-wag         : Whelan and Goldman 2001\n";	
		cerr << "\t\t-jtt         : Jones, Taylor, Thornton 1992\n";	
		cerr << "\t\t-mtrev       : Hadachi and Hasegawa 1996\n";
		cerr << "\t\t-gtr         : General Time Reversible\n";
		cerr << "\t\t-poisson     : Poisson matrix, all relative rates equal to 1 (Felsenstein 1981)\n";
		cerr << "\tunderlying across-site rate variations\n";
		cerr << "\t\t-uni         : uniform rates across sites\n";
		// cerr << "\t\t-gamma       : gamma distributed rates\n";
		// cerr << "\t\t-invgamma    : gamma distributed rates + one category of invariable sites\n";
		cerr << "\t\t-dgam <ncat> : discrete gamma. ncat = number of categories (4 by default)\n";
		cerr << "\t\t-cgam        : continuous gamma distribution\n";
		// cerr << "\texamples: pb -ncat 10 -poisson -uni -d mydata.emb mychain\n"; 
		// cerr << "\t          pb -max -gtr -gamma -cgam -T treefile -d mydata.emb mychain\n"; 
		// cerr << "\tdefault : -cat -poisson -gamma -dgam 4 \n";
		cerr << "relaxed clock models:\n";
		cerr << "\t-cl  : strict molecular clock\n";
		cerr << "\t-ln  : log normal (Thorne et al, 1998)\n";
		cerr << "\t-cir : CIR process (Lepage et al, 2007)\n";
		cerr << "\t-wn  : white noise (flexible but non-autocorrelated clock)\n";
		cerr << "\t-ugam: independent gamma multipliers (Drummond et al, 2006)\n";
		cerr << "priors on divergence times (default: uniform):\n";
		cerr << "\t-dir : dirichlet\n";
		cerr << "\t-bd  : birth death\n";
		cerr << "\n";
		cerr << "\t-cal <calibrations> : impose a set of calibrations\n";
		cerr << "\t-rp  <age>          : impose a prior mean age for root\n";
		cerr << "\t                      (if not specified, root age is uniformly distributed)\n";
		cerr << "additional options\n";
		cerr << "\t-x <every> <until>  : saving frequency, and chain length\n";
		cerr << "\t-f                  : forcing checks\n";
		cerr << "\t-s                  : \"saveall\" option. without it, only the trees are saved\n";
		cerr << '\n';
		cerr << "pb <name>\n";
		cerr << "\tstarts an already existing chain\n";

		exit(1);
	}

	cerr << '\n';

	// if qsub mode
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

	// if restarting an already existing chain
	if ((autorestart) || ((initfile == "") && (datafile == "") && (seplist == "")))	{
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
			ofstream os((name + ".log").c_str(), APPEND);
			MCParameters* mParam = chain.mParam;

			if (nchain>1)	{
				if (mParam->Nchain != nchain)	{
					cerr << "error : number of chains should be the same as when the chain was created: " << mParam->Nchain << '\n';
					exit(1);
				}
				mParam->BurnIn = burnin;
				mParam->FinalBeta = cutoff;
			}
			if (x_option)	{
				if (every != mParam->SaveEvery)	{
					cerr << "error : saving frequency should be the same as when the chain was created: " << mParam->SaveEvery << '\n';
					exit(1);
				}
				if ((until != -1) && (mParam->HowManySaved > until))	{
					cerr << "error : number of points saved is greater than proposed upper limit\n";
					exit(1);
				}

				mParam->StopAfter = until;
			}

			os << '\n';
			os << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\n';
			os << "how many saved : " << chain.mParam->HowManySaved << '\n';
			os << "\nchain restarted\n\n";
			os.close();

			cout << "current log likelihood: - " << mParam->GetCurrentState()->logSampling() << '\n';
			cerr << "how many saved : " << chain.mParam->HowManySaved << '\n';
			cout << "\nchain restarted\n\n";
			chain.Start();
			exit(1);
		}
	}

	// below: creating a new chain

	// checking whether a chain of same name already exists
	if ((!force) && (ifstream((path + name + ".param").c_str())))	{
		cerr << "Chain " << name << " seems to already exist\n";
		cerr << "use \"-f\" option to override\n";
		exit(1);
	}
	
	// open the log file: will contain all the details about model specification
	ofstream os((name + ".log").c_str());

	os << "init seed : " << Random::init_seed << '\n';
	os.flush();
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
	if (rr != poisson)	{
		mParam->ZipGTR = zipgtr;
		mParam->ZipPrior = zipprior;
		if (zipprior == 2)	{
			mParam->ZipGTRDP = 0;
		}
	}

	mParam->MutMode = mutmode;
	mParam->SelMode = selmode;
	mParam->GWf = gwf;
	mParam->StartDP = nstart;
	if (rrprior)	{
		mParam->RRPrior = GammaDistributed;
	}
	else	{
		mParam->RRPrior = Exponential;
	}
 
	if (tempering)	{
		os << "simulated tempering\n";
		os << "from " << initbeta << " every " << betastep << "\n";
		os << "burnin : " << burnin << '\n';	
		mParam->Tempering = 1;
		mParam->MSMode = Thermo;
		mParam->BurnIn = burnin;
		mParam->InitBeta = initbeta;
		mParam->BetaStep = betastep;
	}
	if (! fixtopo)	{
		if (clock && ! normalapprox)	{
			cerr << "error: cannot use relaxed clock models under free topology\n";
			exit(1);
		}
		if (outgroup != "")	{
			cerr << "error: cannot specifiy an outgroup under free topology\n";
		}
	}
	if (rec != NoRecoding)	{
		deleteconstant = 1;
	}
	if (deleteconstant)	{
		mParam->DeleteConstant = Yes;
	}
	if (nchain > 1)	{
		mParam->Nchain = nchain;
		os << nchain << "chains in parallel\n";
		os << '\n';
	}

	if ((datafile != "") ||(seplist != ""))	{
		if (normalapprox == 0)	{
			cerr << "read data\n";
			cerr.flush();
			mParam->ReadDataFromFile(datafile);
			cerr << "read data ok\n";
			cerr.flush();
			os << "data file : " << datafile << '\n';
			os << "number of taxa : " << mParam->Ntaxa << '\n';
			os << "number of sites: " << mParam->Nsite << '\n';
			/*
			for (int i=0; i<mParam->Ntaxa; i++)	{
				cout << "taxa " << i+1 << '\t' << mParam->SpeciesNames[i] << '\n';
				cout.flush();
			}
			cout << '\n';
			*/
		}
		else if (normalapprox == 1)	{
			os << "normal approximation\n";
			if (seplist != "")	{
				mParam->ReadSepCov(seplist);
				os << "separate analysis\n";
				os << "list of variance covariance matrices : " << seplist << '\n';
			}
			if (datafile != "")	{
				mParam->ReadCov(datafile);
				os << "variance covariance matrix: " << datafile << '\n'; 
			}
		}
	}

	if (prior)	{
		cerr << "prior\n";
		mParam->Beta0 = 0;
	}
	if (contdatafile != "")	{
		mParam->ReadContDataFromFile(contdatafile);
	}

	if (genebl)	{
		os << "gene branch length multipliers\n";
		mParam->GeneBLMode = Yes;
		mParam->GeneBLMultiplier = Yes;
		mParam->SeparateModelSwitch = 1;
		mParam->ReadPartition(genepartition);
		os << "separate model. partition defined in file " << genepartition << '\n';
	}

	if (lengthprior == 0)	{
		os << "prior on lengths: exponential\n";
		mParam->LengthPrior = Exponential;
	}
	if (lengthprior == 1)	{
		os << "prior on lengths: gamma\n";
		mParam->LengthPrior = GammaDistributed;
	}

	if (fixmeanlength)	{
		os << "prior mean length is fixed : " << meanlength << '\n';
	}
	else	{
		os << "prior mean length is random exponential of mean 0.1\n";
	}
	if (mbl)	{
		mParam->MBL = mbl;
	}

	mParam->PopEff = popeff;

	if (rec)	{
		mParam->Recoding = rec;
		mParam->RecodingFile = recfile;
		mParam->LoadRecoding();
		os << "data recoding: ";
		if (rec == Dayhoff6)	{
			os << "dayhoff 6 states\n";
		}
		else if (rec == Dayhoff4)	{
			os << "dayhoff 4 states\n";
		}
		else if (rec == HP)	{
			os << "hydrophobic/polar recoding\n";
		}
		else	{
			os << "custom, as defined in " << recfile << "\n";
		}
	}	
	if (deleteconstant)	{
		os << "constant sites deleted\n";
	}
	os << '\n';

	if (withpartial)	{
		mParam->WithPartial = 1;
		os << "fast computation / high memory use \n";
	}
	else	{
		os << "slow computation / low memory use\n";
	}

	if (fixtopo)	{
		if (treefile == "")	{
			cerr << "error: no topology file was specified\n";
			exit(1);
		}
		os << "fixed tree topology: " << treefile << '\n';
	}
	if (outgroup != "")	{
		os << "outgroup : " << outgroup << '\n';
	}
	os << '\n';
	
	if ((rr != poisson) && (! forcecat))	{
		ncat = 1;
	}
	if (ncat == -2)	{
		os << "fixed profiles taken from " << statfix << '\n';
		if (freeweight)	{
			os << "free weights\n";
		}
		else	{
			ncat = -3;
			os << "fixed weights\n";
		}
	}
	else if (ncat == 0)	{
		os << "cat model: mixture of a free number of profiles\n";
		if (rr != poisson)	{
			// statcenter = 0;
			mParam->AlphaPrior = Exponential;
		}
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
		if (empfreq == 0)	{
			os << "free equilibrium frequencies (inferred from data)\n";
		}
		else if (empfreq == 1)	{
			os << "empirical equilibrium frequencies\n";
		}
		else if (empfreq == 2)	{
			os << "site-specific equilibrium frequencies\n";
		}
		else	{
			os << "fixed equilibrium frequences, from published matrix\n";
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

	os << "relative exchange rates: ";
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
	else if (rr == waghssp)	{
		os << "wag_hssp (Le and Gascuel 2008)\n";
	}
	else if (rr == lg)	{
		os << "lg (Le and Gascuel 2008)\n";
	}
	else if (rr == mtrev)	{
		os << "mtrev\n";
	}
	else if (rr == custom)	{
		os << "custom, from file: " << rrfile << '\n';
	}
		
	if (nh)	{
		mParam->NH = nh;
		mParam->NHNcatMax = nnh;
		mParam->NHPrior = nhprior;
		os << "non-homogeneous model (foster-like)\n";
		os << "prior on distorters: ";
		if (nhprior == 0)	{
			os << "dirichlet\n";
		}
		else if (nhprior == 1)	{
			os << "diagonal gaussian\n";
		}
		else	{
			os << "generalised multivariate gaussian\n";
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

	if (covarion == 1)	{
		os << "covarion model\n";
		mParam->HeteroMode = Covarion;
		if (external)	{
			mParam->ExternalCov = Yes;
		}
		else	{
			mParam->ExternalCov = No;
		}
	}

	if (clock==1)	{
		os << "strict molecular clock\n";
	}
	else if (clock == 2)	{
		os << "lognormal autocorrelated relaxed clock\n";
	}
	else if (clock == 3)	{
		os << "cir autocorrelated relaxed clock\n";
	}
	else if (clock == 4)	{
		os << "flexible, non-autocorrelated clock (white noise process)\n";
	}
	else if (clock == 4)	{
		os << "flexible, non-autocorrelated clock (uncorrelated gamma)\n";
	}
	if (statfix != "")	{
		mParam->ReadStatFix(statfix);
	}
	if (clock)	{
		if (genebl)	{
			cerr << "-gbl option cannot be used in a clock relaxation context\n";
			exit(1);
		}
		if (calibration != "")	{
			mParam->ReadCalibration(calibration);
			os << "calibrations : " << calibration << "\n";
			if (timeprior != 0)	{
				if (chi == -1)	{
					cerr << "error: calibrations are valid only\n";
					cerr << "\t1; under a uniform prior for divergence times\n";
					cerr << "\t2; under a birth death prior for divergence times, but with fixed hyperparameters (-bd <chi1> <chi2>)\n";
					cerr << "\t2; under a Dirichlet prior for divergence times, but with fixed hyperparameter (-dir <chi>)\n";

					exit(1);
				}
			}
		}
		if (meanscale != 0)	{
			if (calibration == "")	{
				cerr << "error : no point specifying a root age in the absence of any calibration\n";
				exit(1);
			}
			if (meanscale <0)	{
				mParam->ScalePrior = PowerLaw;
				os << "prior on root age: powerlaw\n";
			}
			else	{
				mParam->ScalePrior = Exponential;
				mParam->MeanScale = meanscale;
				os << "prior on root age: exponential of mean " << meanscale << "\n";
			}
		}
		else	{
			os << "prior on root age: uniform\n";
		}
		os << "prior for divergence times : ";
		if (timeprior == 1)	{
			mParam->TimePrior = BD;
			mParam->Chi = chi;
			mParam->Chi2 = chi2;
			os << "birth death\n";
		}
		else if (timeprior ==2)	{
			mParam->TimePrior = Dirich;
			mParam->Chi = chi;
			os << "dirichlet\n";
		}
		else	{
			os << "uniform\n";
		}
	}
	if (outgroup != "")	{
		mParam->ReadOutGroup(outgroup);
	}
	os << '\n';

	if (rr == custom)	{
		mParam->ReadCustomRR(rrfile);
	}

	
	// everything is all set. Now call the Init() function
	// to finish the initialisation

	if (qmod)	{
		mParam->Qmode = Yes;
		os << "CAT GTR N model\n";
	}
	if (initfile != "")	{
		mParam->InitFromFile(initfile);
	}
	else	{
		if (clock)	{
			mParam->InitClock(clock, separate);
		}
		if (! normalapprox)	{
			mParam->Init(rr,empfreq, ncat, statcenter, ras, discrate, fixmeanlength, meanlength);
		}
	}

	if (! mParam->DataOK)	{
		cerr << "error: should specify a datafile\n";
		exit(1);
	}


	if (rec)	{
		mParam->WriteDataToFile(mParam->DataFileSpec + ".recoded");
	}

	// if a tree has been specified, set it as initial tree
	if (treefile != "")	{
		ifstream is((path + treefile).c_str());
		mParam->GetCurrentState()->SetTree(is);
		if (fixtopo)	{
			mParam->FixTopo = 1;
		}
	}

	// if no initfile, make the array of MCMC updates
	if (initfile == "")	{
		mParam->InitMove();
	}

	
	cerr << "creating chain\n";
	cerr.flush();
	// create a chain ...
	Chain chain(mParam,name);

	mParam->IncrementalDP = 0;
	mParam->IncrementalRateDP = 0;
	mParam->EmpiricalDP = 0;

	mParam->ZipSub = 1;

	cerr << "rep\talpha\tnmode\tgamma\tlength\n";
	os << "rep\talpha\tnmode\tgamma\tlength\n";
	mParam->SumOverModes = No;
	for (int rep=0; rep<nrep; rep++)	{
		ostringstream s;
		s << name << rep;
		mParam->GetCurrentState()->Reinitialise();
		mParam->Update();
		mParam->GetCurrentState()->mLogPrior = 0;
		mParam->GetCurrentState()->mLogSampling = 0;
		mParam->GetCurrentState()->mLogPosterior = 0;
		chain.SavePoint();
		mParam->GetCurrentState()->SimulateData();
		mParam->DataFileSpec = s.str() + ".phy";
		mParam->WriteDataToFile(mParam->DataFileSpec);
		PhyloBayes* pb = mParam->GetCurrentState();
		os << rep << '\t' << pb->alpha << '\t' << pb->Nmode << '\t' << pb->gamma << '\t' << pb->GetLength() << '\t' << pb->GetStationaryEntropy() << '\n';
		cerr << rep << '\t' << pb->alpha << '\t' << pb->Nmode << '\t' << pb->gamma << '\t' << pb->MeanLength << '\t' << pb->GetLength() << '\t' << pb->GetStationaryEntropy() << '\n';
	}
	os.close();
}


