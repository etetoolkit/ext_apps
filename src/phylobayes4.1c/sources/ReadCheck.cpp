
#include "phylo.h"


void Cred95(int b, double p, int& infb, int& supb)	{

	double threshold = 0.05;

	double logmult[b+1];
	double logp = log(p);
	double logq = log(1 - p);

	double tot = 0;
	double max = 0;
	for (int i=0; i<=b; i++)	{
		logmult[i] = i * logp + (b-i) * logq;
		for (int j=0; j<i; j++)	{
			logmult[i] += log(b-j) - log(j+1);
		}
			
		if ((i == 0) || (max < logmult[i]))	{
			max = logmult[i];
		}
	}
	for (int i=0; i<b; i++)	{
		tot += exp(logmult[i] - max);
	}
	double mult[b+1];
	for (int i=0; i<=b; i++)	{
		mult[i] = exp(logmult[i] - max) / tot;
	}

	double total = 0;
	infb = 0;
	while (total < threshold)	{
		total += mult[infb];
		infb++;
	}

	total = 0;
	supb = b;
	while (total < threshold)	{
		total += mult[supb];
		supb--;
	}
	cerr << p << '\t' << b << '\t' << infb << '\t' << supb << '\t' << tot << '\n';
}


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
	int topo = 0;

	int clustermodes = 0;
	int fromsample = 0;

 	double mindist = 0.03;
 	int minsize = 10;

	int nrep = -1;

	double cutoff = 0.5;

	int verbose = 0;

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
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-topo")	{
				topo = 1;
			}
			else if (s == "-s")	{
				fromsample = 1;
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
			else if (s == "-v")	{
				verbose = 1;
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
		cerr << "readcheck -nrep <nrep> [-topo -x <burnin> <every> <until>] <basename> \n";
		cerr << '\n';
		exit(1);
	}

	if (SampleName == "")	{
		SampleName = ChainName + "_sample";
	}


	try	{
		string name = ChainName + ".chain";
		MCParameters* mParam = new MCParameters();
		ifstream* true_is = 0;
		if (fromsample || (! ifstream(name.c_str())))	{
			name = ChainName + "_sample.sample";
			fromsample = 1;
			if (! ifstream(name.c_str()))	{
				cerr << "error : does not recognise chain\n";
				exit(1);
			}
			true_is = new ifstream(name.c_str());
			if (nrep == -1)	{
				*true_is >> *mParam;
				int tmp;
				*true_is >> tmp >> tmp >> nrep >> tmp;
			}
		}
		else	{
			true_is = new ifstream(name.c_str());
			ifstream Param_is((ChainName+".param").c_str());
			if (! Param_is)	{
				cerr << "?? non existing chain : " << ChainName << '\n';
				exit(1);
			}
			Param_is >> *mParam;
			if (nrep == -1)	{
				nrep = mParam->HowManySaved;
			}
		}
		mParam->Update();
		PhyloBayes* truepb = mParam->GetCurrentState();
		cerr << nrep << '\n';

		int* tothisto = new int[ncat];
		int* truehisto = new int[ncat];
		for (int i=0; i<ncat; i++)	{
			tothisto[i] = 0;
			truehisto[i] = 0;
		}
		
		double alpha[nrep];
		double gamma[nrep];
		double length[nrep];
		double statent[nrep];
		double rrent[nrep];
		cout << "\t\t\t\talpha\tnmode\tgamma\tlength\tstatent\trrent\n";
		for (int rep=0; rep<nrep; rep++)	{
		cerr << "get true pb\n";
		cerr.flush();
			*true_is >> *truepb;
		cerr << "ok\n";
			if (verbose)	{
				cerr << "rep : " << rep << "\n";
			}
			ostringstream s;
			if (fromsample)	{
				s << ChainName << "_sample_" << rep << ".ali";
			}
			else	{
				s << "check" << ChainName << rep << ".phy";
			}

			cerr << s.str() << "\t";
			Sample sample(s.str(),burnin,every,until,Path);
			cerr << "check\n";
			sample.Check(truepb,verbose, topo, ncat, truehisto, tothisto, alpha[rep],gamma[rep],length[rep],statent[rep],rrent[rep]);
			cerr << "ok\n";
		}	
		ofstream os((ChainName + ".summary").c_str());
		os << '\n';
		os << "#ext\tnmode\tgamma\tlength\tstatent\trrent\n";
		os << '\n';
		double x = 100.0 / ncat / 2;
		for (int i=1; i<ncat; i++)	{
			int alphaprop = 0;
			int gammaprop = 0;
			int lengthprop = 0;
			int statentprop = 0;
			int rrentprop = 0;
			for (int rep=0; rep<nrep; rep++)	{
				if ((alpha[rep] < (i*x)) || (alpha[rep] > (100-i*x))) alphaprop++;
				if ((gamma[rep] < (i*x)) || (gamma[rep] > (100-i*x))) gammaprop++;
				if ((length[rep] < (i*x)) || (length[rep] > (100-i*x))) lengthprop++;
				if ((statent[rep] < (i*x)) || (statent[rep] > (100-i*x))) statentprop++;
				if ((rrent[rep] < (i*x)) || (rrent[rep] > (100-i*x))) rrentprop++;
			}
			int supb, infb;
			Cred95(nrep, 2*i*x / 100,infb,supb);
			os <<  (int) 2*i*x << '\t' << alphaprop << '\t' << gammaprop << '\t' << lengthprop << '\t' << statentprop << '\t' << rrentprop << '\t' << '[' << infb << ';' << supb << ']' << '\n';
		}
		os << '\n';

		if (topo)	{
			ofstream os((ChainName + ".topo.summary").c_str());
			os << "#topo: \n";
			for (int i=0; i<ncat; i++)	{
				os << (int) (2*i+1)*x << '\t';
				if (tothisto[i])	{
					os << (int) ((100.0 * truehisto[i]) / tothisto[i]);
				}
				else	{
					os << "-";
				}
				int supb, infb;
				Cred95(tothisto[i], (2*i+1)*x / 100,infb,supb);
				os << '\t' << truehisto[i] << '\t' << tothisto[i] << '\t' << "[" << infb << "," << supb << "]" << '\n';
			}
		}
		os.close();
		string cat = "cat " + ChainName + ".summary";
		string cat2 = "cat " + ChainName + ".topo.summary";
		system(cat.c_str());
		system(cat2.c_str());
	}
	catch(...)	{
		exit(1);
	}
}
