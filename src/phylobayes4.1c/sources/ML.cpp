#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	int ncycle = 10;
	int burnin = 100;
	int every = 10;
	int until = 1000;
	int mlsize = -1;

	int keep = 0;

	string datafile = "";
	string initfile = "";
	string treefile = "";
	string name = "";
	string path = "./";
	string prepath = "";
	string directory = "";

	RRMode rr = gtr;
	RASMode ras = gam;
	int discrate = 0; // 0 : continuous 
	int empfreq = 0; // 1 : global freq, 2: site specific freq
	int ncat = 1; // 0 : cat, -1 max, otherwise, fixed number of modes
	int statcenter = 1;

	int qmode = 0;
		// 0 : serial submission
		// 1 : parallel submission
		// 2 : qsub submission
		// 3 : qprep submission


	// read arguments
	try	{
		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-t")	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-i")	{
				i++;
				initfile = argv[i];
			}

			else if (s == "-keep")	{
				keep = 1;
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
			}
			else if (s == "-gamma")	{
				ras = gam;
			}
			else if (s == "-invgamma")	{
				ras = invgam;
			}
			else if (s == "-discrate")	{
				discrate =  4;
				i++;
				string temp = argv[i];
				if (IsInt(temp))	{
					discrate = Int(temp);
				}
				else	{
					i--;
				}
			}
			else if (s == "-cat")	{
				ncat = 0;
			}
			else if (s == "-max")	{
				ncat = -1;
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
			}
			else if (s == "-statfree")	{
				statcenter = 1;
			}
			else if (s == "-statflat")	{
				statcenter = 0;
			}
			else if (s == "-mlsize")	{
				i++;
				mlsize = atoi(argv[i]);
			}

			else if (s == "-ncycle")	{
				i++;
				ncycle = atoi(argv[i]);
			}
			else if (s == "-p")	{
				i++;
				path = argv[i];
			}
			else if ((s == "-qsub") || (s == "-qprep"))	{
				if (s == "-qsub")	{
					qmode = 2;
				}
				if (s == "-qprep")	{
					qmode = 3;
				}
				i++;
				prepath = argv[i];
				i++;
				directory = argv[i];
			} 
			else if (s == "-bg")	{
				qmode = 1;
			}

			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
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
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			cerr << "give it a name ! \n";
			throw(0);
		}
	}
	catch(...)	{
		cerr << "ml : wrong command\n";
		exit(1);
	}

	if (mlsize == -1)	{
		mlsize = until;
	}

	ofstream* qsub_os = 0;
	if (qmode > 1)	{
		qsub_os = new ofstream((name + ".qsub").c_str());
	}

	if ((treefile == "") && (initfile == "") && (datafile == ""))	{
		if (qmode>1)	{
			path = prepath + directory;
		}
		Chain chain(name, 0, path);
		chain.ML(ncycle, burnin, every, until, mlsize, keep);
		exit(1);
	}

	int Ntree = 1;
	Tree* tree = 0;

	if (treefile == "")	{ // assumes only one tree; but should be specified by the initfile
		if (initfile == "")	{
			cerr << "error: should specify a tree\n";
			exit(1);
		}
	}
	else	{
		ifstream is(treefile.c_str());
		is >> Ntree;
		tree = new Tree[Ntree];
		for (int n=0; n<Ntree; n++)	{
			tree[n].ReadFromStream(is, 1);
		}
	}


	MCParameters* mParam = new MCParameters();
	if (datafile != "")	{
		mParam->ReadDataFromFile(datafile);
	}

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

	if (! Ntree)	{
		string currentname = name;
		if (qmode == 0)	{
			Chain chain(mParam, currentname, path);
			chain.ML(ncycle, burnin, every, until, mlsize, keep);
		}
		else if (qmode == 1)	{
			Chain chain(mParam, currentname, path);
			ostringstream launch;
			launch << "../ml ";
			if (path != "./")	{
				launch << "-p " << path << ' ';
			}
			launch << "-x " << burnin << ' ' << every << ' ' << until;
			launch << " -ncycle " << ncycle << " -mlsize " << mlsize << ' ';
			if (keep)	{
				launch << " -keep ";
			}
			launch << name << " &";
			cerr << launch.str() << '\n';
			system(launch.str().c_str());
		}
		else	{
			path = prepath + directory;
			Chain chain(mParam, currentname, path);
			ofstream os((currentname + ".launch").c_str());
			os << prepath << "ml ";
			os << " -x " << burnin << ' ' << every << ' ' << until;
			os << " -ncycle " << ncycle;
			os << " -mlsize " << mlsize;
			if (keep)	{
				os << " -keep ";
			}
			os  << " -p "  << path << ' ' << currentname << '\n';
			(*qsub_os) << "qsub -o /baie/home/usertest/out -e /baie/home/usertest/err" << currentname << ".launch\n";
		}
	}
	else	{

		for (int n=0; n<Ntree; n++)	{

			// prune the tree of all supernumerary species, and init the parameter
			Tree inittree(&tree[n]);
			TaxaParameters taxaparam(&inittree);

			for (int i=0; i<taxaparam.Ntaxa; i++)	{
				string species = taxaparam.GetSpeciesName(i);
				int found = 0;
				int j=0;
				while ( (!found) && (j<mParam->Ntaxa) )	{
					found = (species == mParam->SpeciesNames[j]);
					j++;
				}
				if (! found)	{
					// cerr << species << " not found : eliminate...";
					// cerr.flush();
					inittree.Eliminate(species);
					// cerr << "ok\n";
					// cerr.flush();
				}
			}
			inittree.Simplify();
			// inittree.Phylip(cerr,0,0,1,0);
			// cerr << "\n";
			// cerr.flush();

			mParam->GetCurrentState()->SetTree(inittree);
			mParam->GetCurrentState()->SetBranchLengths(Random);

			ostringstream s;
			s << name << "_" << n;
			string currentname = s.str();

			if (qmode == 0)	{
				Chain chain(mParam, currentname, path);
				chain.ML(ncycle, burnin, every, until, mlsize, keep);
			}
			else if (qmode == 1)	{
				Chain chain(mParam, currentname, path);
				ostringstream launch;
				launch << "../ml ";
				if (path != "./")	{
					launch << "-p path ";
				}
				launch << "-x " << burnin << ' ' << every << ' ' << until;
				launch << " -ncycle " << ncycle << " -mlsize " << mlsize << ' ';
				if (keep)	{
					launch << "-keep ";
				}
				launch << name << " &";
				cerr << launch.str() << '\n';
				cerr.flush();
				system(launch.str().c_str());
			}
			else	{
				path = prepath + directory;
				Chain chain(mParam, currentname, path);
				ofstream os((currentname + ".launch").c_str());
				os << prepath << "ml ";
				os << " -x " << burnin << ' ' << every << ' ' << until;
				os << " -ncycle " << ncycle;
				os << " -mlsize " << mlsize;
				if (keep)	{
					os << " -keep ";
				}
				os  << " -p "  << path << ' ' << currentname << '\n';
				(*qsub_os) << "qsub -e " << path << currentname << ".err -o  " << path << currentname << ".out " << currentname << ".launch\n";
			}
		}
	}

	if (qmode>1)	{
		qsub_os->close();
		string chmod = "chmod +x " + name + ".qsub";
		system(chmod.c_str());
	}
	if (qmode == 2)	{
		string submit = "./" + name + ".qsub";
		system(submit.c_str());
	}
}

