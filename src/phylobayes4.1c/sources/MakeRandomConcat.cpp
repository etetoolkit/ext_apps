#include "phylo.h"

int main(int argc, char* argv[])	{

	// initialise random 
	Random::Random();

	if (argc == 1)	{
		cerr << "makerandomconcat list number output\n";
		exit(1);
	}

	cerr << "reading list\n";
	ifstream is(argv[1]);
	int nrep = atoi(argv[2]);
	string infile = argv[1];
	string output = argv[3];

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

	int member[N];
	int ngene[nrep];
	for (int rep=0; rep<nrep; rep++)	{
		ngene[rep] = 0;
	}
	
	for (int i=0; i<N; i++)	{
		member[i] = (int) (Random::Uniform() * nrep);
		ngene[member[i]] ++;
	}

	int total = 0;
	for (int rep=0; rep<nrep; rep++)	{
		total += ngene[rep];
		cerr << rep << '\t' << ngene[rep] << '\n';
	}
	
	for (int rep=0; rep<nrep; rep++)	{
		ostringstream s;
		s << output << "_" << rep;
		string streamname = s.str();
		ofstream os((streamname + ".name").c_str());
		os <<  ngene[rep] << '\n';
		for (int i=0; i<N; i++)	{
			if (member[i] == rep)	{
				os << name[i] << '\n';
			}
		}
		os.close();
		string concat = "../makeprotconcat " + streamname + ".name " + streamname;
		system(concat.c_str()); 
	}
}

	
