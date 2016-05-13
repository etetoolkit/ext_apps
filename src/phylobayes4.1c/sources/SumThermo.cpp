#include "phylo.h"

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	double total = 0;
	double totalerror = 0;
	double totaldiscreteerror = 0;
	while (! is.eof())	{
		string tmp;
		is >> tmp;
		double temp = 0;
		is >> temp;
		total += temp;
		is >> temp;
		totalerror += temp;
		is >> temp;
		totaldiscreteerror += temp;
	}
	is.close();
	/*
	ostringstream s;
	s << "echo " << "--- >> " << argv[1] << ';';
	s << "echo " << total << "   " << totalerror << " >> " << argv[1];
	system(s.str().c_str());
	*/
	cerr << '\n';
	cerr << "total : " << total << '\t' << totalerror << '\t' << totaldiscreteerror << '\n';
	cerr << '\n';
	ofstream os(((string) argv[1] + ".total").c_str());
	os << "total          : " << total << '\n';
	os << "total error    : " << totalerror << '\n';
	os << "discrete error : " << totaldiscreteerror << '\n';

}
