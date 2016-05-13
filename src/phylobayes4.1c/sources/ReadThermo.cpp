#include "phylo.h"

int main(int argc, char* argv[])	{

	// processing command

	string ChainName = "";
	double min = 0;
	double max = 1;
	int thick = 0;

	if (argc == 1)	{
		cerr << "readthermo -loo -p <chainname>\n";
		exit(1);
	}
	int i = 1;
	while (i < argc)	{
		string s = argv[i];

		if (s == "-th")	{
			thick = 1;
		}
		else if (s == "-min")	{
			i++;
			min = atof(argv[i]);
		}
		else if (s == "-max")	{
			i++;
			max = atof(argv[i]);
		}
		else if (ChainName == "")	{
			ChainName = argv[i];
		}
		else	{
			cerr << "error in command\n";
			exit(1);
		}
		i++;
	}
	

	Chain* chain = new Chain(ChainName,0,"");
	double support = chain->ReadThermo(min,max,thick);
	cout << '\n';
	cout << "support : " << support << '\n';
	cout << '\n';

}

