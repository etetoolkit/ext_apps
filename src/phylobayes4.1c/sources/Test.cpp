#include "phylo.h"

int main(int argc, char* argv[])	{


	string filename = argv[1];

	MCParameters* mParam = new MCParameters;
	mParam->ReadDataFromFile(filename);
	mParam->Recoding = HP;
	mParam->LoadRecoding();
	mParam->RecodeData();

	string outfile = argv[2];
	mParam->WriteDataToFile(outfile);
}
