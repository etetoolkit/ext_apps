#include <sstream>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])	{

	double total = 0;
	double oldx = 0;
	double oldy = 0;
	bool started = false;
	ifstream is(argv[1]);
	ofstream os(argv[2]);
	while (! is.eof())	{
		double x,y,offset;
		is >> x >> y >> offset;
		if (!started)	{
			started = true;
		}
		else	{
			total += 0.5 * (y + oldy) * (x - oldx);
		}
		oldx = x;
		oldy = y;
		os << x << '\t' << log(x) / log(10.0) << '\t' << total + offset << '\n';
	}
}


