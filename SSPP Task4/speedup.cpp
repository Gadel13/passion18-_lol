#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char**  argv)
{
	ifstream IN(argv[1], ios_base::in);
	ofstream OUT ("speedup.txt", ios_base::out | ios_base::trunc);

	int n;
	double ti1;

	double ti;


	IN >> n;
	IN >> ti;

	if (n != 1)
	{
		cout << "ERROR SPEEDRUN PROG first 1 cpu!!!" << endl;
		return 1;
	}

	ti1 = ti;

	while(!IN.eof())
	{
		OUT << n << ' ' << ti1/(double)ti << " " << ti1/((double)ti * n) << endl;
		IN >> n;
		IN >> ti;

	}


	IN.close();
	OUT.close();
	return 0;
}