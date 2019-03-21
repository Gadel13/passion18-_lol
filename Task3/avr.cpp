#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>

using namespace std;

int main(int argc, char** argv)//filename
{
	ifstream IN;
	IN.open(argv[1], ios::in);
	int n;
	double a;
	double rez = 0;



	for(int i = 0; i<60; i++)
	{
		IN >> n >> a;

		rez += a;

	}

	rez = rez / 60;

	IN.close();

	cout << rez << endl;

	return 0;

}