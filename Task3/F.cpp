#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>

using namespace std;

int main(int argc, char** argv) // n filename_in_1_noize filename_in_2_ideal filename_out
{
	unsigned n;
	n = atoi(argv[1]);
	unsigned int size = pow(2,n);

	ifstream A,B;
	A.open(argv[2], ios::binary);
	B.open(argv[3], ios::binary);

	complex <double> rez = 0;
	double ar,ai,br,bi;

	vector< complex<double> > a(size), b(size);

	for(int i = 0; i < size; i++)
	{
		A.read((char*)&ar, sizeof(double) );
		A.read((char*)&ai, sizeof(double) );

		B.read((char*)&br, sizeof(double) );
		B.read((char*)&bi, sizeof(double) );

		a[i].real(ar);
		a[i].imag(ai);

		b[i].real(br);
		b[i].imag(-bi);

		rez = rez + a[i]*b[i];
	}

	A.close();
	B.close();

	rez = abs(rez);
	rez = rez*rez;



	ofstream OUT;
	OUT.open(argv[4], ios::out | ios::app);
	OUT << n << "  " << 1 - rez.real() << endl;
	OUT.close();

	return 0;
}