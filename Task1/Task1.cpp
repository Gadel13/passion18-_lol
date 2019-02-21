#include <iostream>
#include <complex>
#include <vector>
#include <omp.h>
#include <math.h>
#include <cstdio>
#include <ctime>
#include <cstdlib>

using namespace std;

void transform(vector< complex<double> > &a, vector< complex<double> > &b, vector< vector< complex<double> > > &U, unsigned n, unsigned k)
{
	unsigned int size = pow(2,n), i, ik = pow(2,n-k);

	#pragma omp parallel for private(i) shared(a,b,U,ik) 
		for(i = 0; i < size; i++)
			if(ik & i)
				b[i] = U[1][0]*a[ i^ik ] + U[1][1]*a[i];
			else
				b[i] = U[0][0]*a[i] + U[0][1]*a[i | ik];

}

int main(int argc, char const *argv[]) // n, k
{
	unsigned n,k;
	n = atoi(argv[1]);
	k = atoi(argv[2]);

	unsigned int size = pow(2,n);

	vector< complex<double> > a(size), b(size);
	vector< vector<complex<double> > > U(2);
	U[0].reserve(2);
	U[1].reserve(2);
	U[0][0] = 1.0 / sqrt(2);
	U[0][1] = U[0][0];
	U[1][0] = U[0][0];
	U[1][1] = -U[0][0];

	unsigned seed;
	double time = -omp_get_wtime();
	unsigned int i;
	double norma = 0;

	#pragma omp parallel private(seed, i) shared(a, size, time)
	{
		seed = (long)(-1 * time) + omp_get_thread_num();

		#pragma omp for reduction(+: norma)
			for(i = 0; i < size; i++)
			{
				a[i].real( (double)rand_r(&seed) );
				a[i].imag( (double)rand_r(&seed) );
				norma += abs(a[i] * a[i]);
			}
	}

	//cout << norma << " NORMA!!" << endl;
	double del = sqrt(norma);
	//cout << del << "DEL!!!" << endl;

	#pragma omp parallel for 
	{
		for(i = 0; i < size; i++)
			a[i] /= del;
	}

	transform(a,b,U,n,k);

	time += omp_get_wtime();
	cout << "TIME: " << time << endl;
	for(unsigned int i = 0; i < size; i++)
    {
    	cout << a[i] << " " << b[i] << endl;
    }



	return 0;
}