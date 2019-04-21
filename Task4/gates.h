#ifndef gates
#define gates

#include <mpi.h>
#include <omp.h>
#include <complex>
#include <vector>
#include <cmath>
using namespace std;

void NOT(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k);
void H(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k);
void nH(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n);
void CNOT(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned control, unsigned k);
void Rw(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k, double fi);
void CRw(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned control, unsigned k, double fi);

#endif
