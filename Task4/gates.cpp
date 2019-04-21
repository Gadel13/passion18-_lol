#include "gates.h"
#include<iostream>

using namespace std;

void NOT(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k) {
  cout << "NOT GATE" << endl;

  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  unsigned size = pow(2, n), ik;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;

  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  H[0][0] = 0;
  H[0][1] = 1;
  H[1][0] = 1;
  H[1][1] = 0;

  ik = pow(2, n-k);

  int flag = 1;  //   0 - cur proc, 1 - other proc
  complex<double> a_tmp;

  int addr = (myid*local_size)^ik;
  int recv_addr = addr/local_size;
  if (myid == recv_addr) { {
    flag = 0;
  } } else {
    tmp_buf = new double[2*local_size];
    tmp_buf1 = new double[2*local_size];
    unsigned i;

    #pragma omp parallel
    {
      #pragma omp for private(i)
        for (i = 0; i < local_size; i++) {
            tmp_buf1[2*i] = a[i].real();
            tmp_buf1[2*i + 1] = a[i].imag();
        }
      }

      MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
      delete[] tmp_buf1;
  }

  unsigned i;
  #pragma omp parallel
  {
    #pragma omp for private(i, a_tmp)
      for (i = 0; i < local_size; i++) {
        if (!flag) {
          a_tmp = a[((i + myid*local_size)^ik) % local_size];
        } else {
          a_tmp.real(tmp_buf[2*i]);
          a_tmp.imag(tmp_buf[2*i+1]);
        }

        if (ik & (i + myid*local_size))
          b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
        else
          b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
      }
  }

  if (flag)
    delete[] tmp_buf;
}

void H(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k) {
  cout << "H GATE" << endl;

  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned size = pow(2, n), ik;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;


  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  H[0][0] = 1.0 / sqrt(2);
  H[0][1] = H[0][0];
  H[1][0] = H[0][0];
  H[1][1] = -H[0][0];

  ik = pow(2, n-k);

  int flag = 1;  //  0 - cur proc, 1 - other proc
  complex<double> a_tmp;

  int addr = (myid*local_size)^ik;
  int recv_addr = addr/local_size;
  if (myid == recv_addr) { {
    flag = 0;
  } } else {
    tmp_buf = new double[2*local_size];
    tmp_buf1 = new double[2*local_size];
    unsigned i;

    #pragma omp parallel
    {
      #pragma omp for private(i)
        for (i = 0; i < local_size; i++) {
          tmp_buf1[2*i] = a[i].real();
          tmp_buf1[2*i + 1] = a[i].imag();
        }
    }

    MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
    delete[] tmp_buf1;
  }

  unsigned i;
  #pragma omp parallel
  {
    #pragma omp for private(i, a_tmp)
      for (i = 0; i < local_size; i++) {
        if (!flag) {
          a_tmp = a[((i + myid*local_size)^ik) % local_size];
        } else {
          a_tmp.real(tmp_buf[2*i]);
          a_tmp.imag(tmp_buf[2*i+1]);
        }

        if (ik & (i + myid*local_size))
          b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
        else
          b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
      }
  }

  if (flag)
    delete[] tmp_buf;
}

void nH(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n) {
  cout << "nH GATE" << endl;

  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned size = pow(2, n), ik;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;


  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  H[0][0] = 1.0 / sqrt(2);
  H[0][1] = H[0][0];
  H[1][0] = H[0][0];
  H[1][1] = -H[0][0];

  for (unsigned k = 1; k <= n; k++) {
    MPI_Barrier(MPI_COMM_WORLD);
    ik = pow(2, n-k);

    int flag = 1;  //  0 - cur proc, 1 - other proc
    complex<double> a_tmp;

    int addr = (myid*local_size)^ik;
    int recv_addr = addr/local_size;
    if (myid == recv_addr) {
      flag = 0;
    } else {
      tmp_buf = new double[2*local_size];
      tmp_buf1 = new double[2*local_size];
      unsigned i;

      #pragma omp parallel
      {
        #pragma omp for private(i)
          for (i = 0; i < local_size; i++) {
            tmp_buf1[2*i] = a[i].real();
            tmp_buf1[2*i + 1] = a[i].imag();
          }
      }

      MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
      delete[] tmp_buf1;
    }

      unsigned i;
      #pragma omp parallel
      {
        #pragma omp for private(i, a_tmp)
        for (i = 0; i < local_size; i++) {
          if (!flag) {
           a_tmp = a[((i + myid*local_size)^ik) % local_size];
          } else {
            a_tmp.real(tmp_buf[2*i]);
            a_tmp.imag(tmp_buf[2*i+1]);
          }

          if (ik & (i + myid*local_size))
            b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
          else
            b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
        }
      }

      if (flag)
        delete[] tmp_buf;

      #pragma omp parallel
      {
        #pragma omp for private(i)
         for (unsigned i = 0; i < local_size; i++)
          a[i] = b[i];
      }
    }
}

void CNOT(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned control, unsigned k) {
  cout << "CNOT GATE" << endl;
  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned size = pow(2, n), ik, icontrol;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;

  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  //  NOT
  H[0][0] = 0;
  H[0][1] = 1;
  H[1][0] = 1;
  H[1][1] = 0;

  icontrol = pow(2, n-control);
  ik = pow(2, n-k);

  int flag = 1;  //  0 - cur proc, 1 - other proc
  complex<double> a_tmp;

  int addr = (myid*local_size)^ik;
  int recv_addr = addr/local_size;
  if (myid == recv_addr) {
    flag = 0;
  } else {
    tmp_buf = new double[2*local_size];
    tmp_buf1 = new double[2*local_size];
    unsigned i;

    #pragma omp parallel
    {
      #pragma omp for private(i)
        for (i = 0; i < local_size; i++) {
          tmp_buf1[2*i] = a[i].real();
          tmp_buf1[2*i + 1] = a[i].imag();
        }
    }

    MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
    delete[] tmp_buf1;
  }

  unsigned i;
  #pragma omp parallel
  {
    #pragma omp for private(i, a_tmp)
      for (i = 0; i < local_size; i++) {
        if (!flag) {
          a_tmp = a[((i + myid*local_size)^ik) % local_size];
        } else {
          a_tmp.real(tmp_buf[2*i]);
          a_tmp.imag(tmp_buf[2*i+1]);
        }

        if (icontrol & (i + myid*local_size)) {  //   if control == 1 do NOT
          if (ik & (i + myid*local_size)) {
            b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
          } else {
            b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
          }
        } else {
          b[i] = a[i];
        }
      }
  }

  if (flag)
    delete[] tmp_buf;
}

void Rw(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned k, double fi) {
  cout << "Rw GATE" << endl;
  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned size = pow(2, n), ik;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;


  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  H[0][0] = 1;
  H[0][1] = 0;
  H[1][0] = 0;
  H[1][1].real(cos(fi));
  H[1][1].imag(sin(fi));

  ik = pow(2, n-k);

  int flag = 1;  //  0 - cur proc, 1 - other proc
  complex<double> a_tmp;

  int addr = (myid*local_size)^ik;
  int recv_addr = addr/local_size;
  if (myid == recv_addr) {
    flag = 0;
  } else {
    tmp_buf = new double[2*local_size];
    tmp_buf1 = new double[2*local_size];
    unsigned i;

    #pragma omp parallel
    {
      #pragma omp for private(i)
        for (i = 0; i < local_size; i++) {
          tmp_buf1[2*i] = a[i].real();
          tmp_buf1[2*i + 1] = a[i].imag();
        }
    }

      MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
      delete[] tmp_buf1;
    }

  unsigned i;
  #pragma omp parallel
  {
    #pragma omp for private(i, a_tmp)
    for (i = 0; i < local_size; i++) {
      if (!flag) {
        a_tmp = a[((i + myid*local_size)^ik) % local_size];
      } else {
          a_tmp.real(tmp_buf[2*i]);
          a_tmp.imag(tmp_buf[2*i+1]);
        }

        if (ik & (i + myid*local_size))
          b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
        else
          b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
      }
  }

  if (flag)
    delete[] tmp_buf;
}

void CRw(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n, unsigned control, unsigned k, double fi) {
  cout << "CRw GATE" << endl;
  MPI_Status status;
  int numprocs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  unsigned size = pow(2, n), ik, icontrol;
  unsigned local_size = size/numprocs;
  double *tmp_buf;
  double *tmp_buf1;


  vector< vector<complex<double> > > H(2);
  H[0].reserve(2);
  H[1].reserve(2);

  H[0][0] = 1;
  H[0][1] = 0;
  H[1][0] = 0;
  H[1][1].real(cos(fi));
  H[1][1].imag(sin(fi));

  ik = pow(2, n-k);
  icontrol = pow(2, n-control);

  int flag = 1;  //  0 - cur proc, 1 - other proc
  complex<double> a_tmp;

  int addr = (myid*local_size)^ik;
  int recv_addr = addr/local_size;
  if (myid == recv_addr) {
    flag = 0;
  } else {
    tmp_buf = new double[2*local_size];
    tmp_buf1 = new double[2*local_size];
    unsigned i;

    #pragma omp parallel
    {
      #pragma omp for private(i)
        for (i = 0; i < local_size; i++) {
          tmp_buf1[2*i] = a[i].real();
          tmp_buf1[2*i + 1] = a[i].imag();
        }
    }

    MPI_Sendrecv(tmp_buf1, 2*local_size, MPI_DOUBLE, recv_addr, 123, tmp_buf, 2*local_size, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
    delete[] tmp_buf1;
  }

  unsigned i;
  #pragma omp parallel
  {
    #pragma omp for private(i, a_tmp)
      for (i = 0; i < local_size; i++) {
        if (!flag) {
          a_tmp = a[((i + myid*local_size)^ik) % local_size];
        } else {
          a_tmp.real(tmp_buf[2*i]);
          a_tmp.imag(tmp_buf[2*i+1]);
        }

        if (icontrol & (i + myid*local_size)) {  // if control == 1 do Rw
          if (ik & (i + myid*local_size)) {
            b[i] = H[1][0]*a_tmp + H[1][1]*a[i];
          } else {
            b[i] = H[0][0]*a[i] + H[0][1]*a_tmp;
          }
        } else {
          b[i] = a[i];
        }
      }
  }

  if (flag)
    delete[] tmp_buf;
}