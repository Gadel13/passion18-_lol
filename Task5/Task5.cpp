#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "mpi.h"
#include "omp.h"
#include "string.h"
#include "math.h"
#include "gates.h"

using namespace std;

void fileread(vector< complex<double> > &a, char* filename, unsigned n, int myid, unsigned local_size, double &norma) {
     MPI_File IN;
     MPI_Status status;

     double *tmp_buf;
     tmp_buf = new double[local_size*2];

     MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &IN);
     MPI_File_seek(IN, myid*local_size*2*sizeof(double), MPI_SEEK_SET);
     MPI_File_read_all(IN, tmp_buf, local_size*2, MPI_DOUBLE, &status);
     unsigned i;
     #pragma omp parallel
     {
          #pragma omp for private(i) reduction(+: norma)
               for (i = 0; i < local_size; i++) {
                    a[i].real(tmp_buf[2*i]);
                    a[i].imag(tmp_buf[2*i + 1]);
                    norma += abs(a[i] * a[i]);
               }
     }

     delete[] tmp_buf;
     MPI_File_close(&IN);
}

void QFT(vector< complex<double> > &a, vector< complex<double> > &b, unsigned n) {
    for ( unsigned i = 1; i <= n; i++) {
        H(a, b, n, i);
        a = b;
        unsigned m = 2;
        for (unsigned j = i + 1; j <= n; j++) {
            double fi_m = (2 * M_PI) / (pow(2, m));
            CRw(a, b, n, j, i, fi_m);
            a = b;
            m++;
        }
    }
    a = b;
    MPI_Status status;
    int numprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for (unsigned i = 0; i < a.size(); i++) {
        int my_addr = i + myid*a.size();
        int want_addr = 0;
        for(unsigned j = 0; j < n; j++) {
            unsigned mask = pow(2, j);
            if ( (my_addr & mask) != 0 )
                want_addr = want_addr^((int)pow(2, n-1-j));
        }
        int recv_addr = want_addr/a.size();

        if ( myid != recv_addr ) {
            double atmp[2], btmp[2];
            atmp[0] = a[i].real();
            atmp[1] = a[i].imag();

            MPI_Sendrecv(atmp, 2, MPI_DOUBLE, recv_addr, 123, btmp, 2, MPI_DOUBLE, recv_addr, 123, MPI_COMM_WORLD, &status);
            b[i].real(btmp[0]);
            b[i].imag(btmp[1]);
        } else
            b[i] = a[want_addr%a.size()];
    }
}

// n, <1 - rand, 0 - file>, <"if read_file "infile.dat" else "rand_data_file.dat" >, <outfile.dat>
int main(int argc, char** argv) {
    int numprocs, myid;
    unsigned n, enter;
    n = atoi(argv[1]);
    enter = atoi(argv[2]);
    double max_time = 0;
    double time = 0;

    unsigned size = pow(2, n);

    unsigned seed;
    double norma = 0;

    int rc;
    rc = MPI_Init(&argc, &argv);
    if (rc) {
       cout << "ERROR! MPI_Init " << endl;
       MPI_Abort(MPI_COMM_WORLD, rc);
      }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    unsigned local_size;
    local_size = (unsigned)size/(unsigned)numprocs;

    vector< complex<double> > a(local_size), b(local_size), c;

    if ( enter ) {
        double time_rand;
        if ( myid == 0 )
            time_rand = MPI_Wtime();

        MPI_Bcast(&time_rand, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        unsigned i;
        #pragma omp parallel private(seed, i)
        {
            seed = (int64_t)time_rand + myid*omp_get_num_threads() + omp_get_thread_num();

            #pragma omp for reduction(+: norma)
            for (i = 0; i < local_size; i++) {
            // a[i] = (double)(i + myid*local_size) + 1;
            a[i].real((double)rand_r(&seed));
            a[i].imag((double)rand_r(&seed));
            norma += abs(a[i] * a[i]);
            }
        }

        if ( strcmp(argv[3], "0") != 0 ) {
            if ( myid == 0 ) {
                ifstream file;
                file.open(argv[3], ios::trunc | ios::binary | ios::out);
                file.close();
            }

            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Status status;

            MPI_File OUT;
            MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &OUT);
            MPI_File_seek(OUT, myid*local_size*2*sizeof(double), MPI_SEEK_SET);
            double *tmp_buf;
            tmp_buf = new double[2*local_size];

            #pragma omp parallel
            {
                #pragma omp for private(i)reduction(+: norma)
                    for (i = 0; i < local_size; i++) {
                        tmp_buf[2*i] = a[i].real();
                        tmp_buf[2*i + 1] = a[i].imag();
                }
            }

            MPI_File_write_all(OUT, tmp_buf, local_size*2, MPI_DOUBLE, &status);

            delete[] tmp_buf;

            MPI_File_close(&OUT);
        }
    } else {
        fileread(a, argv[3], n, myid, local_size, norma);
    }

    double del = 0.0;

    MPI_Allreduce(&norma, &del, (int)1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    del = sqrt(del);

    unsigned i;
    #pragma omp parallel
    {
        #pragma omp for private(i)
            for (unsigned i = 0; i < local_size; i++) {
                a[i] /= del;
            }
    }

    c = a;


    time = MPI_Wtime();
    QFT(a, b, n);
    time = MPI_Wtime() - time;

    MPI_Reduce(&time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    if (myid == 0) {
        ifstream file;
        file.open(argv[4], ios::trunc | ios::binary | ios::out);
        file.close();
    }

    MPI_Status status;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File OUT;
    MPI_File_open(MPI_COMM_WORLD, argv[4], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &OUT);
    MPI_File_seek(OUT, myid*local_size*2*sizeof(double), MPI_SEEK_SET);
    double *tmp_buf;
    tmp_buf = new double[2*local_size];

    #pragma omp parallel
    {
        #pragma omp for private(i)
            for (i = 0; i < local_size; i++) {
                tmp_buf[2*i] = b[i].real();
                tmp_buf[2*i + 1] = b[i].imag();
            }
    }

    MPI_File_write_all(OUT, tmp_buf, local_size*2, MPI_DOUBLE, &status);

    delete[] tmp_buf;

    MPI_File_close(&OUT);

    for (int i = 0; i < numprocs; i++) {
        if (i == myid) {
            for (unsigned j = 0; j < local_size; j++)
                cout << c[j] << " --> " << b[j] << endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0) {
        cout << "QFT TIME: " << max_time << endl;
    }


    MPI_Finalize();
    return 0;
}
