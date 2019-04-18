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

// n, <1 - rand, 0 - file>, <"if read_file "infile.dat" else "rand_data_file.dat" >, <outfile.dat>
int main(int argc, char** argv) {
    int numprocs, myid;
    unsigned n, enter;
    n = atoi(argv[1]);
    enter = atoi(argv[2]);

    unsigned size = pow(2, n);

    unsigned seed;
    double norma = 0;

    unsigned contr = 0, k = 0;


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

    vector< complex<double> > a(local_size), b(local_size);

    if ( enter ) {
        double time_rand;
        if ( myid == 0 )
            time_rand = MPI_Wtime();

        MPI_Bcast(&time_rand, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        unsigned i;
        #pragma omp parallel private(seed, i)
        {
            seed = (long)time_rand + myid*omp_get_num_threads() + omp_get_thread_num();

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
    } else
    {
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


    int type;  // H - 1, nH - 2, NOT - 3, CNOT - 4, Rw - 5, CRw - 6
    double fi = 0;  // for Rw and CRw

    if (myid == 0) {
        char gate[10];
        cout << "ENTER NAME OF GATE (H/nH/NOT/CNOT/Rw/CRw)" << endl;
        cin >> gate;
        if (strcmp(gate, "H") == 0) {
            type = 1;
            cout << "ENTER number of qbit" << endl;
            cin >> k;
        } else
            if (strcmp(gate, "nH") == 0) {
                type = 2;
            } else
                if (strcmp(gate, "NOT") == 0) {
                    type = 3;
                    cout << "ENTER number of qbit" << endl;
                    cin >> k;
                } else
                    if (strcmp(gate, "CNOT") == 0) {
                        type = 4;
                        cout << "ENTER number of contr qbit" << endl;
                        cin >> contr;
                        cout << "ENTER number of qbit" << endl;
                        cin >> k;
                    } else
                        if (strcmp(gate, "Rw") == 0) {
                            type = 5;
                            cout << "ENTER number of qbit" << endl;
                            cin >> k;
                            cout << "ENTER fi" << endl;
                            cin >> fi;
                        } else
                            if (strcmp(gate, "CRw") == 0) {
                                type = 6;
                                cout << "ENTER number of contr qbit" << endl;
                                cin >> contr;
                                cout << "ENTER number of qbit" << endl;
                                cin >> k;
                                cout << "ENTER fi" << endl;
                                cin >> fi;
                            }
    }

    MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&contr, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fi, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (type == 1)
        H(a, b, n, k);
    else
    {
        if (type == 2)
            nH(a, b, n);
        else
            if (type == 3)
                NOT(a, b, n, k);
            else
                if (type == 4)
                    CNOT(a, b, n, contr, k);
                else
                    if (type == 5)
                        Rw(a, b, n, k, fi);
                    else
                        if (type == 6)
                            CRw(a, b, n, contr, k, fi);
    }


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

    // for (int i = 0; i < numprocs; i++) {
    //     if (i == myid) {
    //         for (int j = 0; j < local_size; j++)
    //             cout << a[j] << " --> " << b[j] << endl;
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    // }


    MPI_Finalize();
    return 0;
}