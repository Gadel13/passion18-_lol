#include <mpi.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <cstdio>
#include <cstdlib>
	
using namespace std;

void fileread(vector< complex<double> > &a, char* filename, unsigned n, unsigned k,int myid, unsigned int local_size, double &norma)
{
	MPI_File IN;
	MPI_Status status;
	unsigned int size = pow(2,n);

	double *tmp_buf;
	tmp_buf = new double[local_size*2];

	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &IN);
	MPI_File_seek(IN,myid*local_size*2*sizeof(double), MPI_SEEK_SET);
	MPI_File_read(IN, tmp_buf, local_size*2, MPI_DOUBLE, &status);

	for(int i = 0; i < local_size; i++)
	{
		a[i].real(tmp_buf[2*i]);
		a[i].imag(tmp_buf[2*i + 1]);
		norma += abs(a[i] * a[i]);
	}

	delete[] tmp_buf;
	
	MPI_File_close(&IN);
}

void transform(vector< complex<double> > &a, vector< complex<double> > &b, vector< vector< complex<double> > > &U, unsigned n, unsigned k, int myid, int numprocs)
{
	MPI_Status status;
	unsigned int size = pow(2,n), ik = pow(2,n-k);
	unsigned int local_size = size/numprocs;

	int flag = 1; // 0 - cur proc, 1 - other proc
	double *tmp_buf;
	complex<double> a_tmp;

	int addr = ( myid*local_size )^ik;
	int recv_addr = addr/local_size;
	if (myid == recv_addr)
		flag = 0;
	else
	{
		tmp_buf = new double[2*local_size];
		for(int i = 0; i < local_size; i++)
		{
			tmp_buf[2*i] = a[i].real();
			tmp_buf[2*i + 1] = a[i].imag();
		}
		MPI_Sendrecv_replace(tmp_buf,2*local_size,MPI_DOUBLE,recv_addr,123,recv_addr,123,MPI_COMM_WORLD,&status);
	}

	for(int i = 0; i < local_size; i++)
	{

		if (!flag)
		{
			a_tmp = a[( ( i + myid*local_size )^ik ) % local_size];
		}
		else
		{
			a_tmp.real(tmp_buf[2*i]);
			a_tmp.imag(tmp_buf[2*i+1]);
		}

		if( ik & ( i + myid*local_size ) )
			b[i] = U[1][0]*a_tmp + U[1][1]*a[i];
		else
			b[i] = U[0][0]*a[i] + U[0][1]*a_tmp;
	}

	if (flag)
		delete[] tmp_buf;

}

int main(int argc, char** argv) // n, k, <file - 0, rand - 1>, <if file filename.txt>
{
	
	int numprocs, myid;
	unsigned n,k, enter;
	n = atoi(argv[1]);
	k = atoi(argv[2]);
	enter = atoi(argv[3]);

	unsigned int size = pow(2,n);
	unsigned int local_size ;

	vector< complex<double> > a(size), b(size);
	vector< vector<complex<double> > > U(2);
	U[0].reserve(2);
	U[1].reserve(2);
	U[0][0] = 1.0 / sqrt(2);
	U[0][1] = U[0][0];
	U[1][0] = U[0][0];
	U[1][1] = -U[0][0];

	unsigned seed;
	double norma = 0;

	double time = 0.0;
	double max_time = 0.0;

	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	}

  	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  	local_size = size/numprocs;

  	time = -MPI_Wtime();

  	if (enter)
  	{
  		double time_rand;
  		if(myid == 0)
  		{
  			time_rand = MPI_Wtime();
  			MPI_Bcast(&time_rand, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  		}
  		seed = (long)time_rand + myid;
  		for(int i = 0; i < local_size; i++)
		{
  			//a[i] = (double)(i + myid*local_size);
			a[i].real( (double)rand_r(&seed) );
			a[i].imag( (double)rand_r(&seed) );
			norma += abs(a[i] * a[i]);
		}

  	}else
  	{
  		fileread(a,argv[4],n,k,myid,local_size,norma);
  	}


  	double del = 0.0;

  	//MPI_Barrier(MPI_COMM_WORLD);
  	MPI_Allreduce(&norma, &del, (int)1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  	del = sqrt(del);


  	for(int i = 0; i < local_size; i++)
		a[i] /= del;


	transform(a,b,U,n,k,myid,numprocs);


	time += MPI_Wtime();
	MPI_Status status;
	
	MPI_File OUT;
	MPI_File_open(MPI_COMM_WORLD, "REZ.dat", MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &OUT);
	MPI_File_seek(OUT,myid*local_size*2*sizeof(double), MPI_SEEK_SET);
	double *tmp_buf;
	tmp_buf = new double[2*local_size];
	for(int i = 0; i < local_size; i++)
	{
		tmp_buf[2*i] = b[i].real();
		tmp_buf[2*i + 1] = b[i].imag();
	}


	MPI_File_write(OUT,tmp_buf,local_size*2, MPI_DOUBLE,&status);

	// for(int i = 0; i < local_size; i++)
	// 	cout << a[i] << " --> " << b[i] << endl;

	delete[] tmp_buf;

	MPI_File_close(&OUT);


	MPI_Reduce(&time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	//MPI_Barrier(MPI_COMM_WORLD);

	if(myid == 0)
		cout << "TIME:" << max_time << endl;




    MPI_Finalize();
	return 0;
}