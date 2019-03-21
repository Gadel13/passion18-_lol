#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <string.h>
	
#define eps 0.01

using namespace std;

double normal_dis_gen()
{
	double S = 0;
	for(int i = 0; i < 12; ++i)
	{	
		S += (double)rand()/RAND_MAX;
	}
	return S-6;
}


void fileread(vector< complex<double> > &a, char* filename, unsigned n,int myid, unsigned int local_size, double &norma)
{
	MPI_File IN;
	MPI_Status status;
	unsigned int size = pow(2,n);

	double *tmp_buf;
	tmp_buf = new double[local_size*2];

	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &IN);
	MPI_File_seek(IN,myid*local_size*2*sizeof(double), MPI_SEEK_SET);
	MPI_File_read(IN, tmp_buf, local_size*2, MPI_DOUBLE, &status);
	int i;
	#pragma omp parallel 
	{
		#pragma omp for private(i) reduction(+: norma)
			for(i = 0; i < local_size; i++)
			{
				a[i].real(tmp_buf[2*i]);
				a[i].imag(tmp_buf[2*i + 1]);
				norma += abs(a[i] * a[i]);
			}
	}

	delete[] tmp_buf;
	
	MPI_File_close(&IN);
}

void transform(vector< complex<double> > &a, vector< complex<double> > &b, vector< vector< complex<double> > > &U, unsigned n, int myid, int numprocs, unsigned noise)
{
	MPI_Status status;
	unsigned int size = pow(2,n), ik;
	unsigned int local_size = size/numprocs;
	double *tmp_buf;

	double teta;
	vector< vector<complex<double> > > H(2);
	H[0].reserve(2);
	H[1].reserve(2);

	srand(time(NULL));
	for(int k = 1; k <= n; k++)
	{
		if(noise == 1)
		{
			if(myid == 0)
			teta = eps*normal_dis_gen();
			MPI_Bcast(&teta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}else
			teta = 0;


		H[0][0] = U[0][0]*cos(teta) + U[0][1]*(-sin(teta));
		H[0][1] = U[0][0]*sin(teta) + U[0][1]*cos(teta);
		H[1][0] = U[1][0]*cos(teta) + U[1][1]*(-sin(teta));
		H[1][1] = U[1][0]*sin(teta) + U[1][1]*cos(teta);

		ik = pow(2,n-k);

		int flag = 1; // 0 - cur proc, 1 - other proc
		complex<double> a_tmp;

		int addr = ( myid*local_size )^ik;
		int recv_addr = addr/local_size;
		if (myid == recv_addr)
			flag = 0;
		else
		{
			tmp_buf = new double[2*local_size];
			int i;

			#pragma omp parallel 
			{
				#pragma omp for private(i)
					for(i = 0; i < local_size; i++)
					{
						tmp_buf[2*i] = a[i].real();
						tmp_buf[2*i + 1] = a[i].imag();
					}
			}

			MPI_Sendrecv_replace(tmp_buf,2*local_size,MPI_DOUBLE,recv_addr,123,recv_addr,123,MPI_COMM_WORLD,&status);
		}

		int i;
		#pragma omp parallel
		{
			#pragma omp for private(i, a_tmp)
			for(i = 0; i < local_size; i++)
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
				for(int i = 0; i < local_size; i++)
					a[i] = b[i];
		}
	}

}

int main(int argc, char** argv) // n, <file - 0, rand - 1>, <0 - ideal and (1-F), 1 - noise only>, <if file filename.txt read, if rand filename out gener vec, 0 - no file>,
{
	
	int numprocs, myid, numthread;
	unsigned n, enter, noise;
	n = atoi(argv[1]);
	enter = atoi(argv[2]);
	noise = atoi(argv[3]);

	unsigned int size = pow(2,n);
	unsigned int local_size ;

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

	vector< complex<double> > a(local_size), b(local_size);


  	if (enter)
  	{
  		double time_rand;
  		if(myid == 0)
  			time_rand = MPI_Wtime();

  		MPI_Bcast(&time_rand, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  		int i;
  		#pragma omp parallel private(seed,i)  	
  		{
		  	seed = (long)time_rand + myid*omp_get_num_threads() + omp_get_thread_num();

		  	#pragma omp for reduction(+: norma)
		  		for(i = 0; i < local_size; i++)
				{
		  			//a[i] = (double)(i + myid*local_size);
					a[i].real( (double)rand_r(&seed) );
					a[i].imag( (double)rand_r(&seed) );
					norma += abs(a[i] * a[i]);
				}
		}

		if( strcmp(argv[4], "0") != 0 && argc == 5) 
		{
			if (myid == 0)
  			{
  				ifstream file;
  				file.open(argv[4],ios::trunc | ios::binary | ios::out);
  				file.close();
  			}

			MPI_Status status;
	
			MPI_File OUT;
			MPI_File_open(MPI_COMM_WORLD, argv[4], MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &OUT);
			MPI_File_seek(OUT,myid*local_size*2*sizeof(double), MPI_SEEK_SET);
			double *tmp_buf;
			tmp_buf = new double[2*local_size];

			#pragma omp parallel 
			{
				#pragma omp for private(i)reduction(+: norma)
					for(i = 0; i < local_size; i++)
					{
						tmp_buf[2*i] = a[i].real();
						tmp_buf[2*i + 1] = a[i].imag();
					}
			}

			MPI_File_write(OUT,tmp_buf,local_size*2, MPI_DOUBLE,&status);

			delete[] tmp_buf;

			MPI_File_close(&OUT);
		}

  	}else
  	{
  		fileread(a,argv[4],n,myid,local_size,norma);
  	}

  	double del = 0.0;

  	MPI_Allreduce(&norma, &del, (int)1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  	del = sqrt(del);



	// for(int i = 0; i < local_size; i++)
	// 	cout << a[i] << endl;

	int i;
	#pragma omp parallel 
	{
		#pragma omp for private(i)
	  		for(int i = 0; i < local_size; i++)
	  		{
	  			numthread = omp_get_num_threads();
				a[i] /= del;
	  		}
	}

  	time = -MPI_Wtime();

	transform(a,b,U,n,myid,numprocs,noise);


	time += MPI_Wtime();

	if (myid == 0)
	{
		ifstream file;
		file.open("REZ.dat",ios::trunc | ios::binary | ios::out);
		file.close();
	}

	MPI_Status status;
	
	MPI_File OUT;
	if(noise)
		MPI_File_open(MPI_COMM_WORLD, "REZ.dat", MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &OUT);
	else
		MPI_File_open(MPI_COMM_WORLD, "REZideal.dat", MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL, &OUT);
	MPI_File_seek(OUT,myid*local_size*2*sizeof(double), MPI_SEEK_SET);
	double *tmp_buf;
	tmp_buf = new double[2*local_size];

	#pragma omp parallel
	{
		#pragma omp for private(i) 
			for(i = 0; i < local_size; i++)
			{
				tmp_buf[2*i] = b[i].real();
				tmp_buf[2*i + 1] = b[i].imag();
			}
	}

	MPI_File_write(OUT,tmp_buf,local_size*2, MPI_DOUBLE,&status);

	// for(int i = 0; i < local_size; i++)
	// 	cout << " --> " << b[i] << endl;

	delete[] tmp_buf;

	MPI_File_close(&OUT);


	MPI_Reduce(&time,&max_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);


	if(myid == 0)
	{
		ofstream file;
  		file.open("data.txt",ios::out | ios::app);
  		file << n << endl;
  		if (enter != 0)
  		{
  			file << "random" << endl;
  			if(argc == 5 && strcmp(argv[4], "0") != 0)
  				file << "Rand_rez : " << argv[4] << endl;
  		}else
  		{
  			file << "Read from: " << argv[4] << endl;
  		}

  		if (noise)
  			file << "with noise" << endl;
  		else
  			file << "ideal " << endl;
  		file << "numprocs : " << numprocs << endl;
  		file << "Threads : " << numthread << endl;
  		file << "TIME : " << max_time << endl;

  		file << "--------------------" << endl;
  		file.close();
	}




    MPI_Finalize();
	return 0;
}