#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char** argv) // A*b = c    ./main <namefile A> <namefile b> <namefile c> 
{
	int numprocs;
	int myid;
	int mode = 0; // 0 n >= m , 1 n < m
	double* C;

	double maxtime = 0.0;
  	double sumtime = 0.0;

  	double time =0.0;


	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	} 

  	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

	if (myid == 0)
	{
		if (argc != 4)
		{
			cout << "ERROR! argc != 4 " << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}


	}


	ifstream A(argv[1], ios_base::binary);

	if(!A.is_open())
	{	
		cout << "ERROR file not found: " << argv[1] << " " << endl;
		return 1;
	}

	ifstream B(argv[2], ios_base::binary);

	if(!B.is_open())
	{	
		cout << "ERROR file not found: " << argv[2] << " " << endl;
		return 1;
	}


	char type1,type2;

	A.read((char*)&type1,sizeof(char));
	B.read((char*)&type2,sizeof(char));

	int n1,m1,n2,m2;

	A.read((char*)&n1,sizeof(int));
	A.read((char*)&m1,sizeof(int));

	B.read((char*)&n2,sizeof(int));
	B.read((char*)&m2,sizeof(int));

	if(m2 != 1)
	{
		cout << "matrix b is not stlb" << endl;
		return 1;
	}

	if (m1 != n2)
	{
		cout << "ERROR check matrix sizes" << endl;
		return 1;
	}

	if ( n1 >= m1 )
		mode = 0;
	else
		mode = 1;


	if(myid == 0)
	{
		C = new double[n1];
		for(int i = 0; i < n1; i++)
			C[i] = 0;
	}


	if(mode == 0)
	{
		int block = n1/numprocs;
		int osta = n1%numprocs;
		if (block == 0)
			osta = n1;

		int recvcount[numprocs];
		for(int i = 0; i < numprocs; i++)
		{
			if (i == numprocs-1)
				recvcount[i] = block + osta;
			else
				recvcount[i] = block;
		}

		int displs[numprocs];
		for(int i = 0; i < numprocs;i++)
			displs[i] = i*block;


		if ( myid != numprocs-1)
		{
			double** bufA;
			bufA = new double*[block];
			for (int i = 0; i < block; i++)
				bufA[i] = new double[m1];

			A.seekg(block * m1 * myid * sizeof(double),ios_base::cur);


			for (int i = 0; i < block; i++)
				for (int j = 0; j < m1; j++)
					A.read((char*)&bufA[i][j],sizeof(double));

			double* bufB = new double[n2];
			for (int i = 0; i < n2; i++)
				B.read((char*)&bufB[i],sizeof(double));


			double* bufC = new double[block];

			for (int i = 0; i < block; i++)
				bufC[i] = 0;

			time = MPI_Wtime();
			for(int i = 0; i < block; i ++)
				for(int j = 0; j < m1; j++)
					bufC[i] += bufA[i][j]* bufB[j];


			time = MPI_Wtime() - time;
			MPI_Gatherv(bufC, block, MPI_DOUBLE, C, recvcount,displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			for (int i = 0; i < block; i++)
				delete[] bufA[i];
			delete[] bufA;

			delete[] bufB;
			delete[] bufC;

		}
		else
		{
			double** bufA;
			bufA = new double*[block+osta];
			for (int i = 0; i < block+osta; i++)
				bufA[i] = new double[m1];

			A.seekg(block * m1 * myid * sizeof(double),ios_base::cur);

			for (int i = 0; i < block+osta; i++)
				for (int j = 0; j < m1; j++)
					A.read((char*)&bufA[i][j],sizeof(double));

			double* bufB = new double[n2];
			for (int i = 0; i < n2; i++)
				B.read((char*)&bufB[i],sizeof(double));


			double* bufC = new double[block+osta];

			for (int i = 0; i < block+osta; i++)
				bufC[i] = 0;

			time = MPI_Wtime();

			for(int i = 0; i < block+osta; i++)
				for(int j = 0; j < m1; j++)
					bufC[i] += bufA[i][j]* bufB[j];


			time = MPI_Wtime() - time;
			MPI_Gatherv(bufC, block+osta, MPI_DOUBLE, C, recvcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			for (int i = 0; i < block+osta; i++)
				delete[] bufA[i];
			delete[] bufA;

			delete[] bufB;
			delete[] bufC;

		}
	}
	else
	{
		int block = m1/numprocs;
		int osta = m1%numprocs;
		if (block == 0)
			osta = m1;

		if( myid != numprocs-1)
		{
			double** bufA;
			bufA = new double*[n1];
			for (int i = 0; i < n1; i++)
				bufA[i] = new double[block];

			for (int i = 0; i < n1; i++)
			{
				A.seekg(block * myid * sizeof(double),ios_base::cur);

				for (int j = 0; j < block; j++)
					A.read((char*)&bufA[i][j],sizeof(double));
				
				A.seekg((m1 - block * myid - block) * sizeof(double),ios_base::cur);
			}


			double* bufB = new double[block];

			B.seekg(block * myid * sizeof(double),ios_base::cur);
			for (int i = 0; i < block; i++)
				B.read((char*)&bufB[i],sizeof(double));

			double* bufC = new double[n1];

			for (int i = 0; i < n1; i++)
				bufC[i] = 0;

			time = MPI_Wtime();

			for(int i = 0; i < n1; i++)
				for(int j = 0; j < block; j++)
					bufC[i] += bufA[i][j]* bufB[j];

			time = MPI_Wtime() - time;
			MPI_Reduce(bufC, C, n1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			for (int i = 0; i < n1; i++)
				delete[] bufA[i];
			delete[] bufA;

			delete[] bufB;
			delete[] bufC;

		}
		else
		{
			double** bufA;
			bufA = new double*[n1];
			for (int i = 0; i < n1; i++)
				bufA[i] = new double[block + osta];

			for (int i = 0; i < n1; i++)
			{
				A.seekg(block * myid * sizeof(double),ios_base::cur);

				for (int j = 0; j < block + osta; j++)
					A.read((char*)&bufA[i][j],sizeof(double));
			}


			double* bufB = new double[block + osta];

			B.seekg(block * myid * sizeof(double),ios_base::cur);
			for (int i = 0; i < block + osta; i++)
				B.read((char*)&bufB[i],sizeof(double));

			double* bufC = new double[n1];

			for (int i = 0; i < n1; i++)
				bufC[i] = 0;

			time = MPI_Wtime();

			for(int i = 0; i < n1; i++)
				for(int j = 0; j < block + osta; j++)
					bufC[i] += bufA[i][j]* bufB[j];

			time = MPI_Wtime() - time;
			MPI_Reduce(bufC, C, n1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

			for (int i = 0; i < n1; i++)
				delete[] bufA[i];
			delete[] bufA;

			delete[] bufB;
			delete[] bufC;

		}


		A.close();
		B.close();

	}


	if (myid == 0)
	{

		ofstream rez(argv[3],ios::binary);
		char c = 'd';
		rez.write((char*)&c, sizeof(char) );
		
		rez.write((char*)&n1, sizeof(int) );
		rez.write((char*)&m2, sizeof(int) );
		for (int i = 0; i < n1; i++)
			rez.write((char*)&C[i], sizeof(double) );

		rez.close();

		ofstream TT("TT.txt", ios_base::out|ios::app);
		ofstream MT("MT.txt", ios_base::out|ios::app);
		MT << numprocs << " " << maxtime << endl;
		TT << numprocs << " " << sumtime << endl;
		TT.close();
		MT.close();
	}


	MPI_Finalize();
	return 0;

}