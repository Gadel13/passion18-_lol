#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
	int a;
	int b = atoi(argv[2]);
	bool* array;
	array = new bool[b-1];
	int numprocs;
	int myid;
	bool* rez;
	rez = new bool[b-1];
	double startwtime = 0.0;
	double time =0.0;
  	double endwtime;

	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	} 

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

	double maxtime = 0.0;

	if( myid == 0)
	{
		if (argc != 4)
		{
			cout << "ERROR! argc != 4 " << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		a = atoi(argv[1]);

		for(int i = 0; i < b - 1; i++)
		{
			array[i] = 1;
			rez[i] = 0;
		}

		startwtime = MPI_Wtime();

		for(int i = 0; ( i <= ( (int)sqrt(b) - 2 ) ); i++ )
		{
			if (array[i] == 0)
				continue;
			for(int j = i+i+2; j < ( b - 1 ) ; j += i+2 )
			{
				array[j] = 0;
			}
		}

		if( ( (int)sqrt(b) + 1 ) > a)
			a = (int)sqrt(b) + 1;

	}

	MPI_Bcast(&array[0], b - 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int block = (b-a+1)/numprocs;
	time = MPI_Wtime();

	if ( myid != numprocs-1)
	{
		for(int i =  a-2 + block*myid; ( i < ( a-2 + block*(myid+1) ) ); i++ )
		{
			if (array[i] == 0)
				continue;
			for(int j = i+i+2; ( j < ( a-2 + block*(myid+1) ) ); j += i+2 )
			{
				array[j] = 0;
			}
		}
	}else
	{
		for(int i =  a-2 + block*myid; ( i < ( a-2 + (b-a+1)%numprocs + block*(myid+1) ) ); i++ )
		{
			if (array[i] == 0)
				continue;
			for(int j = i+i+2; ( j < ( a-2 + (b-a+1)%numprocs + block*(myid+1) ) ); j += i+2 )
			{
				array[j] = 0;
			}
		}
	}

	time = MPI_Wtime() - time;
	MPI_Reduce(&array[0], &rez[0], b-1, MPI_CHAR, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	endwtime = MPI_Wtime();
	ofstream OUT(argv[3], ios_base::out);
	if(myid == 0)
	{
		int kol = 0;
		for(int i = atoi(argv[1]) - 2; i< (atoi(argv[2]) - 1 ); i++)
		{
			if(rez[i] != 0)
			{
				kol++;
				OUT << i+2 << " ";
			}
		}
		cout<< "KOL: " << kol << endl;
		ofstream TT("TT.txt", ios_base::out|ios::app);
		ofstream MT("MT.txt", ios_base::out|ios::app);
		MT << numprocs << " " << maxtime << endl;
		TT << numprocs << " " << (endwtime-startwtime) << endl;
		TT.close();
		MT.close();
		OUT.close();
	}

	delete[] array;
	delete[] rez;

	MPI_Finalize();
	return 0;


}