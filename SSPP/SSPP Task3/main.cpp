#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	bool* array;
	array = new bool[(int)sqrt(b) - 1];
	int numprocs;
	int myid;
	double maxtime = 0.0;
  	double sumtime = 0.0;
	double time =0.0;
  	double sumtime0 = 0.0;
  	int kol = 0;
  	int sumkol = 0;

	int rc;
	if (rc= MPI_Init(&argc, &argv)) 
 	{ 
  		cout << "ERROR! MPI_Init " << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
  	} 

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
	MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

	ofstream OUT(argv[3], ios_base::out|ios::app);


	if( myid == 0)
	{
		if (argc != 4)
		{
			cout << "ERROR! argc != 4 " << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (int i = 0; i < ( (int)sqrt(b)-1 ); i++)
				array[i] = 1;

		sumtime = MPI_Wtime();
		for (int i = 0; i < ( (int)sqrt(b)-1 ); i++)
		{
			if(array[i] == 0)
				continue;
			for(int j = i + i+2; j < ( (int)sqrt(b)-1 ); j+= i+2)
				array[j] = 0;
		}


		if(((int)sqrt(b)+1) > a)
		{

			for (int i = a - 2; i < ((int)sqrt(b) - 1); i++ )
			{
				if(array[i] != 0)
				{
					kol++;
					OUT << i + 2 << " ";
				}
			}
				a = (int)sqrt(b) +1;
		}

		sumtime = MPI_Wtime() - sumtime;

	}
	

	MPI_Bcast(&array[0], (int)sqrt(b) - 1, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int block = (b-a+1)/numprocs;
	int osta = (b-a+1)%numprocs;
	if (block == 0)
		osta == b - a + 1;



	if ( myid != numprocs-1 )
	{
		bool* tmp;
		tmp = new bool[block];
		for(int i = 0; i < block; i++)
			tmp[i] = 1;
		time = MPI_Wtime();
		for(int i = 0; i < block; i++)
		{

			if (tmp[i] == 0)
				continue;

			for(int j = 0; j < ((int)sqrt(b) -1); j++)
			{
				if(array[j] != 0)
					if((a + myid*block + i)%(j+2) == 0)
					{
						tmp[i] = 0;
						break;
					}
			}
		}


		time = MPI_Wtime() - time;

		for(int i = 0; i < block; i++)
		{
			if(tmp[i] != 0)
			{
				OUT << a + myid*block + i << " ";
				kol++;
			}
		}

		delete[] tmp;
		
	}else
	{
		bool* tmp;
		tmp = new bool[block + osta];
		for(int i = 0; i < (block+ osta); i++)
		{
			tmp[i] = 1;
		}
		time = MPI_Wtime();
		for(int i = 0; i < (block+ osta); i++)
		{
			if (tmp[i] == 0)
				continue;

			for(int j = 0; j < ((int)sqrt(b) -1); j++)
			{
				if(array[j] != 0)
				{
					if((a + myid*block + i)%(j+2) == 0)
					{
						tmp[i] = 0;
						break;
					}
				}
			}
		}



		time = MPI_Wtime() - time;

		for(int i = 0; i < (block+ osta); i++)
		{
			if(tmp[i] != 0)
			{
				OUT << a + myid*block + i << " ";
				kol++;
			}
		}
		delete[] tmp;

	}

	MPI_Reduce(&kol, &sumkol, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if(myid == 0)
	{
		cout<< "KOL: " << sumkol << endl;
		ofstream TT("TT.txt", ios_base::out|ios::app);
		ofstream MT("MT.txt", ios_base::out|ios::app);
		MT << numprocs << " " << maxtime << endl;
		TT << numprocs << " " << sumtime << endl;
		TT.close();
		MT.close();
	}

	OUT.close();
	delete[] array;

	MPI_Finalize();
	return 0;


}