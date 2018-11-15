#include <iostream>
#include <pthread.h>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

typedef struct
{
	int* kolthr;
	bool* arraythr;
	double* timeth;
	int arrayend;
	int athr;
	int bthr;
	char* outfile;
} pthrData;

void* threadFunc(void* thread_data)
{
	pthrData* data = (pthrData*) thread_data;
	struct timespec start,stop;

	bool* tmp;
	tmp = new bool[data->bthr - data->athr + 1];
	for(int i = data->athr; i <= data->bthr; i++)
	{
		tmp[i - data->athr] = 1;
	}
	clock_gettime(CLOCK_REALTIME, &start);
	for(int i = data->athr; i <= data->bthr; i++)
	{

		if (tmp[i - data->athr] == 0)
			continue;

		for(int j = 0; j < data->arrayend; j++)
		{
			if(data->arraythr[j] != 0)
				if((i)%(j+2) == 0)
				{
					tmp[i - data->athr] = 0;
					break;
				}
		}
	}
	clock_gettime(CLOCK_REALTIME, &stop);

	*(data->timeth) += (double)(stop.tv_sec-start.tv_sec + 1e-9 * (stop.tv_nsec - start.tv_nsec));

	ofstream OUT(data->outfile, ios_base::out|ios::app);

	for(int i = data->athr; i <= data->bthr; i++)
	{
		if(tmp[i - data->athr] != 0)
		{
			OUT << i << endl;
			*(data->kolthr) += 1;
		}
	}

	OUT.close();

	delete[] tmp;
}


int main(int argc, char** argv) // ./main a b rez.txt <num threads>
{
	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	bool* array;
	array = new bool[(int)sqrt(b) - 1];
	int numprocs = atoi(argv[4]);
	int kolforthr[numprocs];
	double time[numprocs];
	double maxtime = 0.0;
	double sumtime = 0.0;

	struct timespec start,stop;

	for(int i = 0; i < numprocs; i++)
	{
		kolforthr[i] = 0;
		time[i] = 0.0;
	}

	int kol = 0;

	if (argc != 5)
	{
		cout << "ERROR! argc != 5 " << endl;
		return 1;
	}

	for (int i = 0; i < ( (int)sqrt(b)-1 ); i++)
			array[i] = 1;


	clock_gettime(CLOCK_REALTIME, &start);

	for (int i = 0; i < ( (int)sqrt(b)-1 ); i++)
	{
		if(array[i] == 0)
			continue;
		for(int j = i + i+2; j < ( (int)sqrt(b)-1 ); j+= i+2)
			array[j] = 0;
	}

	clock_gettime(CLOCK_REALTIME, &stop);
	

	if(((int)sqrt(b)+1) > a)
	{
		ofstream OUT(argv[3], ios_base::out|ios::app);

		for (int i = a - 2; i < ((int)sqrt(b) - 1); i++ )
		{
			if(array[i] != 0)
			{
				kol++;
				OUT << i + 2 << endl;
			}
		}

		a = (int)sqrt(b) +1;

		OUT.close();
	}


	time[0] += (double)(stop.tv_sec-start.tv_sec + 1e-9 * (stop.tv_nsec - start.tv_nsec));

	pthread_t* threads = (pthread_t*) malloc(numprocs * sizeof(pthread_t));
	pthrData* threadData = (pthrData*) malloc(numprocs * sizeof(pthrData));

	int block = (b-a+1)/numprocs;
	int osta = (b-a+1)%numprocs;
	if (block == 0)
		osta = b - a + 1;

	for (int i = 0; i < numprocs; i++)
	{
		threadData[i].kolthr = &kolforthr[i];
		threadData[i].arraythr = array;
		threadData[i].arrayend = (int)sqrt(b) - 1;
		threadData[i].outfile = argv[3];
		threadData[i].timeth = &time[i];

		threadData[i].athr = a + i*block;
		if (i == (numprocs -1) )
			threadData[i].bthr = a + (i+1)*block + osta - 1;
		else
			threadData[i].bthr = a + (i+1)*block - 1;


		pthread_create(&(threads[i]), NULL, threadFunc, &threadData[i]);

	}

	for(int i = 0; i < numprocs; i++)
		pthread_join(threads[i], NULL);

	for(int i = 0; i < numprocs; i++)
	{
		sumtime += time[i];
		kol += kolforthr[i];

		if( time[i] > maxtime)
			maxtime = time[i];
	}

	cout << "KOL: " << kol << endl;

	ofstream TT("TT.txt", ios_base::out|ios::app);
	ofstream MT("MT.txt", ios_base::out|ios::app);
	MT << numprocs << " " << maxtime << endl;
	TT << numprocs << " " << sumtime << endl;
	TT.close();
	MT.close();

	free(threads);
	free(threadData);

	delete[] array;

	return 0;


}