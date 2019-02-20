#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <papi.h>

using namespace std;

#define standart 32
#define optim 52
#define NUM_EVENTS 3


void block_mult(float** &A, float** &B, float** &C, int na, int mab, int mb,  const int sizeofblock, int type)
{	
	int Events[NUM_EVENTS] = { PAPI_L1_DCM, PAPI_L2_DCM, PAPI_TOT_CYC};
	long long values[NUM_EVENTS];

	clock_t start,stop;

	switch(type)
	{
		case 0: //ijk
		{	
			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
				cout << "ERRRRR";
			start = clock();
			for(int i = 0; i < (na/sizeofblock) + ( ( na % sizeofblock == 0 ) ? 0 : 1 ); i++)
				for(int j = 0; j < ( mb/sizeofblock ) + ( ( mb % sizeofblock == 0 ) ? 0 : 1 ); j++)
					for(int k = 0; k < ( mab/sizeofblock ) + ( ( mab % sizeofblock == 0 ) ? 0 : 1 ); k++)
						for(int a = i * sizeofblock; (a < (i+1)*sizeofblock) && (a < na); a++ )
							for(int b = j *sizeofblock; (b < (j+1)*sizeofblock) && (b < mb); b++ )
								for(int c = k * sizeofblock; (c < (k+1)*sizeofblock) && (c < mab); c++ )
									C[a][b] += A[a][c] * B[c][b];
			stop = clock();
			if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
				cout << "ERRRRR";
			break;
		}
		case 1: //ikj
		{	
			if (PAPI_start_counters(Events, NUM_EVENTS) != PAPI_OK)
				cout << "ERRRRR";
			start = clock();
			for(int i = 0; i < (na/sizeofblock) + ( ( na % sizeofblock == 0 ) ? 0 : 1 ); i++)
				for(int k = 0; k < ( mab/sizeofblock ) + ( ( mab % sizeofblock == 0 ) ? 0 : 1 ); k++)
					for(int j = 0; j < ( mb/sizeofblock ) + ( ( mb % sizeofblock == 0 ) ? 0 : 1 ); j++)
						for(int a = i * sizeofblock; (a < (i+1)*sizeofblock) && (a < na); a++ )
							for(int c = k * sizeofblock; (c < (k+1)*sizeofblock) && (c < mab); c++ )
								for(int b = j *sizeofblock; (b < (j+1)*sizeofblock) && (b < mb); b++ )
									C[a][b] += A[a][c] * B[c][b];
			stop = clock();
			if (PAPI_stop_counters(values, NUM_EVENTS) != PAPI_OK)
				cout << "ERRRRR";
			break;
		}
	}

	if(sizeofblock == standart)
	{
		if(type == 0)
		{
			ofstream DATA32ijk("DATA32ijk.txt", ios::in|ios::app);
			DATA32ijk << na << ' ' << type << ' ' << ' ' << sizeofblock << ' ' << (double)(stop-start)/CLOCKS_PER_SEC << ' ' <<values[0] << ' ' << values[1] << ' ' <<  values[2] << endl ;
			DATA32ijk.close();
		}else
		{
			ofstream DATA32ikj("DATA32ikj.txt", ios::in|ios::app);
			DATA32ikj << na << ' ' << type << ' ' << ' ' << sizeofblock << ' ' << (double)(stop-start)/CLOCKS_PER_SEC << ' ' <<values[0] << ' ' << values[1] << ' ' <<  values[2] << endl ;
			DATA32ikj.close();
		}
	}else
	{
		ofstream DATA50ikj("DATA50ikj.txt", ios::in|ios::app);
		DATA50ikj << na << ' ' << type << ' ' << ' ' << sizeofblock << ' ' << (double)(stop-start)/CLOCKS_PER_SEC << ' ' <<values[0] << ' ' << values[1] << ' ' <<  values[2] << endl ;
		DATA50ikj.close();

	}

}

int main (int argc, char** argv) // /.main A.dat B.dat C.dat <mode> <block>
								// mode: 0-ijk 1-ikj ; block: 0-standart, 1-optimal
{

	if (argc != 6)
	{
		cout << "ERROR enter correct data" <<endl;
		return 1;
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

	if ( m1 != n2 )
	{

		cout << "ERROR incorrect size of matrix A and B" << endl;
		return 1;
	}

	if (type1 != 'f' || type2 != 'f')
	{
		cout << "ERROR incorrect type of matrix A or/and B" << endl;
		return 1;
	}

	float **mA;
	mA = new float*[n1];
	for (int i = 0; i < n1; i++)
		mA[i] = new float[m1];
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < m1; j++)
			A.read((char*)&mA[i][j],sizeof(float));

	float **mB;
	mB = new float*[n2];
	for (int i = 0; i < n2; i++)
		mB[i] = new float[m2];
	for (int i = 0; i < n2; i++)
		for (int j = 0; j < m2; j++)
			B.read((char*)&mB[i][j],sizeof(float));

	float **mC;
	mC = new float*[n1];
	for (int i = 0; i < n1; i++)
		mC[i] = new float[m2];
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < m2; j++)
			mC[i][j] = 0;


	if(atoi(argv[5]))
		block_mult(mA, mB, mC, n1, m1, m2, optim, atoi(argv[4]) );
	else
		block_mult(mA, mB, mC, n1, m1, m2, standart, atoi(argv[4]) );

	ofstream rez(argv[3],ios::binary);
	char c = 'f';
	rez.write((char*)&c, sizeof(char) );
	
	rez.write((char*)&n1, sizeof(int) );
	rez.write((char*)&m2, sizeof(int) );

	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < m2; j++)
			rez.write((char*)&mC[i][j], sizeof(float) );
	}

	rez.close();
	for (int i = 0; i < n1; i++)
		delete[] mA[i];
	delete[] mA;

	for (int i = 0; i < n1; i++)
		delete[] mB[i];
	delete[] mB;

	for (int i = 0; i < n1; i++)
		delete[] mC[i];
	delete[] mC;



	
	A.close();
	B.close();

	return 0;
	
}