#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
using namespace std;

#define multmatr(x1, x2, y1, y2, z1, z2,mat1,mat2,mat3) for(x1 = 0; x1<x2; x1++) for(y1 = 0; y1<y2; y1++) for(z1 = 0 ; z1<z2; z1++) mat3[i][j] += mat1[i][k] * mat2[k][j]; 

int main (int argc, char** argv) // /.main A.dat B.dat C.dat <mode>
								// mode: 0-ijk 1-ikj 2-kij 3-jik 4-jki 5-kji
{

	if (argc != 5)
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

	if (type1 == type2 && type1 == 'd')
	{
		double **mA;
		mA = new double*[n1];
		for (int i = 0; i < n1; i++)
			mA[i] = new double[m1];
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m1; j++)
				A.read((char*)&mA[i][j],sizeof(double));

		double **mB;
		mB = new double*[n2];
		for (int i = 0; i < n2; i++)
			mB[i] = new double[m2];
		for (int i = 0; i < n2; i++)
			for (int j = 0; j < m2; j++)
				B.read((char*)&mB[i][j],sizeof(double));

		double **mC;
		mC = new double*[n1];
		for (int i = 0; i < n1; i++)
			mC[i] = new double[m2];
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m2; j++)
				mC[i][j] = 0;

		ofstream TIME("DATA.txt", ios::in|ios::app);
		int i,j,k;
		clock_t start,stop;

		switch ( atoi(argv[4]) )
		{
			case 0:
			{
				start = clock();
				multmatr(i, n1, j, m2, k, m1,mA,mB,mC);
				stop = clock();
				TIME << "0 ";
				break;
			}
			case 1:
			{
				start = clock();
				multmatr(i, n1, k, m1, j, m2,mA,mB,mC);
				stop = clock();
				TIME << "1 ";
				break;
			}
			case 2:
			{
				start = clock();
				multmatr(k, m1, i, n1, j, m2,mA,mB,mC);
				stop = clock();
				TIME << "2 ";
				break;
			}
			case 3:
			{
				start = clock();
				multmatr(j, m2, i, n1, k, m1,mA,mB,mC);
				stop = clock();
				TIME << "3 ";
				break;
			}
			case 4:
			{
				start = clock();
				multmatr(j, m2, k, m1, i, n1,mA,mB,mC);
				stop = clock();
				TIME << "4 ";
				break;
			}
			case 5:
			{
				start = clock();
				multmatr(k, m1, j, m2, i, n1,mA,mB,mC);
				stop = clock();
				TIME << "5 ";
				break;
			}
		}


		TIME << (double)(stop-start)/CLOCKS_PER_SEC << "\n";

		TIME.close();



		ofstream rez(argv[3],ios::binary);
		char c = 'd';
		rez.write((char*)&c, sizeof(char) );
		
		rez.write((char*)&n1, sizeof(int) );
		rez.write((char*)&m2, sizeof(int) );

		for (i = 0; i < n1; i++)
		{
			for (j = 0; j < m2; j++)
				rez.write((char*)&mC[i][j], sizeof(double) );
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

	} else
		if (type1 == type2 && type1 == 'f')
	{
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

		ofstream TIME("DATA.txt", ios::in|ios::app);
		int i,j,k;
		clock_t start,stop;

		switch ( atoi(argv[4]) )
		{
			case 0:
			{
				start = clock();
				multmatr(i, n1, j, m2, k, m1,mA,mB,mC);
				stop = clock();
				TIME << "0 ";
				break;
			}
			case 1:
			{
				start = clock();
				multmatr(i, n1, k, m1, j, m2,mA,mB,mC);
				stop = clock();
				TIME << "1 ";
				break;
			}
			case 2:
			{
				start = clock();
				multmatr(k, m1, i, n1, j, m2,mA,mB,mC);
				stop = clock();
				TIME << "2 ";
				break;
			}
			case 3:
			{
				start = clock();
				multmatr(j, m2, i, n1, k, m1,mA,mB,mC);
				stop = clock();
				TIME << "3-jik ";
				break;
			}
			case 4:
			{
				start = clock();
				multmatr(j, m2, k, m1, i, n1,mA,mB,mC);
				stop = clock();
				TIME << "4 ";
				break;
			}
			case 5:
			{
				start = clock();
				multmatr(k, m1, j, m2, i, n1,mA,mB,mC);
				stop = clock();
				TIME << "5 ";
				break;
			}
		}


		TIME << (double)(stop-start)/CLOCKS_PER_SEC << "\n";

		TIME.close();



		ofstream rez(argv[3],ios::binary);
		char c = 'f';
		rez.write((char*)&c, sizeof(char) );
		
		rez.write((char*)&n1, sizeof(int) );
		rez.write((char*)&m2, sizeof(int) );

		for (i = 0; i < n1; i++)
		{
			for (j = 0; j < m2; j++)
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
	

	} 

	
	A.close();
	B.close();


	return 0;
}