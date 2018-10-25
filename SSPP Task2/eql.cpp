#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

#define EPS 1e-11

using namespace std;

int main(int argc, char** argv) // ./eql <file> <true file>
{

	if (argc != 3)
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
	

	if ( m1 != m2 || n1 != n2 )
	{

		cout << "ERROR incorrect size of matrix A or B" << endl;
		return 1;
	}

	if (type1 == type2 && type1 == 'd')
	{
		double a,b;
	
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m1; j++)
			{
				A.read((char*)&a,sizeof(double));
				B.read((char*)&b,sizeof(double));


				if ( std::abs(std::abs(a) - std::abs(b) ) > EPS )
				{
					cout << "FILES " << argv[1] << " and " << argv[2] <<" NOT EQL!!!" << endl;
					A.close();
					B.close();
					return 1;
				}
			}
	} else
	 if (type1 == type2 && type1 == 'f')
	{

		float a,b;
	
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m1; j++)
			{
				A.read((char*)&a,sizeof(float));
				B.read((char*)&b,sizeof(float));
				if (  std::abs(std::abs(a) - std::abs(b) ) > EPS )
				{
					cout << "FILES " << argv[1] << " and " << argv[2] <<" NOT EQL!!!" << endl;
					A.close();
					B.close();
					return 1;
				}
			}
	} else 
	{
		cout << "ERROR type1 != type2" << endl;
		A.close();
		B.close();
		return 1;
	}



	return 0;
}