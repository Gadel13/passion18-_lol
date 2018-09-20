#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

#define EPS 0.0001

using namespace std;

int main(int argc, char** argv) // ./eql <file> <true file>
{

	if (argc != 3)
	{
		cout << "ERROR enter correct data" <<endl;
		return 1;
	}

	ifstream A(argv[1], ios_base::in);

	if(!A.is_open())
	{	
		cout << "ERROR file not found: " << argv[1] << " " << endl;
		return 1;
	}

	ifstream B(argv[2], ios_base::in);

	if(!B.is_open())
	{	
		cout << "ERROR file not found: " << argv[2] << " " << endl;
		return 1;
	}

	char type1,type2;

	A >> type1;
	B >> type2;

	int n1,m1,n2,m2;

	A >> n1 >> m1;
	B >> n2 >> m2;

	if ( m1 != m2 || n1 != n2 )
	{

		cout << "ERROR incorrect size of matrix A or B" << endl;
		return 1;
	}


	double a,b;
	
		for (int i = 0; i < n1; i++)
			for (int j = 0; j < m1; j++)
			{
				A >> a;
				B >> b;
				if ( abs(abs(a) - abs(b) ) <= EPS )
				{
					cout << "FILES " << argv[1] << " and " << argv[2] <<" NOT EQL!!!" << endl;
					A.close();
					B.close();
					return 1;
				}
			}



	return 0;
}