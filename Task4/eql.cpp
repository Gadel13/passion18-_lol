#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>



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

	double a,b;

	A.read((char*)&a,sizeof(double));
	B.read((char*)&b,sizeof(double));
	while( !( A.eof() ) || !( B.eof() ) )
	{
		if( abs(a - b) > 0.01 )
			return(1);
		A.read((char*)&a,sizeof(double));
		B.read((char*)&b,sizeof(double));
	}

		if( A.eof()  &&  B.eof())
		{
			return(0);
		}else
			return(1);

}