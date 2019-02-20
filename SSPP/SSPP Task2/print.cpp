#include <iostream>
#include <fstream>
using namespace std;

int main (int argc, char** argv)
{
	if (argc != 2)
	{
		cout << "ERROR enter correct data" <<endl;
		return 1;
	}

	ifstream pr(argv[1], ios::binary);

	if(!pr.is_open())
	{	
		cout << "ERROR file not found" << endl;
		return 1;
	}

	char c;

	pr.read((char*)&c, sizeof(char) );
	//pr >> c;
	cout << c << " ";

	int n,m;

	pr.read((char*)&n, sizeof(int) );
	//pr >> n;
	cout << n << " ";

	pr.read((char*)&m, sizeof(int) );
	//pr >> m;
	cout << m << " " << endl;

	double k;

	for (int i = 1; i <= n; i++)
	{
		for(int j = 1; j <= m; j++)
		{
			pr.read((char*)&k, sizeof(double) );
			//pr >> k;
			cout << k << " " ;
		}
		cout << endl;
	}



	pr.close();
	return 0;

}