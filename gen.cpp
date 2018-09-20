#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
//#include <cstdlib>
#include <string.h>

using namespace std;
 
int main(int argc, char** argv) // ./get <type> <n> <m> <name.dat>
{
	if (argc != 5)
	{
		cout << "ERROR enter correct data (./get <type> <n> <m> <name.dat>)" << endl;
		return 1;
	}

	ofstream rez(argv[4]);
	rez << argv[1] << " " << atoi(argv[2]) << " " << atoi(argv[3]) << "\n";

	srand( time(0) );
	for (int i = 1; i <= atoi(argv[2]); i++)
	{
		for (int j = 1; j <= atoi(argv[3]); j++ )
		{
			if ( strcmp( argv[1] , "d" ) == 0 )
			{
				double dd = -5000;
				int t = rand()%10000;
				dd += (double)rand()/(double)RAND_MAX + (double)(t);
				rez << dd << " " ;
			}
			else
			{
				float ff = -2500;
				int t = rand()%5000;
				ff += (float)rand()/(float)RAND_MAX + (float)(t);
				rez << ff << " " ;
			}

		}
		rez << "\n";
	}

	rez.close();
	return 0;
}