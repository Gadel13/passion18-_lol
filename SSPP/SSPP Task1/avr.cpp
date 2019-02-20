#include <iostream>
#include <fstream>
using namespace std;
 int main(int argc, char** argv)
 {
 	double time[6];
 	for (int i = 0; i<6; i++)
 		time[i] = 0;

 	ifstream IN(argv[1]);
 	int n;
 	double t;
 	for (int i = 1; i <= 60; i++)
 	{
 		IN >> n;
 		IN >> t;
 		time[n] += t;
 	}
 	IN.close();


 	ofstream OUT(argv[1]);
 	for (int i = 0; i<6; i++)
 	{
 		time[i] /= 10;
 		OUT << i << " " << time[i] << endl;
 	}
 	OUT.close();
 	return 0;

 }

