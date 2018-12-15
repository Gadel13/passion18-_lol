#include <ctime>
#include <iostream>
#include <fstream>




using namespace std;

int main() {
	srand(time(0));
	int a[512];
	for(int i = 0 ; i < 512; ++i) {
		a[i] = i;
	}
	for (int i = 0 ; i < 512; ++i) {
		int k = rand() % 512;
		int temp = a[i];
		a[i] = a[k];
		a[k] = temp;
	}
	ofstream fout("my.map");
	for (int i = 0 ; i < 512; ++i) {
		fout << (a[i] / 8 / 8 ) % 8 << ' ' << (a[i] / 8) % 8 << ' ' << a[i] % 8 << ' ' << 0 << endl;
	}
	fout.close();
	return 0;
}