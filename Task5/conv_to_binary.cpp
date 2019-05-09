#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>

using namespace std;

int main(int argc, char** argv) {  // ./conv_to_binary <file> <bin file>
    if (argc != 3) {
        cout << "ERROR enter correct data" <<endl;
        return 1;
    }

    ifstream A;
    A.open(argv[1], ios::in);

    if (!A.is_open()) {
        cout << "ERROR file not found: " << argv[1] << " " << endl;
        return 1;
    }

    ofstream B;
    B.open(argv[2], ios::binary | ios::out | ios::trunc);

    double a, b;

    while ( !(A.eof()) ) {
        A >> a >> b;
        B.write((char*)&a, sizeof(double));
        B.write((char*)&b, sizeof(double));
    }

    A.close();
    B.close();

    return 0;
}