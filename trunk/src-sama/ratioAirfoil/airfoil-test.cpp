#include<iostream>
#include<cstdlib>
#include <fstream>
#include "../tools/utilities.h"
using namespace std;

#define DIM 24
int main()
{
	double bounds[DIM][2];
	ifstream bFile("minmax24.dat");
	for(int i = 0; i < DIM; i++) {
		bFile >> bounds[i][0] >> bounds[i][1];
		cout << bounds[i][0] << " " << bounds[i][1] << endl;
	}
	bFile.close();


	double x[DIM];
	ofstream oFile("ex-alpha24.dat");	
	for(int i = 0; i < DIM; i++) {
		x[i] = (bounds[i][0] + bounds[i][1])/2;
		oFile << x[i] << endl;
	}
	oFile.close();

	long long startT = getTime();
	system("./ratioairfoilObjIn 0.5 2.0");
	cout << "Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;

	ifstream rFile("pdResult.dat");
	double a, b, c;
	rFile >> a >> b >> c;
	rFile.close();
	cout << "Fitness : " << c << endl;

	cout << "End" << endl;
	
	return 0;
}
