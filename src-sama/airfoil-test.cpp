#include<iostream>
#include<cstdlib>
#include<unistd.h>
#include <fstream>
#include <string.h>
#include "tools/utilities.h"
using namespace std;

#define DIM 24
int test_ratio()
{
	chdir("./ratioAirfoil");

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

	chdir("..");
	cout << "End" << endl;
	
	return 0;
}

int test_inverse()
{
	chdir("./airfoil_inverse");

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
	system("./lana 0.5 2.0");
	cout << "Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;

	ifstream rFile("result.dat");
	double a, b, c;
	rFile >> c;
	rFile.close();
	cout << "Fitness : " << c << endl;

	chdir("..");
	cout << "End" << endl;
	
	return 0;
}

int main( int argc, char** argv ) 
{
	if (argc <= 1) {
		cout << "Missing parameter.." << endl;
		exit(0);
	}
	cout << "`inverse' for Airfoil Inverse problem" << endl;
	cout << "`ratio' for Airfoil ratio problem" << endl;
	if ( !strcmp(argv[1], "ratio") ) 
		test_ratio();
	else
		test_inverse();

	return 0;
}
