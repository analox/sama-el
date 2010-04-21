/**
 * 
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#ifndef UniformDistFunc_H
#define UniformDistFunc_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "DistributionFunc.h"
using namespace std;

class UniformDistFunc: public DistributionFunc{
public:
	vector<double> lBound;
	vector<double> uBound;
	int nDim;
	double hyperVol;
	double dValue;
	
	/**
	 * 
	 * @param l Lower bound
	 * @param u Upper bound
	 */
	UniformDistFunc(vector<double> l, vector<double> u);
	
	void setBounds(vector<double> l, vector<double> u);
	
	virtual double getDensityVal(vector<double> x);
	string toString();

};
#endif

