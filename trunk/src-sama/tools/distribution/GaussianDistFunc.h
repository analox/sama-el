/**
 * Gaussian distribution
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#ifndef GaussianDistFunc_H
#define GaussianDistFunc_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include "DistributionFunc.h"
#include "../Matrix.cpp"
using namespace std;

class GaussianDistFunc: public DistributionFunc {
public:
	vector<double> mean;
	_Matrix<double> covMat;
	_Matrix<double> inverseCovMat;
	double rsDet;
	double rsPi;
	/**
	 * 
	 * @param m Mean
	 * @param cov Covariance matrix
	 */
	GaussianDistFunc(vector<double> m, _Matrix<double> cov);
	void setMean(vector<double> a);
	vector<double> & getMean();
	_Matrix<double> & getCovMat() {return this->covMat; }
	void setCovMat(_Matrix<double> mat);
	virtual	double getDensityVal(vector<double> x);

	string toString();
};
#endif
