#ifndef WeightCal_H
#define WeightCal_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <limits>
#include "DistributionFunc.h"
#include "tools/distribution/UniformDistFunc.h"
#include "tools/distribution/GaussianDistFunc.h"
#include "../Matrix.cpp"
using namespace std;
/**
 * Provide weights to samples according to given distribution
 * @author lemi0005
 *
 */
class WeightCal {
public:
	/**
	 * 
	 * @param p
	 * @param array
	 * @return Vector of weighted RELEVANT samples, otherwise vector of 0s for irrelevant samples
	 */
	static vector<double> weightCal(DistributionFunc* p, _Matrix<double> array);

};
#endif
