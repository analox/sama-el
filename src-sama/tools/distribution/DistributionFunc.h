/**
 * Probability Density Functions
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#ifndef DistributionFunc_H
#define DistributionFunc_H
#include <vector>
using namespace std;
class DistributionFunc {
public:
	virtual double getDensityVal(vector<double> x) 
	{
		cout << "This is called" << endl;
		return 0;
	}
};
#endif
