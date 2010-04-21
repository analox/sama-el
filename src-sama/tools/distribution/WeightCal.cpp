
#include "WeightCal.h"
	/**
	 * 
	 * @param p
	 * @param array
	 * @return Vector of weighted RELEVANT samples, otherwise vector of 0s for irrelevant samples
	 */
	vector<double> WeightCal::weightCal(DistributionFunc* p, _Matrix<double> array)
	{
		double ESP = numeric_limits<double>::min();

		int nSamples = array.m;
		int nDim = array.n;
		vector<double> w(nSamples);		

		double sum = 0;
		
		// for each samples
		for (int i = 0; i < nSamples; i++)
		{
			// initialize 
			w[i] = 0;
			
			// copy values to temp
			vector<double> temp(nDim);
			for (int k = 0; k < nDim; k++)
				temp[k] = array.at(i, k);
	
			w[i] = p->getDensityVal(temp);
			//cout << "Weights " << i << ": " << w[i] << endl;
			// for normalization
			sum += w[i];
		}
		if (sum < nDim*ESP) {
			cout << "Irrelevant samples..." << endl;
			for (int i = 0; i < nSamples; i++) w[i] = 0;
		}
		else {
			for (int i = 0; i < nSamples; i++) w[i] = w[i]/sum;
		}
		
		//for (int i = 0; i < nSamples; i++) cout << w[i] << ", ";
		//cout << endl;
	
		return w;
	}
	
	/* UNIT TEST */
	int _WeightCal_main() 
	{
	// gaussian 
	/*	*/
	double _mean[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	vector<double> mean(_mean, _mean+12);
	double _sigma[] = {3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0,
     			   0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3
							};
	_Matrix<double> sigMat(_sigma, 12, 12);
	cout << "Create Gaussian distribution " << endl;
	GaussianDistFunc dist(mean, sigMat);
	cout << "here";
	fflush(stdout);

		
	// uniform 
	/*
		double _l[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		vector<double> l(_l, _l+12);
		double _u[] = {3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3};
		vector<double> u(_u, _u + 12);
		UniformDistFunc dist(l, u);
	*/
		// input
		double _array[] = {
				1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02,
				2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 
				1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 
				1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5
				};
		_Matrix<double> array(_array, 4, 12);
		double fitness[] = {3, 1, 1, 2};
		
		// weight
		vector<double> weight;
		double sum;
		cout<< "Distribution: " << dist.toString() << endl;
		weight = WeightCal::weightCal(&dist, array);
		sum = 0;
		for (int i = 0; i < weight.size(); i++) {
			cout << weight[i] << ", ";
			sum += weight[i] * fitness[i];
		}
		cout << endl;
		cout << "Average fitness = " << sum << endl;
		return 0;
	}
