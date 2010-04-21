/**
 * Gaussian distribution
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#include "GaussianDistFunc.h"

	/**
	 * 
	 * @param m Mean
	 * @param cov Covariance matrix
	 */
	GaussianDistFunc::GaussianDistFunc(vector<double> m, _Matrix<double> cov)
	{
		double PI = acos(-1.00);
		rsPi = sqrt(2*PI);
		this->mean = m;
		this->covMat = cov;
		
		this->inverseCovMat = _Matrix<double>(cov.m, cov.n);
		_Matrix<double>::inverse(covMat, inverseCovMat);
		this->rsDet = sqrt(covMat.det());
	}
	void GaussianDistFunc::setMean(vector<double> a)
	{
		this->mean = a;
	}
	vector<double> & GaussianDistFunc::getMean()
	{
		return this->mean;
	}

	void GaussianDistFunc::setCovMat(_Matrix<double> mat)
	{
		this->covMat = mat; 
		
		this->inverseCovMat = _Matrix<double>(mat.m, mat.n);
		_Matrix<double>::inverse(covMat, inverseCovMat);
		this->rsDet = sqrt(covMat.det());
	
	}
	
	double GaussianDistFunc::getDensityVal(vector<double> x) 
	{
		// TODO Auto-generated method stub
		double res = 0;
		if (x.size() != mean.size()) return -1;
		vector<double> dist2Mean(x.size());
		for (int i = 0; i < x.size(); i++)
			dist2Mean[i] = x[i] - mean[i];
		_Matrix<double> d2Mat(dist2Mean, 1, dist2Mean.size());
		_Matrix<double> d2MatTran(dist2Mean.size(),1);
	 	_Matrix<double>::transpose(d2Mat, d2MatTran);

		_Matrix<double> temp(1,1);
		_Matrix<double> iTemp(dist2Mean.size(),1);
		_Matrix<double>::multiply(inverseCovMat, d2MatTran, iTemp);
		_Matrix<double>::multiply(d2Mat, iTemp, temp);
		
		double kernelVal = temp.at(0, 0);
		double upper = exp(-kernelVal/2);
		res = upper;
		for (int i = 0; i < x.size(); i++)
			res = res/ rsPi;
		res = res/rsDet;
		return res;
	}
	string GaussianDistFunc::toString() {
		ostringstream res;
		res << "Gassian distribution - \n";
		res << "[Mean]: " ;
		for (int i = 0; i < mean.size(); i++) 
			res << mean[i] << ",";
		res << "\n[Cov]: \n";
		res << covMat;	
		return res.str();
	}
/*
	static public void main(String [] args)
	{
		double [] mean = {1, 1};
		double [] [] sigma = {{0.25 ,0.3},{0.3 , 1}};
		_Matrix sigMat = new _Matrix(sigma);
		double [] x = {1.01, 1.02};
		GaussianDistFunc func = new GaussianDistFunc(mean, sigMat);
		System.out.println("PDF: " + func.getDensityVal(x));
	}
*/

