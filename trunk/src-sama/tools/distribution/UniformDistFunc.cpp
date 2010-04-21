/**
 * 
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
#include "UniformDistFunc.h"

	UniformDistFunc::UniformDistFunc(vector<double> l, vector<double> u)
	{
		setBounds(l,u);
	}
	
	void UniformDistFunc::setBounds(vector<double> l, vector<double> u)
	{
		if (l.size() != u.size())
		{
			cout << "Dimension mismatched..." << endl;
			return;
		}
		this->lBound = l;
		this->uBound = u;
		this->nDim = lBound.size();
		hyperVol = 1;
		for (int i = 0; i < nDim; i++)
		{
			double temp = uBound[i] - lBound[i];
			if (temp > 0) hyperVol *= temp;
			if (temp < 0)
			{
				cout <<"Dimension mismatched..." << endl;
				return;
			}	
		}
		if (hyperVol != 0) dValue = 1/hyperVol;	
	}

	 double UniformDistFunc::getDensityVal(vector<double> x) {
		// TODO Auto-generated method stub
		bool inRange = true;
		for (int i = 0; i < nDim; i++)
			if ((x[i] < lBound[i]) || (x[i] > uBound[i]))
				inRange = false;
		if (!inRange) {
	//		cout << "Out of bound" << endl;
			return 0;
		}
		return dValue;
	}
	 string UniformDistFunc::toString() {
		ostringstream res;
		res << "Uniform distribution - \n";
		for (int i = 0; i < lBound.size(); i++) {
			res << "[" << lBound[i] << "," << uBound[i] << "]x";
		}
		return res.str();
	}

	int _UniformDistFunc_main()
	{
		
		double a[] = {1, 1};
		vector<double>  l(a, a+2);
		double b[] = {1, 3};
		vector<double>  u(b, b+2);
		double c[] = {1, 2};
		vector<double>  x(c, c+2);
		UniformDistFunc* func = new UniformDistFunc(l, u);
		cout << "PDF: " << func->getDensityVal(x) << endl;
		delete func;
		return 0;	
	}



