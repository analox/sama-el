
#ifndef _GaussKernel_H_
#define _GaussKernel_H_


#include <ReClaM/KernelFunction.h>
#include <LinAlg/linalg.h>


//! \brief Guassian Kernel with independent scaling of every axis
//!
//! \f$ k(x, z) = exp(- (x-z)^T D^T D (x-z) ) \f$ with diagonal matrix D
//!
//! In case computation speed matters, it is equivalent
//! to linearly transform the input data according to
//! the diagonal matrix D before applying the standard
//! #RBFKenerl with \f$ \gamma = 1 \f$.
class DiagGaussKernel : public KernelFunction
{
public:
	//! Constructor
	DiagGaussKernel(int dim, double sigma = 1.0)
	{
		// initialize D with the unit matrix
		double gamma = 0.5 / (sigma * sigma);
		parameter.resize(dim);
		parameter = gamma;
	}

	//! Destructor
	~DiagGaussKernel()
	{
	}


	//! Evaluates the kernel function.
	double eval(const double* x1, const double* x2, unsigned int dim)
	{
		int i;
		double a;
		double dist2 = 0.0;
		for (i=0; i<dim; i++)
		{
			a = parameter(i) * (x1[i] - x2[i]);
			dist2 += a * a;
		}
		return exp(-dist2);
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters.
	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		derivative.resize(dim, false);
		int i;
		double a, b;
		double dist2 = 0.0;
		for (i=0; i<dim; i++)
		{
			a = x1[i] - x2[i];
			b = parameter(i) * a;
			dist2 += b * b;
			derivative(i) = -2.0 * a * b;
		}
		double ret = exp(-dist2);
		for (i=0; i<dim; i++) derivative(i) *= ret;
		return ret;
	}

	//! compute the single RBF coefficient gamma, which is
	//! in this context best defined as the determinant of
	//! the matrix D, squared.
	double computeGamma()
	{
		double det = 1.0;
		int p, pc = getParameterDimension();
		for (p=0; p<pc; p++) det *= parameter(p);
		return det * det;
	}
};


//! \brief General Guassian Kernel
//!
//! \f$ k(x, z) = exp(- (x-z)^T M^T M (x-z) ) \f$ with symmetric matrix M
//!
//! In case computation speed matters, it is equivalent
//! to linearly transform the input data according to
//! the matrix M before applying the standard
//! #RBFKenerl with \f$ \gamma = 1 \f$.
//! This matrix can be determined calling #getTransformation.
class GeneralGaussKernel : public KernelFunction
{
public:
	//! Constructor
	GeneralGaussKernel(int dim, double sigma = 1.0)
	{
		// initialize M with the unit matrix
		double gamma = 0.5 / (sigma * sigma);
		parameter.resize(dim * (dim+1) / 2);
		parameter = 0.0;
		int i;
		for (i=0; i<dim; i++) parameter(i*(i+1)/2 + i) = gamma;
	}

	//! Destructor
	~GeneralGaussKernel()
	{
	}


	//! Evaluates the kernel function.
	double eval(const double* x1, const double* x2, unsigned int dim)
	{
		Array<double> d(dim);
		int i, j, k;
		double a;
		double dist2 = 0.0;
		for (i=0; i<dim; i++) d(i) = x1[i] - x2[i];
		for (i=0; i<dim; i++)
		{
			k = i*(i+1)/2;
			a = 0.0;
			for (j=0; j<dim; j++)
			{
				a += d(j) * parameter(k);
				k++;
				if (j >= i) k += j;
			}
			dist2 += a * a;
		}
		return exp(-dist2);
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters.
	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		int size = dim * (dim+1) / 2;
		derivative.resize(size, false);
		Array<double> d(dim);
		Array<double> t(dim);
		int i, j, k;
		double a;
		double dist2 = 0.0;
		double ret;
		for (i=0; i<dim; i++) d(i) = x1[i] - x2[i];
		for (i=0; i<dim; i++)
		{
			k = i*(i+1)/2;
			a = 0.0;
			for (j=0; j<dim; j++)
			{
				a += d(j) * parameter(k);
				k++;
				if (j >= i) k += j;
			}
			t(i) = a;
			dist2 += a * a;
		}
		ret = exp(-dist2);
		k = 0;
		for (i=0; i<dim; i++)
		{
			for (j=0; j<i; j++)
			{
				derivative(k) = -2.0 * (d(i) * t(j) + d(j) * t(i)) * ret;
				k++;
			}
			derivative(k) = -2.0 * (d(i) * t(i)) * ret;
			k++;
		}
		return ret;
	}

	//! This method fills in the quadratic matrix trans
	//! with the linear transformation the kernel applies
	//! to the input data before applying the standard
	//! Gaussian kernel \f$ k(x, z) = exp(-(x-z)^2) \f$.
	void getTransformation(Array<double>& trans)
	{
		int dim = (int)floor(sqrt(2.0 * getParameterDimension()));
		int i, j, k = 0;
		trans.resize(dim, dim, false);
		for (i=0; i<dim; i++)
		{
			for (j=0; j<=i; j++)
			{
				trans(i, j) = parameter(k);
				trans(j, i) = parameter(k);
				k++;
			}
		}
	}

	//! compute the single RBF coefficient gamma, which is
	//! in this context best defined as the determinant of
	//! the matrix M, squared.
	double computeGamma()
	{
		double det;
		int i, j, k = 0;
		int pc = getParameterDimension();

		Array2D<double> M(pc, pc);
		Array2D<double> v(pc, pc);
		Array<double> d(pc);

		for (i=0; i<pc; i++)
		{
			for (j=0; j<=i; j++)
			{
				M(i, j) = parameter(k);
				M(j, i) = parameter(k);
				k++;
			}
		}
		det = detsymm(M, v, d);

		return det * det;
	}
};


#endif
