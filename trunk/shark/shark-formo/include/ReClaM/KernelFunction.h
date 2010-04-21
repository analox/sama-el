//===========================================================================
/*!
 *  \file KernelFunction.h
 *
 *  \brief This file contains the definition of a kernel
 *         function as well as basic examples and related stuff.
 *
 *  \author  T. Glasmachers
 *  \date    2005
*/

#ifndef _KernelFunction_H_
#define _KernelFunction_H_


#include <math.h>
#include <ReClaM/Model.h>


//! forward declaration
class C_SVM;


/*!
 *
 *  \brief Definition of a kernel function as a ReClaM model
 *
 *  Please note that the interpretation of a kernel function
 *  as a ReClaM Model involves some complications. We have
 *  to cope with the problem that the Model usually produces
 *  one output per input, while the KernelFunction produces
 *  one output for a pair of inputs. For simplicity, the
 *  input data are required to be matrices composed of exactly
 *  two columns for KernelFunction Models.
 *  For this reason, the KernelFunction provides an additional
 *  interface through its virtual members #eval and evalDerivative
 *  which should be overriden. The inherited members #model
 *  and #modelDerivative redirect to this interface.
 *
 */
class KernelFunction : public Model
{
public:
	//! Constructor
	KernelFunction()
	{
	}

	//! Destructor
	virtual ~KernelFunction()
	{
	}


	//! Evaluates the kernel function.
	virtual double eval(const double* x1, const double* x2, unsigned int dim) = 0;

	//! Evaluates the kernel function on a const object.
	double eval(const double* x1, const double* x2, unsigned int dim) const
	{
		KernelFunction* pT = const_cast<KernelFunction*>(this);
		return pT->eval(x1, x2, dim);
	}

	//! Evaluates the kernel function.
	double eval(const Array<double>& x1, const Array<double>& x2)
	{
		SIZE_CHECK(x1.ndim() == 1);
		SIZE_CHECK(x2.ndim() == 1);
		RANGE_CHECK(x1.nelem() == x2.nelem());

		return eval(x1.elemvec(), x2.elemvec(), x1.nelem());
	}

	//! Evaluates the kernel function on a const object.
	double eval(const Array<double>& x1, const Array<double>& x2) const
	{
		KernelFunction* pT = const_cast<KernelFunction*>(this);
		return pT->eval(x1, x2);
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters.
	virtual double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		int i, ic = parameter.dim(0);
		double temp;
		double ret = eval(x1, x2, dim);
		derivative.resize(ic, false);
		for (i=0; i<ic; i++)
		{
			temp = parameter(i);
			parameter(i) += epsilon;
			derivative(i) = (eval(x1, x2, dim) - ret) / epsilon;
			parameter(i) = temp;
		}
		return ret;
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters
	//! on a const object.
	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative) const
	{
		KernelFunction* pT = const_cast<KernelFunction*>(this);
		return pT->evalDerivative(x1, x2, dim, derivative);
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters.
	double evalDerivative(const Array<double>& x1, const Array<double>& x2, Array<double>& derivative)
	{
		SIZE_CHECK(x1.ndim() == 1);
		SIZE_CHECK(x2.ndim() == 1);
		SIZE_CHECK(x1.nelem() == x2.nelem());

		return evalDerivative(x1.elemvec(), x2.elemvec(), x1.nelem(), derivative);
	}

	//! Evaluates the kernel function and computes
	//! its derivatives w.r.t. the kernel parameters
	//! on a const object.
	double evalDerivative(const Array<double>& x1, const Array<double>& x2, Array<double>& derivative) const
	{
		KernelFunction* pT = const_cast<KernelFunction*>(this);
		return pT->evalDerivative(x1, x2, derivative);
	}

	//! The Model behaviour of the KernelFunction is to interpret
	//! the input as a matrix of two vectors and to compute the
	//! kernel value on them.
	void model(const Array<double>& input, Array<double> &output)
	{
		output.resize(1, false);
		output(0) = eval(input[0], input[1]);
	}

	//! Same as #model, additionally the derivatives w.r.t. all kernel
	//! parameters are computed.
	void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		evalDerivative(input[0], input[1], derivative);
	}

	//! Same as #model, additionally the derivatives w.r.t. all kernel
	//! parameters are computed.
	void modelDerivative(const Array<double>& input, Array<double> &output, Array<double>& derivative)
	{
		output.resize(1, false);
		output(0) = evalDerivative(input[0], input[1], derivative);
	}

	friend class C_SVM;
};


//! \brief Linear Kernel, parameter free
class LinearKernel : public KernelFunction
{
public:
	LinearKernel()
	{
		parameter.resize(0, false);
	}

	~LinearKernel()
	{
	}


	double eval(const double* x1, const double* x2, unsigned int dim) const
	{
		unsigned int i;
		double ret = 0.0;
		for (i=0; i<dim; i++) ret += x1[i] * x2[i];
		return ret;
	}

	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		derivative.resize(0, false);
		return KernelFunction::eval(x1, x2, dim);
	}
};


//! \brief Polynomial Kernel
class PolynomialKernel : public KernelFunction
{
public:
	PolynomialKernel(int degree, double offset)
	{
		parameter.resize(2, false);
		parameter(0) = degree;
		parameter(1) = offset;
	}

	~PolynomialKernel()
	{
	}


	double eval(const double* x1, const double* x2, unsigned int dim)
	{
		unsigned int i;
		double dot = 0.0;
		for (i=0; i<dim; i++) dot += x1[i] * x2[i];
		return pow(dot + parameter(0), parameter(1));
	}

	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		// TODO (later...)
		// Only fill in the derivative w.r.t. the offset?
		return 0.0;
	}
};


/*!
 *
 *  \brief Definition of the RBF Gaussian kernel
 *
 *  A special but very important type of kernel is the Gaussian
 *  normal distribution density kernel
 *  \f[
 *  exp(-\gamma \|x_1 - x_2\|^2)
 *  \f]
 *  It has a single parameter \f$\gamma > 0\f$ controlling the kernel
 *  width \f$\sigma = \sqrt{\gamma / 2}\f$.
 */
class RBFKernel : public KernelFunction
{
public:
	//! Constructor
	RBFKernel(double gamma)
	{
		parameter.resize(1, false);
		parameter(0) = gamma;
	}

	//! Destructor
	~RBFKernel()
	{
	}


	//! Computes the Gaussian kernel
	double eval(const double* x1, const double* x2, unsigned int dim)
	{
		double tmp, dist2 = 0.0;
		unsigned int i;
		for (i=0; i<dim; i++)
		{
			tmp = x1[i] - x2[i];
			dist2 += tmp * tmp;
		}
		return exp(-parameter(0) * dist2);
	}

	//! Computes the Gaussian kernel and the derivative w.r.t \f$\gamma\f$.
	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
	{
		derivative.resize(1, false);
		double tmp, dist2 = 0.0;
		unsigned int i;
		for (i=0; i<dim; i++)
		{
			tmp = x1[i] - x2[i];
			dist2 += tmp * tmp;
		}
		double ret = exp(-parameter(0) * dist2);
		derivative(0) = -dist2 * ret;
		return ret;
	}
};


//
// ********************
// ***  DECREPATED  ***
// ********************
//
// This class is not available any more
// because there is no systematic way to
// incorporate class specific penalties
// at the kernel level.
// Please use the norm2-parameter of the
// C_SVM class instead!
//
// /*!
//  *
//  *  \brief Kernel modification used for 2-norm C-SVM training
//  *
//  *  \par The Kernel2norm class modifies an existing kernel object
//  *  in the following way:
//  *  \f[
//  *       k(x, z) = k_{base}(x, z) + \frac{1}{C} \delta_{x, z}
//  *  \f]
//  *  where \f$ \delta \f$ is the Kronecker delta which is 1 if
//  *  x and z are equal and 0 otherwise.
//  *  This kernel function allows for the training of a 2-norm
//  *  C-SVM where C becomes a kernel parameter and the training
//  *  is carried out in the hard margin case, which means
//  *  \f$ C = \infty \f$.
//  *
//  *  \par The modified version of the kernel should only be used
//  *  during training. For prediction the original function should
//  *  be used. In practice, this will not make any difference, as
//  *  the two version only differ on finitely many points, that is
//  *  on the training set.
//  *
// */
// class Kernel2norm : public KernelFunction
// {
// public:
// //! Constructor
// 	Kernel2norm(KernelFunction* kernel, double C)
// 	{
// 		baseKernel = kernel;
// 
// 		int k, kc = kernel->getParameterDimension();
// 		parameter.resize(kc + 1, false);
// 		for (k=0; k<kc; k++) parameter(k) = kernel->getParameter(k);
// 		parameter(kc) = C;
// 		factor_plus = 1.0;
// 		factor_minus = 1.0;
// 	}
// 
// 	Kernel2norm(KernelFunction* kernel, double Cplus, double Cminus)
// 	{
// 		baseKernel = kernel;
// 		double C = sqrt(Cplus * Cminus);
// 
// 		int k, kc = kernel->getParameterDimension();
// 		parameter.resize(kc + 1, false);
// 		for (k=0; k<kc; k++) parameter(k) = kernel->getParameter(k);
// 		parameter(kc) = C;
// 		factor_plus = Cplus / C;
// 		factor_minus = Cminus / C;
// 	}
// 
// 	//! Destructor
// 	~Kernel2norm()
// 	{
// 	}
// 
// 
// 	//! Computes the 2norm kernel
// 	double eval(const double* x1, const double* x2, unsigned int dim)
// 	{
// 		int kc = baseKernel->getParameterDimension();
// 		double ret = baseKernel->eval(x1, x2, dim);
// 		unsigned int d;
// 		for (d=0; d<dim; d++) if (x1[d] != x2[d]) break;
// 		if (d == dim)
// 		{
// 			ret += 1.0 / parameter(kc);
// 		}
// 		return ret;
// 	}
// 
// 	//! Computes the 2norm kernel and its derivative
// 	double evalDerivative(const double* x1, const double* x2, unsigned int dim, Array<double>& derivative)
// 	{
// 		int k, kc = baseKernel->getParameterDimension();
// 		derivative.resize(kc + 1, false);
// 		Array<double> baseDerivative;
// 		double ret = baseKernel->evalDerivative(x1, x2, dim, baseDerivative);
// 		for (k=0; k<kc; k++) derivative(k) = baseDerivative(k);
// 		unsigned int d;
// 		for (d=0; d<dim; d++) if (x1[d] != x2[d]) break;
// 		if (d == dim)
// 		{
// 			double f = 1.0 / parameter(kc);
// 			ret += f;
// 			derivative(kc) = -f*f;
// 		}
// 		else derivative(kc) = 0.0;
// 		return ret;
// 	}
// 
// 	void setParameter(unsigned int index, double value)
// 	{
// 		if (index == baseKernel->getParameterDimension())
// 			parameter(index) = fabs(value);
// 		else
// 		{
// 			parameter(index) = value;
// 			baseKernel->setParameter(index, value);
// 		}
// 	}
// 
// protected:
// 	KernelFunction* baseKernel;
// 	double factor_plus;
// 	double factor_minus;
// };


#endif
