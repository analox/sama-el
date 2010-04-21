/*!
 *  \brief Implementation of the negative Kernel Target
 *  Alignment (KTA) as proposed by Nello Cristianini.
 *
 *  \par Kernel Target Alignment measures how good a kernel
 *  fits a binary classification training set. It is invariant
 *  under kernel rescaling. To turn it into an error function
 *  (i.e., to make minimization meaningful) the negative
 *  KTA is implemented.
 *
 *  \par The KTA has two important properties: It is
 *  differentiable and independent of the actual classifier.
 *
 *  \par The KTA as originally proposed by Nello Cristianini
 *  is not properly arranged for unbalanced datasets. Thus,
 *  the #negativeBKTA (balanced kernel target alignment) class
 *  implements an invariant version, which is preferable for
 *  unbalanced datasets.
 *
 *  \par The classes negativeKTA and negativeBKTA accept two
 *  kinds of models: KernelFunction derived classes and C_SVM
 *  with 2-norm slack penalty.
 *  In the first case, only the kernel parameters are considered
 *  in the derivative computation. In the second case,
 *  in addition the derivatives w.r.t. the parameters \f$ C_+ \f$
 *  and \f$ C_- \f$ are computed.
*/


#ifndef _KTA_H_
#define _KTA_H_


#include <ReClaM/ErrorFunction.h>
#include <ReClaM/KernelFunction.h>
#include <math.h>


//! Implementation of the negative Kernel Target Alignment (KTA)
//! as proposed by Nello Cristianini
class negativeKTA : public ErrorFunction
{
public:
	//! Constructor
	negativeKTA()
	{
	}

	//! Destructor
	~negativeKTA()
	{
	}


	//! Computes the negative Kernel Target Alignment between
	//! the target and the kernel function output on the input.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		// check the model type
		C_SVM* pCSVM = dynamic_cast<C_SVM*>(&model);
		KernelFunction* pKernel = dynamic_cast<KernelFunction*>(&model);
		bool norm2 = false;
		if (pCSVM != NULL)
		{
			if (! pCSVM->is2norm()) throw "[negativeKTA::error] C_SVM does not use 2-norm slack penalty.";
			pKernel = pCSVM->getSVM()->getKernel();
			norm2 = true;
		}
		else if (pKernel == NULL) throw "[negativeKTA::error] model is not a valid C_SVM or KernelFunction.";

		int i, j, T = input.dim(0);
		double k;
		double A = 0.0;
		double S = 0.0;
		for (i=0; i<T; i++)
		{
			for (j=0; j<i; j++)
			{
				k = pKernel->eval(input[i], input[j]);
				A += 2.0 * target(i) * target(j) * k;
				S += 2.0 * k * k;
			}
			k = pKernel->eval(input[i], input[i]);
			if (norm2) k += (target(i) > 0.0) ? 1.0 / pCSVM->get_Cplus() : 1.0 / pCSVM->get_Cminus();
			A += k;
			S += k * k;
		}
		return -A / (T * sqrt(S));
	}

	//! Computes the negative Kernel Target Alignment between
	//! the target and the kernel function output on the input.
	//! The partial derivatives of the negative KTA w.r.t. the
	//! model parameters are returned in the derivative parameter.
	double errorDerivative(Model& model, const Array<double>& input, const Array<double>& target, Array<double>& derivative)
	{
		// check the model type
		C_SVM* pCSVM = dynamic_cast<C_SVM*>(&model);
		KernelFunction* pKernel = dynamic_cast<KernelFunction*>(&model);
		bool norm2 = false;
		if (pCSVM != NULL)
		{
			if (! pCSVM->is2norm()) throw "[negativeKTA::errorDerivative] C_SVM does not use 2-norm slack penalty.";
			pKernel = pCSVM->getSVM()->getKernel();
			norm2 = true;
		}
		else if (pKernel == NULL) throw "[negativeKTA::errorDerivative] model is not a valid C_SVM or KernelFunction.";

		int i, j, T = input.dim(0);
		int d, D = model.getParameterDimension();
		double k;
		double A = 0.0;
		double S = 0.0;
		double yy;
		Array<double> a(D);
		Array<double> b(D);
		Array<double> der(D);
		derivative.resize(D, false);

		a = 0.0;
		b = 0.0;
		double aCplus = 0.0;
		double aCminus = 0.0;
		double bCplus = 0.0;
		double bCminus = 0.0;
		derivative = 0.0;

		for (i=0; i<T; i++)
		{
			for (j=0; j<i; j++)
			{
				k = pKernel->evalDerivative(input[i], input[j], der);
				yy = target(i) * target(j);
				for (d=0; d<D; d++)
				{
					a(d) += 2.0 * yy * der(d);
					b(d) += 2.0 * k * der(d);
				}
				A += 2.0 * yy * k;
				S += 2.0 * k * k;
			}
			k = pKernel->evalDerivative(input[i], input[i], der);
			if (norm2)
			{
				if (target(i) > 0.0)
				{
					double invC = 1.0 / pCSVM->get_Cplus();
					k += invC;
					aCplus -= invC * invC;
					bCplus -= k * invC * invC;
				}
				else
				{
					double invC = 1.0 / pCSVM->get_Cminus();
					k += invC;
					aCminus -= invC * invC;
					bCminus -= k * invC * invC;
				}
			}
			for (d=0; d<D; d++)
			{
				a(d) += der(d);
				b(d) += k * der(d);
			}
			A += k;
			S += k * k;
		}

		double N = T * sqrt(S);
		for (d=0; d<D; d++)
		{
			derivative(d + 1) = (A * b(d) / S - a(d)) / N;
		}
		if (norm2)
		{
			double cr = pCSVM->getCRatio();
			derivative(0) = (A * (bCplus + cr * bCminus) / S - (aCplus + cr * aCminus)) / N;
		}
		else derivative(0) = 0.0;

		return -A / N;
	}
};


//! \brief Balanced version of the #negativeKTA.
//!
//! \par The Balanced Kernel Target Alignment measure is a variant
//! of the Kernel Target Alignment.
//! This version of the measure is invariant under the fractions of
//! positive and negative examples.
//!
//! \par It is computed as follows:
//! Let p be the number of positive examples and let n be the number
//! of negative examples. We define an inner product between kernel
//! (Gram) matrices M and N as
//! \f[
//!     \langle M, N \rangle := \sum_{i, j} \lambda_{ij} M_{ij} N_{ij}
//! \f]
//! with positive coeffitients
//! \f$ \lambda_{ij} = n/p \f$ if \f$ y_i = y_j = +1 \f$,
//! \f$ \lambda_{ij} = 1 \f$ if \f$ y_i \not= y_j \f$ and
//! \f$ \lambda_{ij} = p/n \f$ if \f$ y_i = y_j = -1 \f$.
//!
//! Then we apply the usual definition of the kernel target alignment
//! \f[
//!     \hat A = \frac{\langle yy^T, K \rangle}{\|yy^T\| \cdot \|K\|}
//! \f]
//! where the norm \f$ \| M \| = \sqrt{\langle M, M \rangle} \f$ is
//! defined according to the inner product defined above.
//!
//! \par To my knowledge, the balancing modification has not been
//! proposed in the literature so far.
//!
//! \author: T. Glasmachers
//!
//! \date: 2006
//!
class negativeBKTA : public ErrorFunction
{
	public:
	//! Constructor
	negativeBKTA()
	{
	}

	//! Destructor
	~negativeBKTA()
	{
	}


	//! Computes the negative Balanced Kernel Target Alignment between
	//! the target and the kernel function output on the input.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		// check the model type
		C_SVM* pCSVM = dynamic_cast<C_SVM*>(&model);
		KernelFunction* pKernel = dynamic_cast<KernelFunction*>(&model);
		bool norm2 = false;
		if (pCSVM != NULL)
		{
			if (! pCSVM->is2norm()) throw "[negativeBKTA::error] C_SVM does not use 2-norm slack penalty.";
			pKernel = pCSVM->getSVM()->getKernel();
			norm2 = true;
		}
		else if (pKernel == NULL) throw "[negativeBKTA::error] model is not a valid C_SVM or KernelFunction.";

		int i, j, T = input.dim(0);

		int n = 0;
		int p = 0;
		double mult_pp, mult_nn;
		for (i=0; i<T; i++) if (target(i) > 0.0) p++; else n++;
		mult_pp = ((double)n) / ((double)p);
		mult_nn = ((double)p) / ((double)n);

		double k, temp;
		double A = 0.0;
		double S = 0.0;
		for (i=0; i<T; i++)
		{
			for (j=0; j<i; j++)
			{
				k = pKernel->eval(input[i], input[j]);
				if (target(i) > 0.0 && target(j) > 0.0)
				{
					temp = 2.0 * mult_pp * k;
					A += temp;
					S += temp * k;
				}
				else if (target(i) < 0.0 && target(j) < 0.0)
				{
					temp = 2.0 * mult_nn * k;
					A += temp;
					S += temp * k;
				}
				else
				{
					temp = 2.0 * k;
					A -= temp;
					S += temp * k;
				}
			}
			k = pKernel->eval(input[i], input[i]);
			if (norm2) k += (target(i) > 0.0) ? 1.0 / pCSVM->get_Cplus() : 1.0 / pCSVM->get_Cminus();
			if (target(i) > 0.0)
			{
				temp = mult_pp * k;
				A += temp;
				S += temp * k;
			}
			else
			{
				temp = mult_nn * k;
				A += temp;
				S += temp * k;
			}
		}
		return -0.5 * A / (sqrt(S * ((double)p) * ((double)n)));
	}

	//! Computes the negative Balanced Kernel Target Alignment between
	//! the target and the kernel function output on the input.
	//! The partial derivatives of the negative BKTA w.r.t. the
	//! model parameters are returned in the derivative parameter.
	double errorDerivative(Model& model, const Array<double>& input, const Array<double>& target, Array<double>& derivative)
	{
		// check the model type
		C_SVM* pCSVM = dynamic_cast<C_SVM*>(&model);
		KernelFunction* pKernel = dynamic_cast<KernelFunction*>(&model);
		bool norm2 = false;
		if (pCSVM != NULL)
		{
			if (! pCSVM->is2norm()) throw "[negativeBKTA::errorDerivative] C_SVM does not use 2-norm slack penalty.";
			pKernel = pCSVM->getSVM()->getKernel();
			norm2 = true;
		}
		else if (pKernel == NULL) throw "[negativeBKTA::errorDerivative] model is not a valid C_SVM or KernelFunction.";

		int i, j, T = input.dim(0);

		int n = 0;
		int p = 0;
		double mult, amult, mult_pp, mult_nn;
		for (i=0; i<T; i++) if (target(i) > 0.0) p++; else n++;
		mult_pp = ((double)n) / ((double)p);
		mult_nn = ((double)p) / ((double)n);

		int d, D = model.getParameterDimension();
		int kc = pKernel->getParameterDimension();
		double k;
		double A = 0.0;
		double S = 0.0;
		Array<double> a(kc);
		Array<double> b(kc);
		Array<double> der(kc);
		derivative.resize(D, false);

		a = 0.0;
		b = 0.0;
		double aCplus = 0.0;
		double aCminus = 0.0;
		double bCplus = 0.0;
		double bCminus = 0.0;
		derivative = 0.0;

		for (i=0; i<T; i++)
		{
			for (j=0; j<i; j++)
			{
				k = pKernel->evalDerivative(input[i], input[j], der);
				if (target(i) > 0.0 && target(j) > 0.0) mult = amult = 2.0 * mult_pp;
				else if (target(i) < 0.0 && target(j) < 0.0) mult = amult = 2.0 * mult_nn;
				else { mult = -2.0; amult = 2.0; }

				for (d=0; d<kc; d++)
				{
					a(d) += mult * der(d);
					b(d) += amult * k * der(d);
				}
				A += mult * k;
				S += amult * k * k;
			}
			k = pKernel->evalDerivative(input[i], input[i], der);
			if (norm2)
			{
				if (target(i) > 0.0)
				{
					double invC = 1.0 / pCSVM->get_Cplus();
					k += invC;
					aCplus -= mult_pp * invC * invC;
					bCplus -= mult_pp * k * invC * invC;
				}
				else
				{
					double invC = 1.0 / pCSVM->get_Cminus();
					k += invC;
					aCminus -= mult_nn * invC * invC;
					bCminus -= mult_nn * k * invC * invC;
				}
			}
			if (target(i) > 0.0) mult = mult_pp;
			else mult = mult_nn;
			for (d=0; d<kc; d++)
			{
				a(d) += mult * der(d);
				b(d) += mult * k * der(d);
			}
			A += mult * k;
			S += mult * k * k;
		}

		double N = 2.0 * sqrt(((double)n) * ((double)p) * S);
		if (norm2)
		{
			double cr = pCSVM->getCRatio();
			derivative(0) = (A * (bCplus + bCminus / cr) / S - (aCplus + aCminus / cr)) / N;
			for (d=0; d<kc; d++)
			{
				derivative(d + 1) = (A * b(d) / S - a(d)) / N;
			}
		}
		else
		{
			for (d=0; d<kc; d++)
			{
				derivative(d) = (A * b(d) / S - a(d)) / N;
			}
		}

		return -A / N;
	}
};


#endif
