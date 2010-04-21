
#ifndef _RadiusMargin_H_
#define _RadiusMargin_H_


#include <ReClaM/ErrorFunction.h>
#include <ReClaM/svm.h>
#include <ReClaM/C_Solver.h>
#include <math.h>
#include <fstream>


using namespace std;


//!
//! \brief Squared Radius-Margin-Quotient
//!
//! \author T. Glasmachers
//! \date 2006
//!
//! \par
//! The Radius-Margin-Quotient \f$ R^/\gamma^2 \f$ is a well
//! known quantity in statistical learning theory which can be
//! turned into a generalization bound after normalization.
//! Its minimization has been proposed for SVM model selection.
//!
//! \par
//! The Radius-Margin-Quotient depends on the kernel parameters
//! as well as the SVM regularization parameter C. However, as
//! the derivative of the \f$ \alpha \f$ vector w.r.t. C is hard
//! to compute (in fact, it involves the inverse of the matrix
//! \f$ Q_{ij} = y_i y_j k(x_i, x_j) \f$ which is considered
//! to be too large to fit into memory), the #errorDerivative
//! member only computes the derivative w.r.t the kernel
//! parameters. However, if the 2-norm slack penalty formuation
//! is used, C can be interpreted as a kernel matrix parameter
//! and the derivative can be computed easily.
//!
class RadiusMargin : public ErrorFunction
{
public:
	//! Constructor
	RadiusMargin()
	{
	}

	//! Destructor
	~RadiusMargin()
	{
	}


	//! Computes the radius margin quotient.
	double error(Model& model, const Array<double>& input, const Array<double>& target)
	{
		SVM* pSVM = NULL;
		double Cplus, Cminus;
		bool norm2;

		// check the model type
		C_SVM* csvm = dynamic_cast<C_SVM*>(&model);
		if (csvm == NULL)
			throw "[RadiusMargin::error] model is not a valid C_SVM.";

		Array<double> tmp;
		csvm->model(input[0], tmp);			// copy the parameters into the kernel
		Cplus = csvm->get_Cplus();
		Cminus = csvm->get_Cminus();
		norm2 = csvm->is2norm();
		pSVM = csvm->getSVM();

		// compute the coefficients alpha and beta
		Array<double> alpha;
		Array<double> beta;
		solveProblems(pSVM, Cplus, Cminus, input, target, alpha, beta, norm2);

		KernelFunction* kernel = pSVM->getKernel();
		unsigned int i, j, examples = input.dim(0);
		double R2, w2;
		double kv;
		double bi, bj;

		w2 = 0.0;
		R2 = 0.0;
		for (j=0; j<examples; j++)
		{
			bj = beta(j);
			for (i=0; i<examples; i++)
			{
				bi = beta(i);
				kv = kernel->eval(input[i], input[j]);
				if (i == j && norm2)
				{
					if (target(i) > 0.0) kv += 1.0 / csvm->get_Cplus();
					else kv += 1.0 / csvm->get_Cminus();
				}

				R2 -= bi * bj * kv;
				if (i == j)
				{
					w2 += alpha(i);
					R2 += bi * kv;
				}
			}
		}
		return R2 * w2;
	}

	//! Computes the radius margin quotient and its derivatives
	//! as proposed by Chapelle, Vapnik, Bousquet and Mukherjee in
	//! "Choosing Multiple Parameters for Support Vector Machines",
	//! sections 6.1 and 6.2. Please note that in the 1-norm slack
	//! penalty case the derivative w.r.t. C is returned as zero,
	//! which is actually not the case.
	double errorDerivative(Model& model, const Array<double>& input, const Array<double>& target, Array<double>& derivative)
	{
		SVM* pSVM = NULL;
		double Cplus, Cminus;
		bool norm2;

		// check the model type
		C_SVM* csvm = dynamic_cast<C_SVM*>(&model);
		if (csvm == NULL)
			throw "[RadiusMargin::errorDerivative] model is not a valid C_SVM or C_SVM2.";

		Array<double> tmp;
		csvm->model(input[0], tmp);			// copy the parameters into the kernel
		Cplus = csvm->get_Cplus();
		Cminus = csvm->get_Cminus();
		norm2 = csvm->is2norm();
		pSVM = csvm->getSVM();

		// compute the coefficients alpha and beta
		Array<double> alpha;
		Array<double> beta;
		solveProblems(pSVM, Cplus, Cminus, input, target, alpha, beta, norm2);

		KernelFunction* kernel = pSVM->getKernel();
		unsigned int i, j, examples = input.dim(0);
		unsigned int k, kc = kernel->getParameterDimension();
		double R2, w2;
		double kv;
		double ayi, ayj, bi, bj;
		Array<double> dw2(kc);
		Array<double> dR2(kc);
		Array<double> dkv(kc);
		derivative.resize(kc + 1, false);

		w2 = 0.0;
		dw2 = 0.0;
		R2 = 0.0;
		dR2 = 0.0;
		double dR2dCplus = 0.0;
		double dR2dCminus = 0.0;
		double dw2dCplus = 0.0;
		double dw2dCminus = 0.0;
		for (j=0; j<examples; j++)
		{
			ayj = alpha(j) * target(j);
			bj = beta(j);
			for (i=0; i<examples; i++)
			{
				ayi = alpha(i) * target(i);
				bi = beta(i);
				kv = kernel->evalDerivative(input[i], input[j], dkv);
				if (i == j)
				{
					if (norm2)
					{
						double invC;
						if (target(i) > 0.0)
						{
							invC = 1.0 / csvm->get_Cplus();
							kv += invC;
							dR2dCplus -= bi * invC * invC;
							dw2dCplus += (ayi * ayi) * invC * invC;
						}
						else
						{
							invC = 1.0 / csvm->get_Cminus();
							kv += invC;
							dR2dCminus -= (bi - bi * bi) * invC * invC;
							dw2dCminus += (ayi * ayi) * invC * invC;
						}
					}
					w2 += alpha(i);
					R2 += bi * kv;
					for (k=0; k<kc; k++) dR2(k) += bi * dkv(k);
				}
				R2 -= bi * bj * kv;
				for (k=0; k<kc; k++)
				{
					dw2(k) -= ayi * ayj * dkv(k);
					dR2(k) -= bi * bj * dkv(k);
				}
			}
		}

		if (norm2)
		{
			double cr = csvm->getCRatio();
			derivative(0) = R2 * (dw2dCplus + cr * dw2dCminus) + w2 * (dR2dCplus + cr * dR2dCminus);
		}
		else
		{
			derivative(0) = 0.0;
		}
		for (k=0; k<kc; k++) derivative(k + 1) = R2 * dw2(k) + w2 * dR2(k);

		return R2 * w2;
	}

protected:
	//! Helper function calling the quadratic problem solver
	//! to obtain the coefficients \alpha and \beta.
	void solveProblems(SVM* pSVM, double Cplus, double Cminus, const Array<double>& input, const Array<double>& target, Array<double>& alpha, Array<double>& beta, bool norm2)
	{
		unsigned int i, examples = input.dim(0);
		KernelFunction* kernel = pSVM->getKernel();
		Array<double> ones(examples);
		Array<double> inf(examples);
		Array<double> C(examples);
		Array<double> diag(examples);
		double tmp;
		C_Solver* pSolver;

		alpha.resize(examples, false);
		beta.resize(examples, false);
		alpha = 0.0;
		beta = 1.0 / examples;
		ones = 1.0;
		inf = 1e100;
		for (i=0; i<examples; i++)
		{
			C(i) = (target(i) > 0.0) ? Cplus : Cminus;
			diag(i) = kernel->eval(input[i], input[i]);
		}

		// solve the first problem to obtain the alpha coefficients
		pSolver = new C_Solver(kernel, input, target, ones, target, C, alpha, tmp, false, 100, 0.001, "", true, norm2);
		delete pSolver;

		// solve the second problem to obtain the beta coefficients
		if (norm2)
			pSolver = new C_Solver(kernel, input, diag, ones, ones, C, beta, tmp, false, 100, 0.001, "", true, norm2);
		else
			pSolver = new C_Solver(kernel, input, diag, ones, ones, inf, beta, tmp, false, 100, 0.001, "", true, norm2);
		delete pSolver;
	}
};


#endif
