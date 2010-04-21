
#ifndef _svm_H_
#define _svm_H_


#include <ReClaM/Model.h>
#include <ReClaM/Optimizer.h>
#include <ReClaM/KernelFunction.h>
#include <ReClaM/C_Solver.h>


//!
//! \brief Support Vector Machine (SVM) as a ReClaM #Model
//!
//! \author T. Glasmachers
//! \date 2005
//!
//! \par The SVM class provides a Support Vector Machine as a
//! parametric #Model, that is, it computes a linear expansion
//! of the kernel function with fixed training examples as one
//! component. It can be viewed as a parametrized family of maps
//! from the input space to the reals.
//!
//! \par The parameter array of the SVM class defines the
//! affine linear solution in feature space, usually described
//! as a vector \f$ \alpha \f$ and a real valued offset b.
//! Note that different SVM training procedures impose
//! constraints on the possible values these parameters are
//! allowed to take.
//!
//! \par In ReClaM, the SVM as a model is used for prediction
//! only. That is, it does not impose any training scheme and
//! could in theory be training using any error measure and any
//! optimizer. In practice, one wants to apply standard SVM
//! training schemes to efficiently find the SVM solution.
//! Therefore, the SVM should be trained with the special
//! optimizer derived class #C_SVM. This class implements the
//! so-called C-SVM as a training scheme and uses a quadratic
//! program solver to obtain the optimal solution.
//!
class SVM : public Model
{
public:
	//! Constructor
	//!
	//! \param  pKernel      kernel function to use for training and prediction
	//! \param  bSignOutput  true if the SVM should output binary labels, false if it should output real values function evaluations
	SVM(KernelFunction* pKernel, bool bSignOutput = true);

	//! Constructor
	//!
	//! \param  pKernel      kernel function to use for training and prediction
	//! \param  input        training data points
	//! \param  target       training data labels
	//! \param  bSignOutput  true if the SVM should output binary labels, false if it should output real values function evaluations
	SVM(KernelFunction* pKernel, const Array<double>& input, const Array<double>& target, bool bSignOutput = true);

	//! Destructor
	~SVM();

	//! \par
	//! As the SVM can be constructed without training data
	//! points, although these are needed for the computation
	//! of the model, this member makes the training data
	//! known to the SVM. This method is usually called by
	//! the #SVM_Optimizer class, but in case a model was
	//! loaded from a file it can be necessary to invoke the
	//! function manually.
	//!
	//! \par
	//! A side effects, the method eventually rescales the
	//! parameter vector. In any case, all parameters are reset.
	//! The number of examples and the input space dimension
	//! are overwritten.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	//! \param  input        training data points
	//! \param  target       training data labels
	void setTrainingData(const Array<double>& input, const Array<double>& target);

	//! compute the SVM prediction on data
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	void model(const Array<double>& input, Array<double> &output);

	//! \par
	//! The modelDerivative member computes the derivative
	//! of the SVM function w.r.t. its parameters. Although
	//! this information is probably never used, it is easy
	//! to compute :)
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	void modelDerivative(const Array<double>& input, Array<double>& derivative);

	//! \par
	//! The modelDerivative member computes the derivative
	//! of the SVM function w.r.t. its parameters. Although
	//! this information is probably never used, it is easy
	//! to compute :)
	//! In addition, the model output is computed.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative);

	//! \par
	//! retrieve one of the coeffitients from the solution
	//! vector usually referred to as \f$ \alpha \f$
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	inline double get_alpha(int index)
	{
		return parameter(index);
	}

	//! \par
	//! retrieve the solution offset, usually referred to
	//! as \f$ b \f$
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	inline double get_b()
	{
		return parameter(examples);
	}

	//! return the kernel function object
	inline KernelFunction* getKernel()
	{
		return kernel;
	}

	//! return the training data points
	inline const Array<double>& get_points()
	{
		return *x;
	}

	//! return the training data labels
	inline const Array<double>& get_labels()
	{
		return *y;
	}

	//! return the number of training examples
	inline unsigned int get_examples()
	{
		return examples;
	}

	//! return the input space dimension
	inline unsigned int get_dimension()
	{
		return dimension;
	}

	//! \par
	//! Load the complete SVM model including
	//! kernel parameters, alpha, b and the
	//! support vectors.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	bool LoadSVMModel(std::istream& is);

	//! \par
	//! Save the complete SVM model including
	//! kernel parameters, alpha, b and the
	//! support vectors.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	bool SaveSVMModel(std::ostream& os);

protected:
	//! kernel function
	KernelFunction* kernel;

	//! true if x and y are allocated by LoadModel
	bool bOwnMemory;

	//! training data points
	const Array<double>* x;

	//! training data labels
	const Array<double>* y;

	//! true if the SVM outputs the binary label \f$ \pm 1 \f$ only,
	//! false if the SVM outputs the value of the linear feature space function.
	bool signOutput;

	//! number of training data points #x and labels #y
	unsigned int examples;

	//! input space dimension of the training data points #x
	unsigned int dimension;
};


//!
//! \brief Meta Model for SVM training
//!
//! \author T. Glasmachers
//! \date 2006
//!
//! \par
//! The C-SVM is a training scheme for Support Vector Machines.
//! It defines a positive constant C controlling the solution
//! complexity.
//!
//! \par
//! In the 1-norm SVM formulation, this constant
//! bounds the \f$ alpha_i \f$ parameters of the SVM from above,
//! limiting their contribution to the solution.
//! For the 2-norm SVM formulation, the \f$ alpha_i \f$ are not
//! bounded from above. Instead, the kernel matrix is reglarized
//! by adding \f$ 1 / C \f$ to the diagonal entries.
//!
//! \par
//! The C_SVM class is a meta model which is based on an #SVM
//! model. Its parameter vector is composed of the SVM hyper-
//! parameters, that is the complexity parameter C and the
//! parameters of a #KernelFunction objects.
//!
//! \par
//! For unbalanced data, it has proven fertile to consider class
//! specific complexity parameters \f$ C_+ \f$ and \f$ C_- \f$
//! instead of a single constant C. These parameters are usually
//! couples by the equation \f$ \ell_- C_+ = \ell_+ C_- \f$, where
//! $\ell_{\pm}$ are the class magnitudes. Therefore, the C_SVM
//! model introduces only the \f$ C_+ \f$ parameter as a model
//! parameter, whereas \f$ C_- = \ell_ / \ell_+ C_+ \f$ can be
//! computed from this parameter. As only the \f$ C_+ \f$ parameter
//! is subject to optimization, all derivatives have to be computed
//! w.r.t \f$ C_+ \f$, even if \f$ C_- \f$ is used for the
//! computation.
//!
//! \par
//! The parameter \f$ C_+ \f$ can be stored in a way allowing for
//! unconstrained optimization, that is, the exponential function
//! is used to compute the value from the model parameter. This
//! parameterization is even better suited for optimization.
//!
class C_SVM : public Model
{
public:
	//! Constructor
	//!
	//! \param  pSVM     Pointer to the SVM to be optimized.
	//! \param  Cplus    initial value of \f$ C_+ \f$
	//! \param  Cminus   initial value of \f$ C_- \f$
	//! \param  norm2    true if 2-norm slack penalty is to be used
	//! \param  unconst  true if the parameters are to be represented as log(C). This allows for unconstrained optimization.
	C_SVM(SVM* pSVM, double Cplus, double Cminus, bool norm2 = false, bool unconst = false);

	//! Descructor
	~C_SVM();


	//! return the unerlying #SVM model
	inline SVM* getSVM()
	{
		return svm;
	}

	//! Copy the parameters into the underlying #KernelFunction object
	//! and let the #SVM object compute the model prediction.
	void model(const Array<double>& input, Array<double>& output);

	//! overloaded version of Model::setParameter
	void setParameter(unsigned int index, double value);

	//! return the parameter C for positive class examples
	double get_Cplus();

	//! return the parameter C for negative class examples
	double get_Cminus();

	//! return true if the 2-norm slack penalty is in use
	inline bool is2norm()
	{
		return norm2penalty;
	}

	//! return the quotient \f$ C_- / C_+ \f$
	inline double getCRatio()
	{
		return C_ratio;
	}

protected:
	//! pointer to the underlying #SVM model
	SVM* svm;

	//! true if the 2-norm slack penalty is to be used
	bool norm2penalty;

	//! Complexity constant used for positive training examples
	double C_plus;

	//! Complexity constant used for negative training examples
	double C_minus;

	//! fraction \f$ C_- / C_+ \f$
	double C_ratio;

	//! if true the C-parameters are computed via \f$ C = exp(\tilde C) \f$ .
	bool exponential;
};


//!
//! \brief Optimizer for SVM training by quadratic programming
//!
//! \author T. Glasmachers
//! \date 2006
//!
//! Although the SVM_Optimizer fits into the ReClaM concept
//! of an #Optimizer, it should not be used the usual way,
//! for three reasons:
//! <ul>
//!   <li>It does not depend on the error function object
//!     provided as a parameter to #optimize.</li>
//!   <li>It does not return an error value. The return is zero
//!     in any case.</li>
//!   <li>It is not an iterative optimizer. The first call to
//!     #optimize already finds the optimal solution.</li>
//! </ul>
//! Nevertheless, the SVM_Optimizer is one of the most important
//! classes in the SVM context, as it calls the #C_Solver class
//! to solve the SVM dual problem, computing the #SVM
//! coefficients \f$ alpha \f$ and b.
//!
class SVM_Optimizer : public Optimizer
{
public:
	//! Constructor
	SVM_Optimizer();

	//! Destructor
	~SVM_Optimizer();


	//! \brief Default initialization
	//!
	//! \par
	//! This is the default #Model initialization method.
	//! If the model parameter is not a valid #C_SVM meta
	//! model, the method will throw an exception. This
	//! method is only provided for compatibility with the
	//! #Model interface. If possible the specialized init
	//! method should be used.
	//! 
	//! \param  model    Meta model containing complexity constant and kernel
	void init(Model& model);

	//! \brief Initialization of the SVM training optimizer.
	//!
	//! \par
	//! This initialization method is specifically designed
	//! for the needs of SVM training. It should be used instead
	//! of the default initialization method whenever possible.
	//!
	//! \param  model    The model parameter is passed to the default init method.
	//! \param  verbose  If true, the optimizer outputs status information to the standard output.
	//! \param  mb       Megabytes of main memory used as a cache for the kernel matrix.
	void init(C_SVM& model, bool verbose, unsigned int mb = 100);

	//! \brief Default #Optimizer interface
	//!
	//! \par
	//! The optimize member uses the C_Solver class to train the
	//! Support Vector Machine. If the model parameter supplied
	//! is not a valid #SVM object, the method checks for a
	//! #C_SVM object to determine the unserlying #SVM.
	//! Otherwise it throws an exception.
	//!
	//! \par
	//! Note that no #ErrorFunction object is needed for SVM
	//! training. If no ErrorFunction object is available,
	//! the static dummy reference #dummyError can be used.
	//!
	//! \param   model   The #SVM model to optimize
	//! \param   error   The error object is not used at all. The #dummyError can be used.
	//! \param   input   Training input used for the optimization
	//! \param   target  Training labels used for the optimization
	//! \return  As there is no error evaluation, the function returns 0 in any case.
	double optimize(Model& model, ErrorFunction& error, const Array<double>& input, const Array<double>& target);

	//! \brief Trains the #SVM with the given dataset.
	//!
	//! \par
	//! This version of the optimize member is well suited for
	//! #SVM training, but it does not fit into the ReClaM concept.
	//! Thus, a less specialized version of optimize is available,
	//! which overrides the corresponding #Optimizer member.
	//!
	//! \param   model   The #SVM model to optimize
	//! \param   input   Training input used for the optimization
	//! \param   target  Training labels used for the optimization
	void optimize(SVM& model, const Array<double>& input, const Array<double>& target);

	//! The static member dummyError is provided for calls to the
	//! #optimize method in case there is no #ErrorFunction object
	//! available. Do not use or call any members of dummyError,
	//! as it is an invalid reference. It is only provided for
	//! making syntactically correct calls to the #optimize method.
	static ErrorFunction& dummyError;

	//! Return the #C_Solver used by the previous call to #optimize.
	//! The #C_Solver object may be of use for the calling procedure
	//! as it proviedes some solution statistics as the number of
	//! iterations used for training and the final objective
	//! function value.
	inline C_Solver* get_Solver()
	{
		return solver;
	}

protected:
	//! quadratic program solver used in the last #optimize call
	C_Solver* solver;

	//! upper bound for the variables corresponding to positive class examples
	double Cplus;

	//! upper bound for the variables corresponding to negative class examples
	double Cminus;

	//! true if the 2-norm slack penalty is to be used
	bool norm2penalty;

	//! should the quadratic program solver output its progress?
	bool printInfo;

	//! kernel cache size in megabytes used by the quadratic program solver
	unsigned int cacheMB;
};


#endif
