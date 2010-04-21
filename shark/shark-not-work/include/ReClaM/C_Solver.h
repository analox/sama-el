
#ifndef _C_Solver_H_
#define _C_Solver_H_


#include <ReClaM/Model.h>
#include <ReClaM/KernelFunction.h>
#include <Array/Array.h>


//!
//! \brief Quadratic Program Solver for the C-SVM.
//!
//! \author T. Glasmachers
//! \date 2005
//!
//! \par
//! The training procedure iteratively solves the C-SVM dual problem:
//! \f[
//!   maximize f(\alpha) = l^T \alpha - \frac{1}{2} \alpha^T diag(d) K diag(d) \alpha
//! \f]
//! s.t. \f$ \sum_{i=1}^\ell e_i \alpha_i = 0 \f$
//! and \f$ 0 \leq \alpha_i \leq u_i \,\forall\, i \f$
//!
//! \par
//! The matrix K is the positive definite Gram matrix of the kernel function
//! evaluated on the data. The matrix diag(d) K diag(d) is refered to as Q.
//! The linear part l is an arbitrary vector. The vector d is usually nonzero
//! in every component. The equality constraint normal e may only consist of
//! entries +1 or -1, while the box constraint upper edge vector u has to be
//! non-negative in every component.
//!
//! \par
//! The solver uses a SMO-like decomposition algorithm with MVP, HMG or
//! LIBSVM-2.8 working set selection algorithm, all providing guaranteed
//! convergence to the optimum.
//! The HMG algorithm is known to be fast especially on large problems,
//! while the LIBSVM-2.8 algorithm is fastest if the matrix Q fits into the
//! kernel cache. Therefore, the solver's default working set selection
//! strategy is a hybrid mixture of these two algorithms.
//!
class C_Solver
{
public:
	//! The constructor solves the quadratic program given as follows:
	//!
	//! \par
	//! maximize \f$ f(\alpha) = l^T \alpha - \frac{1}{2} \alpha^T D K D \alpha \f$ <br>
	//! such that \f$ e^T \alpha = const \f$ <br>
	//! and \f$ 0 \leq \alpha_i \leq u_i \,\forall i \f$ <br>
	//!
	//! \par
	//! If the norm2 flag is set, two changes are made to the
	//! quadratic problem in order to solve the 2-norm SVM formulation:
	//! The upper box constraint is removed and its value is instead
	//! misused to regularize the kernel matrix K by 
	//! \f[
	//!     K \leftarrow K + diag(1 / u_i) \f] \enspace.
	//! \f]
	//!
	//! \param pKernel          kernel function used for the computation of K
	//! \param data             data to evaluate the kernel function
	//! \param linearPart       linear part l of the target function
	//! \param diagModificator  vector; diagonal of the matrix D
	//! \param equalityNormal   vector e
	//! \param inequalityUpper  vector u of upper bounds
	//! \param solutionVector   input: initial feasible vector \f$ \alpha \f$; output: solution \f$ \alpha^* \f$
	//! \param solutionOffset   output: affine offset b of the C-SVM solution
	//! \param verbose          turn status output to stdout on or off (off by default)
	//! \param cacheMB          megabytes of main memory to use as a kernel cache (100 by default)
	//! \param eps              solution accuracy, in terms of the maximum KKT violation (0.001 by default)
	//! \param strategy         working set selection strategy to use; values: "" (default), "MVP", "HMG", "LIBSVM28"
	//! \param shrinking        turn shrinking on or off (off by default)
	//! \param norm2            use 2-norm slack penalties; ignore inequalityUpper
	C_Solver(KernelFunction* pKernel,
			const Array<double>& data,
			const Array<double>& linearPart,
			const Array<double>& diagModification,
			const Array<double>& equalityNormal,
			const Array<double>& inequalityUpper,
			Array<double>& solutionVector,
			double& solutionOffset,
			bool verbose = false,
			int cacheMB = 100,
			double eps = 0.001,
			const char* strategy = "",
			bool shrinking = true,
			bool norm2 = false);

	//! Descructor
	virtual ~C_Solver();


	//! Return the number of iterations used to solve the problem
	inline unsigned int iterations()
	{
		return iter;
	}

	//! Compute the objective value \f$ f(\alpha^*) \f$,
	//! see the description of the constructor.
	double getObjectiveValue();

protected:
	//! Kernel function defining the kernel Gram matrix
	KernelFunction* kernel;

	//! Array of data vectors for kernel evaluations
	Array<const double*> x;

	//! This vector modifies the quadratic part of the target function
	Array<double> diagMod;

	//! linear part of the target function
	Array<double> linear;

	//! Normal vector of the equality constraint
	Array<double> eqNormal;

	//! Box constraint upper ends. The lower ends are always 0.
	Array<double> box;

	//! Solution candidate
	Array<double> alpha;

	//! Number of variables
	unsigned int examples;

	//! Input space dimension
	unsigned int dimension;

	//! Number of currently active variables
	unsigned int active;

	//! permutation of the variables alpha, gradient, etc.
	Array<unsigned int> permutation;

	//! diagonal entries of Q
	//! The diagonal array is of fixed size and not subject to shrinking.
	Array<double> diagonal;

	//! gradient of the objective function
	//! The gradient array is of fixed size and not subject to shrinking.
	Array<double> gradient;

	//! indicator of the first decomposition iteration
	bool bFirst;

	//! first component of the previous working set
	unsigned int old_i;

	//! second component of the previous working set
	unsigned int old_j;

	//! stopping condition - solution accuracy
	double epsilon;

	//! This flag indicated that 2-norm slack penalty is to be used.
	bool norm2penalty;

	//! \brief Select the most violatig pair (MVP)
	//!
	//! \return true if the solution is already sufficiently optimal
	//!  \param i  first working set component
	//!  \param j  second working set component
	bool MVP(unsigned int& i, unsigned int& j);

	//! \brief Select a working set according to the hybrid maximum gain (HMG) algorithm
	//!
	//! \return true if the solution is already sufficiently optimal
	//!  \param i  first working set component
	//!  \param j  second working set component
	bool HMG(unsigned int& i, unsigned int& j);

	//! \brief Select a working set according to the second order algorithm of LIBSVM 2.8
	//!
	//! \return true if the solution is already sufficiently optimal
	//!  \param i  first working set component
	//!  \param j  second working set component
	bool Libsvm28(unsigned int& i, unsigned int& j);

	//! \brief Select a working set
	//!
	//! \par
	//! This member is implemented as a wrapper to HMG.
	//! \return true if the solution is already sufficiently optimal
	//!  \param i  first working set component
	//!  \param j  second working set component
	virtual bool SelectWorkingSet(unsigned int& i, unsigned int& j);

	//! Choose a suitable working set algorithm
	void SelectWSS();

	//! Shrink the problem
	void Shrink();

	//! Active all variables
	void Unshrink(bool complete = false);

	//! true if the problem has already been unshrinked
	bool bUnshrinked;

	//! exchange two variables via the permutation
	void Exchange(unsigned int i, unsigned int j);

	//! Kernel cache
	unsigned int cacheSize;						// current cache size in floats
	unsigned int cacheMaxSize;					// maximum cache size in floats
	struct tCacheEntry							// data held for every cache entry
	{
		float* data;							// float array containing a matrix row
		int length;								// length of this matrix row
		int older;								// next older entry
		int newer;								// next newer entry
	};
	std::vector<tCacheEntry> cacheEntry;		// cache entry description
	float* cacheTemp;							// single kernel row
	int cacheNewest;							// index of the newest entry
	int cacheOldest;							// index of the oldest entry

	//! append the entry to the ordered list
	void cacheAppend(int var);

	//! remove the entry from the ordered list
	void cacheRemove(int var);

	//! add an entry to the cache and append
	//! it to the ordered list
	void cacheAdd(int var, unsigned int length);

	//! remove an entry from the cache and the ordered list
	void cacheDelete(int var);

	//! resize a cache entry
	void cacheResize(int var, unsigned int newlength);

	//! completely clear the cache
	void cacheClear();

	//! \brief Return the active subset of a kernel matrix row
	//!
	//! \par
	//! This method returns an array of float with at least
	//! the entries in the interval [begin, end[ filled in.
	//! If #temp is set to true, the computed values are not
	//! stored in the cache.
	//! \param k      index identifying the matrix row
	//! \param begin  first entry to be filled in
	//! \param end    last entry to be filled in +1
	//! \param temp   are the return values temporary or should they be cached?
	float* Q_row(unsigned int k, int begin, int end, bool temp = false);

	//! solution statistics: number of iterations
	unsigned int iter;

	//! should the solver print its status to the standard output?
	bool printInfo;

	//! working set selection strategy to follow
	const char* WSS_Strategy;

	//! pointer to the currently used working set selection algorithm
	bool (C_Solver::*currentWSS)(unsigned int&, unsigned int&);
};


#endif
