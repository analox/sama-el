
#include <math.h>
#include <ReClaM/C_Solver.h>

#include <iostream>
#include <iomanip>


using namespace std;


////////////////////////////////////////////////////////////////////////////////


C_Solver::C_Solver(KernelFunction* pKernel,
				   const Array<double>& data,
				   const Array<double>& diagModification,
				   const Array<double>& linearPart,
				   const Array<double>& equalityNormal,
				   const Array<double>& inequalityUpper,
				   Array<double>& solutionVector,
				   double& solutionOffset,
				   bool verbose,
				   int cacheMB,
				   double eps,
				   const char* strategy,
				   bool shrinking,
				   bool norm2)
{
	examples = data.dim(0);
	dimension = data.dim(1);

	kernel = pKernel;
	diagMod = diagModification;
	linear = linearPart;

	eqNormal = equalityNormal;
	box = inequalityUpper;

	alpha = solutionVector;

	printInfo = verbose;
	epsilon = eps;
	WSS_Strategy = strategy;

	norm2penalty = norm2;

	unsigned int a, i, j;
	float* qi;
	float* qj;

	double lowerBound = -1e100;
	double upperBound = 1e100;
	double sum = 0.0;
	int freeVars = 0;
	double value;

	// prepare the cache
	cacheMaxSize = 1048576 * cacheMB / sizeof(float);
	cacheSize = 0;
	if (cacheMaxSize < 2 * examples) throw "[C_Solver::C_Solver] invalid cache size";
	cacheTemp = (float*)malloc(examples * sizeof(float));
	if (cacheTemp == NULL) throw "[C_Solver::C_Solver] out of memory error";
	cacheNewest = -1;
	cacheOldest = -1;
	cacheEntry.resize(examples);
	for (i=0; i<examples; i++)
	{
		cacheEntry[i].data = NULL;
		cacheEntry[i].length = 0;
		cacheEntry[i].older = -2;
		cacheEntry[i].newer = -2;
	}

	// prepare the solver internal variables
	active = examples;
	x.resize(examples, false);
	diagonal.resize(examples, false);
	permutation.resize(examples, false);
	gradient.resize(examples, false);
	for (i=0; i<examples; i++)
	{
		if (box(i) <= 0.0) throw "[C_Solver::C_Solver] The feasible region is empty.";
		x(i) = data[i].elemvec();
		if (norm2)
			diagonal(i) = diagMod(i) * diagMod(i) * (kernel->eval(x(i), x(i), dimension) + 1.0 / box(i));
		else
			diagonal(i) = diagMod(i) * diagMod(i) * kernel->eval(x(i), x(i), dimension);
		permutation(i) = i;
	}
	for (i=0; i<examples; i++)
	{
		value = linear(i);
		if (alpha(i) != 0.0)
		{
			qi = Q_row(i, 0, examples);
			for (a=0; a<examples; a++) value -= qi[a] * alpha(a);
		}
		gradient(i) = value;
	}

	bFirst = true;
	bUnshrinked = false;

	unsigned int shrinkCounter = (active < 1000) ? active : 1000;

	SelectWSS();

	// decomposition loop
	if (printInfo) cout << "{" << flush;
	iter = 0;
	while (true)
	{
		// select a working set and check for optimality
		if (SelectWorkingSet(i, j))
		{
			if (printInfo) cout << "+" << flush;
			// the solution may be optimal - check again
			if (SelectWorkingSet(i, j))
			{
				// seems to be optimal
				if (printInfo) cout << "*" << flush;

				if (! shrinking) break;

				// do costly unshrinking
				Unshrink();
				shrinkCounter = 1;

				// check again on the whole problem
				if (SelectWorkingSet(i, j)) break;
			}
		}

		double ee = eqNormal(i) * eqNormal(j);
		double ai = alpha(i);
		double aj = alpha(j);
		double Ci = box(i);
		double Cj = box(j);

		// get the rows of Q corresponding to the working set
		qi = Q_row(i, 0, active);
		qj = Q_row(j, 0, active);

		// update alpha, that is, solve the sub-problem defined by i and j
		double nominator = gradient(i) - ee * gradient(j);
		double denominator = diagonal(i) + diagonal(j) - 2.0 * ee * qi[j];
		double mu = nominator / denominator;
		if (mu < -ai) mu = -ai;
		else if (! norm2penalty && mu > Ci - ai) mu = Ci - ai;
		if (ee * mu > aj) mu = ee * aj;
		else if (! norm2penalty && ee * mu < aj - Cj) mu = ee * (aj - Cj);
		alpha(i) += mu;
		alpha(j) -= ee * mu;

		// update the gradient
		for (a=0; a<active; a++) gradient(a) -= mu * (qi[a] - ee * qj[a]);

		shrinkCounter--;
		if (shrinkCounter == 0)
		{
			// shrink the problem
			if (shrinking) Shrink();

			shrinkCounter = (active < 1000) ? active : 1000;
			// shrinkCounter = active;
		}

		iter++;
		if (printInfo)
		{
			if ((iter & 1023) == 0) cout << "." << flush;
		}
	}
	if (printInfo) cout << endl << "} #iterations=" << iter << endl;

	Unshrink();

	// compute the offset b
	for (a=0; a<active; a++)
	{
		value = eqNormal(a) * gradient(a);

		if (eqNormal(a) > 0.0)
		{
			if (alpha(a) == 0.0)
			{
				if (value < upperBound) upperBound = value;
			}
			else if (! norm2penalty && alpha(a) == box(a))
			{
				if (value > lowerBound) lowerBound = value;
			}
			else
			{
				sum += value;
				freeVars++;
			}
		}
		else
		{
			if (alpha(a) == 0.0)
			{
				if (value > lowerBound) lowerBound = value;
			}
			else if (! norm2penalty && alpha(a) == box(a))
			{
				if (value < upperBound) upperBound = value;
			}
			else
			{
				sum += value;
				freeVars++;
			}
		}

		if (freeVars > 0)
			solutionOffset = sum / freeVars;						// stabilized exact value
		else
			solutionOffset = 0.5 * (upperBound + lowerBound);		// best estimate
	}

	// return alpha
	for (i=0; i<examples; i++)
	{
		solutionVector(permutation(i)) = alpha(i);
	}
}

C_Solver::~C_Solver()
{
	cacheClear();
	free(cacheTemp);
}


double C_Solver::getObjectiveValue()
{
	unsigned int i;
	double objective;

	// Compute the dual objective value
	Unshrink(true);
	objective = 0.0;
	for (i=0; i<examples; i++)
	{
		objective += (1.0 + gradient(i)) * alpha(i);
	}
	objective *= 0.5;
	return objective;
}

bool C_Solver::MVP(unsigned int& i, unsigned int& j)
{
	double violation_i = -1e100;
	double violation_j = -1e100;
	unsigned int a;

	for (a=0; a<active; a++)
	{
		if (eqNormal(a) > 0.0)
		{
			if (norm2penalty || alpha(a) < box(a))
			{
				if (gradient(a) > violation_i)
				{
					violation_i = gradient(a);
					i = a;
				}
			}
			if (alpha(a) > 0.0)
			{
				if (-gradient(a) > violation_j)
				{
					violation_j = -gradient(a);
					j = a;
				}
			}
		}
		else
		{
			if (norm2penalty || alpha(a) < box(a))
			{
				if (gradient(a) > violation_j)
				{
					violation_j = gradient(a);
					j = a;
				}
			}
			if (alpha(a) > 0.0)
			{
				if (-gradient(a) > violation_i)
				{
					violation_i = -gradient(a);
					i = a;
				}
			}
		}
	}

	// MVP stopping condition
	return (violation_i + violation_j < epsilon);
}

bool C_Solver::HMG(unsigned int& i, unsigned int& j)
{
	if (bFirst)
	{
		bFirst = false;
		i = 0;
		for (j=1; j<active; j++) if (eqNormal(j) != eqNormal(i)) return false;
		return true;
	}

	// check the corner condition
	if (norm2penalty)
	{
		if ((alpha(old_i) <= 1e-8) && (alpha(old_j) <= 1e-8))
		{
			if (printInfo) cout << "^" << flush;
			return MVP(i, j);
		}
	}
	else
	{
		double C_i = box(old_i);
		double C_j = box(old_j);
		double eps_i = 1e-8 * C_i;
		double eps_j = 1e-8 * C_j;
		if ((alpha(old_i) <= eps_i || alpha(old_i) >= C_i - eps_i)
				&& ((alpha(old_j) <= eps_j || alpha(old_j) >= C_j - eps_j)))
		{
			if (printInfo) cout << "^" << flush;
			return MVP(i, j);
		}
	}

	// generic situation: use the MG selection
	unsigned int a;
	double aa, ab;					// alpha values
	double da, db;					// diagonal entries of Q
	double ga, gb;					// gradient in coordinates a and b
	double gain;
	double Ca, Cb;
	double denominator;
	float* q;
	double mu_max, mu_star;

	double best = 0.0;
	double mu_best = 0.0;

	// try combinations with b = old_i
	q = Q_row(old_i, 0, active);
	ab = alpha(old_i);
	db = diagonal(old_i);
	Cb = box(old_i);
	gb = gradient(old_i);
	for (a=0; a<active; a++)
	{
		if (a == old_i || a == old_j) continue;

		aa = alpha(a);
		da = diagonal(a);
		Ca = box(a);
		ga = gradient(a);

		if (eqNormal(a) * eqNormal(old_i) > 0.0)
		{
			denominator = (da + db - 2.0*q[a]);
			mu_max = (ga - gb) / denominator;
			mu_star = mu_max;

			if (mu_star < -aa) mu_star = -aa;
			else if (! norm2penalty && mu_star > Ca - aa) mu_star = Ca - aa;
			if (mu_star > ab) mu_star = ab;
			else if (! norm2penalty && mu_star < ab - Cb) mu_star = ab - Cb;

			gain = mu_star * (2.0 * mu_max - mu_star) * denominator;
		}
		else
		{
			denominator = (da + db + 2.0*q[a]);
			mu_max = (ga + gb) / denominator;
			mu_star = mu_max;

			if (mu_star < -aa) mu_star = -aa;
			else if (! norm2penalty && mu_star > Ca - aa) mu_star = Ca - aa;
			if (mu_star < -ab) mu_star = -ab;
			else if (! norm2penalty && mu_star > Cb - ab) mu_star = Cb - ab;

			gain = mu_star * (2.0 * mu_max - mu_star) * denominator;
		}

		// select the largest gain
		if (gain > best)
		{
			best = gain;
			mu_best = mu_star;
			i = a;
			j = old_i;
		}
	}

	// try combinations with old_j
	q = Q_row(old_j, 0, active);
	ab = alpha(old_j);
	db = diagonal(old_j);
	Cb = box(old_j);
	gb = gradient(old_j);
	for (a=0; a<active; a++)
	{
		if (a == old_i || a == old_j) continue;

		aa = alpha(a);
		da = diagonal(a);
		Ca = box(a);
		ga = gradient(a);

		if (eqNormal(a) * eqNormal(old_j) > 0.0)
		{
			denominator = (da + db - 2.0*q[a]);
			mu_max = (ga - gb) / denominator;
			mu_star = mu_max;

			if (mu_star < -aa) mu_star = -aa;
			else if (! norm2penalty && mu_star > Ca - aa) mu_star = Ca - aa;
			if (mu_star > ab) mu_star = ab;
			else if (! norm2penalty && mu_star < ab - Cb) mu_star = ab - Cb;

			gain = mu_star * (2.0 * mu_max - mu_star) * denominator;
		}
		else
		{
			denominator = (da + db + 2.0*q[a]);
			mu_max = (ga + gb) / denominator;
			mu_star = mu_max;

			if (mu_star < -aa) mu_star = -aa;
			else if (! norm2penalty && mu_star > Ca - aa) mu_star = Ca - aa;
			if (mu_star < -ab) mu_star = -ab;
			else if (! norm2penalty && mu_star > Cb - ab) mu_star = Cb - ab;

			gain = mu_star * (2.0 * mu_max - mu_star) * denominator;
		}

		// select the largest gain
		if (gain > best)
		{
			best = gain;
			mu_best = mu_star;
			i = a;
			j = old_j;
		}
	}

	// stopping condition
	return (fabs(mu_best) < epsilon);
}

bool C_Solver::Libsvm28(unsigned int& i, unsigned int& j)
{
	i=0;
	j=1;

	double violation_i = -1e100;
	double violation_j = -1e100;
	unsigned int a;

	// find the first index of the MVP
	for (a=0; a<active; a++)
	{
		if (eqNormal(a) > 0.0)
		{
			if (norm2penalty || alpha(a) < box(a))
			{
				if (gradient(a) > violation_i)
				{
					violation_i = gradient(a);
					i = a;
				}
			}
		}
		else
		{
			if (alpha(a) > 0.0)
			{
				if (-gradient(a) > violation_i)
				{
					violation_i = -gradient(a);
					i = a;
				}
			}
		}
	}

	// find the second index using second order information
	float* q = Q_row(i, 0, active);
	double best = 0.0;
	for (a=0; a<active; a++)
	{
		if (eqNormal(a) > 0.0)
		{
			if (alpha(a) > 0.0)
			{
				if (-gradient(a) > violation_j) violation_j = -gradient(a);

				double grad_diff = violation_i - gradient(a);
				if (grad_diff > 0)
				{
					double quad_coef = diagonal(i) + diagonal(a) - 2.0 * eqNormal(i) * q[a];
					double obj_diff = (grad_diff * grad_diff) / quad_coef;

					if (obj_diff > best)
					{
						best = obj_diff;
						j = a;
					}
				}
			}
		}
		else
		{
			if (norm2penalty || alpha(a) < box(a))
			{
				if (gradient(a) > violation_j) violation_j = gradient(a);

				double grad_diff = violation_i + gradient(a);
				if (grad_diff > 0)
				{
					double quad_coef = diagonal(i) + diagonal(a) + 2.0 * eqNormal(i) * q[a];
					double obj_diff = (grad_diff * grad_diff) / quad_coef;

					if (obj_diff > best)
					{
						best = obj_diff;
						j = a;
					}
				}
			}
		}
	}

	// MVP stopping condition
	return (violation_i + violation_j < epsilon);
}

bool C_Solver::SelectWorkingSet(unsigned int& i, unsigned int& j)
{
	// dynamic working set selection call
	bool ret = (this->*(this->currentWSS))(i, j);

	// bool ret = Libsvm28(i, j);
	// bool ret = HMG(i, j);
	old_i = i;
	old_j = j;
	return ret;
}

void C_Solver::SelectWSS()
{
	if (strcmp(WSS_Strategy, "MVP") == 0)
	{
		// most violating pair, used e.g. in LIBSVM 2.71
		currentWSS = &C_Solver::MVP;
	}
	else if (strcmp(WSS_Strategy, "HMG") == 0)
	{
		// hybrid maximum gain, suitable for large problems
		currentWSS = &C_Solver::HMG;
	}
	else if (strcmp(WSS_Strategy, "LIBSVM28") == 0)
	{
		// LIBSVM 2.8 second order algorithm
		currentWSS = &C_Solver::Libsvm28;
	}
	else
	{
		// default strategy:
		// use HMG as long as the problem does not fit into the cache,
		// use the LIBSVM 2.8 algorithm afterwards
		if (active * active > cacheMaxSize)
			currentWSS = &C_Solver::HMG;
		else
			currentWSS = &C_Solver::Libsvm28;
	}
}

void C_Solver::Shrink()
{
	std::vector<unsigned int> shrinked;
	unsigned int a;

	double g1 = 1e100;
	double g2 = 1e100;
	double v, g;

	// find largest KKT violations
	for (a=0; a<active; a++)
	{
		g = gradient(a);
		v = alpha(a);
		if (eqNormal(a) > 0.0)
		{
			if (v > 0.0 && g < g1) g1 = g;
			if ((norm2penalty || v < box(a)) && -g < g2) g2 = -g;
		}
		else
		{
			if (v > 0.0 && g < g2) g2 = g;
			if ((norm2penalty || v < box(a)) && -g < g1) g1 = -g;
		}
	}

	if (! bUnshrinked && -g1 - g2 < 10.0 * epsilon)
	{
		// unshrink the problem at this accuracy level
		if (printInfo) cout << "#" << flush;
		Unshrink();
		bUnshrinked = true;
		SelectWSS();
		return;
	}

	// identify the variables to shrink
	for (a=0; a<active; a++)
	{
		if (a == old_i) continue;
		if (a == old_j) continue;
		g = gradient(a);
		v = alpha(a);
		if (eqNormal(a) > 0.0)
		{
			if (v == 0.0)
			{
				if (g > g1) continue;
			}
			else if (! norm2penalty && v == box(a))
			{
				if (g < -g2) continue;
			}
			else continue;
		}
		else
		{
			if (v == 0.0)
			{
				if (g > g2) continue;
			}
			else if (! norm2penalty && v == box(a))
			{
				if (g < -g1) continue;
			}
			else continue;
		}

		// In this moment no feasible step including this variable
		// can improve the objective. Thus deactivate the variable.
		shrinked.push_back(a);
		if (cacheEntry[a].length > 0) cacheDelete(a);
	}

	int s, sc = shrinked.size();
	if (sc == 0)
	{
		return;
	}
	unsigned int new_active = active - sc;

	// exchange variables such that shrinked variables
	// are moved to the ends of the lists.
	unsigned int k, high = active;
	for (s=sc-1; s>=0; s--)
	{
		k = shrinked[s];
		high--;

		// exchange the variables "k" and "high"
		Exchange(k, high);
	}

	// shrink the cache entries
	for (a=0; a<examples; a++)
	{
		if (cacheEntry[a].length > (int)new_active) cacheResize(a, new_active);
	}

	active = new_active;

	SelectWSS();
}

void C_Solver::Unshrink(bool complete)
{
	if (printInfo) cout << "[" << flush;
	if (active == examples)
	{
		if (printInfo) cout << "]" << flush;
		return;
	}

	unsigned int i, a;
	float* q;
	double v, g;
	double g1 = 1e100;
	double g2 = 1e100;

	// compute the inactive gradient components (quadratic time complexity)
	unsigned int ac = active;
	active = examples;
	for (a=ac; a<examples; a++) gradient(a) = 1.0;
	for (i=0; i<examples; i++)
	{
		v = alpha(i);
		if (v == 0.0) continue;

		q = Q_row(i, ac, examples, true);
		for (a=ac; a<examples; a++) gradient(a) -= q[a] * v;
	}

	if (complete)
	{
		active = examples;
		return;
	}

	// find largest KKT violations
	for (a=0; a<active; a++)
	{
		g = gradient(a);
		v = alpha(a);
		if (eqNormal(a) > 0.0)
		{
			if (v > 0.0 && g < g1) g1 = g;
			if ((norm2penalty || v < box(a)) && -g < g2) g2 = -g;
		}
		else
		{
			if (v > 0.0 && g < g2) g2 = g;
			if ((norm2penalty || v < box(a)) && -g < g1) g1 = -g;
		}
	}

	// identify the variables to activate
	for (a=active; a<examples; a++)
	{
		if (a == old_i) continue;
		if (a == old_j) continue;
		g = gradient(a);
		v = alpha(a);
		if (eqNormal(a) > 0.0)
		{
			if (v == 0.0)
			{
				if (g <= g1) continue;
			}
			else if (! norm2penalty && v == box(a))
			{
				if (g >= -g2) continue;
			}
		}
		else
		{
			if (v == 0.0)
			{
				if (g <= g2) continue;
			}
			else if (! norm2penalty && v == box(a))
			{
				if (g >= -g1) continue;
			}
		}
		Exchange(active, a);
		active++;
	}

	if (printInfo) cout << active << "]" << flush;
}

#define XCHG_A(t, a, i, j) {t temp; temp = a(i); a(i) = a(j); a(j) = temp;}
#define XCHG_V(t, a, i, j) {t temp; temp = a[i]; a[i] = a[j]; a[j] = temp;}
// (i <= j) is required!
void C_Solver::Exchange(unsigned int i, unsigned int j)
{
	if (i == j) return;

	int t;

	// check the previous working set
	if (old_i == i) old_i = j;
	else if (old_i == j) old_i = i;

	if (old_j == i) old_j = j;
	else if (old_j == j) old_j = i;

	// exchange entries in the simple lists
	XCHG_A(const double*, x, i, j);
	XCHG_A(double, diagMod, i, j);
	XCHG_A(double, eqNormal, i, j);
	XCHG_A(double, box, i, j);
	XCHG_A(double, linear, i, j);
	XCHG_A(double, alpha, i, j);
	XCHG_A(unsigned int, permutation, i, j);
	XCHG_A(double, diagonal, i, j);
	XCHG_A(double, gradient, i, j);

	// update the ordered cache list predecessors and successors
	t = cacheEntry[i].older;
	if (t != -2)
	{
		if (t == -1) cacheOldest = j;
		else cacheEntry[t].newer = j;
		t = cacheEntry[i].newer;
		if (t == -1) cacheNewest = j;
		else cacheEntry[t].older = j;
	}
	t = cacheEntry[j].older;
	if (cacheEntry[j].older != -2)
	{
		if (t == -1) cacheOldest = i;
		else cacheEntry[t].newer = i;
		t = cacheEntry[j].newer;
		if (t == -1) cacheNewest = i;
		else cacheEntry[t].older = i;
	}

	// exchange the cache entries
	XCHG_V(tCacheEntry, cacheEntry, i, j);

	// exchange all cache row entries
	unsigned int k, l;
	for (k=0; k<examples; k++)
	{
		l = cacheEntry[k].length;
		if (j < l)
		{
			XCHG_V(float, cacheEntry[k].data, i, j);
		}
		else if (i < l)
		{
			// only one element is available from the cache
			if (norm2penalty && i == k)
				cacheEntry[k].data[k] = diagMod(k) * diagMod(k) * (kernel->eval(x(k), x(k), dimension) + 1.0 / box(k));
			else
				cacheEntry[k].data[i] = diagMod(i) * diagMod(k) * kernel->eval(x(i), x(k), dimension);
		}
	}
}

void C_Solver::cacheAppend(int var)
{
	if (cacheNewest == -1)
	{
		cacheNewest = var;
		cacheOldest = var;
		cacheEntry[var].older = -1;
		cacheEntry[var].newer = -1;
	}
	else
	{
		cacheEntry[cacheNewest].newer = var;
		cacheEntry[var].older = cacheNewest;
		cacheEntry[var].newer = -1;
		cacheNewest = var;
	}
}

void C_Solver::cacheRemove(int var)
{
	if (cacheEntry[var].older == -1)
		cacheOldest = cacheEntry[var].newer;
	else
		cacheEntry[cacheEntry[var].older].newer = cacheEntry[var].newer;

	if (cacheEntry[var].newer == -1)
		cacheNewest = cacheEntry[var].older;
	else
		cacheEntry[cacheEntry[var].newer].older = cacheEntry[var].older;

	cacheEntry[var].older = -2;
	cacheEntry[var].newer = -2;
}

void C_Solver::cacheAdd(int var, unsigned int length)
{
	cacheEntry[var].length = length;
	cacheEntry[var].data = (float*)(void*)malloc(length * sizeof(float));
	if (cacheEntry[var].data == NULL) throw "[C_Solver::cacheAppend] out of memory error";
	cacheSize += length;

	cacheAppend(var);
}

void C_Solver::cacheDelete(int var)
{
	free(cacheEntry[var].data);
	cacheSize -= cacheEntry[var].length;

	cacheEntry[var].data = NULL;
	cacheEntry[var].length = 0;

	cacheRemove(var);
}

void C_Solver::cacheResize(int var, unsigned int newlength)
{
	int diff = newlength - cacheEntry[var].length;
	if (diff == 0) return;
	cacheSize += diff;
	cacheEntry[var].length = newlength;
	cacheEntry[var].data = (float*)(void*)realloc((void*)cacheEntry[var].data, newlength * sizeof(float));
	if (cacheEntry[var].data == NULL) throw "[C_Solver::cacheResize] out of memory error";
}

void C_Solver::cacheClear()
{
	unsigned int e, ec = cacheEntry.size();
	for (e=0; e<ec; e++)
	{
		if (cacheEntry[e].data != NULL) free(cacheEntry[e].data);
		cacheEntry[e].data = NULL;
		cacheEntry[e].length = 0;
		cacheEntry[e].older = -1;
		cacheEntry[e].newer = -1;
	}
	cacheOldest = -1;
	cacheNewest = -1;
	cacheSize = 0;
}

float* C_Solver::Q_row(unsigned int k, int begin, int end, bool temp)
{
	if (temp)
	{
		// return temporary data
		if (cacheEntry[k].length > begin)
		{
			memcpy(cacheTemp + begin, cacheEntry[k].data + begin, cacheEntry[k].length - begin);
			begin = cacheEntry[k].length;
		}
		const double* x_k = x(k);
		double d_k = diagMod(k);
		int a;
		for (a=begin; a<end; a++) cacheTemp[a] = d_k * diagMod(a) * kernel->eval(x_k, x(a), dimension);
		if (norm2penalty && (int)k >= begin && (int)k < end)
			cacheTemp[k] = d_k * d_k * (kernel->eval(x_k, x_k, dimension) + 1.0 / box(k));
		return cacheTemp;
	}
	else
	{
		// the data will be stored in the cache
		int l = cacheEntry[k].length;
		while (cacheSize + end - l > cacheMaxSize)
		{
			if (cacheOldest == (int)k)
			{
				cacheRemove(k);
				cacheAppend(k);
			}
			cacheDelete(cacheOldest);
		}
		if (l == 0)
		{
			cacheAdd(k, end);
		}
		else
		{
			cacheResize(k, end);
			if ((int)k != cacheNewest)
			{
				cacheRemove(k);
				cacheAppend(k);
			}
		}

		// compute remaining entries
		if (l < end)
		{
			const double* x_k = x(k);
			double d_k = diagMod(k);
			int a;
			float* p = cacheEntry[k].data + l;
			for (a=l; a<end; a++)
			{
				*p = d_k * diagMod(a) * kernel->eval(x_k, x(a), dimension);
				p++;
			}
			if (norm2penalty && (int)k >= l && (int)k < end)
				cacheEntry[k].data[k] = d_k * d_k * (kernel->eval(x_k, x_k, dimension) + 1.0 / box(k));
		}

		return cacheEntry[k].data;
	}
}
