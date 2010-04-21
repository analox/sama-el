
#ifndef _optRprop_H_
#define _optRprop_H_


#include <ReClaM/Optimizer.h>


//===========================================================================
/*!
 *  \brief This class offers methods for the usage of the improved
 *         Resilient-Backpropagation-algorithm with weight-backtracking.
 *
 *  The Rprop algorithm is an improvement of the algorithms with adaptive
 *  learning rates (as the Adaptive Backpropagation algorithm by Silva
 *  and Ameida, please see AdpBP.h for a description of the 
 *  working of such an algorithm), that uses increments for the update
 *  of the model parameters, that are independant from the absolute partial
 *  derivatives. This makes sense, because large flat regions
 *  in the search space (plateaus) lead to small absolute partial
 *  derivatives and so the increments are chosen small, but the increments
 *  should be large to skip the plateau. In contrast, the absolute partial
 *  derivatives are very large at the "slopes" of very "narrow canyons",
 *  which leads to large increments that will skip the minimum lying
 *  at the bottom of the canyon, but it would make more sense to
 *  chose small increments to hit the minimum. <br>
 *  So, the Rprop algorithm only uses the signs of the partial derivatives
 *  and not the absolute values to adapt the parameters. <br>
 *  Instead of individual learning rates, it uses the parameter
 *  \f$\Delta_i^{(t)}\f$ for model parameter \f$w_i,\ i = 1, \dots, n\f$ in 
 *  iteration \f$t\f$, where the parameter will be adapted before the 
 *  change of the model parameters: <br>
 *
 *  \f$
 *  \Delta_i^{(t)} = \Bigg\{ 
 *  \begin{array}{ll}
 *  min( \eta^+ \cdot \Delta_i^{(t-1)}, \Delta_{max} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} > 0 \\
 *  max( \eta^- \cdot \Delta_i^{(t-1)}, \Delta_{min} ), & \mbox{if\ } 
 *  \frac{\partial E^{(t-1)}}{\partial w_i} \cdot
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \\
 *  \Delta_i^{(t-1)}, & \mbox{otherwise}
 *  \end{array}
 *  \f$ 
 *
 *  The parameters \f$\eta^+ > 1\f$ and \f$0 < \eta^- < 1\f$ control
 *  the speed of the adaptation. To stabilize the increments, they are
 *  restricted to the interval \f$[\Delta_{min}, \Delta_{max}]\f$. <br>
 *  After the adaptation of the \f$\Delta_i\f$ the update for the 
 *  model parameters will be calculated as
 *
 *  \f$
 *  \Delta w_i^{(t)} := - \mbox{sign} 
 *  \left( \frac{\partial E^{(t)}}{\partial w_i}\right) \cdot \Delta_i^{(t)}
 *  \f$
 *        
 *  Furthermore, weight-backtracking will take place to increase the
 *  stability of the method. In contrast to the original Rprop algorithm
 *  with weight-backtracking (see RpropPlus) this weight-backtracking 
 *  is improved by additionally taken the error of the last iteration
 *  \f$t - 1\f$ into account. <br>
 *  The idea of this modification is, that a change of the sign of the
 *  partial derivation \f$\frac{\partial E}{\partial w_i}\f$
 *  only states, that a minimum was skipped and not, whether this step
 *  lead to an approach to the minimum or not. <br>
 *  By using the old error value the improved weight-backtracking only
 *  undoes changes, when the error has increased and only the parameters
 *  \f$w_i\f$ are reset to the old values, where a sign change of
 *  \f$\frac{\partial E}{\partial w_i}\f$ has taken place. <br>
 *  So the new weight-backtracking rule is: <br>
 *
 *  \f$
 *  \mbox{if\ } \frac{\partial E^{(t-1)}}{\partial w_i} \cdot 
 *  \frac{\partial E^{(t)}}{\partial w_i} < 0 \mbox{\ then} \{
 *  \f$
 *
 *  \f$  
 *  \begin{array}{lll}
 *   \Delta w_i^{(t)} = \bigg\{ &
 *   - \Delta w_i^{(t-1)}, & \mbox{if\ } E^{(t)} > E^{(t - 1)} \\
 *   & 0, & otherwise \\
 *  \frac{\partial E^{(t)}}{\partial w_i} := 0 
 *  \end{array}
 *  \f$
 *
 *  \f$\}\f$
 *
 *  , where the assignment of zero to the partial derivative of the error
 *  leads to a freezing of the increment in the next iteration. <br>
 *
 *  This modification of the weight backtracking leads to a better
 *  optimization on artifical, paraboloidal error surfaces. <br>
 *
 *  For further information about the algorithm, please refer to: <br>
 *  
 *  Christian Igel and Michael H&uuml;sken, <br>
 *  "Empirical Evaluation of the Improved Rprop Learning Algorithm". <br>
 *  In Neurocomputing Journal, 2002, in press <br> 
 *
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *
 */
class IRpropPlus : public Optimizer
{
	public:

	void init(Model& model)
	{
		initRprop(model);

		int i, ic = model.getParameterDimension();
		for (i=0; i<ic; i++) delta(i) = 0.01 * model.getParameter(i);
	}

//===========================================================================
/*!
 *  \brief Prepares the Rprop algorithm for the currently used model.
 *
 *  Internal variables of the class instance are initialized and memory
 *  for used structures adapted to the used model. An initial value
 *  for the parameter \f$\Delta\f$ is assigned to all parameters of the model.
 *
 *  \param mode    Model to be optimized
 *  \param _delta0 Initial value for the parameter \f$\Delta\f$. 
 *                 The default value is "0.01".
 *  \param np      The increase factor \f$\eta^+\f$, by default set to "1.2".
 *  \param nm      The decrease factor \f$\eta^-\f$, by default set to "0.5".
 *  \param dMax    The upper limit of the increments \f$\Delta w_i^{(t)}\f$.
 *                 The default value is "50".
 *  \param dMin    The lower limit of the increments \f$\Delta w_i^{(t)}\f$.
 *                 The default value is "1e-6".
 *  \return none
 *
 *  \author  C. Igel, M. H&uuml;sken
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
void initRprop(Model& model, double _delta0 = 0.01, double np = 1.2, double nm = 0.5, double dMax = 50.0, double dMin = 1e-6)
{
	deltaw.resize(model.getParameterDimension(), false);
	delta.resize(model.getParameterDimension(), false);
	dedwOld.resize(model.getParameterDimension(), false);
	delta  = _delta0; 
	deltaw = 0; 
	dedwOld  = 0; 
#ifndef __SOLARIS__
	oldE = std::numeric_limits< double >::max( );
#else
	oldE = 1e309;
#endif
	increaseFactor = np;
	decreaseFactor = nm;
	upperStepLimit = dMax;
	lowerStepLimit = dMin;
}


//===========================================================================
/*!
 *  \brief Performs a run of the Rprop algorithm.
 *
 *  The error \f$E^{(t)}\f$ of the used model for the current iteration
 *  \f$t\f$ is calculated and the values of the model parameters \f$w_i,\ 
 *  i = 1, \dots, n\f$ and the parameters \f$\Delta_i\f$ are 
 *  adapted depending on this error.
 *
 *  \param errorfunction   The errorfunction to minimize.
 *  \param model           The model to evaluate the error function on.
 *  \param input           The input patterns used for the training of the model.
 *  \param target          The target values for to the input patterns.
 *  \return none
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *
 */  
double optimize(Model& model, ErrorFunction& errorfunction, const Array<double>& input, const Array<double>& target)
{
	int p, pc = model.getParameterDimension();

	double param;
	Array<double> dedw(pc);

	double newE = errorfunction.errorDerivative(model, input, target, dedw);

	for (p=0; p<pc; p++)
	{
		param = model.getParameter(p);

		if (dedw(p) * dedwOld(p) > 0.0)
		{
			delta(p) = min(upperStepLimit, increaseFactor * delta(p));
			deltaw(p) = delta(p) * -sgn(dedw(p));
			model.setParameter(p, param + deltaw(p));
			dedwOld(p) = dedw(p);
		}
		else if (dedw(p) * dedwOld(p) < 0.0)
		{
			delta(p) = max(lowerStepLimit, decreaseFactor * delta(p));
			if(oldE < newE)
			{
				model.setParameter(p, param - deltaw(p));
			}
			dedwOld(p) = 0.0;
		}
		else
		{
			deltaw(p) = delta(p) * -sgn(dedw(p));
			model.setParameter(p, param + deltaw(p));
			dedwOld(p) = dedw(p);
		}

		if (! model.isFeasible())
		{
			model.setParameter(p, param);
			delta(p) = max(lowerStepLimit, decreaseFactor * delta(p));
		}
	}

	oldE = newE;
	return newE; 
}


protected:
//===========================================================================
/*!
 *  \brief Determines the sign of "x".
 *
 *  \param x The value of which the sign shall be determined.
 *  \return "-1", if the sign of \em x is negative, "0" otherwise.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
	int sgn(double x) 
	{ 
		if (x > 0) return 1; if (x < 0) return -1; return 0; 
	}

	//! Used to replace the standard definition.
	double min(double x, double y) { return x < y ? x : y; }

	//! Used to replace the standard definition.
	double max(double x, double y) { return x > y ? x : y; }

	//! The error of the last iteration.
	double oldE;

	//! The final update values for all model parameters.
	Array<double> deltaw;     

	//! The last error gradient.
	Array<double> dedwOld;

	//! The absolute update values (increment) for all model parameters.
	Array<double> delta;

	//! The increase factor \f$\eta^+\f$, by default set to "1.2".
	double increaseFactor;

	//! The decrease factor \f$\eta^-\f$, by default set to "0.5".
	double decreaseFactor;

	//! The upper limit of the increments \f$\Delta w_i^{(t)}\f$. The default value is "50".
	double upperStepLimit;

	//! The lower limit of the increments \f$\Delta w_i^{(t)}\f$. The default value is "1e-6".
	double lowerStepLimit;
};


#endif
