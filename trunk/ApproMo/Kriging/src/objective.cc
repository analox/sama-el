/* objective.cc
 * simple quadratic equation in four dimensions with rigged starting point
 */

#include "objective.h" 
#include "math.h"
#include <cstdio>
#include <iostream>

long n=4;

/* I defined my stoppingStepLength to be tolerance, which I must define here.
 * Note that the value chosen is machine/compiler dependent as approximately
 *   the square root of machine epsilon
 */
double tolerance = 10e-8; //sqrt(fabs(3.0 * (4.0/3.0 -1.0) -1.0));

/* Initial trial step length, reflecting the refinement of the search lattice
 */
//double initialStep = .25;

/* Returns initial step length for the search.
 * This function is not technically necessary, as long as initialStepLength 
 * is set to a real number in the header file.
 */
//double getInitialStep()
//{
//  return initialStep;
//}

/* The user's function to be minimized.
 * vars is the number of variables in the problem.
 * double * x is a pointer to an array of doubles,
 *   which represent the point to be evaluated
 * double & f will be the value of the function
 *   at that point upon return
 * flag will be 1 if the function can be successfully
 *   evaluated, 0 otherwise
 * Note that pattern searches have enough flexibility
 *   that they need not necessarily be able to get the
 *   value of all possible points
 */
void fcn(long vars, Vector<double> &x, double & f, bool & flag, void* nothing)
{
  if(vars!=n) {flag=false; return;}
  double tone = 0.0;
  for(long i = 0; i < vars; i++)
    {
      tone += x.begin()[i] * x.begin()[i];
    }
  f = tone / vars;
  flag = true;
  return;
}

/* initpt returns the user's initial guess
 * vars is the number of variables of the problem
 * double * x is a pointer (NULL),
 *   which will point to a newly allocated array
 *   of doubles that represent the initial point
 */
void initpt(long vars, double *&x)
{
  double grace = 1;
  x = new double[vars];
  for(long i = 0; i < vars; i++)
    {
      x[i] = grace;
    }

}








