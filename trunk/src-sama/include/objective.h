//objective.h
//this example file declares the functions and variable n as required
//by the pattern search implementations
//simple definitions for certain parameters may be given here as well.

#if !defined _userfile_
#define _userfile_

#include "vec.h"

using namespace std;

//double getInitialStep();

void fcn(long, Vector<double> &, double &, bool &, void *);

void initpt(long, double *&);

extern long n;
extern double tolerance;
extern bool newInitialPoint;
extern bool joke;   // just to see what happens---pls  No prob.  Links fine--remove

/* Initial trial step length, reflecting the refinement of the search lattice
 */
//#define initialStepLength getInitialStep()

/*  The search will halt if the trial step length is less than 
 *  or equal to the stoppingStepLength after a reduction of the 
 *  trial step length.
 */ 
//#define stoppingStepLength tolerance

/*  The maxCalls variable is not a hard limit on the number of function
 *  evaluations that may be performed by the search.  The search should
 *  always halt at a point from which none of the trial steps can find 
 *  improvement to be consistent with the formal theory regarding pattern
 *  searches.
 */
//#define maxCalls 1000000

#endif






