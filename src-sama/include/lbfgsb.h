#ifndef __LBFGSB__
#define __LBFGSB__

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include <vector.h>
#include <Array/Array.h>


//
// function prototype for evaluator function of the local search
//
void LSEvaluate( vector<double> xStart, double& predicted );

//
// function prototype for local search
//
void findLocalOpt( vector<double> xx, vector<double>& newX, Array<double> boundTR );

#endif

