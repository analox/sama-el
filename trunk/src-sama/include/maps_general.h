/*6/1/99*/

#ifndef _MAPS_GENERAL_
#define _MAPS_GENERAL_

#include "vec.h"  /*Chris Siefert's modified version of R. Pozo's TNT libraries*/
#include "cppmat.h" /*Chris Siefert's modified version of R. Pozo's TNT libraries*/
#include "rngs.h" /*Steve Park's Random Number Generation Libraries*/

#define round(x) ((long) (x + 0.5))
#define sqr(x) ((x)*(x))
#define TOLERANCE 10e-14
#define ERROR -1

/*status flag stuff. */
#define NOTINIT 0
#define READY 1
#define MAPS_DONE 2

#ifndef DEBUG
#define DEBUG 1
#endif

using namespace std;

/*
  debuglevels set at compile-time
  DEBUG = 0 - dead silent
  DEBUG = 1 - soft - errors only
  DEBUG = 2 - verbose - function calls and returns, mostly
  DEBUG = 3 - quite loud
  DEBUG = 4 - as loud as a rock concert.
*/


typedef Vector<double> mp_vector;
typedef Vector<long> ml_vector;
typedef Matrix<double> pt_collect;

struct Bounds {
  /*public:*/
  mp_vector lower;
  mp_vector upper;
};

#endif




