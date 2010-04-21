
/* This is the header for several transcendental mathematical functions.
   The source code was copied from
   Lau, H. T.  Numerical library in C for scientists and engineers
   1995 by CRC Press Inc.
   pgs 539 - 543, 570-576

*/

#ifndef GAMMA_H
#define GAMMA_H

#include <cmath>
#include "maps_general.h"
#include "Dyn_alloc.h"

#ifndef COMPSVD_F
#define COMPSVD_F

using namespace std;

const int MAXCHOLREPS = 1; // Maximum number of times to repeat the 
// cholesky iteration loop before using an SVD

//const char errmsg[40] = "Memory allocation failed.  Exiting";

bool SolveSystem( const Matrix<double>& A, Vector<double>& x, Vector<double>&
                  b);
   /***
       INPUT: Matrix A, vector b, and vector x through which the answer will be returned
       OUTPUT: true if clapack returns OK
       EFFECT: Calls CLAPACK's dgesv function, to solve the linear system Ax=b,
       and stores the solution in x.
   ***/

bool ComputeSVD(const pt_collect &A, pt_collect* U, pt_collect* D,
                  pt_collect* VT);
   /****
   INPUT: Matrix A, to take the SVD of, and pointers to three matrices that will
   contain what is returned by the SVD routine. 
   its main diagonal.
   OUTPUT: true if CLAPACK's dgesvd returns OK.
   EFFECT: Calls CLAPACK's dgesvd routine to compute the Singular Value
   Decomposition(SVD) of A.  Recall that a SVD of A yields U, D, and V',
   s.t. A=U*D* V', where U and V are orthogonal matrices, and D is a diagonal
   matrix with the singular values on its main diagonal.
  
   Borrowed from Chris Siefert's ParamEstimate Class 6/14/99
   ****/

bool ComputeCholesky( const pt_collect & A, pt_collect* U);
 /****
      INPUT: Matrix A, of which to take the Cholesky factorization and a pointer
      to the matrix in which to store U.  A = U`U
      OUTPUT: true if CLAPACK's dpotrf returns OK.
      EFFECT: Calls CLAPACK's dpotrf routine to compute the Cholesky factorization
      of A.  Computes A = U'U where U is upper triangular
 ****/

bool ComputeInverse( const pt_collect & A, pt_collect* INV);
 /****
      INPUT: Matrix A, of which to find the inverse using the Cholesky
      factorization and a pointer to the matrix in which to store the inverse INV.  A INV = I
      OUTPUT: true if CLAPACK's dpotri returns OK.
      EFFECT: Calls CLAPACK's dpotri routine to compute the Cholesky factorization
      of A, and find the inverse from it.
      based roughly on Chris Siefert's ComputeSVD code, modified by Anthony Padula
 ****/

#endif


double recipgamma(double x, double *odd, double *even);
/*  Delivers 1/Gamma(1-x);

    x:  double;
    entry:  this argument should satisfy -0.5 <= x <= 0.5;
       (actually the gamma function is calculated for 1-x, ie. if one wants to
       calculate 1/Gamma(1), one has to set x to 0);

    odd:  double *;
    exit:   the odd part of 1/Gamma(1-x) divided by (2x);
    ie. (1/Gamma(1-x) - 1/Gamma(1+x)) / (2x);

    even:  double *;
    exit:   the even part of 1/Gamma(1-x) divided by (2x);
    ie. (1/Gamma(1-x) + 1/Gamma(1+x)) / (2x);

*/

double gamma(double x);
/* Delivers the value of the gamma function at x.

   x: double;
   entry:  the argument; if one of the following three conditions is fulfilled
   the overflow will occur:

   1. the argument is too large;
   2. the argument is a non-positive integer;
   3. the argument is too close to a large ( in absolute value) non-positive
   integer.

   functions used: recipgamma, loggamma.
*/

double loggamma( double x);
/* delivers the value of the natural logarithm of the gamma function at x;

   x: double;
   entry:  this argument must be positive

*/

void bessk01( double x, double *k0, double *k1);
/* computes the value of the modified Bessel function of the third kind of order
   0 and 1, ie. K_0(x) and K_1(x)

   x: double;
   entry:  the argument of the Bessel functions; x>0;

   k0: double *;
   exit:   k0 has the value of the modified Bessel function of the third kind of
   order zero with argument x;

   k1: double *;
   exit:   k1 has the value of the modified Bessel function of the third kind of
   order one with argument x;
*/

void bessk( double x, int n, double k[]);
/* Generates an array of modified Bessel functions of the third kind of order j,
   K_j(x).

   Calls bessk01 to generate the order zero and order one terms, then computes
   higher order terms by recursion:

   K_{j+1}(x) = K_{j-1}(x)+(2*j/x)*K_{j}(x)


   x: double;
   entry:  the argument of the Bessel functions; x>0

   n: int;
   entry:  the upper bound of the indices of array k; n>=0;

   k: double k[0:n];
   exit:   k[j] is the value of the modified Bessel function of the third kind
   of order j with argument x, j=0,...,n.
*/

void nonexpbessk01( double x, double *k0, double * k1);
/* computes the value of the modified Bessel function of the third kind of order
   0 and 1, ie. K_0(x) and K_1(x), multiplied by e^x.

   x: double;
   entry:  the argument of the Bessel functions; x>0;

   k0: double *;
   exit:   k0 has the value of the modified Bessel function of the third kind of
   order zero with argument x, multiplied by e^x;

   k1: double *;
   exit:   k1 has the value of the modified Bessel function of the third kind of
   order one with argument x, multiplied by e^x;
*/
     
double MachineEpsilon();
// Computes and returns Machine Epsilon for this platform
#endif
