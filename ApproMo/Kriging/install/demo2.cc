/* This is the demo main program for the Krigifier and the Kriging approximation classes.

 We instruct the krigifier to do Latin Hypercube sampling.
 It creates a random function, using the parameters specified on the function call. 
 It chooses to use a user-defined trend, described in TrendFunction().  It then 
 evaluates the function at several sites, and uses these values to
 create an approximation, choosing to use the Gaussian isotropic correlation function 
 and passing in the known value of theta.  
 The approximation is then evaluated systematically and the results are output to the file 
 "results2.out" while simple statistics on the data are output to the screen

 Differences from "demo.cc":
 1. krigifier parameters passed through function call, rather than being read from a file
 2. krigifier given specific rng seed and stream to use
 3. krigifier uses user-defined trend 10+3 ||x||^{2} + ||x||^{5}
 4. approximation given theta, instead of doing estimation
 5. call analyzevalues() instead of scanfunction() in order to get more information
 6. Latin Hypercubes are used instead of Uniform Random sampling in the krigifier 

 Anthony D. Padula    adpadu@wm.edu
 5/31/00
*/

/* INSTRUCTIONS FOR USE
 To compile this, you need to have downloaded the files :

 krigify.cc, krigify.h, demo.cc rngs.c, rngs.h, rvgs.c, rvgs.h,
 krig.cc, krig.h, gamma.cc, gamma.h, ParamEstimate.cc, ParamEstimate.h,
 PatternSearch.cc, PatternSearch.h, CompassSearch.cc, CompassSearch.h,
 approx.cc, approx.h, Makefile

 You will also need the libraries CLAPACK, BLAS, and F2C, which can be obtained
 from http://www.netlib.org/
 You should modify the macros at the top of the Makefile to suit your compiler and libraries

 Then, on your command line you run

% make demo2
% ./demo2

*/

/* WHAT TO DO WITH THE OUTPUT 
Once you have run this program, you can plot the approximation using gnuplot.  Supposing
your CWD contains results.out, in gnuplot you should run the commands:

  set parametric
  splot "results2.out" with lines

This will produce a nice 3-d graph of the approximation  

Alternatively, you can compile and run plotpoints by

% g++ plotpoints.cc -o plotpoints
% ./plotpoints results2.out results2.ps

This produces a graph in postscript format, which can then be viewed or printed.
*/


#include "krig.h"
#include "krigify.h"
#include "rngs.h"
#include <cmath>

using namespace std;

/* Fills X with uniform random points between the bounds.  Thus
   lower[j] <= X[i][j] <= upper[j]   for all 0<=j<p
*/
void randPoints( Matrix<double> & X, Vector<double>& lower, Vector<double> & upper);


/* We will use this trend instead of the default quadratic */
void TrendFunction( const Vector<double>& x, double& result)
{
  double nx;
  Vector<double> temp;
  temp = x;
  nx = temp.l2norm();
  result = 10+ (3 * nx*nx) + pow(nx,5);
}

int main( void ) {

  // The class object for the approximation
  krigapprox MyApproximation; 
  
  // The class object for the random function (i.e. krigifier function)
  randfunc MyRandomFunction( 98567303, // Seed for the random number generator
                             42);      // Stream for the random number generator 

  int m = 100;                 // The number of initial design sites used to create the approximation
  long p = 2;                       // The dimension of the space

  Matrix<double> X(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y(m);         // Vector to store the function values at the initial design sites 

  Vector<double> lower(2), upper(2); // Vectors to store the upper and lower bounds of the space 

  Vector<double> theta(2);     // Vector for passing theta information to the approximation
  Vector<double> alpha(1);

  Matrix<double> beta2(2,2);        // Going to fill up with zeros, to have stuff to pass
  Vector<double> beta1(2);
  beta2 = 0;
  beta1 = 0;

  lower = 0; // set lower bounds to zero
  upper = 1; // set upper bounds to one

  theta[0] = 50; // set theta = 50
  alpha = theta[1] = 2;  // set alpha = 2

  ofstream outf;

  outf.open( "results2.out");  // open an output file 

  MyRandomFunction.chooseInitialDesign(2); // Latin Hypercubes

  MyRandomFunction.setTrend(TrendFunction); // Use this trend

  // Initialize and create the random function
  MyRandomFunction.setvalues(p,      // p 
			     135,    // n
			     0.0,     // unused because of user-defined trend
			     beta1, // more zeros
			     beta2,    // even more zeros
			     beta1, // yet again more zeros
			     alpha, // alpha
			     theta, // theta
			     1000,     // sigma2
			     lower,
			     upper);  

  X.newsize(m,p);  // Make X an n by p matrix, so that each row contains one point 
 
  randPoints(X, lower, upper); // Randomly generate a set of initial design sites

  // Now evaluate the random function at each of these sites
  for( int i = 0; i < m; i++)
    y[i] = MyRandomFunction.evalf(X.row(i));


  MyApproximation.setvalues(X,y,1, (&(theta[0])));  // Create the approximation

  /* Now evaluate the approximation at each point on a 100 by 100 grid in the space, and
     output these results to the file.  Results are printed in a format:

     Point Value

     So, for example, if we evaluate the function at the point x = (0.4, 0.6), and get a 
     function value y = f(x) = 23.34, then the following line would be output:

     0.4 0.6 23.34
  */

  analyzevalues( MyApproximation, outf, 100, lower, upper );

  outf.close();
  return 0;
}  // end main

void randPoints(Matrix<double> & X, Vector<double>& lower, Vector<double>& upper)
{
  long n,p;
  double temp;
  n = X.num_rows();
  p = X.num_cols();

  if( (p != lower.dim()) || (p != upper.dim()))
    cout << "randPoints::ERROR: dimension of bounds not equal to num_cols of Matrix" << endl;
  
  SelectStream(84);
  //PutSeed(7583863); // Use this seed
  PutSeed(-1);     // Get seed from system time clock
  for( long i = 0; i < n; i++)
    for( long j = 0; j < p; j++) {
      temp = lower[j];
      X[i][j] = temp + Random()*(upper[j] - temp);
    } // end for

  return;
} // end randPoints
