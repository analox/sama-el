/* demo5.cc
   This is the demo main program for the Krigifier and the Kriging approximation classes
 It creates a random function, using the parameters specified in "krig.dat",
 evaluates it at several sites, and uses these values to
 create an approximation.  The approximation is then evaluated systematically and the results
 are output to the file "results5.out"

 This is indentical to demo.cc except that it uses the krigapproxWithGLS class that has the
 extra quadratic trend.

 Anthony D. Padula    adpadu@wm.edu
 5/31/00
*/

/* INSTRUCTIONS FOR USE
 To compile this, you need to have downloaded the files :

 krigify.cc, krigify.h, demo5.cc rngs.c, rngs.h, rvgs.c, rvgs.h,
 krig.cc, krig.h, gamma.cc, gamma.h, ParamEstimate.cc, ParamEstimate.h,
 PatternSearch.cc, PatternSearch.h, CompassSearch.cc, CompassSearch.h,
 approx.cc, approx.h, Makefile

 gls.cc, gls.h

 You will also need the libraries CLAPACK, BLAS, and F2C, which can be obtained
 from http://www.netlib.org/
 You should modify the macros at the top of the Makefile to suit your compiler and libraries

 Then, on your command line you run

% make demo5
% ./demo5

*/

/* WHAT TO DO WITH THE OUTPUT 
Once you have run this program, you can plot the approximation using gnuplot.  Supposing
your CWD contains results5.out, in gnuplot you should run the commmands:

  set parametric
  splot "results5.out" with lines

This will produce a nice 3-d graph of the approximation  


Alternatively, you can compile and run plotpoints by

% g++ plotpoints.cc -o plotpoints
% plotpoints results5.out results5.ps

This produces a graph in postscript format, which can then be viewed or printed.
*/


#include "krig.h"
#include "krigify.h"
#include "rngs.h"
#include "gls.h"

using namespace std;

void randPoints( Matrix<double> & X, Vector<double>& lower, Vector<double> & upper);
/* Fills X with uniform random points between the bounds.  Thus
   lower[j] <= X[i][j] <= upper[j]   for all 0<=j<p
*/



int main( void ) {

  // The class object for the approximation
  krigapproxWithGLS MyApproximation; 
  
  // The class object for the random function (i.e. krigifier function)
  randfunc MyRandomFunction;

  int m = 200;                 // The number of initial design sites used to create the approximation
  int p;                       // The dimension of the space ( Given in "krig.dat" )
  Matrix<double> X(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y(m);         // Vector to store the function values at the initial design sites 

  Vector<double> lower, upper; // Vectors to store the upper and lower bounds of the space 

  ofstream outf;

  outf.open( "results5.out");  // open an output file 

  MyRandomFunction.setvalues("krig.dat");  // Initialize and create the random function

  cout << "Done setting values" << endl;
   /* The data will be read in the format:
   p = long integer
   n = long integer
   beta0 = double
   beta1 = p doubling-point values separated by spaces
   beta2 = p*p doubling-point values read along the rows
      Ex. the Matlab matrix [1 2; 3 4] would be input
           1 2 3 4
   x0 = p doubling-point values separated by spaces
   alpha = double                    
   theta = double                  
   sigma2 = double
   lower = p floating-point values separated by spaces
   upper = p floating-point values separated by spaces
   */

  p = MyRandomFunction.getp();  // Fetch the dimension of the space
  X.newsize(m,p);  // Make X an n by p matrix, so that each row contains one point 
  lower = MyRandomFunction.getlower();  // Fetch the lower bounds
  upper = MyRandomFunction.getupper();  // Fetch the upper bounds

  randPoints(X, lower, upper); // Randomly generate a set of initial design sites

  // Now evaluate the random function at each of these sites
  for( int i = 0; i < m; i++)
    y[i] = MyRandomFunction.evalf(X.row(i));

  cout << "creating approximation..." << endl;

  double theta[2];
  theta[0] = MyRandomFunction.gettheta();
  theta[1] = MyRandomFunction.getalpha();
  MyApproximation.setvalues(X,y,1,theta );  // Create the approximation

  /* Now evaluate the approximation at each point on a 100 by 100 grid in the space, and
     output these results to the file.  Results are printed in a format:

     Point Value

     So, for example, if we evaluate the function at the point x = (0.4, 0.6), and get a 
     function value y = f(x) = 23.34, then the following line would be output:

     0.4 0.6 23.34
  */

  cout << "Done creating the approximation.  Now scanning" << endl;

  scanfunction( MyApproximation, outf, 100, lower, upper );

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
  PutSeed(-1); // Use this seed
  //PutSeed(-1);     // Get seed from system time clock
  for( long i = 0; i < n; i++)
    for( long j = 0; j < p; j++) {
      temp = lower[j];
      X[i][j] = temp + Random()*(upper[j] - temp);
    } // end for

  return;
} // end randPoints


