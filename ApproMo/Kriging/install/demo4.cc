/* demo4.cc
   This is the demo main program for the Krigifier and the Kriging approximation classes
 It creates a random function, using the parameters specified in "krig.dat",
 evaluates it at several sites, and uses these values to
 create an approximation.  The approximation is then evaluated systematically and the results
 are output to the file "results4.out"

 This demo is very similar to the one in demo.cc, except that it also evaluates the
 derivatives of the kriging approximation instead of the approximation function itself.
 The program keeps track of the smallest point it finds and displays this location
 along with the derivative at this location.

 Anthony D. Padula    adpadu@wm.edu
 6/7/00
*/

/* INSTRUCTIONS FOR USE
 To compile this, you need to have downloaded the files :

 krigify.cc, krigify.h, demo4.cc rngs.c, rngs.h, rvgs.c, rvgs.h,
 krig.cc, krig.h, gamma.cc, gamma.h, ParamEstimate.cc, ParamEstimate.h,
 PatternSearch.cc, PatternSearch.h, CompassSearch.cc, CompassSearch.h,
 approx.cc, approx.h, Makefile

 You will also need the libraries CLAPACK, BLAS, and F2C, which can be obtained
 from http://www.netlib.org/
 You should modify the macros at the top of the Makefile to suit your compiler and libraries

 Then, on your command line you run

% make demo4
% ./demo4

*/

/* WHAT TO DO WITH THE OUTPUT 
Once you have run this program, you can plot the approximation using gnuplot.  Supposing
your CWD contains results.out, in gnuplot you should run the commmands:

  set parametric
  splot "results4.out" with lines

This will produce a nice 3-d graph of the approximation  


Alternatively, you can compile and run plotpoints by

% g++ plotpoints.cc -o plotpoints
% ./plotpoints results.out results.ps

This produces a graph in postscript format, which can then be viewed or printed.
*/


#include "krig.h"
#include "krigify.h"
#include "rngs.h"


using namespace std;

void randPoints( Matrix<double> & X, Vector<double>& lower, Vector<double> & upper);
/* Fills X with uniform random points between the bounds.  Thus
   lower[j] <= X[i][j] <= upper[j]   for all 0<=j<p
*/



int main( void ) {

  // The class object for the approximation
  krigapprox MyApproximation; 
  
  // The class object for the random function (i.e. krigifier function)
  randfunc MyRandomFunction;

  int m = 100;                 // The number of initial design sites used to create the approximation
  int p;                       // The dimension of the space ( Given in "krig.dat" )
  Matrix<double> X(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y(m);         // Vector to store the function values at the initial design sites 

  Vector<double> lower, upper; // Vectors to store the upper and lower bounds of the space 

  ofstream outf;

  outf.open( "results4.out");  // open an output file 

  MyRandomFunction.setvalues("krig.dat");  // Initialize and create the random function

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

  // Create the approximation
  MyApproximation.setvalues(X,   // design sites
			    y,   // corresponding function values y[i] = f(X.row(i))
			    1);  // choice of correlation function

  /* Now evaluate the approximation at each point on a 100 by 100 grid in the space, and
     output these results to the file.  Results are printed in a format:

     Point Value

     So, for example, if we evaluate the function at the point x = (0.4, 0.6), and get a 
     function value y = f(x) = 23.34, then the following line would be output:

     0.4 0.6 23.34
  */


  Vector<double> derivsmallest(p);
  double valuesmallest, result;
  Vector<double> current(p), smallest(p);
  Vector<double> increment(p);
  long j;
  double temp = 1.0/100;

  smallest = HUGE_VAL;
  derivsmallest = HUGE_VAL;
  valuesmallest = HUGE_VAL;
  current = lower;
  increment = upper-lower;
  increment = increment * temp;

  while( current[0] < (upper[0]+increment[0]/2.0) ) {
    for( int i = 0; i < p; i++ ) {
      outf << current[i] << ' ';
    } // end for
    result = MyApproximation.evalf(current);
    outf << result << endl;
    if( result < valuesmallest) {
      derivsmallest = MyApproximation.evalderiv(current);
      smallest = current;
      valuesmallest = result;
    }
      
    current[p-1] += increment[p-1];
    for( j = p-1; (j > 0)&&(current[j] > (upper[j]+increment[j]/2.0)); j--) {
      current[j] = lower[j];
      current[j-1] += increment[j];
      if( j-1 == 0)
	outf << endl;
    } // end for
  } // end while   

  cout << "Smallest value " << valuesmallest 
       << " found at (" << smallest[0] << ", " << smallest[1] 
       << ")\nwith a derivative of ("<< derivsmallest[0] << ", " 
       << derivsmallest[1] << ")"<< endl;
 
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
  //PutSeed(74728294); // Use this seed
  PutSeed(-1);     // Get seed from system time clock
  for( long i = 0; i < n; i++)
    for( long j = 0; j < p; j++) {
      temp = lower[j];
      X[i][j] = temp + Random()*(upper[j] - temp);
    } // end for

  return;
} // end randPoints


