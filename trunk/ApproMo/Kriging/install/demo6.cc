/* demo6.cc

   This is an extremely simplistic demo.  It estimates the quadratic function
 using the glsapprox class, 
and then evaluates the estimate at a number of sites, outputing the results
to "results6.out"  

   Compile this program by typing on the command line

% make demo6

   Run the program by typing

% ./demo6

   Users can then run 

% make plotpoints
% ./plotpoints results6.out results6.ps
% gv results6.ps

   to view a plot of the generated function.


Anthony Padula
adpadu@wm.edu
6/3/2000
*/

#include "rngs.h"
#include "gls.h"

using namespace std;

// just a 2d quadratic function to estimate
double quadfunc( const Vector<double>& point ) {
  double p0, p1;
  p0 = point[0];
  p1 = point[1];
  return 5.0 - 10.0*p0 - 10.0*p1 + 10.0*p1*p1 + 10.0*p0*p0;
}

void randPoints( Matrix<double>& X, Vector<double>& lower, Vector<double>& upper);
/* Fills X with uniform random points between the bounds.  Thus
   lower[j] <= X[i][j] <= upper[j] for all 0<=j<p and all i
*/

int main( void ) {

  int n = 10;                    // number of points to give the estimate

  long p = 2;                       // the dimension of the space

  Vector<double> lower(p), upper(p);   // bounds on the space

  Vector<double> x(p);             // point at which to evaluate 
                                // the estimated function 

  double y;                     // y = f(x)

  Matrix<double> Points(n,p);        // input to the gls estimate
  Vector<double> Values(n);        // more input

  ofstream outf;


  outf.open("results6.out");

  glsapprox TheEstimate;

  lower[0] = 0; lower[1] = 0;
  upper[0] = 1; upper[1] = 1;  // work on the unit square

  randPoints(Points,lower,upper); // Randomly generate a set of points

  for( int i = 0; i < n; i++) 
    Values[i] = quadfunc( Points.row(i));

  TheEstimate.setvalues(Points, Values);

  x = (lower + upper)*0.5;  // put x in the center of the space

  y = TheEstimate.evalf(x);

  cout << "We find that f(" << x[0];
  for(int i = 1; i < p; i++)
    cout << ", " << x[i];
  cout << ") = " << y << endl;

  cout << "Scanning function and outputing results to \"results6.out\"" << endl;
		
  scanfunction(TheEstimate,       // scan this function, 
               outf,              // output results to this stream,
               100,              // divide each axis into this many points
	                         // and evaluate on the resulting grid,
	       lower,
	       upper);           // between these bounds

  return 0;
}

void randPoints( Matrix<double>& X, Vector<double>& lower, Vector<double>& upper)
{
  long n,p;
  double temp;
  n = X.num_rows();
  p = X.num_cols();

  if( (p != lower.dim()) || (p != upper.dim())) {
    cout << "randPoints::ERROR: dimension of bounds not equal to num_cols of Matrix" << endl;

    exit(1);
  }

  SelectStream(84);
  PutSeed(-1); // use this seed
  for(long i = 0; i < n; i++)
    for( long j = 0; j<p; j++) {
      temp = lower[j];
      X[i][j] = temp+Random()*(upper[j] - temp);
    } // end for

  return;
}
