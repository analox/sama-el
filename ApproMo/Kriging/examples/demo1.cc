/* demo.cc
   This is the demo main program for the Krigifier and the Kriging approximation classes
 It creates a random function, using the parameters specified in "krig.dat",
 evaluates it at several sites, and uses these values to
 create an approximation.  The approximation is then evaluated systematically and the results
 are output to the file "results.out"

 Anthony D. Padula    adpadu@wm.edu
 5/31/00
*/

/* INSTRUCTIONS FOR USE
 To compile this, you need to have downloaded the files :

 krigify.cc, krigify.h, demo5.cc rngs.c, rngs.h, rvgs.c, rvgs.h,
 krig.cc, krig.h, gamma.cc, gamma.h, ParamEstimate.cc, ParamEstimate.h,
 PatternSearch.cc, PatternSearch.h, CompassSearch.cc, CompassSearch.h,
 approx.cc, approx.h, Makefile


 You will also need the libraries CLAPACK, BLAS, and F2C, which can be obtained
 from http://www.netlib.org/
 You should modify the macros at the top of the Makefile to suit your compiler and libraries

 Then, on your command line you run

% make demo
% ./demo

*/

/* WHAT TO DO WITH THE OUTPUT 
Once you have run this program, you can plot the approximation using gnuplot.  Supposing
your CWD contains results.out, in gnuplot you should run the commmands:

  set parametric
  splot "results.out" with lines

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
  //randfunc MyRandomFunction(960587311, 42); // choose the seed for the RNG
  randfunc MyRandomFunction;   // or take a seed from the system clock

  int m = 100;       // The number of initial design sites used to create the approximation
  int p;                       // The dimension of the space ( Given in "krig.dat" )
  Matrix<double> X(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y(m);         // Vector to store the function values at the initial design sites 

  Vector<double> lower, upper; // Vectors to store the upper and lower bounds of the space 

  ofstream outf;

  outf.open( "results.out");  // open an output file 

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
  cout << "p: " << p << endl; // smenzel
  cout << "m: " << m << endl; // smenzel
  X.newsize(m,p);  // Make X an n by p matrix, so that each row contains one point 
  lower = MyRandomFunction.getlower();  // Fetch the lower bounds
  upper = MyRandomFunction.getupper();  // Fetch the upper bounds

  cout << "lower: " << lower << endl; // smenzel
  cout << "upper: " << upper << endl; // smenzel

  randPoints(X, lower, upper); // Randomly generate a set of initial design sites

  // Now evaluate the random function at each of these sites
  for( int i = 0; i < m; i++) {
    y[i] = MyRandomFunction.evalf(X.row(i));
    cout << "X.row(" << i << "): " << X.row(i) << ", y[" << i << "]: " << y[i] << endl;
  }

  // smenzel
  m = 75; p = 2;
  X.newsize(m,p);

  Vector < double > y1(m);

  FILE* fin; 
  char Inputfilename[100];
  strcpy(Inputfilename,((string)"values1" + ".dat").c_str());
  fin = fopen(Inputfilename,"rt");
  double tmpx, tmpy, tmpz;
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpx);
    X[a][0L] = tmpx;
  }
  fclose(fin);
  strcpy(Inputfilename,((string)"values2" + ".dat").c_str());
  fin = fopen(Inputfilename,"rt");
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpy);
    X[a][1L] = tmpy;
  }
  fclose(fin);
  strcpy(Inputfilename,((string)"values3" + ".dat").c_str());
  fin = fopen(Inputfilename,"rt");
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpz);
    y1[a] = tmpz;
    cout << X[a][0L] << " " << X[a][1L] << " " << y1[a] << endl;
  }
  fclose(fin);

  //  X[0L][0L] = 1.0; X[0L][1L] = 1.0;
  //  X[1L][0L] = 1.0; X[1L][1L] = 2.0;
  //  X[2L][0L] = 3.0; X[2L][1L] = 3.0;
  //  X[3L][0L] = 1.0; X[3L][1L] = 4.0;
  //  X[4L][0L] = 4.0; X[4L][1L] = 5.0;
  //  X[5L][0L] = 2.0; X[5L][1L] = 1.0;
  //  X[6L][0L] = 1.0; X[6L][1L] = 2.0;
  //  X[7L][0L] = 3.0; X[7L][1L] = 3.0;
  //  X[8L][0L] = 5.0; X[8L][1L] = 4.0;
  //  X[9L][0L] = 4.0; X[9L][1L] = 5.0;
  //  X[10L][0L] = 3.0; X[10L][1L] = 1.0;
  //  X[11L][0L] = 3.0; X[11L][1L] = 2.0;
  //  X[12L][0L] = 3.0; X[12L][1L] = 3.0;
  //  X[13L][0L] = 1.0; X[13L][1L] = 4.0;
  //  X[14L][0L] = 1.0; X[14L][1L] = 5.0;


//  y1[0] = 10.0;
//  y1[1] = 30.0;
//  y1[2] = -10.0;
//  y1[3] = 0.0;
//  y1[4] = 4.0;
// y1[5] = 0.0;
//  y1[6] = 15.0;
  //y1[7] = -20.0;
//  y1[8] = -10.0;
//  y1[9] = 20.0;
//  y1[10] = 40.0;
//  y1[11] = 250.0;
//  y1[12] = -15.0;
//  y1[13] = 0.0;
//  y1[14] = 230.0;
  
  lower = 0;
  upper = 100;

  cout << "creating approximation..." << endl;

  // Use these commands to fix theta in the approximation
  //double theta[2];
  //theta[0] = MyRandomFunction.gettheta();
  //theta[1] = MyRandomFunction.getalpha();
  //MyApproximation.setvalues(X,y,1,theta );  // Create the approximation

  // Otherwise, do this to estimate theta automatically
  double theta[2];
  theta[0] = 0.02;
  theta[1] = 0.02;
  MyApproximation.setvalues(X,y1,1);  // Create the approximation

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
  //PutSeed(86746386); // Use this seed
  PutSeed(-1);     // Get seed from system time clock
  for( long i = 0; i < n; i++)
    for( long j = 0; j < p; j++) {
      temp = lower[j];
      X[i][j] = temp + Random()*(upper[j] - temp);
    } // end for

  return;
} // end randPoints


