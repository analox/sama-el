/* demo3.cc

   This is an extremely simplistic demo.  It creates an randfunc object, 
randomly generates the function using parameter values from "krig.dat",
and then evaluates the function at a number of sites, outputing the results
to "results3.out"  

   Compile this program by typing on the command line

% make demo3

   Run the program by typing

% ./demo3

   Users can then run 

% make plotpoints
% plotpoints results3.out results3.ps
% gv results3.ps

   to view a plot of the generated function.


Anthony Padula
adpadu@wm.edu
6/3/2000
*/

#include "krigify.h"

using namespace std;

int main( void ) {

  Vector<double> lower, upper;   // bounds on the space

  long p;                       // the dimension of the space

  Vector<double> x;             // point at which to evaluate the random function 
  double y;                     // y = f(x)

  ofstream outf;


  outf.open("results3.out");

  randfunc MyRandomFunction( 960587311, // seed for the random number generator
			     42);       // stream for the rng

  MyRandomFunction.setvalues("krig.dat");

  p = MyRandomFunction.getp();           // fetch dimension from function
  lower = MyRandomFunction.getlower();
  upper = MyRandomFunction.getupper();
  x.newsize(p);  // resize x appropriately

  x = (lower + upper)*0.5;  // put x in the center of the space

  y = MyRandomFunction.evalf(x);

  cout << "We find that f(" << x[0];
  for(int i = 1; i < p; i++)
    cout << ", " << x[i];
  cout << ") = " << y << endl;

  cout << "Scanning function and outputing results to \"results3.out\"" << endl;
			       
  scanfunction(MyRandomFunction,  // scan this function, 
               outf,              // output results to this stream,
               100);              // divide each axis into this many points
	       // and evaluate on the resulting grid,

  return 0;
}




