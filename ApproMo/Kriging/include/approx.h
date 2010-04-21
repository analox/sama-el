#ifndef APPROX_H
#define APPROX_H

#include <cmath>
#include "maps_general.h"
#include <iostream>

using namespace std;

/*
#ifndef BOOL_C
#define BOOL_C
int TRUE = 1;
int FALSE = 0;
#endif
*/
#ifndef PI_C
#define PI_C
double PI = 3.1415926535897932384;
#endif

 int DEFAULTPOINTS = 25;


/* Choose the appropriate machine epsilon for your system, or declare one with
   the appropriate value.  If you don't know machine epsilon, it's easy to find:

   macheps = 1;
   do {
      macheps = macheps/2.0;
   } while ( 1+macheps > 1);
   macheps = macheps*2.0;

   Or else just simply choose the MACHEPS = 1.0 option and the program will run
   automatically to find the correct machine epsilon.
*/

class approx { // abstract base class
public:

  approx();
  virtual ~approx();
  
  // evalf returns the value of the function at the point passed in 
  virtual double evalf(const Vector<double>& point) = 0;
  // point is a p-dim vector at which f is to be evaluated
  
  virtual int setvalues() = 0;  // This is an I/O function to interactively set the
                             // necessary parameters

  virtual int setvalues(char* filename) = 0;
   /* Reads the parameter values in order from a file and runs approximate() to
      build the approximation
      Returns 1 if successful, and 0 if insufficient data in the file.
      Format is the same as the above.
   */

  // functions to get copies of the data members
  long getp();
  long getn();
  
  Matrix<double> getPoints();
  Vector<double> getValues();

  virtual void storeApprox( ostream& outf ) = 0;
/****
     Prints the data of the approximation to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing functions so that they might be re-created later.
****/
  
protected:
  
  void notseterror();
  virtual void deleteall() = 0;  // deallocates any dynamically allocated members and sets
  // them to NULL

  int virtual approximate() = 0;  // creates the approximation using the current values
  // Returns TRUE if Successful, FALSE otherwise.

  bool ValuesSet;
  long p;        // dimension of the space
  long n;        // number of sites
  
  Matrix<double> X;  // n p-dimensional points (resized as necessary)
  Vector<double> y;  // function values for the points in X	
}; // end class approx


void scanfunction( approx& Func, ostream& outf, long pointsperaxis,
                   Vector<double>& lower, Vector<double>& upper);
/**** 
      INPUT: The krigapprox object to be sampled and an ofstream to output the
      values, as well as the number of points per axis to be sampled.  Thus, the
      total number of sample points is approximately (pointsperaxis)^p.  Should
      work for any value of p.
      
      Bound is a 2 x p matrix, with the lower Bound in the first row and the
      upper Bound in the second row.

      EFFECT: Outputs to the stream values in columns by axis.  For example:

      2.3 3.4 13.695
      
      would correspond to the point (2.3, 3.4, 13.695) on the Cartesian
      coordinate system.

****/

void scanfunction( approx& Func, ostream& outf, Vector<double>& lower, Vector<double>& upper);
// A default version which uses 25 pointsperaxis.  You can modify this default
// by changing the 25 above.

double reportsample( approx& Func, Vector<double>& x, ostream& outf);
/****
     INPUT: The krigapprox object to be sampled, the point at which to evaluate
     the function, and the output stream to which it should report. 
     OUTPUT: the value of the function at x - result of Func.evalf(x)
     EFFECT: Outputs to the stream the elements of x, seperated by spaces, and
     then the function value at x, followed by an endline. For example:

     1.0 2.0 3.0 4.0 5.0 6.0 21.0

     would be written if x == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0) and f(x) == 21.0

****/

void scanplane( approx& Func, ostream& outf, long pointsperaxis,
                Vector<double>& lower, Vector<double>& upper );
/**** 
      A simplified, but faster version of scanfunction().
      ONLY works when p=2 !!!
      

      INPUT: The krigapprox object to be sampled and an ofstream to output the
      values.

      Bound is a 2 x p matrix, with the lower Bound in the first row and the
      upper Bound in the second row.
      
      EFFECT: Outputs to the stream values in columns by axis.  For example:

      2.3 3.4 13.695
      
      would correspond to the point (2.3, 3.4, 13.695) on the Cartesian
      coordinate system. 
****/

void analyzevalues( approx& Func, ostream& outf, Vector<double>& lower, Vector<double>& upper );
/****
     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Uses a default value of 25 points per axis.
     
     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.

     
****/

void analyzevalues( approx& Func, ostream& outf, long pointsperaxis,
                    Vector<double>& lower, Vector<double>& upper );
/****
     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Divides each axis into pointsperaxis evenly spaced points, and samples at
     each permutation of values.
     
     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.
****/

#endif



