/* rbf.h
written by Anthony Padula	adpadu@maila.wm.edu
This file contains the class definition for a radial basis
function approximator.  When given a matrix of points and a vector of the
function values at those points, will create a radial basis function
approximation which can then be evaluated.
This procedure is based on equations from a number of references.
*/

/* Instructions
   1. create a rbfapprox object with the statement:
        rbfapprox TestRBF;
   
   2. run setvalues one of three ways:
     a. TestRBF.setvalues() will prompt the user for the necessary information.
     However, this is only recommended for beginning users with small problems.
     
     b. TestRBF.setvalues( filename ) will read in the same information as (a.) from the
     file.  This is much quicker and is recommended for most problems.

     c. TestRBF.setvalues( Points, Values, deg, choice, c, beta) takes the
     necessary info as parameters.  This is useful if you have created the
     information at runtime.  See the function description below for details on
     the parameters.

   3. evaluate the function as often as desired.  This is done using:
        TestRBF.evalf( x);
      where x is a Vector<double> of length p which contains the point at which
      to evaluate.
      
     Last Modified 7/21/99
*/
#ifndef RBF_H
#define RBF_H

#include "approx.h"
#include "gamma.h"

using namespace std;

const int DISPLAY = 0;  // DISPLAY controls to amount of output
// =0   Only print error messages and prompts for data
// =1   Displays b and a as well
// =2   Displays right-A*x, the residual vector as well

class rbfapprox: public approx {
public:
  rbfapprox(); 
  rbfapprox( rbfapprox& ToBeCopied );
  rbfapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
             Vector<double>& Values,  // The vector of function values, so
                                      // f(Points[i])=Values[i]
             long deg,    // The size of the polynomial basis
             // deg < 0  no polynomial basis
             // deg = 0  constant 1 only
             // 0<deg<=p use the first deg monomials and a constant 1 
             int choice, // The flag for which RBF to use (see chooseRBF() below)
             double c,   // c and beta are the parameters for the RBF 
             double beta ); // You may not need both for all choices, but always
                            // pass in something at least
  
  ~rbfapprox(); // Destructor

  // evalf returns the value of the function at the point passed in 
  double evalf(const Vector<double>& point);
  // point is a p-dim vector at which f is to be evaluated
  
  int setvalues();  // This is an I/O function to interactively set the
                             // necessary parameters
   /* The data will be read in the format:
   p = long integer  // the dimension of the space
   n = long integer  // the number of points to be interpolated
   degree = long     // The size of the polynomial basis
                     // deg < 0  no polynomial basis
                     // deg = 0  constant 1 only
                     // 0<deg<=p use the first deg monomials and a constant 1 
   choice = integer  // the choice of radial basis function (see chooseRBF() below)
   c = double        // Will only prompt for c and/or beta if needed
   beta = double
   Points = a matrix of n x p+1 doubles, such that each row is one point follwed
   by the function value at that point.  For example, if f( (2.3, 4.5) ) = 8.2,
   then you would input:
   2.3 4.5 8.2
   This is the same format used by reportsample() and scanfunction() below
   */

  int setvalues(char* filename);
   /* Reads the parameter values in order from a file and runs approximate() to
      build the approximation
      Returns 1 if successful, and 0 if insufficient data in the file.
      Format is the same as the above.
      Ex.
      2
      5
      -1
      3
      1.0
      -1.0
      1.5 1.7 20.6
      0.3 2.6 563.7
      3.6 2.9 27.04
      1.0 0.0 100.3
      3.87 4.96 0.06

      NOTE:  Will only read in c and/or beta if necessary, depending on the
      choice of basis function.
   */

  int setvalues(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
                          Vector<double>& Values,  // The vector of function
                                                   // values, so f(Points[i])=Values[i]
                          long deg,  // The size of the polynomial basis
                          // deg < 0  no polynomial basis
                          // deg = 0  constant 1 only
                          // 0<deg<=p use the first deg monomials and a constant 1 
                          int choice, // The flag for which RBF to use (see chooseRBF() below)
                          double c,   // c and beta are the parameters for the RBF 
                          double beta ); // You may not need both for all choices, but always
                                         // pass in something at least           
  /* This function works like the special constructor, allowing the data to be
     passed in as parameters instead of being read from I/O.
     Stores all data and calls approximate() to create the approximation
  */
  void chooseRBF( int choice, double c, double beta, long deg );
  /*  Allows the user to choose among the possible families of radial basis
      functions.  If the user desires a basis function not available as
      one of the choices, then use of the chooseRBF()
      function below is recommended.

      Currently, the options are:

      0  =>  User specified function (See below)

      1 => Thin-plate splines:  g(r) = r^beta                 beta>0 and beta not even

      2 => Gaussian:            g(r) = exp(-c*(r^2))          c>0

      3 => Multiquadrics        g(r) = (c^2 + r^2)^(beta/2)   beta>-d and beta not even,
                                                           d > beta/2 

      4 => Thin-plate splines   g(r) = (-1)^(beta/2+1)*r^(beta)*log(r)    beta even

      5 => Sobolew splines      g(r) = (2*PI^(c))*K(nu,2*PI*r)*(r^nu)/gamma(c)
                                       nu = c-d/2, c > d/2

      6 => Thin-plate splines   g(r) = (cr)^2 * log(cr)        reference suggests c=1

      7 => Madych & Nelson II pg 221 g(r) = -2*sqrt(pi*(1+r^2))
     

      where K_{nu}(x) is the degree nu modified Bessel function of the third order. 

      This function automatically constructs the new approximation.
      setvalues() must have been run prior to using chooseRBF().
  */
  
  void chooseRBF( void (*RBFtoUse)(double r, double& /*output*/ result ), long deg );
    /*This is an alternative to the previous function.  It allows the user to
      pass in a correlation function of their own.  The RBF must accept the
      double r = ||u-v||, and return the value of the RBF using the second
      parameter.

      This sets the choice to option 0, and choosing another choice later will
      use a different correlation function without overwriting the pointer.
      Thus, you may switch back and forth between a hard-wired function and a
      user specified function by indicating a new choice as directed above.

      NOTE: NOT EVERY ARITHMETIC FUNCTION WORKS AS A RADIAL BASIS FUNCTION
      Be sure to choose deg appropriately for your function.  If you are unsure,
      use deg = p to be safe.
    */ 
    
  // functions to get copies of the data members

  double getc();
  double getbeta();

  void storeApprox( ostream& outf );
/****
     r
     Prints the data of the rbf to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing RBF's so that they might be re-created later.
****/
  
protected:

  void deleteall(){} // no dynamic variables to be deleted
  int approximate();  // creates the approximation using the current values
  // Returns TRUE if Successful, FALSE otherwise.
  
  double polybasis( Vector<double> x, long k);
  double rbfvalue( const double r, double& result);
  /****
       INPUT: r = norm(u-v)
       OUTPUT: The value of the RBF at r. 
       EFFECT: calculates the value of the RBF indicated by the member choice at r.
  ****/  

  int rbffamily; // flag to choose the radial basis function
  long degree;    // size of the polynomial basis
  
  void (*phi)( double r, double& result);

  double c_;
  double beta_;  //
  Vector<double> b;
  Vector<double> a;
}; // end class rbfapprox

#endif

