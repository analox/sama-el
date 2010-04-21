/* gls.h  

   This is the header file for the generalized linear least squares class,
   that will be used to implement a quadratic trend in the MAPS algorithm.

   Currently just estimates a quadratic trend, solving the problem

   y = X*beta + e

   to minimize e.  We choose inner product <u,v> = u Sigma^{-1} v, with Sigma = sigma^{2} I.

   The solution is then \hat{beta} = ( X^{T} X )^{-1} X^{T} y.


   Anthony Padula  6/10/2000
*/

#ifndef GLS_H
#define GLS_H

#include "approx.h"
#include "gamma.h"
#include "krig.h"
#include <cmath>

using namespace std;

const int DISPLAYG = 1; // Sets the amount of output
// 0 - dead silent
// 1 - nothing but error messages
// 2 - also outputs estimates of mu, theta, and tolerance
// 4 - also LOTS of messages updating it's current position
// 5 - also outputs most important matrices.  BE WARNED: Matrices are often HUGE


class glsapprox: public approx {
public:
  glsapprox(); 
  glsapprox( glsapprox& ToBeCopied );
  glsapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
              Vector<double>& Values);  // The vector of function values, so
                                       // f(Points[i])=Values[i]  

  glsapprox(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
	      Vector<double>& Values, // The vector of function
                                      // values, so f(Points[i])=Values[i]
	      Matrix<double>& InvCorr); // The inverse of the correlation matrix 

  virtual ~glsapprox(); // Destructor

  // evalf returns the value of the function at the point passed in 
  // point is a p-dim vector at which f is to be evaluated
  virtual double evalf(const Vector<double>& point);


  // evalderiv returns the vector value of the gradient at the point
  virtual Vector<double> * evalderiv( const Vector<double> & point);

  /* virtual int setvalues();  

     This is an I/O function to interactively set the necessary parameters
     The data will be read in the format:
     p = long integer  // the dimension of the space
     n = long integer  // the number of points to be interpolated
     Points = a matrix of n x p+1 doubles, such that each row is one point follwed
       by the function value at that point.  For example, if f( (2.3, 4.5) ) = 8.2,
       then you would input:
       2.3 4.5 8.2

     This is the same format used by reportsample() and scanfunction()
  */
  virtual int setvalues();  

  
  /* virtual int setvalues(char* filename);

     Reads the parameter values in order from a file and runs approximate() to
     build the approximation
     Returns 1 if successful, and 0 if insufficient data in the file.
     Format is the same as the above.
     Ex.
     2
     5
     1.5 1.7 20.6
     0.3 2.6 563.7
     3.6 2.9 27.04
     1.0 0.0 100.3
     3.87 4.96 0.06
  */
  virtual int setvalues(char* filename);
  

  /* 
     This function works like the special constructor, allowing the data to be
     passed in as parameters instead of being read from I/O.
     Stores all data and calls approximate() to create the approximation
  */
  virtual int setvalues(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
                          Vector<double>& Values);// The vector of function
                                                   // values, so f(Points[i])=Values[i]

  virtual int setvalues(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
                          Vector<double>& Values, // The vector of function
                                                   // values, so f(Points[i])=Values[i]
			  Matrix<double>& InvCorr); // The inverse of the correlation matrix 


  // Pass in the appropriate correlation matrix R so that 
  // <u,v> = u^{T} R^{-1} v  
  // is the inner product used
  virtual int setInverseCorrelation( Matrix<double>& InvCorr );


  // functions to get copies of the data members
  
  // Note that this returns pointer to an object which must be
  // deleted by user to avoid memory leak.
  Vector<double> getbeta(); 
    
  // Returns the vector of residuals e = y - X*beta 
  Vector<double> getresiduals();

  long getk();

  /****
     Prints the data of the approximation to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing functions so that they might be re-created later.
  ****/
  void storeApprox( ostream& outf );

  
protected:

  // deallocates any dynamically allocated members and sets
  // them to NULL
  void deleteall(); 
  
  // creates the approximation using the current values
  // Returns TRUE if Successful, FALSE otherwise.
  int approximate(); 

  Vector<double> beta;  // used by evalf
  Matrix<double>* Rinv;
 
  // p, n, X, and y inherited from approx

  long k;   // k = 1 + p + p*(p+1)/2
}; // end class glsapprox


class krigapproxWithGLS: public krigapprox {
 public:
  ~krigapproxWithGLS(){ deleteall(); }

  // Uses the new gls trend
  virtual double trend( const Vector<double> & point);

  // The derivative of the GLS trend
  virtual Vector<double> * derivtrend( const Vector<double>& point);

 protected:
  // creates the approximation using the current values
  // Will estimate theta if necessary
  // Returns TRUE if Successful, FALSE otherwise.
  virtual int approximate(); 

  glsapprox GLStrend;
}; // end class krigapproxWithGLS
#endif

