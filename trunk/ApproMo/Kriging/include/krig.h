/* krig.h
written by Anthony Padula	adpadu@maila.wm.edu
This file contains the class definition for a kriging approximator.
When given a matrix of points and a vector of the
function values at those points, will create a kriging
approximation which can then be evaluated.
This procedure is based on equations from a number of references.

 *  Permission to use, copy, modify, and distribute this software  
 *  for any purpose without fee is hereby granted, provided that   
 *  this entire notice is included in all copies of any software   
 *  which is or includes a copy or modification of this software   
 *  and in all copies of the supporting documentation for such     
 *  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    
 *  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, THE AUTHOR    
 *  OFFERS NO REPRESENTATION OR WARRANTY OF ANY KIND                   
 *  CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    
 *  FITNESS FOR ANY PARTICULAR PURPOSE.   
*/

/* Instructions
   1. create a krigapprox object with the statement:
        krigapprox TestKRIG;
   
   2. run setvalues one of three ways:
     a. TestKRIG.setvalues() will prompt the user for the necessary information.
     However, this is only recommended for beginning users with small problems.
     
     b. TestKRIG.setvalues( filename ) will read in the same information as (a.) from the
     file.  This is much quicker and is recommended for most problems.

     c. TestKRIG.setvalues( Points, Values, deg, choice, c, beta) takes the
     necessary info as parameters.  This is useful if you have created the
     information at runtime.  See the function description below for details on
     the parameters.

   3. evaluate the function as often as desired.  This is done using:
        TestKRIG.evalf( x);
      where x is a Vector<double> of length p which contains the point at which
      to evaluate.


      NOTE: At this time, parameter estimation only works for correlation
      functions 1 and 2.  Any other correlation function must be given
      parameters by the user.  However, estimating mu always works, and will be
      done in all cases.
      
     Last Modified 7/22/99
*/
#ifndef KRIGA_H
#define KRIGA_H

#include "approx.h"
#include "gamma.h"
#include "ParamEstimate.h"
#include <math.h>

using namespace std;

const int ALPHA = 2;  // Currently used to select alpha, the parameter of the
                      // Gaussian correlation function.

const double DELTA = 0.05;  // Currently used to select the starting delta for
                            // the parameter estimation.

const int DISPLAYK = 1; // Sets the amount of output
// 0 - dead silent
// 1 - nothing but error messages
// 2 - also outputs estimates of mu, theta, and tolerance
// 4 - also LOTS of messages updating it's current position
// 5 - also outputs most important matrices.  BE WARNED: Matrices are often HUGE


class krigapprox: public approx {
public:
  krigapprox(); 
  krigapprox( krigapprox& ToBeCopied );
  krigapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
              Vector<double>& Values,  // The vector of function values, so
                                       // f(Points[i])=Values[i] 
              int choice); // The flag for which correlation function to use
                          // (see chooseCorrelation() below)
         
  krigapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
              Vector<double>& Values,  // The vector of function values, so
                                       // f(Points[i])=Values[i] 
              int choice, // The flag for which correlation function to use
                          // (see chooseCorrelation() below)
              double* theta); // The known value(s) of the parameter(s)
  // theta will either be a single number or an array of numbers.  If
  // there are multiple parameters,
  // they are expected to be stored in the array in order of increasing
  // dimension, and all theta_i before alpha_i.  Thus, in 3 dimensions, the
  // array might contain [1.0, 2.0, 3.0, 0.5, 4.5, 0] if theta_1 = 1.0, 
  // theta_2 = 2.0, theta_3 = 3.0, alpha_1 = 0.5, alpha_2 = 4.5, alpha_3 = 0.0
  // See the chooseCorrelation() function below for more discussion of parameters

  virtual ~krigapprox(); // Destructor

  // evalf returns the value of the function at the point passed in 
  // point is a p-dim vector at which f is to be evaluated
  virtual double evalf(const Vector<double>& point);


  // evalderiv returns the vector value of the first derivative of the
  // interpolant at the point passed in.
  //
  // Note: in some cases, the derivative at interpolated points
  // is undefined.  The function will then print an error message
  // and return HUGE_VAL as the derivative  
  //
  //  STILL UNDER CONSTRUCTION done for 1,2,3,5
  // tested for 1,2,3
  virtual Vector<double> evalderiv( const Vector<double>& point);

  /* virtual int setvalues();  

     This is an I/O function to interactively set the necessary parameters
     The data will be read in the format:
     p = long integer  // the dimension of the space
     n = long integer  // the number of points to be interpolated
     choice = integer  // the choice of  which correlation function to use
                       // (see chooseCorrelation() below)
     query = y/n       // The program will then ask the user if
                       // the parameters are known.
		       // if 'y' then will prompt for the parameters
		       // if 'n' then will estimate the parameters from the data
     Points = a matrix of n x p+1 doubles, such that each row is one point follwed
       by the function value at that point.  For example, if f( (2.3, 4.5) ) = 8.2,
       then you would input:
       2.3 4.5 8.2

     This is the same format used by reportsample() and scanfunction() below
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
     1
     n
     1.5 1.7 20.6
     0.3 2.6 563.7
     3.6 2.9 27.04
     1.0 0.0 100.3
     3.87 4.96 0.06

     NOTE:  Will only read in parameters if necessary, depending on the 4th input.
  */
  virtual int setvalues(char* filename);
  

  /* This function works like the special constructor, allowing the data to be
     passed in as parameters instead of being read from I/O.
     Stores all data and calls approximate() to create the approximation
  */
  virtual int setvalues(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
                          Vector<double>& Values,  // The vector of function
                                                   // values, so f(Points[i])=Values[i]
                          int choice); // The flag for which correlation function to use
                                       // (see chooseCorrelation() below)                        
 

   /* This function works like the special constructor, allowing the data to be
     passed in as parameters instead of being read from I/O.
     Stores all data and calls approximate() to create the approximation
  */
  virtual int setvalues( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
                         Vector<double>& Values,  // The vector of function values, so
                                                  // f(Points[i])=Values[i] 
                         int choice, // The flag for which correlation function to use
                                     // (see chooseCorrelation() below)
                         double* theta); // The known value(s) of the parameter(s)
  // theta will either be a single number or an array of numbers.  If
  // there are multiple parameters,
  // they are expected to be stored in the array in order of increasing
  // dimension, and all theta_i before alpha_i.  Thus, in 3 dimensions, the
  // array might contain [1.0, 2.0, 3.0, 0.5, 4.5, 0] if theta_1 = 1.0, 
  // theta_2 = 2.0, theta_3 = 3.0, alpha_1 = 0.5, alpha_2 = 4.5, alpha_3 = 0.0
  // See the chooseCorrelation() function below for more discussion of parameters

 
  /* void chooseCorrelationFamily( int choice );
     
     Allows the user to choose among the possible families of correlation
     functions.  Will delete mu and theta and reapproximate if the other values
     have been set.
     
     Currently, the options are:
     
     choice = 0  =>  User defined function passed in previously
     
     choice = 1  =>  Gaussian isotropic r(x,y) = exp( -theta * ||y-x||^alpha )
     
     choice = 2  =>  Gaussian product \prod_{j=1}^{p} exp(-theta_j * |s_j - t_j|^alpha)
     
     choice = 3  =>  Product of Linear  \prod_{j=1}^{p} (1-theta_j * |y_j - x_j|)
     
     choice = 5  =>  Mate'rn correlation function
                 \prod_{j=1}^{p} ((theta_j * |y_j - x_j|^alpha_j) /
		 (2^(alpha_j -1) * gamma(alpha_j))) *
		 K_{alpha_j}( theta_j * |y_j - x_j| )

     where K_{nu}(x) is the degree nu modified Bessel function of the third order. 

     The default choice is Gaussian isotropic, so you really only need to run
     this function in order to choose one of the other families.
     
     Notice that choices 2 and 3 take a vector of theta and a single alpha,
     while choice 4 takes a vector of theta and a vector of alpha.  Your input
     when running set values must be appropriate to the correlation function.
  */
  void chooseCorrelationFamily( int choice );
    

  /* void chooseCorrelationFamily( int choice,   // indicates the correlation function
                                double* theta ); // used for passing known
                                                 // parameter values.  If values
                                                 // unknown, pass in theta=NULL.
  
     Allows the user to choose among the possible families of correlation
     functions.  If the user desires a correlation function not available as
     one of the choices, then use of the chooseCorrelationFamily( void)
     function below is recommended.
     
     Currently, the options are:
     
     choice = 0  =>  User specified function (See below)
      
     choice = 1  =>  Gaussian isotropic r(x,y) = exp( -theta * ||y-x||^alpha )
     
     choice = 2  =>  Gaussian product \prod_{j=1}^{p} exp(-theta_j * |y_j - x_j|^alpha)

     choice = 3  =>  Product of Linear  \prod_{j=1}^{p} ( 1-theta_j*|y_j - x_j|)+
		    
     choice = 5  =>  Mate'rn correlation function
                    \prod_{j=1}^{p} ((theta_j * |y_j - x_j|^alpha_j) /
                    (2^(alpha_j -1) * gamma(alpha_j))) *
                    K_{alpha_j}( theta_j * |y_j - x_j| )

     where K_{nu}(x) is the degree nu modified Bessel function of the third order. 

     The default choice is Gaussian isotropic, so you really only need to run
     this function in order to choose one of the other families.  

     Notice that choices 2 and 3 take a vector of theta and a single alpha,
     while choice 4  and 5 takes a vector of theta and a vector of alpha.
     Also, in choice 5, the values of alpha must be all integers.  Your input
     when running set values must be appropriate to the correlation function.
     
     newrand() MUST be called after changing the correlation function, either
     directly or via setvalues().  Otherwise your results will be nonsense.
     
     Theta > 0!
     
     Warning: Do not use the Mate'rn family unless you really know how to set
     the parameters.  Poor choices of parameters tend to cause problems with
     the SVD.  The author of this code CANNOT guarentee proper functioning if
     you use choice 5!  alpha > 0 and alpha integer.  I suggest alpha=1 for all
     alpha and theta around 3.
  */
  void chooseCorrelationFamily( int choice,  
                                double* theta );
  

  /* void chooseCorrelationFamily( void (*corfunc)(const Vector<double>& x,
                                                const Vector<double>& y,
                                                double& result ) );
     This is an alternative to the previous function.  It allows the user to
     pass in a correlation function of their own.  The correlation function
     must accept two equal length vectors as input, and the double parameter is
     used for output of the result.
     
     This sets the choice to option 0, and choosing another choice later will
     use a different correlation function without overwriting the pointer.
     Thus, you may switch back and forth between a hard-wired function and a
     user specified function by indicating a new choice as directed above.

     NOTE: NOT EVERY ARITHMETIC FUNCTION WORKS AS A CORRELATION FUNCTION
  */ 
  void chooseCorrelationFamily( void (*corfunc)(const Vector<double>& /*input*/ x,
                                                const Vector<double>& /*input*/ y,
                                                double& /*output*/ result ) );


  // functions to get copies of the data members
  
  // Note that this returns pointer to an object which must be
  // deleted by user to avoid memory leak.
  Vector<double>* gettheta(); 
  
  // returns v = V(D^-1)U' y
  Vector<double> getv();  
  

  /****
     Prints the data of the kriging to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing functions so that they might be re-created later.
  ****/
  void storeApprox( ostream& outf );

  
protected:

  /* This returns the value of the derivative of |a-b| with respect to b,
   i.e. 1.0 or -1.0 . */
  double derivabs( double a, double b);
  
  /* The value of the trend at this point */
  virtual double trend( const Vector<double> & point);

  /* The value of the first derivative of the trend at point */
  virtual Vector<double>* derivtrend( const Vector<double> & point);
 
  // deallocates any dynamically allocated members and sets
  // them to NULL
  void deleteall(); 

  // returns the proper number of parameters
  int thetasize();  

  // returns the proper size of mu
  int musize();

  int estimateTheta();
  void estimateMu(const Matrix<double> & R);
  
  // creates the approximation using the current values
  // Will estimate theta if necessary
  // Returns TRUE if Successful, FALSE otherwise.
  virtual int approximate(); 

  /****
       virtual double correlation( const Vector<double>& x, const Vector<double>& y);
       
       INPUT: Two points x and y in the space
       OUTPUT: The value of the correlation between them.  Currently, this uses the
       Gaussian correlation function, but plans are to later allow the user
       to select one of several correlation functions at the time of input
       EFFECT: Calls l2norm() on the difference (ie. y - x ) of the vectors and then
       calculates the Gaussian correlation.
  ****/  
   virtual double correlation( const Vector<double>& x, const Vector<double>& y);

  int corfamily; // flag to choose the correlation function
  
  void (*corfunction)( const Vector<double>& x, const Vector<double>& y, double& result);

  double* theta_;
  double* mu_;

  Vector<double> v;  // used by evalf
}; // end class krigapprox

#endif
