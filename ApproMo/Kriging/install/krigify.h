/* krigify.h
written by Anthony Padula	adpadu@maila.wm.edu
This file contains the function definitions for
krigify and evalf to produce random functions
with a specified trend, as if they were being
drawn from a second-order Gaussian function.
This procedure is described in
	Trosset Michael W.
	The Krigifier: A Procedure for Generating Pseudorandom
	Nonlinear Objective Functions for Computational
	Experimentation

This is a modified version of krigify.h that has been upgraded to attempt to use
the Cholesky factorization if possible.  If the Cholesky fails, it will default
back to the SVD.

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
   1. Declare an randfunc object
   2. run chooseCorrelationFamily(), if desired.  Otherwise, defaults
      to using the Gaussian isotropic family.
   3. run setvalues() to set the parameter values in order to define
      the random function.  krigify() will automatically be called in this procedure.
   4. call objectname.evalf(x) to evaluate the random function at a given
      point x in the p-space.
   5. repeat step 4 as many times as desired to produce a sample
   6. a new random function can be defined simply by calling
      objectname.newrand(), which generates a new X and v using the same
      parameters, or objectname.newrandsameX(), which only generates a
      new v. 
      Alternatively, objectname.setvalues() can also be called and fed new
      parameters.

      The text file krigifier_instructions.txt contains more thorough
      instructions and explanations.
      
     Note:  It is possible to create the same function on repeated runs of the same
     program.  If you use the same set of parameters, seed, and rng stream, you
     will generate the same function.  Thus, in order to get a variety of 
     functions, it is IMPORTANT to use the special constructor and feed in
     varying seeds and streams.  This is easy to do.  For example:

     randfunc TestFunction( 2306, 14)

     uses the seed 2306, and stream 14.  The seed can be any positive, long integer,
     while the stream must be a positive integer between 0 and 255.


     For a description of the necessary inputs, see the header for the setvalues()
     function below.
     
     I have included a non-member function scanfunction() that takes an ofstream
     and a randfunc object and outputs a thorough sample of data points from
     the function which can be used for graphing.  scanfunction is only designed
     for 2-D functions.

     The class defaults to using a Gaussian isotropic correlation function, but
     other choices are available through the chooseCorrelationFamily(int choice)
     function.  If you change the correlation family after inputing parameters,
     you will need to run newrand() again.  See the description below for more information


     Last Modified 6/23/99
*/
#ifndef KRIG_H
#define KRIG_H

#include <cmath>
#include "maps_general.h"
#include "rvgs.h"
#include <iostream>
#include "gamma.h"

#ifndef BOOL_C
#define BOOL_C
const int TRUE = 1;
const int FALSE = 0;
#endif

using namespace std;

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
#ifndef MACHEPS_C
#define MACHEPS_C
//const double MACHEPS = 2.22045e-16;     // For Suns
const double MACHEPS = 1.0842e-19;        // For Pentiums
//const double MACHEPS = ????             // Fill in your own here

//const double MACHEPS = 1.0;             // If you don't know macheps, choose
                                          // this one, and it will be set
                                          // properly automatically
#endif

class randfunc {
public:
  // Automatically takes seed from system clock.
  // Defaults to stream 42.
  randfunc();

  //  Note:  will copy seed along with other data, but will not re-seed the rng.  
  randfunc( randfunc& ToBeCopied );

  // Use to feed a seed to rng when creating a new class object
  randfunc( long seedtoplant, int stream );

  // Destructor
  virtual ~randfunc(); 

   /* evalf returns the value of the function specified by the parameters
   of the evalf at the point x

   x: A p-dim vector at which f is to be evaluated
   */
  virtual double evalf(const Vector<double>& x);

  // evalderiv returns the vector value of the first derivative of the
  // interpolant at the point passed in.
  //
  // Note: in some cases, the derivative at interpolated points
  // is undefined.  The function will then print an error message
  // and return HUGE_VAL as the derivative  
  //
  // tested for correlation functions 1,2,3
  virtual Vector<double> evalderiv( const Vector<double>& point);

  // Calls krigify to get new values for X and v
  // returns TRUE if successful
  //         FALSE if fails due to setvalues() not having
  //               been called previously
  int newrand();  

  // Calls krigify to get new values for v
  // returns TRUE if successful
  //         FALSE if fails due to setvalues() not having
  //               been called previously
  int newrandsameX();

  /* virtual void setvalues() 
     This is an I/O function to interactively set the parameter values

     The data will be read in the format:
     p = long integer
     n = long integer
     beta0 = double
     beta1 = p doubling-point values separated by spaces
     beta2 = p*p doubling-point values read along the rows

        Ex. the Matlab matrix [1 2; 3 4] would be input
        1 2 3 4

     x0 = p doubling-point values separated by spaces
     alpha = double                 //  NOTE: Will not request alpha or theta  
     theta = double                 //  if using User specified correlation function
     sigma2 = double
     lower = p floating-point values separated by spaces
     upper = p floating-point values separated by spaces
     
     X and v will be produced by the krigify function
     
     NOTE: Will not request beta's or x0 if using User specified trend
     
  */
  virtual void setvalues(); 

  /*  virtual int setvalues(char* filename);
      
      Reads the parameter values in order from a file and runs
      krigify to get X and v.  Returns 1 if successful, and 0 if insufficient
      data in the file.
     
      Format is the same as the above
      
      Ex.  Trosset's sample data on pg 6 of the aforementioned paper would be
      2
      200
      50.0
      0 0
      100 0 0 100
      0.3 0.4
      1
      50
      100
      0 0
      1 1

      NOTE: If using a User specified correlation function, then alpha and theta
      should be omitted, and will be set to 0 automatically.  Thus, the input
      would look like:

      2
      200
      50.0
      0 0
      100 0 0 100
      0.3 0.4
      100
      0 0
      1 1
      
      NOTE: If using a User specified trend, then all beta's and x0 should be
      omitted, Thus the input would look like:
      
      2
      200
      1
      50
      100
      0 0
      1 1

      or if User specifing both correlation function and trend, then would look
      like:

      2
      200
      100
      0 0
      1 1
  */
  virtual int setvalues(char* filename);


  /*virtual int setvalues(long pin, long nin, double beta0in, Vector<double>& beta1in,
			Matrix<double>& beta2in, Vector<double>& x0in, 
			Vector<double>& alphain, Vector<double>& thetain, 
			double sigma2in, Vector<double>& lowerin,
			Vector<double>& upperin);

    Has the same effect as setvalues() except that information is passed as parameters 
    instead of being read from I/O.
  */
  virtual int setvalues(long pin, long nin, double beta0in, Vector<double>& beta1in,
			Matrix<double>& beta2in, Vector<double>& x0in, 
			Vector<double>& alphain, Vector<double>& thetain, 
			double sigma2in, Vector<double>& lowerin,
			Vector<double>& upperin);


  /*  void chooseCorrelationFamily( int choice );
      
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
  void chooseCorrelationFamily( int choice );  

  /*  void chooseCorrelationFamily( void (*corfunc)(const Vector<double>& x,
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
  

  /* void setTrend( void (*usertrend)( const Vector<double>& x,
                                       double& result));
     This allows the user to specify a non-quadratic trend function.  x is the
     point at which the trend should be evaluated, and result is used to return
     the value of the trend at x

     NOTE: To return to using default quadratic trend, as long as setvalues()
     has been executed to input beta's and x0 at some point, you may call setTrend(NULL)
  */
  void setTrend( void (*usertrend)( const Vector<double>& /*input*/ x,
                                    double& /*output*/ result));

  /* void generateTrend( double tau2, double* beta0 );

     This generates a psuedorandom quadratic trend.  If beta0 == NULL, we use a 
     default value of beta0 = 0.  We choose beta1 = 0.
     We generate Y such that Y[i][j] are Normal(0,tau2) 
     and then set beta2 = transpose(Y) * Y * ( 1.0/(p+1)).  
     x0 is Uniformly chosen from [lower+(upper-lower)/4, lower+(upper-lower)/2] 
     in each dimension

     Our final trend is then beta0 + 0 + (x-x0)^T beta2 (x-x0)

     If tau2 == -1, we use tau2 = 4*sigma2

     This function should be run AFTER the values are set.  This means that you should 
     pass in some useless values for the beta's and x0 in the setvalues data, but these
     values will be later ignored.
  */
  void generateTrend(double, double* );


  /* void chooseInitialDesign(int choice);
     This allows the user to select which initial design (how X is chosen) is used.  
     The default is choice == 1, Uniform Random Sampling

     Currently, the options are:

      choice = 0  =>  User specified function (See below)
      
      choice = 1  =>  Uniform Random Sampling

      choice = 2  =>  Latin Hypercubes
  */
  void chooseInitialDesign(int choice);

  /*  void chooseInitialDesign( void (*userdesign)( const long p,
						    const long n,
				                    const Vector<double>& lower,
				                    const Vector<double>& upper,
					            Matrix<double>& X) )

     This allows the user to specify an alternative initial design.  The
     initial design must accept four inputs: 
     p is the dimension of the space
     n is the number of initial design sites (how many points in the space)
     lower is the vector of lower bounds (rectangular spaces)
     upper is the vector of upper bounds

     X is then the matrix in which the points will be stored.  X is an n by p
     matrix, with one point in each row.
  */
  void chooseInitialDesign( void (*userdesign)( const long /*input*/ p,
						const long /*input*/ n,
				     const Vector<double>& /*input*/ lower,
				     const Vector<double>& /*input*/ upper,
					   Matrix<double>& /*output*/ X) );    
 

  // functions to get copies of the data members
  long getp();
  long getn();
  double getbeta0();
  Vector<double> getbeta1();
  Matrix<double> getbeta2();
  Vector<double> getx0();
  double getalpha();
  double gettheta();
  Matrix<double> getX();
  Vector<double> getv();
  Vector<double> getlower();
  Vector<double> getupper();
  long getSeed();
  int getStream();
  
protected:

  /* This returns the value of the derivative of |a-b| with respect to b,
   i.e. 1.0 or -1.0 . */
  double derivabs( double a, double b);
 
  /* The value of the first derivative of the trend at point */
  Vector<double>* derivtrend( const Vector<double> & point);
 
  // generates a Latin Hypercube design and stores it in X.
  void genLH();

  // Prints out an error message
  void notseterror();

  // deletes all non-null dynamically allocated data members
  void deleteall(); 
  
  // Generates the necessary random values and kriges them to produce
  // an approximate realization
  //
  // flag == 1  => create new X and v
  // flag == 2  => create new v using old X
  void krigify(int flag); 

  /**** virtual double correlation( const Vector<double>& x, const Vector<double>& y);

       INPUT: Two points x and y in the space
       OUTPUT: The value of the correlation between them.  Currently, this uses the
       Gaussian correlation function, but plans are to later allow the user
       to select one of several correlation functions at the time of input
       EFFECT: Calls l2norm() on the difference (ie. y - x ) of the vectors and then
       calculates the Gaussian correlation.
  ****/  
  virtual double correlation( const Vector<double>& x, const Vector<double>& y);
 

  int ValuesSet;
  long p;        // dimension of the space
  long n;        // number of sample sites
  long seed;
  int rngstream;

  // flag to choose correlation family
  // 1 => Gaussian Isotropic  r(x,y) = exp( -theta * ||y-x||^alpha )
  // 2 => Gaussian Product    r(x,y) = \prod_{j=1}^{P} exp(-theta * |s_j - t_j|^alpha)
  // 3 => Product of Linear   r(x,y) = \prod_{j=1}^{p} ( 1-theta*|y_j - x_j|)+
  // 4 => Cubic correlation on the unit cube
  //      r(x,y) = \prod_{j=1}^{p} ( 1 - theta*(y_j - x_j)^2 + alpha*(y_j - x_j)^3 )
  int corfamily;

  // flag to choose the initial design
  // 1 => Uniform Random Sampling
  // 2 => Latin Hypercube
  int designchoice;

  void (*corfunction)( const Vector<double>& x, const Vector<double>& y, double& result);
  void (*trend)(const Vector<double>& x, double& result);
  void (*design)( const long /*input*/ p,
		      const long /*input*/ n,
		      const Vector<double>& /*input*/ lower,
		      const Vector<double>& /*input*/ upper,
		      Matrix<double>& /*output*/ X);

  Vector<double> *thetavec;  // Vectors for storing the parameters of the 
  Vector<double> *alphavec;  // alternate correlation functions. 

  // these four parameters specify the quadratic trend
  // beta0 + (x-x0)^{T} beta1 + (x-x0)^{T} beta2 (x-x0)
  double beta0; 
  Vector<double> *beta1; 
  Matrix<double> *beta2; 
  Vector<double> *x0;    
   
  double alpha;  // specify the correlation function
  double theta;  // specify the correlation function
  double sigma2; // Variance
  
  Vector<double> *lower; // lower bound for each dimension
  Vector<double> *upper; // upper bound for each dimension
  Matrix<double> X;  // n random p-dimensional points (resized as necessary)
  Vector<double> *v;  // function values for the points in X	
}; // end class randfunc


/**** 
      void scanfunction( randfunc& Func, ostream& outf, int pointsperaxis);

      INPUT: The randfunc object to be sampled and an ofstream to output the
      values, as well as the number of points per axis to be sampled.  Thus, the
      total number of sample points is approximately (pointsperaxis)^p.  Should
      work for any value of p.
      EFFECT: Outputs to the stream values in columns by axis.  For example:

      2.3 3.4 13.695
      
      would correspond to the point (2.3, 3.4, 13.695) on the Cartesian
      coordinate system.
****/
void scanfunction( randfunc& Func, ostream& outf, int pointsperaxis);


// A default version which uses 25 pointsperaxis.  You can modify this default
// by changing the 25 above.
void scanfunction( randfunc& Func, ostream& outf);


/****
     double reportsample( randfunc& Func, Vector<double>& x, ostream& outf);

     INPUT: The randfunc object to be sampled, the point at which to evaluate
     the function, and the output stream to which it should report. 
     OUTPUT: the value of the function at x - result of Func.evalf(x)
     EFFECT: Outputs to the stream the elements of x, seperated by spaces, and
     then the function value at x, followed by an endline. For example:

     1.0 2.0 3.0 4.0 5.0 6.0 21.0

     would be written if x == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0) and f(x) == 21.0
****/
double reportsample( randfunc& Func, Vector<double>& x, ostream& outf);


/**** 
      void scanplane( randfunc& Func, ostream& outf, int pointsperaxis);

      A simplified, but faster version of scanfunction().
      ONLY works when p=2 !!!
      

      INPUT: The randfunc object to be sampled and an ofstream to output the
      values.
      EFFECT: Outputs to the stream values in columns by axis.  For example:

      2.3 3.4 13.695
      
      would correspond to the point (2.3, 3.4, 13.695) on the Cartesian
      coordinate system. 
****/
void scanplane( randfunc& Func, ostream& outf, int pointsperaxis);


/****
     void analyzevalues( randfunc& Func, ostream& outf);

     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Uses a default value of 25 points per axis.
     
     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.
****/
void analyzevalues( randfunc& Func, ostream& outf);

/****
     void analyzevalues( randfunc& Func, ostream& outf, int pointsperaxis);

     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Divides each axis into pointsperaxis evenly spaced points, and samples at
     each permutation of values.
     
     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.
****/
void analyzevalues( randfunc& Func, ostream& outf, int pointsperaxis);

#endif






