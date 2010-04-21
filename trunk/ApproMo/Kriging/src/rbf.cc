/*  rbf.cc
 written by Anthony Padula	adpadu@maila.wm.edu
This file contains the function implementations of the setvalues() and evalf()
functions to create and evaluate a radial basis function approximation.
These functions are encapsulated inside the implementation of the rbfapprox
class.

Last Modified 7/21/99
*/
#include "rbf.h"
#include <cassert>

rbfapprox::rbfapprox() { 

   ValuesSet = FALSE;
   rbffamily = 7;  // Default to choice 7 which requires no parameters
   phi = NULL;
   p = 0;
   n = 0;
   c_ = 0;
   beta_ = 0;
   degree = 0;
}

rbfapprox::rbfapprox( rbfapprox& ToBeCopied )   /* Last Modified 7/14/99 */
{
   ValuesSet = ToBeCopied.ValuesSet;
   if( ValuesSet) {
      p = ToBeCopied.p;
      n = ToBeCopied.n;
      c_ = ToBeCopied.c_;
      beta_ = ToBeCopied.beta_;
      rbffamily = ToBeCopied.rbffamily;
      phi = ToBeCopied.phi;
      degree = ToBeCopied.degree;
      
      X.newsize(n,p);
      y.newsize(n);

      X = ToBeCopied.X;
      y = ToBeCopied.y;
   } // end if
} // end Copy constructor

rbfapprox::rbfapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
                      Vector<double>& Values,  // The vector of function values,
                                               // so f(Points[i])=Values[i]
                      long deg,  // The number of extra polynomials to use
                      // deg < 0  no polynomial basis
                      // deg = 0  constant 1 only
                      // 0<deg<=p+1 use the first deg monomials and a constant 1 
                      int choice, // The flag for which RBF to use (see chooseRBF() below)
                      double c,   // c and beta are the parameters for the RBF 
                      double beta) // You may not need both for all choices, but always
{                                   // pass in something at least
  // Last modified 7/14/99
  ValuesSet = TRUE;
  rbffamily = choice;
  c_ = c;
  beta_ = beta; 
  phi = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  degree = deg;
  approximate();
}

rbfapprox::~rbfapprox() {    /* Last Modified 7/14/99 */
  phi = NULL;
} // end destructor

int rbfapprox::setvalues() {  /* Last Modified 7/14/99 */
  long i, j;
  ValuesSet = TRUE;
  cout << "Please input values for:" << endl << "p = ";
  cin >> p;
  cout << " n = ";
  cin >> n;
  cout << " degree = ";
  cin >> degree;
  cout << " choice = ";
  cin >> rbffamily;
  
  switch( rbffamily ) {
  case 2:
  case 3:
  case 5:
  case 6:
    cout << " c = ";
    cin >> c_;
    break;
  default: c_ =0;
  }; // end c switch
  switch( rbffamily ) {
  case 1:
  case 3:
  case 4:
    cout << " beta = ";
    cin >> beta_;
    break;
  default: beta_ = 0;
  }; // end beta switch
  X.newsize(n,p);
  y.newsize(n);
  cout << " Points (by rows) = ";
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      cin >> X[i][j];
    } // end for
    
    cin >> y[i];
  } // end for
   
  cout << "Input received.  Calculating approximation..." << endl;
  if( approximate() == TRUE) {
    cout << "Approximation finished." << endl;
    return TRUE;
  }
  else {
    cout << "Approximation failed." << endl;   
    return FALSE;
  }
} // end setvalues()

int rbfapprox::setvalues(char* filename) {  /* Last Modified 7/14/99 */
  ifstream infile;   
  infile.open( filename, ios::in );
  long i, j;
  
  infile >> p >> n >> degree >>rbffamily;
  
  switch( rbffamily ) {
  case 2:
  case 3:
  case 5:
  case 6:
    infile >> c_;
    break;
  default: c_ = 0;
  }; // end c switch
  switch( rbffamily ) {
  case 1:
  case 3:
  case 4:
    infile >> beta_;
    break;
  default: beta_ = 0;
  }; // end beta switch

  if( infile.eof()) {
    cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
         << endl << " Quitting setvalues()" << endl;
    return FALSE;
  } // end if
  X.newsize(n,p);
  y.newsize(n);
  // Now read in X and y
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      infile >> X[i][j];
    } // end for
    
    infile >> y[i];
  } // end for

  if( infile.eof()) {
    cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
         << endl << " Quitting setvalues()" << endl;
    return FALSE;
  } // end if

   ValuesSet = TRUE;
   if( !approximate() )
      cout << "Approximation failed." << endl;
   infile.close();
   return TRUE;
} // end setvalues(filename)

int rbfapprox::setvalues(Matrix<double>& Points, // The n x p matrix of points, one point to a row
                         Vector<double>& Values,  // The vector of function
                                                  // values, so f(Points[i])=Values[i]
                         long deg,    // The size of the polynomial basis
                         // deg < 0  no polynomial basis
                         // deg = 0  constant 1 only
                         // 0<deg<=p+1 use the first deg monomials and a constant 1 
                         int choice, // The flag for which RBF to use (see chooseRBF() below)
                         double c,   // c and beta are the parameters for the RBF 
                         double beta ) // You may not need both for all choices, but always
                                       // pass in something at least
  /* This function works like the special constructor, allowing the data to be
     passed in as parameters instead of being read from I/O.
     Stores all data and calls approximate() to create the approximation
  */
{
  // Last modified 7/15/99
  ValuesSet = TRUE;
  rbffamily = choice;
  c_ = c;
  beta_ = beta; 
  phi = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  degree = deg;
  return approximate();
}// end setvalues(data)


// functions to access the data members


double rbfapprox::getbeta() {
   if ( !ValuesSet )
        notseterror();
   return beta_;
}

double rbfapprox::getc() {
    if ( !ValuesSet )
        notseterror();
    return c_;
}

int rbfapprox::approximate() { /* Last modified 7/14/99 */
  long s = n;
  Matrix<double> A(1,1);
  Vector<double> diff(p);
  Vector<double> tempv(1);
  Vector<double> right(1);
  double temp;
  long i,j, info;

  if( degree >= 0 ) {
    s += degree+1;
    a.newsize(degree+1);
  }
  else {
    a.newsize(1);
    a[0] = 0;
  }
  
  b.newsize(n);
  A.newsize(s,s);
  tempv.newsize(s);
  right.newsize(s);
  
  for(i=0; i<n; i++) {  // Fill in G
    for(j=0; j<n; j++) {
      diff = X.row(i) - X.row(j);
      rbfvalue(diff.l2norm(), temp);
      A[i][j] = temp;
    } // end for
    right[i] = y[i];
  } // end for

  if (degree >= 0 ) {
    for( i=0; i<n; i++) { // Fill in F and FT
      for( j=n; j<s; j++) {
        temp = polybasis(X.row(i), j-n);
        A[i][j] = temp;
        A[j][i] = temp;
      } // end for
    }// end for
    for( i=n; i<s; i++) { // Fill in rest with zeros
      for( j=n; j<s; j++) {
        A[i][j] = 0;
      } // end for
      right[i] = 0;
    } // end for
  } // end if
  
  // Solve using LAPACK's solvers
    info = SolveSystem( A, tempv, right);
    if ( info != 0 ) {
      cout << "U(" << info <<", " << info << ") in the LU factorization is "
           << "exactly zero, so the solution could not be computed." << endl;
      return FALSE;
    }
    
    // check solution:
    if (DISPLAY > 1 )
      cout << "residual = " << (right - A*tempv) << endl;
    
    for(i = 0; i<n; i++)
      b[i] = tempv[i];
    if( DISPLAY > 0 )
      cout << "b = " << b << endl;
    if (degree >= 0 ) {
      for(j = n; j<s; j++)
        a[j-n] = tempv[j];
      if( DISPLAY > 0 )
        cout << "a = " << a << endl;
    } // end if
    
    return TRUE;
} // end approximate()

double rbfapprox::evalf(const Vector<double>& point) { /* Last modified 7/14/99 */

  double value = 0.0;
  double temp=0.0;
  Vector<double> tempv(p);
  
  if ( !ValuesSet )
    notseterror();

  for(long i=0; i<n; i++) {
    tempv = point - X.row(i);
    rbfvalue(tempv.l2norm() , temp);
    value += b[i]*temp;
  } // end for

  for( long j=0; j<=degree; j++) {
    value += a[j]*polybasis(point, j);
  } // end for
   return value;
} // end evalf()

double rbfapprox::rbfvalue( const double r, double& result )
{
  /* Last modified 7/15/99 */
  double tempd = 0, tempd2 = 0;
  int nu = int(c_ - p/2);
  double k[nu+2];
  switch( rbffamily ) {
  case 1:  // Thin-Plate r^beta
    if ( r == 0 )
      result = 0;
    else
      result = pow(r,beta_);
    break;
  case 2: // Gaussian exp(-c*r^2)
    if( r==0)
      result = 1;
    else
      result = exp( (-c_)*(r*r)); 
    break;
  case 3: // Multiquadric (c^2 + r^2)^(beta/2)
    if( r == 0 )
      result = c_;
    else {
      tempd = (c_*c_ + r*r);
      tempd2 = beta_/2.0;
      if( beta_ == 1 )
        result = sqrt(tempd);
      else if (beta_ == -1)
        result = 1.0/sqrt(tempd);
      else
        result = pow(tempd , tempd2);
    } // end else
    break;
  case 4: // Thin-Plate splines (-1)^(beta/2+1)*r^(beta)*log(r)
    if( r == 0 )
      result = 0;
    else {
      tempd = pow(r, beta_);
      tempd = tempd * log(r);
      result = tempd * pow(-1.0, beta_/2.0+1);
    } // end else
    break;
  case 5: // Sobolew splines (2*PI^(c))*K(nu,2*PI*r)*(r^nu)/gamma(c)
    if( r == 0 )
      result = 0;
    else {
      tempd = 2 * pow(PI, c_);
      bessk(2*PI*r, nu, k);
      tempd = tempd * k[nu];
      tempd = tempd * pow(r, nu);
      result = tempd / gamma(c_);
    } // end else
    break;
  case 6: // Thin-plate splines (cr)^2 * log( cr )
    if( r == 0 )
      result = 0;
    else {
      tempd = (c_*r);
      tempd = tempd * tempd;
      result = tempd * log( (c_*r) );
    } // end else
    break;
  case 7: // Madych & Nelson II pg 221 -2*sqrt(pi*(1+r^2))
    if( r == 0 )
      result = -2*sqrt(PI);
    else {
      tempd = PI * (1 + r*r);
      result = -2 * sqrt(tempd);
    } // end else
    break;
  case 0: phi(r,result);
  }; // end switch
  return result;
} // end rbfvalue()


double rbfapprox::polybasis( Vector<double> x, long k) {
  if ( k < 0)
    return 0;
  else if( k == 0)
    return 1;
  else
    return x[k-1];
} // end polybasis  


void rbfapprox:: chooseRBF( int choice, double c, double beta, long deg ) {
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
  // Last modified 7/15/98
  rbffamily = choice;
  c_ = c;
  beta_ = beta;
  degree = deg;
  if( ValuesSet )
    approximate();
  return;
 }

void rbfapprox::chooseRBF( void (*RBFtoUse)(double r, double& ), long deg ) {
  // Last modified 7/15/99
  rbffamily = 0;
  phi = RBFtoUse;
  degree = deg;
  if( ValuesSet)
    approximate();
  return;
} // end chooseRBF()  


void rbfapprox::storeApprox( ostream& outf ){
/****
     Prints the data of the rbf to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing RBF's so that they might be re-created later.
****/
  long i, j;
  outf << p << endl << n << endl << degree << endl << rbffamily << endl;
  switch( rbffamily ) {
  case 2:
  case 3:
  case 5:
  case 6:
    outf << c_<< endl;
    break;
  }; // end c switch
  switch( rbffamily ) {
  case 1:
  case 3:
  case 4:
    outf << beta_<< endl;
    break;
  } // end beta switch
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      outf << X[i][j] << ' ';
    } // end for
    
    outf << y[i]<< endl;
  } // end for
  return;
} // end storeApprox()







