/*  krig.cc
 written by Anthony Padula	adpadu@maila.wm.edu
This file contains the function implementations of the setvalues() and evalf()
functions to create and evaluate a kriging approximation.
These functions are encapsulated inside the implementation of the krigapprox
class.

Last Modified 6/10/00
Minimal modification 6/17/03 by Anne Shepherd to add sanity 
   check for proper version of cppmat.h to all constructors.

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
#include "krig.h"
#include <cassert>
#include <vector>

krigapprox::krigapprox() { 

#ifndef _changed_f_indexing_to_c_pls
    cerr<<"\nERROR:  using outdated version of cppmat.h---this will";
    cerr<<"\ngive incorrect results.  Use the updated version, which";
    cerr<<"\nchanges fortran-style indexing to C-style in all cases.\n\n";
    exit(1); 
#endif
  ValuesSet = FALSE;
  corfamily = 1;  // Default to Gaussian isotropic
  corfunction = NULL;
  p = 0;
  n = 0;
  theta_ = NULL;
  mu_ = NULL;
}

krigapprox::krigapprox( krigapprox& ToBeCopied )   /* Last Modified 7/14/99 */
{
#ifndef _changed_f_indexing_to_c_pls
    cerr<<"\nERROR:  using outdated version of cppmat.h---this will";
    cerr<<"\ngive incorrect results.  Use the updated version, which";
    cerr<<"\nchanges fortran-style indexing to C-style in all cases.\n\n";
    exit(1); 
#endif
  ValuesSet = ToBeCopied.ValuesSet;
  if( ValuesSet) {
    p = ToBeCopied.p;
    n = ToBeCopied.n;
    corfamily = ToBeCopied.corfamily;
    corfunction = ToBeCopied.corfunction;
    X.newsize(n,p);
    y.newsize(n);
    
    X = ToBeCopied.X;
    y = ToBeCopied.y;
    
    if( ToBeCopied.theta_ != NULL ) {
      new_array(theta_, thetasize());
      for( int i =0; i < thetasize(); i++)
        theta_[i] = ToBeCopied.theta_[i];
    } else
      theta_ = NULL;
    
    if( ToBeCopied.mu_ != NULL ) {
      new_array(mu_, musize());
      for( int j = 0; j < musize(); j++)
        mu_[j] = ToBeCopied.mu_[j];
    } else
      mu_ = NULL;
      
  } // end if
} // end Copy constructor

krigapprox::krigapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
                        Vector<double>& Values,  // The vector of function values,
                                                 // so f(Points[i])=Values[i]
                         int choice) // The flag for which correlation function to use
                          // (see chooseCorrelation() below)
                        {
  // Last modified 7/22/99
#ifndef _changed_f_indexing_to_c_pls
    cerr<<"\nERROR:  using outdated version of cppmat.h---this will";
    cerr<<"\ngive incorrect results.  Use the updated version, which";
    cerr<<"\nchanges fortran-style indexing to C-style in all cases.\n\n";
    exit(1); 
#endif
  ValuesSet = TRUE;
  corfamily = choice;
  corfunction = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  theta_ = NULL;
  mu_ = NULL;
  approximate();
}

        
krigapprox::krigapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
              Vector<double>& Values,  // The vector of function values, so
                                       // f(Points[i])=Values[i] 
              int choice, // The flag for which correlation function to use
                          // (see chooseCorrelation() below)
              double* theta) // The known value(s) of the parameter(s)
  // theta will either be a single number or an array of numbers.  If
  // there are multiple parameters,
  // they are expected to be stored in the array in order of increasing
  // dimension, and all theta_i before alpha_i.  Thus, in 3 dimensions, the
  // array might contain [1.0, 2.0, 3.0, 0.5, 4.5, 0] if theta_1 = 1.0, 
  // theta_2 = 2.0, theta_3 = 3.0, alpha_1 = 0.5, alpha_2 = 4.5, alpha_3 = 0.0
  // See the chooseCorrelation() function below for more discussion of parameters
{
 // Last modified 7/22/99
#ifndef _changed_f_indexing_to_c_pls
    cerr<<"\nERROR:  using outdated version of cppmat.h---this will";
    cerr<<"\ngive incorrect results.  Use the updated version, which";
    cerr<<"\nchanges fortran-style indexing to C-style in all cases.\n\n";
    exit(1); 
#endif
  ValuesSet = TRUE;
  corfamily = choice;
  corfunction = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  if( thetasize() > 0 ) {
    
    new_array( theta_,thetasize());
    for(int i= 0; i<thetasize(); i++)
      theta_[i] = theta[i];
  }
  else
    theta_ = NULL;

  mu_ = NULL;
  approximate();
}


krigapprox::~krigapprox() {    /* Last Modified 7/14/99 */
  if( DISPLAYK > 3 )
    cout << "deleteing" << endl;
  corfunction = NULL;
  
  deleteall();
  if( DISPLAYK > 3 )
    cout << "done deleting" << endl;
} // end destructor

// end of all constructors and destructors


void krigapprox::deleteall() {
  if( theta_ != NULL ) {
    delete [] theta_;
    theta_ = NULL;
  }
  if( mu_ != NULL) {
    delete mu_;
    mu_ = NULL;
  }
} // end deleteall()



// The setvalues() functions

int krigapprox::setvalues() {  /* Last Modified 7/14/99 */
  long i, j;
  char answer = 0;
  deleteall();
  ValuesSet = TRUE;

  cout << "Please input values for:" << endl << "p = ";
  cin >> p;
  cout << " n = ";
  cin >> n;
  cout << " choice = ";
  cin >> corfamily;
  if (corfamily > 0) {
    do {
      cout << "Are the values of the parameter(s) known? (y/n):";
      cin >> answer;
    } while ((answer != 'y')&&(answer != 'n'));
  }
  else // corfamily == 0
    answer = 'n';
  
  if( answer == 'n' )
    theta_ = NULL;
  else {
    new_array( theta_,thetasize());
    switch( corfamily ) {
    case 1: cout << " theta = ";
      cin >> theta_[0];
      cout << " alpha = ";
      cin >> theta_[1];
      break;
    case 2: cout << " alpha = ";
      cin >> theta_[p];
    case 3: cout << " Input the p-vector of theta's: ";
      for( i = 0; i<p; i++)
        cin >> theta_[i];
      break;
    case 4:  
    case 5: cout << " Input the p-vector of theta's: ";
      for( i = 0; i<p; i++)
        cin >> theta_[i];
       cout << " Input the p-vector of alpha's: ";
      for( i = p; i<2*p; i++)
        cin >> theta_[i];
      break;  
    default: theta_=NULL;
    } // end input switch
  } // end else
  X.newsize(n,p);
  y.newsize(n);
  cout << " Points (by rows) = ";
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      cin >> X[i][j];
    } // end for
    
    cin >> y[i];
  } // end for
  if( DISPLAYK > 1 ) 
    cout << "Input received.  Calculating approximation..." << endl;
  
  if( approximate() == TRUE) {
    if( DISPLAYK > 0 )
      cout << "Approximation finished." << endl;
    return TRUE;
  }
  else {
    if( DISPLAYK > 0 )
    cout << "Approximation failed." << endl;   
    return FALSE;
  }
} // end setvalues()
  
int krigapprox::setvalues(char* filename) {  /* Last Modified 7/14/99 */
  ifstream infile;   
  infile.open( filename, ios::in );
  long i, j;
  char answer = 0;

  deleteall();
  infile >> p >> n >> corfamily;
  if (corfamily > 0) {
      infile >> answer;
  }
  else // corfamily == 0
    answer = 'n';
  
  if( answer == 'n' ) {
    theta_ = NULL;
  }
  else if ( answer == 'y') {
    new_array(theta_,thetasize());
    switch( corfamily ) {
    case 1:
      infile >> theta_[0];
      infile >> theta_[1];
      break;
    case 2:
      for( i = 0; i<p; i++)
        infile >> theta_[i];

      infile >> theta_[p];
      break;
    case 3:
      for( i = 0; i<p; i++)
        infile >> theta_[i];
      break;
    case 4:  
    case 5:
      for( i = 0; i<2*p; i++)
        infile >> theta_[i];
      break;  
  default: theta_=NULL;
    }; // end input switch
  } // end if
  else { // answer screwed up
    if( DISPLAYK > 0 )
      cout << "ERROR: INPUT FILE MISTAKE.  LINE 4 should be 'y' or 'n'" << endl;
    ValuesSet = FALSE;
    return FALSE;
  } // end else

  if( infile.eof()) {
    if( DISPLAYK > 0 )
      cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
           << endl << " Quitting setvalues()" << endl;
    return FALSE;
  } // end if
  
  X.newsize(n,p);
  y.newsize(n);
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      infile >> X[i][j];
    } // end for
    
    if( infile.eof()) {
      if( DISPLAYK > 0 )
        cout <<"Error: File terminated prematurely.  Check that all vectors are the correct length"
             << endl << " Quitting setvalues()" << endl;
      return FALSE;
    } // end if
  
    infile >> y[i];
  } // end for
  
   ValuesSet = TRUE;
   if( !approximate() )
     if( DISPLAYK > 0 )
       cout << "ERROR: Approximation failed. Check input file" << endl;
   infile.close();
   return TRUE;
} // end setvalues(filename)

int krigapprox::setvalues( Matrix<double>& Points,// The n x p matrix of points, one point to a row
                        Vector<double>& Values,  // The vector of function values,
                                                 // so f(Points[i])=Values[i]
                         int choice) // The flag for which correlation function to use
                                      // (see chooseCorrelation() below)
 {
  // Last modified 7/22/99
  deleteall();
  ValuesSet = TRUE;
  corfamily = choice;
  corfunction = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  theta_ = NULL;
  mu_ = NULL;
  return approximate();
}

        
int krigapprox::setvalues( Matrix<double>& Points,//The n x p matrix of points, one point to a row
              Vector<double>& Values,  // The vector of function values, so
                                       // f(Points[i])=Values[i] 
              int choice, // The flag for which correlation function to use
                          // (see chooseCorrelation() below)
              double* theta) // The known value(s) of the parameter(s)
  // theta will either be a single number or an array of numbers.  If
  // there are multiple parameters,
  // they are expected to be stored in the array in order of increasing
  // dimension, and all theta_i before alpha_i.  Thus, in 3 dimensions, the
  // array might contain [1.0, 2.0, 3.0, 0.5, 4.5, 0] if theta_1 = 1.0, 
  // theta_2 = 2.0, theta_3 = 3.0, alpha_1 = 0.5, alpha_2 = 4.5, alpha_3 = 0.0
  // See the chooseCorrelation() function below for more discussion of parameters
{
 // Last modified 7/22/99
  deleteall();
  ValuesSet = TRUE;
  corfamily = choice;
  corfunction = NULL;
  X = Points;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  if( thetasize() > 0 ) {
    new_array( theta_,thetasize()); 
    for(int i= 0; i<thetasize(); i++)
      theta_[i] = theta[i];
  }
  else
    theta_ = NULL;
  mu_=NULL;
  return approximate();
}


// end of setvalues() functions


// functions to access the data members

Vector<double>* krigapprox::gettheta() {
  
   if ( !ValuesSet )
     notseterror();
   if( theta_ == NULL)
     return NULL;
   
   Vector<double>* theta = NULL;
   new_array( theta,thetasize());
   for( int i = 0; i < thetasize(); i++)
     (*theta)[i] = theta_[i]; 
   return theta;
}


Vector<double> krigapprox::getv(){
    if ( !ValuesSet )
        notseterror();

    return v;
}


// approximate, evaluate and cohorts


int krigapprox::approximate() { /* Last modified 7/14/99 */
 
  Matrix<double> R(n,n);
  Matrix<double> RINV(n,n);
  Vector<double> r(n);
  Vector<double> randn(n);
  Vector<double> temppointi(p);
  Vector<double> muv(n);
  long i, j;
  double tol;
  bool success;

  if( ValuesSet == FALSE)
    return FALSE;
  
  if (theta_ == NULL)
    estimateTheta();
  
  // Creat v using X and y 
  for( i = 0; i < n; i++)
    R[i][i] = 1;
  if( DISPLAYK > 3 )
    cout << "R = I" << endl;
  for( i=1; i < n; i++) { // Compute the inter-point correlations
    temppointi = X.row(i); 
    for( j=0; j < i; j++) {
      R[i][j] = correlation( temppointi, X.row(j) );
      R[j][i] = R[i][j];
    } // end for
  } // end for

  if( DISPLAYK > 3 )
    cout << "Correlation matrix computed. " << endl;
  if( DISPLAYK > 4 )
    cout << R << endl;

  if (mu_ == NULL)
    estimateMu(R);

  muv = mu_[0];

  success = ComputeInverse(R, (&RINV)); // Calculate the Cholesky Decomposition
  
  int count = 0;
  while( (!success)&&(count < MAXCHOLREPS) ) {
    // add a constant to the diagonal
    // from the Higham book on error analysis
    double TOADD = n*(n*MachineEpsilon())/(1-n*MachineEpsilon());
    
    for( int i = 0; i < n; i++)
      R(i,i) += TOADD;
    
    count++;
    if( DISPLAYK > 1)
      cout << "Looped " << count << endl;
    // Calculate the Cholesky Decomposition
    success = ComputeInverse(R, (&RINV));
    }
  

  if(success)
    {
      v = (RINV * (y-muv));
    } // end if
  else {
    if( DISPLAYK > 0)
      cout << "Cholesky failed:  Correlation matrix may be positive definite."
	   << endl << "Trying SVD...";
    
    Matrix<double> U(n,n), VT(n,n), D(n,n);
    ComputeSVD(R,(&U), (&D), (&VT)); // Calculate the SVD
    if( DISPLAYK > 3 )
      cout << "SVD computed" << endl;
  
    // Set appropriate tolerance
    tol = n*D(0,0)*MachineEpsilon();
    
    if( DISPLAYK > 1 )
      cout << "tol = " << tol << endl;
    for( i = 0; i < n; i++) { 
      if (fabs(D(i,i)) < tol)
	D(i,i) = 0;
      else
	D(i,i) = 1.0/D(i,i);
    } // end for
    if( DISPLAYK > 3 )
      cout << " D^(-1) computed" << endl;

    v = ((U * D * VT) * (y-muv));  

    if( DISPLAYK > 1 )
      cout << "mu = " << mu_[0] << ", theta = " << theta_[0] << endl;

    if( DISPLAYK > 3 )
      cout << "v done."  << endl;
  } // end else for SVD
  
  return TRUE;
} // end approximate()

double krigapprox::evalf(const Vector<double>& point) 
{ /* Last modified 7/22/99 */
   double value = 0;
   Vector<double> r(n);
    
   if ( !ValuesSet )
      notseterror();
  
   for(int i=0; i<n; i++) { 
      r[i] = correlation(X.row(i), point);
   } // end for
   if( DISPLAYK > 4 )
     cout << "r = " << r;
   
   value = trend(point)+ v*r;
   return value;
} // end evalf()

Vector<double>* krigapprox::derivtrend( const Vector<double>& point)
  /* For now, with a constant trend, the derivative is zero */
{  
  Vector<double>* value;
  long p;
  p = point.size();
  new_Vector( value, p);
  (*value) = 0; 
  return value;
}

/* This returns the value of the derivative of |a-b| with respect to b,
   i.e. 1.0 or -1.0 . */

double krigapprox::derivabs( double a, double b ) {
  if( a > b )
    return -1.0;
  else
    return 1.0;
}

Vector<double> krigapprox::evalderiv( const Vector<double>& point) {
  /* Last modified 2/14/00 */
  
   double temp, t2;
   Vector<double> r(n), diff(p);
   Vector<double> Dprime(n);
   Vector<double> valuevec(p), * muprime;
   int alpha;

   if ( !ValuesSet )
      notseterror();

   if( (corfamily == 2)
       ||(corfamily==3)
       ||(corfamily==5)
       ||((corfamily==1)&&(theta_[1] == 1.0)) ) {
     // If necessary, check to see if point is a known data point 
     for( int i  =0; i < n; i++)
       if( point == X.row(i) ) {
	 cout << "Not differentiable at this point.  Returning INF" << endl;
	 valuevec=HUGE_VAL;
	 return valuevec;
       }
   }

   muprime = derivtrend(point);
   for(int k=0; k < p; k++) {
     switch( corfamily ) {
     case 1: /* Gaussian isotropic */
       if( theta_[1] == 2) {  // alpha == 2
	 for(int i=0; i<n; i++) { // create delD/delx_k
	   diff = X.row(i)-point;
	   r[i] = correlation(X.row(i), point);
	   Dprime[i] = 2.0*r[i]*theta_[0]*diff[k];
	 }
       }
       else if(theta_[1]==1) { //alpha == 1     
	 for(int i = 0; i < n; i++) { //create delD/delx_k
	   diff = X.row(i)-point;
	   r[i] = correlation(X.row(i), point);
	   Dprime[i] = theta_[0] * r[i]*diff[k] / (diff.l2norm()) ;
	 }  
       }
       break;
     case 2: /* Gaussian product */
       for( int i =0; i < n; i++) {
	 diff = (X.row(i) - point);
	 r[i] = correlation(X.row(i), point);
	 if( theta_[p] == 1) 
	 Dprime[i] = r[i]*(theta_[k]
			   *(derivabs(X(i,k), point[k])));
	 else if (theta_[p] == 2)
	   Dprime[i] = r[i]*(theta_[k]* 2.0*(fabs(diff[k])))
			     *(derivabs(X(i,k), point[k]));
	 else
	   Dprime[i] = r[i]*(theta_[k]* theta_[p] *(pow(fabs(diff[k]), (theta_[p]-1.0)))
			   *(derivabs(X(i,k), point[k])));
       }
       break;
     case 3: /* Product of linear */
       for( int i =0; i < n; i++) {
	 temp = fabs(X(i,k) - point[k]);
	 r[i] = correlation(X.row(i) , point);
	 Dprime[i] = r[i]*(theta_[k]*derivabs(X(i,k),point[k])) / (1.0-theta_[k]*temp); 
       }
       break;
       /*case 4: // Cubic correlation
       for( int i = 0; i < n; i++) {
	 temp = X[i][k] - point[k];
	 r[i] = correlation(X.row(i), point);
	 Dprime[i]=r[i]*(temp*(2.0*theta_[k] - 3.0*theta_[k+p]*temp) / 
			(temp*temp*(1.0-theta_[k]+theta_[k+p]*temp)));
       }
       break;*/
     case 5: /* Mate'rn */
       for( int i =0; i < n; i++) {
	 t2 = theta_[k];
	 alpha = (int) theta_[k+p];
	 double* Karray = (double*)calloc(alpha+2,sizeof(double));
	 diff = (X.row(i)-point);
	 r[i] = correlation(X.row(i),point);
	 temp = t2*fabs(diff[k]);
	 bessk(temp, alpha, Karray);
	 Dprime[i] = derivabs(point[k], X(i,k))*r[i]*alpha / fabs(diff[k]);
	 Dprime[i] +=( r[i]*( Karray[alpha+1]*temp-alpha*Karray[alpha])
		       *
		       derivabs(X(i,k), point[k])
		       / 
		       ( fabs(diff[k]) * Karray[alpha]));
	 free (Karray);
	 
       }
       break;
     case 0:
       cout << "Cannot find derivative of user defined function.  Returning 0";
     default:
       return 0;
     } // end switch
     valuevec[k] = (*muprime)[k]+ v*Dprime;
   } // end for
   delete muprime;
   return valuevec;
}

double krigapprox::correlation( const Vector<double>& x, const Vector<double>& y )
{
  /* Last modified 7/22/99 */
  long size = 0;
  Vector<double> tempv(1);
  double tempd = 0, tempd2 = 0, tempd3 = 0, result=1;
  double *k = NULL;

  switch( corfamily ) {
  case 1: // Gaussian isotropic
    tempv.newsize(y.size());
    tempv = y-x;
    tempd = tempv.l2norm();
    tempd = pow(tempd,theta_[1]);
    tempd = (-theta_[0])*tempd;
    result = exp(tempd);
    break;
  case 2: // Gaussian product
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd = pow(tempd, theta_[p]);
      tempd = (-theta_[i])*tempd;
      tempd = exp(tempd);
      result = tempd*result;
    } // end for
    break;
  case 3: // Linear Splines
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd = theta_[i]*tempd;
      tempd = 1.0-tempd;
      if( tempd > 0 )
        result = result*tempd;
      else {
        result = 0.0;
        break;
      }
    } // end for
    break;
    /*  case 4: // cubic correlation on the unit cube
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd2 = tempd*tempd;
      tempd = tempd2*tempd;
      tempd2 = theta_[i]*tempd2;
      tempd  = theta_[i+p]*tempd;
      tempd = 1.0-tempd2+tempd;
      result = result*tempd;
    } // end for
    break;*/
  case 5: // Mate'rn flexible family   UNDER CONSTRUCTION
    // NOTE:  alphavec should be all integers
    // Mate'rn tends to be more computationally intensive than other functions
    for( int j = p; j<2*p; j++) {
      if( theta_[j] > size)
        size = (int)theta_[j];
    }
    new_array(k,size+1);
    for( int i = 0; i<p; i++) {
      tempd3 = fabs(y[i] - x[i]);
      tempd = pow(tempd3, (int)theta_[i+p]);
      tempd = tempd * theta_[i];
      tempd2 = theta_[i+p]-1.0;
      tempd2 = pow(2, tempd2);
      tempv[0] = gamma( theta_[i+p] );
      tempd2 = tempd2 * tempv[0];
      tempd = tempd / tempd2;
      tempd3 = tempd3 * theta_[i];
      bessk( tempd3,(int)theta_[i+p], k );
      tempd2 = k[((int)theta_[i+p])];
      result = result * (tempd2 * tempd);
    } // end for
    delete k;
    break;
  case 0: corfunction(x,y,result);
    break;
  default: return 0;
  }; // end switch
  if( DISPLAYK > 4 )
    cout<< "result = " << result << endl;
  return result;
} // end correlation()

double krigapprox::trend( const Vector<double>& point) {
  if( mu_ == NULL )
    return 0.0;
  else
    return mu_[0];
} // end trend()

int krigapprox::thetasize() {
  int ans = 0;
  switch( corfamily) {
  case 0: theta_ = NULL;
    ans = 0;
    break;
  case 1: ans = 2;
    break;
  case 2: ans = p+1;
    break;
  case 3: ans = p;
    break;
  case 4:
  case 5: ans = 2*p;
    break;
  }; // end switch
  return ans;
} // end thetasize()

int krigapprox::musize() {
  return 1;
} // end musize

void krigapprox::estimateMu(const Matrix<double>& R) {
  if( ValuesSet == FALSE )
    notseterror();
  Vector<double> A(n);
  Vector<double> vp(n);
  A = 1;
  SolveSystem( R, vp, A);
  mu_ = new double;
  mu_[0] = (vp*y)/(vp*A);
  
} // end estimateMu()

int krigapprox::estimateTheta() {
  ParamEstimate PE(p);
  int correlate;
  pt_collect* MLE_Matrix = NULL;
  new_Matrix(MLE_Matrix, 1,1);
  mp_vector* vp = NULL;
  new_Vector(vp, n);
  Vector<double> tempv;
  double* sigma2 = new double;
  pt_collect Xp = X;
  mp_vector yp = y;

  mp_vector* theta = NULL;
  new_Vector( theta,thetasize()-1);

  (*theta) = 17;  // just an arbitrary choice
  
  // This was once suggested by trosset, but it behaves funny.
  //(*theta)[0] = -100*log(0.5)/sqr(L);

  (*MLE_Matrix).newsize(n+p, n+p);
  theta_ = NULL;
  mu_ = NULL;
  new_array( theta_, thetasize());
  new_array( mu_, musize());
  switch( corfamily ){
  case 0: return TRUE;
  case 1:
    if(     /*theta_[1] == 1*/ FALSE )  // default to alpha = 2
      correlate = 0;
    else
      correlate = 1;
    theta_[0] = 1.0;
    theta_[1] = 2.0;
    break;
  case 2:
    if( /* theta_[1] == 1*/ FALSE ) // default to alpha = 2
      correlate = 2;
    else
      correlate = 3;
    theta_[p] = 2;
    for(int j = 0; j < p; j++)
      theta_[j] = 1.0;
    break;
  default:
    if( DISPLAYK > 0 )
      cout << "ERROR:  Cannot do estimation for choice " << corfamily
           << ".  Setting parameters to 0." << endl;
    for( int i = 0; i < thetasize(); i++)
      theta_[i] = 0;
    return FALSE;
  }; // end switch
  PE.ResetCorrelationFamily( correlate);
  PE.EstimateConstantTrend( Xp, yp, DELTA, mu_, sigma2, theta, MLE_Matrix, vp );

  for(int k=0; k<(*theta).size(); k++)
    theta_[k] = (*theta)[k];
  switch( corfamily ) {
  case 1: theta_[1] = 2;
    break;
  case 2: theta_[p] = 2;
    break;
  };
  if( MLE_Matrix != NULL )
    delete MLE_Matrix;
  if( vp != NULL )
    delete vp;

  return TRUE;
} // end estimateTheta

void krigapprox::chooseCorrelationFamily( int choice ) {
  /*  Allows the user to choose among the possible families of correlation
      functions.  Will delete mu and theta and reapproximate if the other values
      have been set.

      Currently, the options are:

      choice = 0  =>  User defined function passed in previously

      choice = 1  =>  Gaussian isotropic r(x,y) = exp( -theta * ||y-x||^alpha )

      choice = 2  =>  Gaussian product \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|^alpha)

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

  corfamily = choice;
  if( ValuesSet == TRUE) {
    deleteall();
    approximate();
  }
 }

 void krigapprox::chooseCorrelationFamily( int choice,  // indicates the correlation function
                                         double* theta ) // used for passing known
   // parameter values.  If values unknown, pass in theta=NULL.
  
  /*  Allows the user to choose among the possible families of correlation
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
{
  corfamily = choice;
  deleteall();
  if( theta != NULL ) {
    new_array( theta_, thetasize());
    for( int i = 0; i < thetasize(); i++)
      theta_[i] = theta[i];
  } // end if
  else
    theta_ = NULL;
  if( ValuesSet == TRUE)
    approximate();
}

void krigapprox::chooseCorrelationFamily( void (*corfunc)(const Vector<double>&,
                                                        const Vector<double>&,
                                                        double& ) ) {
  corfamily = 0;
  corfunction = corfunc;
  deleteall();
  if( ValuesSet == TRUE)
    approximate();
} // end chooseCorrelationFamily()


void krigapprox::storeApprox( ostream& outf ){
/****
     Prints the data of the approximation to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing RBF's so that they might be re-created later.
****/
  long i, j;

  outf << p << endl << n << endl<< corfamily << endl;
  if (corfamily > 0) {
    outf << 'y' << endl;
    switch( corfamily ) {
    case 1:
      outf << theta_[0] << ' ' << theta_[1]<< endl;
      break;
    case 2:
      for( i = 0; i<p; i++)
        outf << theta_[i]<< ' ';

      outf << theta_[p] << endl;
      break;
    case 3:
      for( i = 0; i<p; i++)
        outf << theta_[i] << ' ';
      outf << endl;
      break;
    case 4:  
    case 5:
      for( i = 0; i<2*p; i++)
        outf << theta_[i] << ' ';
      outf << endl;
      break;  
    }; // end output switch
  } // end if
  else // corfamily == 0
    outf << 'n' << endl;

  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      outf << X[i][j] << ' ';
    } // end for
    
    outf << y[i] << endl;
  } // end for
  return;
} // end storeApprox()







