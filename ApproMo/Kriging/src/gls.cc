/*  gls.cc
 written by Anthony Padula	adpadu@maila.wm.edu
This file contains the function implementations of the setvalues() and evalf()
functions to create and evaluate a GLS estimate.
These functions are encapsulated inside the implementation of the glsapprox
class.

Last Modified 6/10/2000
*/
#include "gls.h"
#include <cassert>

glsapprox::glsapprox() { 
  k = 0;
  n = 0;
  p = 0;
  ValuesSet = FALSE;
  Rinv = NULL;
}

glsapprox::glsapprox( glsapprox& ToBeCopied )   /* Last Modified 7/14/99 */
{
  ValuesSet = ToBeCopied.ValuesSet;
  if( ValuesSet) {
    k = ToBeCopied.k;
    n = ToBeCopied.n;
    p = ToBeCopied.p;
    X.newsize(n,k);
    y.newsize(n);
    beta.newsize(k);
    Rinv = NULL;

    X = ToBeCopied.X;
    y = ToBeCopied.y;
    beta = ToBeCopied.beta;
      
  } // end if
} // end Copy constructor

glsapprox::glsapprox( Matrix<double>& Points,  // The n x p matrix of points, one point to a row
                      Vector<double>& Values)  // The vector of function values,
                                                 // so f(Points[i])=Values[i]
                       
{
  // Last modified 6/10/2000
  long i, j, g, h;

  ValuesSet = TRUE;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  k = 1+p+(p*(p+1))/2;
  Rinv = NULL;

  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
  }
  cout << "Design matrix X = " << X; 
  approximate();
}

 glsapprox::glsapprox(  Matrix<double>& Points, // The n x p matrix of points, one point to a row
			Vector<double>& Values, // The vector of function
                                                // values, so f(Points[i])=Values[i]
			Matrix<double>& InvCorr) // The correlation matrix 
{
  // Last modified 6/10/2000
  long i, j, g, h;

  ValuesSet = TRUE;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  k = 1+p+(p*(p+1))/2;
  new_Matrix(Rinv,n,n); 
    
  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);

    for( j = 0; j < n; j++)
      (*Rinv)[i][j] = InvCorr[i][j];
  }
  cout << "Design matrix X = " << X; 
  approximate();
}


glsapprox::~glsapprox() {    /* Last Modified 6/10/2000 */
  if( DISPLAYG > 3 )
    cout << "deleteing" << endl;
  
  deleteall();
  if( DISPLAYG > 3 )
    cout << "done deleting" << endl;
} // end destructor

// end of all constructors and destructors


void glsapprox::deleteall() {
  if( Rinv != NULL ) {
    delete Rinv;
    Rinv = NULL;
  }
} // end deleteall()



// The setvalues() functions

int glsapprox::setvalues() {  /* Last Modified 7/14/99 */
  long i, j, g, h;
  Matrix<double> Points;

  cout << "Please input values for:" << endl << "p = ";
  cin >> p;
  cout << " n = ";
  cin >> n;

  k = 1+p+(p*(p+1))/2;
  X.newsize(n,k);
  y.newsize(n);
  beta.newsize(k);
  Points.newsize(n,p);
  Rinv = NULL;
  cout << " Points (by rows) = ";
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      cin >> Points[i][j];
    } // end for
    
    cin >> y[i];
  } // end for
  if( DISPLAYG > 1 ) 
    cout << "Input received.  Calculating approximation..." << endl;
  ValuesSet = TRUE;
  
  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
  }
  //cout << "Design matrix X = " << X;  
 
  if( approximate() == TRUE) {
    if( DISPLAYG > 0 )
      cout << "Approximation finished." << endl;
    return TRUE;
  }
  else {
    if( DISPLAYG > 0 )
    cout << "Approximation failed." << endl;   
    return FALSE;
  }
} // end setvalues()
  
int glsapprox::setvalues(char* filename) {  /* Last Modified 7/14/99 */
  ifstream infile;   
  infile.open( filename, ios::in );
  long i, j, g, h;
  Matrix<double> Points;

  infile >> p >> n;

  k = 1+p+(p*(p+1))/2;
  X.newsize(n,k);
  y.newsize(n);
  beta.newsize(k);
  Points.newsize(n,p);
  Rinv = NULL;
  for( i= 0; i < n; i++) { 
    for( j= 0; j < p; j++) {
      infile >> Points[i][j];
    } // end for
    
    infile >> y[i];
  } // end for
  if( DISPLAYG > 1 ) 
    cout << "Input received.  Calculating approximation..." << endl;
  ValuesSet = TRUE;
  
  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
  }
  //cout << "Design matrix X = " << X;  

  if( !approximate() )
    if( DISPLAYG > 0 )
      cout << "ERROR: Approximation failed. Check input file" << endl;
  infile.close();
  return TRUE;
} // end setvalues(filename)

int glsapprox::setvalues( Matrix<double>& Points,// The n x p matrix of points, one point to a row
                        Vector<double>& Values)  // The vector of function values,
                                                 // so f(Points[i])=Values[i]
{
  // Last modified 6/10/2000
  long i, j, g, h;

  deleteall();
  ValuesSet = TRUE;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  k = 1+p+(p*(p+1))/2;
  Rinv = NULL;
  X.newsize(n,k);
  y.newsize(n);
  beta.newsize(k);
  
  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
  }
  //cout << "Design matrix X = " << X; 
  return approximate(); 
}

int glsapprox::setvalues( Matrix<double>& Points,// The n x p matrix of points, one point to a row
                        Vector<double>& Values,  // The vector of function values,
                                                 // so f(Points[i])=Values[i]
			  Matrix<double>& InvCorr) // The inverse of the correlation matrix
{
  // Last modified 6/10/2000
  long i, j, g, h;

  deleteall();
  ValuesSet = TRUE;
  y = Values;
  n = Values.size();
  p = Points.num_cols();
  k = 1+p+(p*(p+1))/2;
  new_Matrix(Rinv,n,n); 


  /*  X.newsize(n,k);
  y.newsize(n);
  beta.newsize(k);

  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
    for( j = 0; j < n; j++)
      (*Rinv)[i][j] = InvCorr[i][j];
      } */

  //plan to modify to beta0 + (x-x0)^{T} (tau I) (x - x0)
  X.newsize(n,k);
  y.newsize(n);
  beta.newsize(k);

  for( i = 0; i < n; i++) {
    X[i][0] = 1;
    g = p+1;
    for( j = 1; j <=p; j++) {
      X[i][j] = Points[i][j-1];
      
      for(h = j-1;h < p; h++) {
	X[i][g] = Points[i][j-1]*Points[i][h];
	g++;
      }
    }
    assert(g==k);
    for( j = 0; j < n; j++)
      (*Rinv)[i][j] = InvCorr[i][j];
      }

  //cout << "Design matrix X = " << X; 
  return approximate(); 
}
// end of setvalues() functions


int glsapprox::setInverseCorrelation( Matrix<double>& InvCorr )
{
  int num = InvCorr.num_rows();
  if( Rinv == NULL) {
    new_Matrix(Rinv,num,num);
  }

  (*Rinv) = InvCorr;
  return 0;
}


// functions to access the data members

long glsapprox::getk() {
  if( !ValuesSet)
    notseterror();
  return k;
}

Vector<double> glsapprox::getbeta() {
  
   if ( !ValuesSet )
     notseterror();
   
   return beta;
}


Vector<double> glsapprox::getresiduals(){
  Vector<double> resid(n);

  if ( !ValuesSet )
    notseterror();
    
  resid = y - X*beta;

  return resid;
}


// approximate, evaluate and cohorts


int glsapprox::approximate() { /* Last modified 6/10/2000 */
  Matrix<double> tempmat(k,k), tempinv(k,k);
  bool success;
  double tol;

  if( Rinv != NULL ) {
    tempmat = (transpose(X) * (*Rinv)) *  X;
    beta = transpose(X) * ((*Rinv) * y);
  }
  else {
    tempmat = transpose(X) * X;
    beta = transpose(X) * y;
  }

  if( n >= k) 
    success = ComputeInverse( tempmat, &tempinv);
  if( (n < k)||(!success)) {
    Matrix<double> U(k,k), VT(k,k), D(k,k);
    ComputeSVD(tempmat,(&U), (&D), (&VT)); // Calculate the SVD
    
  
    // Set appropriate tolerance
    tol = n*D[0][0]*MachineEpsilon();
    
    for(long i = 0; i < n; i++) { 
      if (fabs(D[i][i]) < tol)
	D[i][i] = 0;
      else
	D[i][i] = 1.0/D[i][i];
    } // end for
    tempinv = U * (D * VT);
  }
    
  beta = tempinv * beta;
  return TRUE;
} // end approximate()

double glsapprox::evalf(const Vector<double>& point) { /* Last modified 6/10/2000 */
   double value = 0;
   Vector<double> tempv(k);  
   long j, g, h; 

   if ( !ValuesSet )
     notseterror();
  
   tempv[0] = 1;
   g = p+1;
   for( j = 1; j <=p; j++) {
     tempv[j] = point[j-1];
     
     for(h = j-1;h < p; h++) {
       tempv[g] = point[j-1]*point[h];
       g++;
     }
   }
   assert(g==k);

   value = tempv * beta;
   
   return value;
} // end evalf()

Vector<double> * glsapprox::evalderiv(const Vector<double>& point) {
  long g, h, j;
  Vector<double> * value;
  new_Vector(value,p);
  (*value) = 0;

  g = p+1;
  for( j = 1; j <=p; j++) {
    (*value)[j-1] += beta[j];
    
    for(h = j-1;h < p; h++) {
      // on square terms, it just gets added twice to get
      // (*value)[j-1] += 2.0*beta[g]*point[j-1]
      
      (*value)[j-1] += (beta[g])*(point[h]);
      (*value)[h] += (beta[g])*(point[j-1]);      
      g++;
    }
  }
  return value;
}

void glsapprox::storeApprox( ostream& outf ){
/****
     Prints the data of the approximation to the output stream in the same format as it
     would be read in from a file.  Thus, this function is appropriate for
     storing RBF's so that they might be re-created later.
****/
  long i, j;

  outf << p << endl << n << endl;
  for( i= 1; i < n; i++) { 
    for( j= 1; j < p; j++) {
      outf << X[i][j] << ' ';
    } // end for
    
    outf << y[i-1] << endl;
  } // end for
  return;
} // end storeApprox()



/****  krigapproxWithGLS  functions ****/

 // Uses the new gls trend
double krigapproxWithGLS::trend( const Vector<double> & point){
  return GLStrend.evalf(point);
}

Vector<double> * krigapproxWithGLS::derivtrend(const Vector<double>& point) {
  return GLStrend.evalderiv(point);
}

// creates the approximation using the current values
// Will estimate theta if necessary
// Returns TRUE if Successful, FALSE otherwise.
int krigapproxWithGLS::approximate() {

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
  if( DISPLAYG > 3 )
    cout << "R = I" << endl;
  for( i=1; i < n; i++) { // Compute the inter-point correlations
    temppointi = X.row(i); 
    for( j=0; j < i; j++) {
      R[i][j] = correlation( temppointi, X.row(j) );
      R[j][i] = R[i][j];
    } // end for
  } // end for

  if( DISPLAYG > 3 )
    cout << "Correlation matrix computed. " << endl;
  if( DISPLAYK > 4 )
    cout << R << endl;

  success = ComputeInverse(R, (&RINV)); // Calculate the Cholesky Decomposition
  
  int count = 0;
  while( (!success)&&(count < MAXCHOLREPS) ) {
    // add a constant to the diagonal
    double TOADD = n*(n*MachineEpsilon())/(1-n*MachineEpsilon());
    
    for( int i = 0; i < n; i++)
      R[i][i] += TOADD;
    
    count++;
    cout << "Looped " << count << endl;
    // Calculate the Cholesky Decomposition
    success = ComputeInverse(R, (&RINV));
    }
  

  if(!success) {
    cout << "Cholesky failed:  Correlation matrix may be positive definite."
	 << endl << "Trying SVD...";
    
    Matrix<double> U(n,n), VT(n,n), D(n,n);
    ComputeSVD(R,(&U), (&D), (&VT)); // Calculate the SVD
    if( DISPLAYK > 3 )
      cout << "SVD computed" << endl;
  
    // Set appropriate tolerance
    tol = n*D[0][0]*MachineEpsilon();
    
    if( DISPLAYG > 1 )
      cout << "tol = " << tol << endl;
    for( i = 0; i < n; i++) { 
      if (fabs(D[i][i]) < tol)
	D[i][i] = 0;
      else
	D[i][i] = 1.0/D[i][i];
    } // end for
    if( DISPLAYG > 3 )
      cout << " D^(-1) computed" << endl;
    RINV = U * (D * VT); 
  
    if( DISPLAYG > 1 )
      cout << "mu = " << mu_[0] << ", theta = " << theta_[0] << endl;
  } // end else for SVD
  
  GLStrend.setvalues(X, y, RINV);   // Generalize Linear Least Squares
  //GLStrend.setvalues(X,y);          // Ordinary Linear Least Squares 
  for( int z = 0; z < n; z++) muv[z] = GLStrend.evalf(X.row(z));
  v = (RINV * (y-muv));  

    if( DISPLAYG > 3 )
      cout << "v done."  << endl;

  return TRUE;
} // end approximate()






