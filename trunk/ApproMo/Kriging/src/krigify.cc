/*  krigify.cc
 written by Anthony Padula	adpadu@maila.wm.edu
This file contains the function implementations for
krigify and evalf to produce random functions
with a specified trend, as if they were being
drawn from a second-order Gaussian function, or other specified correlation.
These functions are encapsulated inside the implementation of the randfunc
class.

This procedure is described in
	Trosset Michael W.
	The Krigifier: A Procedure for Generating Pseudorandom
	Nonlinear Objective Functions for Computational
	Experimentation

Last Modified 6/10/2000

This is a modified version of krigify.cc that has been upgraded to attempt to use
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
#include "krigify.h"
#include <iomanip>
#include <cstdio>

randfunc::randfunc() {

   ValuesSet = FALSE;
   corfamily = 1;  // Default to Gaussian isotropic
   corfunction = NULL;
   trend = NULL;
   designchoice = 1; // Default to Uniform Random Sampling
   design = NULL;
   PlantSeeds(-1); // Seed Steve Park's random number generator
   GetSeed( &seed);
   cout << "Your seed is " << seed << endl;
   rngstream = 42;
   SelectStream(rngstream);  // Default to stream 42
   beta1 = NULL;
   beta2 = NULL;
   v = NULL;
   x0 = NULL;
   lower = NULL;
   upper = NULL;
   thetavec = NULL;
   alphavec = NULL;
}

randfunc::randfunc( randfunc& ToBeCopied )   /* Last Modified 6/23/99 */
{
  deleteall();
  ValuesSet = ToBeCopied.ValuesSet;
  if( ValuesSet) {
    p = ToBeCopied.p;
    n = ToBeCopied.n;
    seed = ToBeCopied.seed;
    rngstream = ToBeCopied.rngstream;
    beta0 = ToBeCopied.beta0;
    corfamily = ToBeCopied.corfamily;
    corfunction = ToBeCopied.corfunction;
    trend = ToBeCopied.trend;
    designchoice = ToBeCopied.designchoice;
    design = ToBeCopied.design;
    
    new_Vector(beta1, p);
    new_Matrix(beta2,p,p);
    X.newsize(n,p);
    new_Vector(x0, p);
    new_Vector(lower,p);
    new_Vector(upper,p);
    new_Vector(v,n);
    beta1->operator=( *(ToBeCopied.beta1));
    beta2->operator=( *(ToBeCopied.beta2));
    x0->operator=( *(ToBeCopied.x0));
    lower->operator=( *(ToBeCopied.lower));
    upper->operator=( *(ToBeCopied.upper));
    X = ToBeCopied.X;
    v->operator=( *(ToBeCopied.v));
    
    alpha = ToBeCopied.alpha;
    theta = ToBeCopied.theta;
    sigma2 = ToBeCopied.sigma2;
    
    if( (corfamily == 2)||(corfamily == 3)||(corfamily == 4)||(corfamily ==
							       5)) {
      if( thetavec != NULL )
	delete thetavec;
      new_Vector(thetavec,p);
      thetavec->operator=( *(ToBeCopied.thetavec));
    }// end if
    
    if( (corfamily == 4)||(corfamily == 5)) {
      if( alphavec != NULL )
	delete alphavec;
      
       new_Vector(alphavec,p);
      alphavec->operator=( *(ToBeCopied.alphavec));
    } // end if
  } // end if 
} // end Copy constructor

randfunc::randfunc(long seedtoplant, int stream) {
   ValuesSet = FALSE;
   corfamily = 1;   // Default to Gaussian isotropic
   corfunction = NULL;
   designchoice = 1; // Default to Uniform Random Sampling
   design = NULL;
   seed = seedtoplant;
   PlantSeeds(seed); // Seed Steve Park's random number generator
  
   if( seed == -1)
     GetSeed( &seed);
   rngstream = stream;
   SelectStream(stream);
   trend = NULL;
   beta1 = NULL;
   beta2 = NULL;
   v = NULL;
   x0 = NULL;
   lower = NULL;
   upper = NULL;
   alphavec = NULL;
   thetavec = NULL;
}

randfunc::~randfunc() {    /* Last Modified 6/23/99 */
  deleteall();
  if( thetavec != NULL ) {
    delete thetavec;
    thetavec = NULL;
  }
  if( alphavec != NULL ) {
    delete alphavec;
    alphavec = NULL;
  }
  corfunction = NULL;
  trend = NULL;
  design = NULL;
} // end destructor

int randfunc::newrand() {
   
   if( ValuesSet == FALSE)
      return FALSE;

   if( v != NULL )
      delete v;

   X.newsize(n,p);
   
   new_Vector(v,n);
   //      cout << "X and v allocated" << endl;
   if (v == NULL) // if allocation failed
	return FALSE; 	

   krigify(1);
   return TRUE;
} // end newrand()

int randfunc::newrandsameX() {

   if( ValuesSet == FALSE)
      return FALSE;

   if( v != NULL )
      delete v;

   new_Vector(v,n);

   krigify(2);
   return TRUE;
} // end newrandsameX()

void randfunc::setvalues() {  /* Last Modified 6/23/99 */
   ValuesSet = TRUE;
   int i, j;
   deleteall();
   cout << "Please input values for:" << endl << "p = ";
   cin >> p;
   cout << " n = ";
   cin >> n;
   if (trend == NULL) {   
     cout << " beta0 = ";
     cin >> beta0;
     cout << "Input the p-vector beta1 = ";

     new_Vector(beta1,p);
     for(i = 0; i < p; i++)
       cin >> (*beta1)[i];
     cout << "Input by rows the p by p Matrix beta2 = ";
     new_Matrix(beta2,p,p);
     for( i = 0; i < p; i++)
       for( j = 0; j < p; j++)
         cin >> (*beta2)(i,j);
     cout << "Input the p-vector x0 = ";
     new_Vector(x0,p); 
     for( i = 0; i < p; i++)
       cin >> (*x0)[i];
   } // end if

   if( (corfamily == 4)||(corfamily == 5)) {
     if( alphavec != NULL )
       delete alphavec;
     
     new_Vector(alphavec,p);
     cout << "input the p-vector alpha = ";
     for( i = 0; i<p; i++)
       cin >> (*alphavec)[i];
   } // end if
   else if((corfamily == 0)||(corfamily==3)) {
     alpha = 0;
   } // end else if
   else {
     cout << " alpha = ";
     cin >> alpha;
   } // end else if

   if( (corfamily == 2)||(corfamily == 3)||(corfamily == 4)||(corfamily == 5)) {
     if( thetavec != NULL )
       delete thetavec;
     
     new_Vector(thetavec,p);
     cout << "input the p-vector theta = ";
     for( i = 0; i<p; i++)
       cin >> (*thetavec)[i];
   }// end if
   else if (corfamily != 0) {
     cout << " theta = ";
     cin >> theta;
   } // end if
   else { // corfamily == 0  => User specified correlation function
     theta = 0;
   } // end else
   
   cout << " sigma2 = ";
   cin >> sigma2;
   cout << "Input the p-vector lower = ";
   new_Vector(lower,p);
   for( i = 0; i < p; i++)
      cin >> (*lower)[i];
   cout << "Input the p-vector upper = ";
   new_Vector(upper,p);
   for( i = 0; i < p; i++) {
     cin >> (*upper)[i];
     if ( (*upper)[i] <= (*lower)[i]) {
       cout << "Upper bound must be greater than lower bound of " << (*lower)[i]
            << "Please re-enter entire vector:";
       i = 0;
     } // end if
   } // end for
   cout << "Input received.  Executing krigify" << endl;
   if( newrand() == TRUE)
      cout << "Krigify finished.  Values set successfully." << endl;
   else
      cout << "Krigify failed.  Memory possibly full." << endl;   
   return;
} // end setvalues()

void IgnoreCommentLines( istream& infile, char commentchar ) {
  char tempchar = commentchar;
  do {
    tempchar = infile.peek();  // look at next char
    if( tempchar == commentchar )
      infile.ignore(255,'\n');
  } while( tempchar == commentchar );
  return;
}


int randfunc::setvalues(char* filename) {  /* Last Modified 6/23/99 */
  ifstream infile;   
  infile.open( filename, ios::in );
  int i, j;

  if( infile.fail() ) {
    cout << "File failed to open " << filename << ".  Returning." << endl;
    return false;
  }
  deleteall();

  /*
  char templine[500];

  templine[0] = '#';
  while(templine[0] == '#')
    fgets(templine, 500, infile);
  sscanf(templine, "%ld", &p);
  templine[0] = '#';
  while(templine[0] == '#')
    fgets(templine, 500, infile);
  sscanf(templine, "%ld", &n);
  if( trend == NULL ) {
    templine[0] = '#';
    while(templine[0] == '#')
      fgets(templine, 500, infile);
    sscanf(templine, "%lf", &beta0);
    templine[0] = '#';
    while(templine[0] == '#')
      fgets(templine, 500, infile);
    new_Vector(beta1,p);
    i=0;
    i = sscanf( templine

    tempchar = '#';
    do {
    infile >> tempchar;
    if( tempchar == '#' )
    infile.ignore(255,'\n');
    } while( tempchar == '#' );
    fseek( infile, -1, SEEK_CUR);
*/  


  IgnoreCommentLines(infile, '#');
  infile >> p;
  infile.ignore(255, '\n');
  IgnoreCommentLines(infile, '#');
  infile >> n;
  infile.ignore(255, '\n');

  if( trend == NULL) {
    IgnoreCommentLines(infile, '#');
    infile >> beta0;
    infile.ignore(255, '\n');
    new_Vector(beta1,p);
    IgnoreCommentLines(infile, '#');
    for(i = 0; i < p; i++)
      infile >> (*beta1)[i];
    infile.ignore(255, '\n');
    new_Matrix(beta2,p,p);
    IgnoreCommentLines(infile, '#');
    for( i = 0; i < p; i++)
      for( j = 0; j < p; j++)
	infile >> (*beta2)(i,j);
    infile.ignore(255, '\n');
    //cout << ", beta2" << endl;
    new_Vector(x0,p);
    IgnoreCommentLines(infile, '#');
    for( i = 0; i < p; i++)
      infile >> (*x0)[i];
    infile.ignore(255, '\n');
  } // end if
  if( infile.eof()) {
    cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
	 << endl << " Quitting setvalues()" << endl;
    return FALSE;
  } // end if
  if( (corfamily == 4)||(corfamily == 5)) {
    if( alphavec != NULL )
      delete alphavec;
    
    new_Vector(alphavec,p);
    IgnoreCommentLines(infile, '#');
    for( i = 0; i<p; i++)
      infile >> (*alphavec)[i];
    infile.ignore(255, '\n');
  } // end if
  else if ((corfamily == 3)||(corfamily==0) ) {
    alpha = 0;
  } // end else if
  else {
    IgnoreCommentLines(infile, '#');
    infile >> alpha;
    infile.ignore(255, '\n');
  } // end else if
  
  if( (corfamily == 2)||(corfamily == 3)||(corfamily == 4)||(corfamily == 5)) {
    if( thetavec != NULL )
      delete thetavec;
    
    new_Vector(thetavec,p);
    IgnoreCommentLines(infile, '#');
    for( i = 0; i<p; i++)
      infile >> (*thetavec)[i];
    infile.ignore(255, '\n');
  }// end if
  else if (corfamily == 0) {
    theta = 0;
  } // end if
  else { // corfamily == 0  => User specified correlation function
    IgnoreCommentLines(infile, '#');
    infile >> theta;
    infile.ignore(255, '\n');
  } // end else

  IgnoreCommentLines(infile, '#');
  infile >> sigma2;
  infile.ignore(255, '\n');
  if( infile.eof()) {
    cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
	 << endl << " Quitting setvalues()" << endl;
    return FALSE;
  } // end if
  new_Vector(lower,p);
  IgnoreCommentLines(infile, '#');
  for( i = 0; i < p; i++) {
    infile >> (*lower)[i];
    if( infile.eof()) {
      cout << "Error: File terminated prematurely.  Check that all vectors are the correct length"
	   << endl << " Quitting setvalues()" << endl;
      return FALSE;
    } // end if
  } // end for
  infile.ignore(255, '\n');
  new_Vector(upper,p);
  IgnoreCommentLines(infile, '#');
  for( i = 0; i < p; i++) {
    infile >> (*upper)[i];
    if( infile.eof()) {
      cout << "Error: File terminated prematurely."
	   <<"  Check that all vectors are the correct length"
	   << endl << " Quitting setvalues()" << endl;
      return FALSE;
    } // end if
    if ( (*upper)[i] <= (*lower)[i]) {
      cout << "Error: Upper bound " << (*upper)[i] << " must be greater than\n lower bound of "
	   <<(*lower)[i] << ".  Returning FALSE.";
      return FALSE;
    } // end if
  } // end for

  infile.ignore(255, '\n');
  ValuesSet = TRUE;
  if( !newrand() )
    cout << "Krigify failed.  Memory possibly full." << endl;
  infile.close();
  return TRUE;
} // end setvalues(filename)

int randfunc::setvalues(long pin, long nin, double beta0in, Vector<double>& beta1in,
			Matrix<double>& beta2in, Vector<double>& x0in, 
			Vector<double>& alphain, Vector<double>& thetain, 
			double sigma2in, Vector<double>& lowerin,
			Vector<double>& upperin)
{
  int i;

  deleteall();
  p = pin;
  n = nin;
  if( trend == NULL) {
    beta0 = beta0in;
    new_Vector(beta1,p);

    (*beta1) = beta1in;
    new_Matrix(beta2,p,p);
    (*beta2) = beta2in;

    new_Vector(x0,p);
    (*x0) = x0in;
  } // end if
  if( (corfamily == 4)||(corfamily == 5)) {
    if( alphavec != NULL )
      delete alphavec;
    
    new_Vector(alphavec,p);
    (*alphavec) = alphain; 
  } // end if
  else if ((corfamily == 3)||(corfamily==0) ) {
    alpha = 0;
  } // end else if
  else {
    alpha = alphain[0];
  } // end else if
  
  if( (corfamily == 2)||(corfamily == 3)||(corfamily == 4)||(corfamily == 5)) {
    if( thetavec != NULL )
      delete thetavec;
    
    new_Vector(thetavec,p);
    (*thetavec) = thetain;
  }// end if
  else if (corfamily == 0) {
    theta = 0;
  } // end if
  else { // corfamily == 0  => User specified correlation function
    theta = thetain[0];
  } // end else
  
  sigma2 = sigma2in;
  new_Vector(lower,p);
  (*lower) = lowerin;
  new_Vector(upper,p);
  (*upper) = upperin;
  for( i=0; i<p; i++)
    if ( (*upper)[i] <= (*lower)[i]) {
      cout << "Error: Upper bound " << (*upper)[i] <<" must be greater than lower bound of "
	   <<(*lower)[i] << ".  Returning FALSE.";
      return FALSE;
    } // end if
  
  ValuesSet = TRUE;
  if( !newrand() )
    cout << "Krigify failed.  Memory possibly full." << endl;
  return TRUE;

} // end setvalues(parameters)

void randfunc::notseterror() {
   cout << "ERROR: Value not set." << endl;
   return;  // FAILS!!!
}


// functions to access the data members

long randfunc::getp() {
   if ( !ValuesSet )
        notseterror();
   return p;
}

long randfunc::getn() {
   if ( !ValuesSet )
        notseterror();
   return n;
}

double randfunc::getbeta0() {
   if ( !ValuesSet )
        notseterror();
   return beta0;
}

Vector<double> randfunc::getbeta1() {
    if ( !ValuesSet )
        notseterror();
    return (*beta1);
}

Matrix<double> randfunc::getbeta2(){
    if ( !ValuesSet )
        notseterror();
    return (*beta2);
}

Vector<double> randfunc::getx0(){
    if ( !ValuesSet )
        notseterror();
    return (*x0);
}

double randfunc::getalpha() {
   if ( !ValuesSet )
        notseterror();
   return alpha;
}
double randfunc::gettheta() {
   if ( !ValuesSet )
        notseterror();
   return theta;
}
Matrix<double> randfunc::getX(){
    if ( !ValuesSet )
        notseterror();
    return (X);
}

Vector<double> randfunc::getv(){
    if ( !ValuesSet )
        notseterror();

    return (*v);
}

Vector<double> randfunc::getlower(){
  if (!ValuesSet )
    notseterror();

  return (*lower);
}

Vector<double> randfunc::getupper(){
  if (!ValuesSet )
    notseterror();

  return (*upper);
}

long randfunc::getSeed() {
  if(!ValuesSet)
    notseterror();
  return seed;
}

int randfunc::getStream(){
  if( !ValuesSet)
    notseterror();
  return rngstream;
}



void randfunc::deleteall() {
  if( beta1 != NULL ) {
     delete  beta1;
     beta1 = NULL;
  }
  if( beta2 != NULL ) {
    delete  beta2;
    beta2 = NULL;
  }
  if( x0 != NULL ) {
    delete  x0;
    x0 = NULL;
  }
  if( lower != NULL) {
    delete lower;
    lower = NULL;
  }
  if( upper != NULL) {
    delete upper;
    upper = NULL;
  }
  if( v != NULL ) {
    delete v;
    v = NULL;
  }
  X.newsize(1,1);
} // end deleteall()

void randfunc::krigify(int flag) { /* Last modified 11/15/99 */

  Matrix<double> R(n,n);
  Matrix<double> L(n,n);
  Vector<double> randn(n);
  Vector<double> temppointi(p), temppointj(p);
  long i, j;
  double sigma = sqrt(sigma2), tol;
  ofstream outr;
  bool success;

  L = 0;
  R = 0;

  switch(flag) {

  case 1: 	/* create X */
    switch(designchoice) {
    case 1: // Uniform random sampling
      for( i=0; i < n; i++)  // Create n random points in the proper domain
	for( j=0; j < p; j++)
	  X[i][j] = ( (*lower)[j] + ((*upper)[j]-(*lower)[j])*(Random()));
      break;
    case 2: // Latin Hypercube
      genLH();
      break;
    case 0: // User defined
      design(p,n,(*lower),(*upper),X);
      break;
    } // end switch

  case 2 :
	/* create v using X */
    for( i = 0; i < n; i++)
      R[i][i] = 1;
  
    for( i=1; i < n; i++) { // Compute the inter-point correlations
      temppointi = X.row(i);      
      for( j=0; j < i; j++) {
        temppointj = X.row(j);
        R[i][j] = correlation( temppointi, temppointj );
        R[j][i] = correlation( temppointj, temppointi );
      } // end for
    } // end for
    
    success = ComputeCholesky(R, (&L)); // Calculate the Cholesky Decomposition

    int count = 0;
   
    while( (!success)&&(count < MAXCHOLREPS)) {
      // add a constant to the diagonal
      double TOADD = n*(n*MachineEpsilon())/(1-n*MachineEpsilon());
      
      for( int c = 0; c < n; c++)
	R(c, c) += TOADD;
      
      count++;
      cout << "Looped " << count << endl;
      // Calculate the Cholesky Decomposition
      success = ComputeCholesky(R, (&L));
      }
    
    if(success)
      {
        for(i =0; i < n; i++ )
          randn[i] = Normal(0,sigma);  // set up vector of n normal random numbers
        (*v) = (L * randn);
      } // end if
    else {
      cout << "Cholesky failed:  Correlation matrix may be positive definite."
           << endl << "Trying SVD...";
      Matrix<double> U(n,n), VT(n,n), D(n,n);
      ComputeSVD(R,(&U), (&D), (&VT)); // Calculate the SVD
      cout << "SVD computed" << endl;


      // Set appropriate tolerance
      if( MACHEPS == 1 ) {
        /*      tol = 1.0;
                do {
                tol = tol/2.0;
                } while( 1+tol > 1);
                tol = tol*2.0;*/
        // Compute macheps using formula in Heath:  |3*(4/3-1)-1|
        tol = 4.0/3.0;
        tol = tol-1.0;
        tol = 3.0*tol;
        tol = tol-1.0;
        tol = fabs(tol);
	cout << "Machine epsilon on this platform is " << tol << endl;
       
        tol = n*D(0,0)*tol;
      } else {
        tol = n*D(0,0)*MACHEPS;
      } // end else
      
      for( i = 0; i < n; i++) { 
        if (D[i][i] < tol)
          D[i][i] = 0;
        else
          D[i][i] = 1/sqrt(D(i,i));
        
        randn[i] = Normal(0,sigma);  // set up vector of n normal random numbers
      } // end for
      //cout << "Random vector created and D^1/2 computed" << endl;
      (*v) = ((U * D * VT) * randn);
      
    } // end else

  } // end switch
	
   return;
} // end krigify()

double randfunc::evalf(const Vector<double>& x) { /* Last modified 6/15/99 */
   double value = 0;
   Vector<double> r(n);
   Vector<double> tempv;

   if ( !ValuesSet )
      notseterror();
  
   for(int i=0; i<n; i++) {
      r[i] = correlation(X.row(i), x);
   } // end for

   if( trend == NULL) {
     tempv = x - (*x0);
     for( int j =0; j < p; j++)
       for( int k=0; k < p; k++)
	 value+= (*beta2)(j,k) * tempv[k] * tempv[j];

     value = beta0 + (*beta1)*tempv + value + (*v)*r;
     //value = (*v)*r;
   }
   else { // Use User specified trend
     trend(x, value);
     value = value + (*v)*r;
   }
   return value;
} // end evalf()


Vector<double>* randfunc::derivtrend( const Vector<double>& point)
{  
  Vector<double>* value;
  Vector<double> tempv;
  long p;
  p = point.size();
  tempv = point-(*x0);
  new_Vector(value,p);
  (*value) = (*beta1) +  ((*beta2)+ transpose(*beta2))*tempv ; 
  return value;
}

/* This returns the value of the derivative of |a-b| with respect to b,
   i.e. 1.0 or -1.0 . */

double randfunc::derivabs( double a, double b ) {
  if( a > b )
    return -1.0;
  else
    return 1.0;
}

Vector<double> randfunc::evalderiv( const Vector<double>& point) {
  /* Last modified 2/14/00 */
  
   double temp, t2;
   Vector<double> r(n), diff(p);
   Vector<double> Dprime(n);
   Vector<double> valuevec(p), * muprime;
   int alphatemp;

   if ( !ValuesSet )
      notseterror();

   if( (corfamily == 2)
       ||(corfamily==3)
       ||(corfamily==5)
       ||((corfamily==1)&&(alpha == 1)) ) {
     // If necessary, check to see if point is a known data point 
     for( int i  =0; i < n; i++)
       if( point == X.row(i) ) {
	 cout << "Not differentiable at this point.  Returning INF" << endl;
	 valuevec=HUGE_VAL;
	 return valuevec;
       }
   }
   //muprime = derivtrend((*x0));
   //cout << "muprime at x0 is " << (*muprime) << endl;
   //delete muprime;
   muprime = derivtrend(point);
   for(int k=0; k < p; k++) {
     switch( corfamily ) {
     case 1: /* Gaussian isotropic */
       if( alpha == 2) {  // alpha == 2
	 for(int i=0; i<n; i++) { // create delD/delx_k
	   diff = X.row(i)-point;
	   r[i] = correlation(X.row(i), point);
	   Dprime[i] = 2.0*r[i]*theta*diff[k];
	 }
       }
       else if(alpha==1) { //alpha == 1     
	 for(int i = 0; i < n; i++) { //create delD/delx_k
	   diff = X.row(i)-point;
	   r[i] = correlation(X.row(i), point);
	   Dprime[i] = theta * r[i]*diff[k] / (diff.l2norm()) ;
	 }  
       }
       break;
     case 2: /* Gaussian product */
       for( int i =0; i < n; i++) {
	 diff = (X.row(i) - point);
	 r[i] = correlation(X.row(i), point);
	 Dprime[i] = r[i]*(-(*thetavec)[k]* 2.0*(pow(fabs(diff[k]), (alpha-1)))
			   *(derivabs(X(i,k), point[k])));
       }
       break;
     case 3: /* Product of linear */
       for( int i =0; i < n; i++) {
	 temp = fabs(X(i,k) - point[k]);
	 r[i] = correlation(X.row(i) , point);
	 Dprime[i] = r[i]*((*thetavec)[k]*derivabs(X(i,k),point[k])) / (1.0-(*thetavec)[k]*temp); 
       }
       break;
       /* case 4: // Cubic correlation 
       for( int i = 0; i < n; i++) {
	 temp = X[i][k] - point[k];
	 r[i] = correlation(X.row(i), point);
	 Dprime[i]=r[i]*(temp*(2.0*(*thetavec)[k] - 3.0*(*alphavec)[k]*temp) / 
			(temp*temp*(1.0-(*thetavec)[k]+(*alphavec)[k]*temp)));
			}
			break; */
     case 5: /* Mate'rn */
       for( int i =0; i < n; i++) {
	 t2 = (*thetavec)[k];
	 alphatemp = (int) (*alphavec)[k];

	 //	 double Karray[alphatemp+2];

	 double* Karray = new double (alphatemp+2);

	 diff = (X.row(i)-point);
	 r[i] = correlation(X.row(i),point);
	 temp = t2*fabs(diff[k]);
	 bessk(temp, alphatemp, Karray);
	 Dprime[i] = derivabs(point[k], X(i,k))*r[i]*alphatemp / fabs(diff[k]);
	 Dprime[i] +=( r[i]*( Karray[alphatemp+1]*temp-alphatemp*Karray[alphatemp])
		       *derivabs(X(i,k), point[k])
		       / 
		       ( fabs(diff[k]) * Karray[alphatemp]));
	 delete Karray;
	 
       }
       break;
     case 0: 
     default:
       cout << "Cannot find derivative of user defined function.  Returning 0";
       return 0;
     } // end switch
     valuevec[k] = (*muprime)[k] + (*v)*Dprime;
   } // end for
   //   cout << "muprime at the point " << point << " is " << (*muprime) << '\n';
   //cout << "v*Dprime at the point is " << (valuevec - (*muprime)) << endl;
   delete muprime;
   return valuevec;
}  // end evalderiv

double randfunc::correlation( const Vector<double>& x, const Vector<double>& y )
{
  /* Last modified 6/21/99 */
  long size = 0;
  Vector<double> tempv(1);
  double tempd = 0, tempd2 = 0, tempd3 = 0, result=1;
  double *k = NULL;

  switch( corfamily ) {
  case 1: {  // Gaussian isotropic
    tempv.newsize(y.size());
    tempv = y-x;
    tempd = tempv.l2norm();
    tempd = pow(tempd,alpha);
    tempd = (-theta)*tempd;
    result = exp(tempd);
  }
    break;
  case 2: // Gaussian product
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd = pow(tempd, alpha);
      tempd = (-(*thetavec)[i])*tempd;
      tempd = exp(tempd);
      result = tempd*result;
    } // end for
    break;
  case 3: // Linear Splines
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd = (*thetavec)[i]*tempd;
      tempd = 1-tempd;
      if( tempd > 0 )
        result = result*tempd;
      else {
        result = 0;
        break;
      }
    } // end for
    break;
    /*case 4: // cubic correlation on the unit cube
    for( int i = 0; i<p; i++) {
      tempd = fabs(y[i]-x[i]);
      tempd2 = tempd*tempd;
      tempd = tempd2*tempd;
      tempd2 = (*thetavec)[i]*tempd2;
      tempd  = (*alphavec)[i]*tempd;
      tempd = 1-tempd2+tempd;
      result = result*tempd;
    } // end for
    break;*/
  case 5: // Mate'rn flexible family
    // NOTE:  alphavec should be all integers
    // Mate'rn tends to be more computationally intensive than other functions
    for( int j = 0; j<p; j++) {
      if( (*alphavec)[j] > size)
        size = (int)(*alphavec)[j];
    }
    new_array(k,size+1);

    for( int i = 0; i<p; i++) {
      tempd3 = fabs(y[i] - x[i]);
      tempd = pow(tempd3, (int)(*alphavec)[i]);
      tempd = tempd * (*thetavec)[i];
      tempd2 = (*alphavec)[i]-1;
      tempd2 = pow(2, tempd2);
      tempv[0] = gamma( (*alphavec)[i] );
      tempd2 = tempd2 * tempv[0];
      tempd = tempd / tempd2;
      tempd3 = tempd3 * (*thetavec)[i];
      bessk( tempd3,(int)(*alphavec)[i], k );
      tempd2 = k[((int)(*alphavec)[i])];
      result = result * (tempd2 * tempd);
    } // end for
    delete k;
    break;
  case 0: corfunction(x,y,result);
    break;
  default: return 0;
  }; // end switch
  return result;
} // end correlation()

double reportsample( randfunc& Func, Vector<double>& x, ostream& outf) {
  long p = x.size();
  double result= 0;
  int i;

  for( i = 0; i < p; i++ ) {
    outf << x[i] << ' ';
  } // end for

  result = Func.evalf(x);
  outf << result << endl;
  return result;
} // end reportsample()

void scanfunction( randfunc& Func, ostream& outf, int pointsperaxis)  // Last Modified 6/16/99
{ 
  long p = Func.getp();
  Vector<double> x(p);
  Vector<double> lower(p);
  Vector<double> upper(p);
  Vector<double> increment(p);
  int j;
  double temp = 1.0/pointsperaxis;

  if( p == 2) {  // use the faster version if possible
    scanplane( Func, outf, pointsperaxis);
    return;
  } // end if
  lower = Func.getlower();
  upper = Func.getupper();
  x = lower;
  increment = upper-lower;
  increment = increment * temp;

  while( x[p-1] < (upper[p-1]+increment[p-1]/2.0) ) {
    reportsample( Func, x, outf );
    x[0] += increment[0];
    for( j = 0; (j < p-1)&&(x[j] > (upper[j]+increment[j]/2.0)); j++) {
      x[j] = lower[j];
      x[j+1] += increment[j];
    } // end for
  } // end while
    
  return;
} // end scanfunction


// A default version which uses 25 pointsperaxis.  You can modify this default
// by changing the 25 above.
void scanfunction( randfunc& Func, ostream& outf) {
  scanfunction(Func, outf, 25);
}


 void scanplane( randfunc& Func, ostream& outf, int pointsperaxis)  // Last Modified 6/16/99
{ 
  double value = 0, stepe, stepn;
  Vector<double> x(2);
  Vector<double> lower(2);
  Vector<double> upper(2);
  
  lower = Func.getlower();
  upper = Func.getupper();

  stepe = (upper[0]-lower[0])/pointsperaxis;
  stepn = (upper[1]-lower[1])/pointsperaxis;

  for( x[0]=lower[0]; x[0]<(upper[0]+stepe/2.0); x[0]+= stepe ) {
    for( x[1]=lower[1]; x[1]<(upper[1]+stepn/2.0); x[1]+=stepn ) {
      outf << x[0] << ' ' << x[1] << ' ';
      value = Func.evalf(x);
      outf << value << endl;
    } // end for
    outf << endl;
  } // end for
  return;
} // end scanplane


 void randfunc::chooseCorrelationFamily( int choice ) {
  /*  Allows the user to choose among the possible families of correlation
      functions.

      Currently, the options are:

      choice = 0  =>  User defined function passed in previously

      choice = 1  =>  Gaussian isotropic r(x,y) = exp( -theta * ||y-x||^alpha )

      choice = 2  =>  Gaussian product \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|^alpha)

      choice = 3  =>  Product of Linear  \prod_{j=1}^{p} (1-theta_j * |y_j - x_j|)

      choice = 4  =>  Cubic correlation on the unit cube
                 \prod_{j=1}^{p} ( 1 - theta_j*(y_j - x_j)^2 + alpha_j*(y_j - x_j)^3 )

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
   
   if( (choice>0)&&(choice < 6) )
     corfamily = choice;
   else if ( (choice == 0)&&(corfunction != NULL) )
     corfamily = choice;
   else
     notseterror();
 }

void randfunc::chooseCorrelationFamily( void (*corfunc)(const Vector<double>&,
                                                        const Vector<double>&,
                                                        double& ) ) {
  corfamily = 0;
  corfunction = corfunc;
} // end chooseCorrelationFamily()

void randfunc::setTrend( void (*usertrend)( const Vector<double>& /*input*/ x,
                                   double& /*output*/ result)) {
  /* This allows the user to specify a non-quadratic trend function.  x is the
     point at which the trend should be evaluated, and result is used to return
     the value of the trend at x.

     Last Modified 6/16/99
  */
   trend = usertrend;
 }


void randfunc::generateTrend(double tau2, double* beta0in) {
  Matrix<double> Y(p+1, p);
  double a, b, tau;

  if(!ValuesSet) {
    notseterror();
    return;
  }

  tau = sqrt(fabs(tau2));
  if( tau == -1 )
    tau = 2.0* sqrt(sigma2);


  if( beta0in == NULL)
    beta0 = 0;
  else
    beta0 = *beta0in;

  if( beta1 == NULL )
    new_Vector(beta1,p);
  *beta1 = 0;

  if( beta2 == NULL)
    new_Matrix(beta2,p,p);

  if( x0 == NULL)
    new_Vector(x0,p);

  for(long i = 0; i < p+1; i++) {
    for(long j = 0; j < p; j++) {
      Y[i][j] = Normal(0, tau);
    }
    if( i < p) {
      a = (*lower)[i] +  (0.25)* ((*upper)[i] - (*lower)[i]);
      b = (*upper)[i] - (0.25)*((*upper)[i] - (*lower)[i]);
      (*x0)[i] = Uniform(a,b);
    }
  }

  (*beta2) = transpose(Y) * Y;
  (*beta2) = (1.0/(p+1.0)) * (*beta2);


  return;
}  

void randfunc::chooseInitialDesign(int choice)
  /* This allows the user to select which initial design (how X is chosen) is used.  
     The default is choice == 1, Uniform Random Sampling

     Currently, the options are:

      choice = 0  =>  User specified function (See below)
      
      choice = 1  =>  Uniform Random Sampling

      choice = 2  =>  Latin Hypercubes
  */
{
  if( (choice > 0)&&(choice < 3)) 
    designchoice = choice;
  else if( (choice == 0)&&(design != NULL) )
    designchoice = choice;
  else
    notseterror();
}

void randfunc::chooseInitialDesign( void (*userdesign)( const long /*input*/ p,
						const long /*input*/ n,
				     const Vector<double>& /*input*/ lower,
				     const Vector<double>& /*input*/ upper,
					   Matrix<double>& /*output*/ X) )    
  /* This allows the user to specify an alternative initial design.  The
     initial design must accept four inputs: 
     p is the dimension of the space
     n is the number of initial design sites (how many points in the space)
     lower is the vector of lower bounds (rectangular spaces)
     upper is the vector of upper bounds

     X is then the matrix in which the points will be stored.  X is an n by p
     matrix, with one point in each row.
  */
{
  designchoice = 0;
  design = userdesign;
}

void randfunc::genLH() // generates a Latin Hypercube design and stores it in X.
  /* This code was modified from C.M. Sieferts LatinHCDesign class 
     6/1/00 */
{
  Vector<long> a(n);/*Partition of a given dimension*/
  long s,temp;
  double int_width;/*width of interval*/
  long i,j; /*counters*/
  
  SelectStream(rngstream+1);  

  for(i=0;i<p;i++) { /*loop on dimension*/
    int_width = ((*upper)[i]-(*lower)[i])/n;    
    for(j=0;j<n;j++) a[j]=j;/*clear array*/

    /*Shuffle the array*/
    for(j=0;j<n;j++) {
      s=Equilikely(j,n-1);
      temp=a[s];
      a[s]=a[j];
      a[j]=temp;
    }/*end for*/
    /*END SHUFFLE*/

    for(j=0;j<n;j++) 
      X[j][i]=(*lower)[i]+Uniform(a[j]*int_width,(a[j]+1)*int_width);
  
  }/*end for dimension*/
  return;
}

void analyzevalues( randfunc& Func, ostream& outf) {
  analyzevalues( Func, outf, 25);
  return;
}
/****
     Works much like scanfunction, except that it reports simple statistics on
     the points sampled to standard I/O, such as maximum value, minimum value,
     average value, and the approximate location of the 25th and 75th
     percentiles.

     Warning:  This function stores all sampled points in memory.  In higher
     dimensions, this may fill up the available memory.  Use at your own risk.
****/


void analyzevalues( randfunc& Func, ostream& outf, int pointsperaxis)
  /*Accepts a randfunc object and two filenames.  Will create output files using
    the filenames.  filename1 will be the result of scanfunction().  filename 2
    will contain an analysis of these results.  While not exact, this analysis
    should provide the user with an idea of the range of values the function
    assumes.

  */
{
  double max, min, avg=0;
  long psize = Func.getp(), i, numpoints = 0, numpoints2 = 0, count = 0;
  Vector<double> current(1), maxvec(1), minvec(1);

  long p =psize;
  Vector<double> x(p);
  Vector<double> lower(p);
  Vector<double> upper(p);
  Vector<double> increment(p);
  int j;
  double temp = 1.0/pointsperaxis;
  Matrix<double> points (1,1);
  
  lower = Func.getlower();
  upper = Func.getupper();
  x = lower;
  increment = upper-lower;
  increment = increment * temp;
  temp = pow((double)pointsperaxis+2.0, (int)p);
  points.newsize((long)temp,p+1);
  
  while(( x[p-1] < (upper[p-1]+increment[p-1]/2.0) )&& (count<(long)temp)) {

    for(int k = 0; k <=p; k++) {
      points[count][k] = x[k];
      outf << x[k] << ' ';
    } // end for
    points[count][p] = Func.evalf(x);
    outf << points[count][p] << endl;
    count++;
    
    x[0] += increment[0];
    for( j = 0; (j < p-1)&&(x[j] > (upper[j]+increment[j]/2.0)); j++) {
      x[j] = lower[j];
      x[j+1] += increment[j];
    } // end for
  } // end while
  temp = count;
  
  // Process results of scanfunction()
  current.newsize(psize+1);
  maxvec.newsize(psize+1);
  minvec.newsize(psize+1);  

  count = 1;
  current = points.row(count);
  count++;
  avg = avg + current[psize];
  numpoints++;
  maxvec = current;
  minvec = current;
  max = maxvec[psize];
  min = minvec[psize];
  
  while( count < temp ) {
    current = points.row(count);
    count++;
    
    avg = avg+current[psize];
    numpoints++;
    if( current[psize] > max ) {
      maxvec = current;
      max = maxvec[psize];
    } else if ( current[psize] < min ) {
       minvec = current;
       min = minvec[psize];
    } // end elseif
  } // end while

  avg = avg/numpoints;

  cout <<  "\nMaximum value of " << max << " found at (" << maxvec[0];
  for( i=1; i<psize; i++)
    cout << ", " << maxvec[i];
  cout << ")" << endl;
  cout <<  "Minimum value of " << min << " found at (" << minvec[0];
  for( i=1; i<psize; i++)
    cout << ", " << minvec[i];
  cout << ")" << endl;
  
  cout << "The average value of the points sampled was " << avg << endl;
  
  count = 1;

  numpoints = 0;
  numpoints2 = 0;
  max = 0;
  min = 0;
  
  current = points.row(count);
  count++;
  if( current[psize] >= avg) {
    max += current[psize];
    numpoints++;
  } else {
    min += current[psize];
    numpoints2++;
  }
    
  while( count < temp) {
    current = points.row(count);
    count++;
    if( current[psize] >= avg) {
      max += current[psize];
      numpoints++;
    } else {
      min += current[psize];
      numpoints2++;
    } // end else
  } // end while

  max = max/numpoints;
  min = min/numpoints2;

  cout << "\n25th percentile is approximately " << min << endl;
  cout << "75th percentile is approximately " << max << endl;

  return;
}  // end analyzevalues()




