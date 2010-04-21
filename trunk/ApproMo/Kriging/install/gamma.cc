
/* This is the source for several transcendental mathematical functions.
   The source code was copied from
   Lau, H. T.  Numerical library in C for scientists and engineers
   1995 by CRC Press Inc.
   pgs 539 - 543

   We also have several gateway functions to CLAPACK for computing SVD's, 
   Cholesky factorizations, finding inverses, and solving systems of equations.
*/

#include "gamma.h"


#ifndef COMPSVD_FD
#define COMPSVD_FD

extern "C" {
  int dgesvd_(char*,char*,long*,long*,double*,long*,double*,double*,long*,
            double*,long*,double*,long*, long*);
  int dgesv_(long *, long*, double*, long*, long*, double*, long*, long*);
  int dpotrf_(char*, long*, double*, long*, long*);
  int dpotri_(char*, long*, double*, long*, long*);
};

bool SolveSystem( const Matrix<double>& A, Vector<double>& x, Vector<double>&
                  b) {
   /***
       INPUT: square matrix A, vector b, and vector x through which the answer will be returned
       OUTPUT: true if clapack returns OK
       EFFECT: Calls CLAPACK's dgesv function, to solve the linear system Ax=b,
       and stores the solution in x.
   ***/
  if(A.num_rows()!=A.num_cols()) return false;

  long N = A.num_rows();
  long NRHS = 1;
  double* MA=(double*)malloc(sizeof(double)*N*N);
  long LDA = N;
  long* IPIV=(long*)malloc(sizeof(long)*N*N);
  double* MB=(double*)malloc(sizeof(double)*N);
  long LDB = N;
  long INFO;

  for(long i=0;i<N;i++) {
    for(long j=0;j<N;j++)
      MA[N*j+i]=A[i][j]; /*This should flip the Matrix into column-major order
                           ala Fortran.*/

    MB[i] = b[i];
  } // end for
   
  dgesv_(&N,&NRHS,MA,&LDA,IPIV,MB,&LDB, &INFO);

  if(INFO<0) {free(MA);free(IPIV);free(MB);return false;}
#if DEBUG>0
  else if (INFO>0) cerr<<"SOLVE : Warning, System is not positive definite\n";
#endif
  
  /*Copy data to all the outgoing matrices*/
  for(long j=0;j<N;j++) {
    x[j] = MB[IPIV[j]];
  }/*end for*/

#if DEBUG>=4
  cout<<"SOLVE :x="<<(x);
#endif  
  /*Free memory!*/
  free(MA);free(MB); free(IPIV);
  return true;  
}/*end SolveSystem() */


bool ComputeSVD(const pt_collect &A, pt_collect* U, pt_collect* D,
                pt_collect* VT)
{
/****
INPUT: Matrix A, to take the SVD of, and pointers to three matrices that will
contain what is returned from dgesvd_.
OUTPUT: true if CLAPACK's dgesvd_ returns OK.
EFFECT: Calls CLAPACK's degsvd routine to compute the Singular Value
Decomposition(SVD) of A.  Recall that a SVD of A yields U, D, and V',
s.t. A=U*D* V', where U and V are orthogonal matrices, and D is a diagonal
matrix with the singular values on its main diagonal.
NOTES: 6/1/99 TESTED OK
****/  
  if(A.num_rows()!=A.num_cols()) return false;
  if(A.num_rows()!= U->num_rows()) return false;
  if(A.num_rows()!= U->num_cols()) return false;
  if(A.num_rows()!= VT->num_rows()) return false;
  if(A.num_rows()!= VT->num_cols()) return false;
  if(A.num_rows()!= D->num_rows()) return false;
  if(A.num_rows()!= D->num_rows()) return false;

  
  long md=A.num_rows();
  long idx, LWORK=md*md+5*md;
  double* MA=(double*)malloc(sizeof(double)*md*md);
  double* MU=(double*)malloc(sizeof(double)*md*md);  
  double* MVT=(double*)malloc(sizeof(double)*md*md);
  double* MS=(double*)malloc(sizeof(double)*md);
  double* WORK=(double*)malloc(sizeof(double)*LWORK);
  char ch='A';

  for(long i=0;i<md;i++)
    for(long j=0;j<md;j++)
      MA[md*j+i]=A[i][j]; /*This should flip the Matrix into column-major order
                            ala Fortran.*/
      
  dgesvd_(&ch,&ch,&md,&md,MA,&md,MS,MU,&md,MVT,&md,WORK,&LWORK,&idx);

  if(idx<0) {free(MA);free(MU);free(MVT);free(MS);free(WORK);return false;}
#if DEBUG>0
  else if (idx>0) cerr<<"SVD : WARNING - SVD has "<<idx<<" unconverged superdiagonal values\n";
#endif
  
  /*Copy data to all the outgoing matrices*/
  for(long j=0;j<md;j++) {/*cols*/
    idx=md*j;/*Partial computation of index... for speed*/
    for(long i=0;i<md;i++) {/*rows*/
      (*U)[i][j]=MU[idx+i];
      (*VT)[i][j]=MVT[idx+i];
      if(i==j) (*D)[i][j]=MS[i];
      else (*D)[i][j]=0.0;
    }/*end for*/
  }/*end for*/

#if DEBUG>=4
  cout<<"SVD :A="<<(A);
  cout<<"SVD :U="<<(*U);
  cout<<"SVD :D="<<(*D);
  cout<<"SVD :VT="<<(*VT);
  cout<<"SVD :U*D*VT="<<((*U)*(*D)*(*VT));
#endif  
  /*Free memory!*/
  free(MA);free(MU);free(MVT);free(MS);free(WORK);
  return true;  
}/*end ComputeSVD*/



bool ComputeCholesky( const pt_collect & A, pt_collect* L)
 /****
      INPUT: Matrix A, which should be positive definite, 
      and of which to take the Cholesky factorization and a pointer
      to the matrix in which to store L.  A = L L'
      OUTPUT: true if CLAPACK's dpotrf returns OK.
      EFFECT: Calls CLAPACK's dpotrf routine to compute the Cholesky factorization
      of A.  Computes A = L L' where L is lower triangular
      based on Chris Siefert's ComputeSVD code, modified by Anthony Padula
      
 ****/
{
 if(A.num_rows()!=A.num_cols()) return false;
  
  long md=A.num_rows(), idx;
  double* MA=(double*)malloc(sizeof(double)*md*md);
  char ch='L';

  for(long i=0;i<md;i++)
    for(long j=0;j<md;j++)
      MA[md*j+i]=A[i][j]; /*This should flip the Matrix into column-major order
                            ala Fortran.*/
      
  dpotrf_(&ch, &md, MA, &md, &idx); 

  if(idx<0) {free(MA);return false;}
  else if (idx>0) {
    if( DEBUG > 0)
      cerr<<"CHOLESKY : WARNING - A is not positive definite at row "
	  <<idx<<endl;
    free(MA);
    return false;
  }
  
  /*Copy data to all the outgoing matrices*/
  for(long j=0;j<md;j++) {/*cols*/
    idx=md*j;/*Partial computation of index... for speed*/
    for(long i=0;i<md;i++) {/*rows*/
      if( i < j )
	(*L)[i][j] = 0;
      else
	(*L)[i][j]=MA[idx+i];
    }/*end for*/
  }/*end for*/

#if DEBUG>=4
  cout<<"CHOLESKY :A="<<(A);
  cout<<"CHOLESKY :L="<<(*L);
  
#endif  
  /*Free memory!*/
  free(MA);
  return true;  
}/*end ComputeCholesky */

bool ComputeInverse( const pt_collect & A, pt_collect* INV)
 /****
      INPUT: Matrix A,  which should be positive definite,
      of which to find the inverse using the Cholesky
      factorization and a pointer to the matrix in which to store the inverse INV.  
      A * INV == I
      OUTPUT: true if CLAPACK's dpotri_ returns OK.
      EFFECT: Calls CLAPACK's dpotrf_ routine to compute the Cholesky factorization
      of A, and then calls doprti_ to find the inverse from it.
      Based roughly on Chris Siefert's ComputeSVD code, modified by Anthony Padula
      6/1/00
 ****/
{
 if(A.num_rows()!=A.num_cols()) return false;
  
  long md=A.num_rows(), idx;
  Matrix<double> LINV(1,1);
  double* MA=(double*)malloc(sizeof(double)*md*md);
  char ch='U';

  for(long i=0;i<md;i++)
    for(long j=0;j<md;j++)
	MA[md*j+i]=A[i][j]; /*This should flip the Matrix into column-major order
			      ala Fortran.*/
      
  dpotrf_(&ch, &md, MA, &md, &idx); 
  if(idx==0)
    dpotri_(&ch, &md, MA, &md, &idx); 

  if(idx<0) {free(MA);return false;}
  else if (idx>0) { 
    cerr<<"INVERSE : WARNING - A is not positive definite at row "
	<<idx<<endl;
    free(MA);
    return false;
  }

  /*Copy data to all the outgoing matrices*/
  for(long j=0;j<md;j++) {/*cols*/
    idx=md*j;/*Partial computation of index... for speed*/
    for(long i=0;i<=j;i++) {/*rows*/
	(*INV)[i][j]=MA[idx+i];
	(*INV)[j][i]=MA[idx+i];
    }/*end for*/
  }/*end for*/

#if DEBUG>=4
  cout<<"INVERSE :A="<<(A);
  cout<<"INVERSE :INV="<<(*INV);

  // check
  Matrix<double> I(1,1);
  I.newsize(md,md);
  I = 0;
  I = (A * (*INV));
  cout << "INVERSE :A * INV = " << I;
#endif  

  /*Free memory!*/
  free(MA);
  return true;  
}/*end ComputeInverse */

#endif


double recipgamma(double x, double *odd, double *even) {
/*  Delivers 1/Gamma(1-x);

    x:  double;
    entry:  this argument should satisfy -0.5 <= x <= 0.5;
       (actually the gamma function is calculated for 1-x, ie. if one wants to
       calculate 1/Gamma(1), one has to set x to 0);

    odd:  double *;
    exit:   the odd part of 1/Gamma(1-x) divided by (2x);
    ie. (1/Gamma(1-x) - 1/Gamma(1+x)) / (2x);

    even:  double *;
    exit:   the even part of 1/Gamma(1-x) divided by (2x);
    ie. (1/Gamma(1-x) + 1/Gamma(1+x)) / (2x);

*/

  int i;
  double alfa, beta, x2, b[13];

  b[1] = -0.283876542276024;  b[2] = -0.076852840844786;
  b[3] =  0.001706305071096;  b[4] =  0.001271927136655;
  b[5] =  0.000076309597586;  b[6] = -0.000004971736704;
  b[7] = -0.000000865920800;  b[8] = -0.000000033126120;
  b[9] =  0.000000001745136;  b[10]=  0.000000000242310;
  b[11]=  0.000000000009161;  b[12]= -0.000000000000170;
  x2=x*x*8.0;
  alfa = -0.000000000000001;
  beta = 0.0;
  for(i=12; i>=2; i -= 2) {
    beta = -(alfa*2.0+beta);
    alfa = -beta*x2-alfa+b[i];
  }
  *even=(beta/2.0+alfa)*x2-alfa+0.921870293650453;
  alfa = -0.000000000000034;
  beta = 0.0;
  for(i=11; i>=1; i -= 2) {
    beta = -(alfa*2.0+beta);
    alfa = -beta*x2-alfa+b[i];
  }
  *odd = (alfa+beta)*2.0;
  return (*odd)*x+(*even);
} // end recipgamma()


double gamma(double x) {
/* Delivers the value of the gamma function at x.

   x: double;
   entry:  the argument; if one of the following three conditions is fulfilled
   the overflow will occur:

   1. the argument is too large;
   2. the argument is a non-positive integer;
   3. the argument is too close to a large ( in absolute value) non-positive
   integer.

   functions used: recipgamma, loggamma.
*/

  int inv;
  double y,s,f=0.0,g,odd,even;
  if ( x < 0.5 ) {
    y = x - floor(x/2.0)*2;
    s = 3.14159265358979;
    if ( y >= 1.0 ) {
      s = -s;
      y = 2.0-y;
    }
    if ( y >= 0.5 ) y=1.0-y;
    inv=1;
    x=1.0-x;
    f=s/sin(3.14159265358979*y);
  } else
    inv = 0;
  if ( x > 22.0)
    g=exp(loggamma(x));
  else {
    s=1.0;
    while (x > 1.5) {
      x = x-1.0;
      s *= x;
    }
    g=s/recipgamma(1.0-x,&odd,&even);
  }
  return (inv ? f/g : g);
} // end gamma()


double loggamma( double x) {
/* delivers the value of the natural logarithm of the gamma function at x;

   x: double;
   entry:  this argument must be positive

*/

  int i;
  double r,x2,y,f,u0,u1,u,z,b[19];

  if (x > 13.0) {
    r = 1.0;
    while ( x <= 22.0) {
      r /= x;
      x += 1.0;
    }
    x2 = -1.0/(x*x);
    r=log(r);
    return log(x)*(x-0.5)-x+r+0.918938533204672+
      (((0.595238095238095e-3*x2+0.793650793650794e-3)*x2
        +0.277777777777778e-2)*x2+0.833333333333333e-1)/x;
  } else {
    f=1.0;
    u0=u1=0.0;
    b[1]  = -0.0761141616704358;  b[2]  = 0.0084323249659328;
    b[3]  = -0.0010794937263286;  b[4]  = 0.0001490074800369;
    b[5]  = -0.0000215123998886;  b[6]  = 0.0000031979329861;
    b[7]  = -0.0000004851693012;  b[8]  = 0.0000000747148782;
    b[9]  = -0.0000000116382967;  b[10] = 0.0000000018294004;
    b[11] = -0.0000000002896918;  b[12] = 0.0000000000461570;
    b[13] = -0.0000000000073928;  b[14] = 0.0000000000011894;
    b[15] = -0.0000000000001921;  b[16] = 0.0000000000000311;
    b[17] = -0.0000000000000051;  b[18] = 0.0000000000000008;
    if(x< 1.0) {
      f=1.0/x;
      x+=1.0;
    } else
      while ( x > 2.0) {
        x -=1.0;
        f *=x;
      }
    f = log(f);
    y=x+x-3.0;
    z=y+y;
    for (i=18; i>=1; i--) {
      u=u0;
      u0=z*u0+b[i]-u1;
      u1=u;
    }
    return (u0*y+0.491415393029387-u1)*(x-1.0)*(x-2.0)+f;
  }
} // end loggamma


void bessk01( double x, double *k0, double *k1) {
/* computes the value of the modified Bessel function of the third kind of order
   0 and 1, ie. K_0(x) and K_1(x)

   x: double;
   entry:  the argument of the Bessel functions; x>0;

   k0: double *;
   exit:   k0 has the value of the modified Bessel function of the third kind of
   order zero with argument x;

   k1: double *;
   exit:   k1 has the value of the modified Bessel function of the third kind of
   order one with argument x;
*/

  if (x <= 1.5) {
    int k;
    double c,d,r,sum0,sum1,t,term,t0,t1;
    sum0=d=log(2.0/x)-0.5772156649015328606;
    sum1 = c = -1.0 - 2.0*d;
    r=term=1.0;
    t=x*x/4.0;
    k=1;
    do {
      term *= t*r*r;
      d += r;
      c -= r;
      r = 1.0/(k+1);
      c -= r;
      t0=term*d;
      t1=term*c*r;
      sum0 += t0;
      sum1 += t1;
      k++;
    } while (fabs(t0/sum0)+fabs(t1/sum1) > 1.0e-15);
    *k0 = sum0;
    *k1 = (1.0+t*sum1)/x;
  } else {
    double expx;
    expx = exp(-x);
    nonexpbessk01(x,k0,k1);
    *k1 *= expx;
    *k0 *= expx;
  }
}

void bessk( double x, int n, double k[] ) {
/* Generates an array of modified Bessel functions of the third kind of order j,
   K_j(x).

   Calls bessk01 to generate the order zero and order one terms, then computes
   higher order terms by recursion:

   K_{j+1}(x) = K_{j-1}(x)+(2*j/x)*K_{j}(x)


   x: double;
   entry:  the argument of the Bessel functions; x>0

   n: int;
   entry:  the upper bound of the indices of array k; n>=0;

   k: double k[0:n];
   exit:   k[j] is the value of the modified Bessel function of the third kind
   of order j with argument x, j=0,...,n.
*/

  int i;
  double k0,k1,k2;

  bessk01(x,&k0,&k1);
  k[0]=k0;
  if( n > 0 ) k[1] = k1;
  x=2.0/x;
  for (i=2; i<=n; i++) {
    k[i]=k2=k0+x*(i-1)*k1;
    k0 = k1;
    k1 = k2;
  }
}
  
void nonexpbessk01( double x, double *k0, double * k1) {
/* computes the value of the modified Bessel function of the third kind of order
   0 and 1, ie. K_0(x) and K_1(x), multiplied by e^x.

   x: double;
   entry:  the argument of the Bessel functions; x>0;

   k0: double *;
   exit:   k0 has the value of the modified Bessel function of the third kind of
   order zero with argument x, multiplied by e^x;

   k1: double *;
   exit:   k1 has the value of the modified Bessel function of the third kind of
   order one with argument x, multiplied by e^x;
*/
  if (x <= 1.5) {
    double expx;
    expx = exp(x);
    bessk01(x,k0,k1);
    *k0 *= expx;
    *k1 *= expx;
  } else if (x <= 5.0) {
    int i,r;
    double t2,s1,s2,term1,term2,sqrtexpr,exph2,x2;
    static double fac[20] = {0.90483741803596, 0.67032004603564,
                             0.40656965974060, 0.20189651799466,
                             0.82084998623899e-1, 0.27323722447293e-1,
                             0.74465830709243e-2, 0.16615572731739e-2,
                             0.30353913807887e-3, 0.45399929762485e-4,
                             0.55595132416500e-5, 0.55739036926944e-6,
                             0.45753387694459e-7, 0.30748798795865e-8,
                             0.16918979226151e-9, 0.76218651945127e-11,
                             0.28111852987891e-12, 0.84890440338729e-14,
                             0.2098791048793e-15, 0.42483542552916e-17};
    s1=0.5;
    s2=0.0;
    r=0;
    x2=x+x;
    exph2=1.0/sqrt(5.0*x);
    for( i=0; i<=19; i++) {
      r += 1;
      t2 = r*r/10.0;
      sqrtexpr=sqrt(t2/x2+1.0);
      term1=fac[i]/sqrtexpr;
      term2=fac[i]*sqrtexpr*t2;
      s1 += term1;
      s2 += term2;
    }
    *k0 = exph2*s1;
    *k1 = exph2*s2*2.0;
    } else {
      int r,i;
      double br,br1,br2,cr,cr1,cr2,ermin1, erplus1, er,f0,f1,expx,y,y2;
      static float dr[14]={0.27545e-15, -0.172697e-14, 0.1136042e-13,
                           -0.7883236e-13, 0.58081063e-12, -0.457993633e-11,
                           0.3904375576e-10, -0.36454717921e-9,
                           0.379299645568e-8, -0.450473376411e-7,
                           0.63257510850049e-6, -0.11106685196665e-4,
                           0.26953261276272e-3, -0.11310504646928e-1};
      y=10.0/x-1.0;
      y2 = y+y;
      r=30;
      br1=br2=cr1=cr2=erplus1=er=0.0;
      for (i=0; i<=13; i++) {
        r -=2;
        br=y2*br1-br2+dr[i];
        cr=cr1*y2-cr2+er;
        ermin1=r*dr[i]+erplus1;
        erplus1=er;
        er = ermin1;
        br2=br1;
        br1=br;
        cr2=cr1;
        cr1=cr;
      }
      f0=y*br1-br2+0.9884081742308258;
      f1=y*cr1-cr2+er/2.0;
      expx=sqrt(1.5707963267949/x);
      *k0 = f0 *= expx;
      *k1 = (1.0+0.5/x)*f0+(10.0/x/x)*expx*f1;
    }
}

// End Bessel and Gamma function definitions


double MachineEpsilon() {
  // Uses formula in  Heath:  |3*(4/3-1)-1|
    double tol = 4.0/3.0;
    tol = tol-1.0;
    tol = 3.0*tol;
    tol = tol-1.0;
    tol = fabs(tol);
    return tol;
}
