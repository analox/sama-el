/****
ParamEstimate Class - 6/27/00 cmsief
This wonderfully useful class should take care of all of your parameter
estimation needs.     
****/

#ifndef _MAPS_PARAM_
#define _MAPS_PARAM_

#include "maps_general.h"
#include "PatternSearch.h"
#include "f2c.h"

using namespace std;

extern "C" {
  int dgesvd_(char *, char *, long*, long*, double *, long *, double*,
              double *, long *, double *, long *, double *, long *, long *);
  int dpotrf_(char*, long*, double*, long*, long*);
  int dpotri_(char*, long*, double*, long*, long*);
};
/*CLAPACK routines*/

#define START_THETA_STEP (0.1 * GetThetaUpperBound(curr_delta))
#define STOP_THETA_STEP (0.000001 * GetThetaUpperBound(curr_delta))

/*correlation family numbers*/
/*There are to be used to pick your correlation function type.
 See the private member functions for more details.*/
#define COR_EXP_ISOTROPIC 0
#define COR_GAUSS_ISOTROPIC 1
#define COR_EXP_PRODUCT 2
#define COR_GAUSS_PRODUCT 3
#define COR_CUBIC_ISOTROPIC 4

class ParamEstimate {

public:
  ParamEstimate();
  ParamEstimate(ParamEstimate &PE);
  ParamEstimate(long dim);
  /*Prep work and call SetDefault. This is the default.*/
  ParamEstimate(int correlate, long dim, PatternSearch &PS1);
  /*WARNING: REQUIRES A FULLY INITIALIZED PATTERNSEARCH*/
  ~ParamEstimate();
  ParamEstimate& operator=(ParamEstimate &PE);
  

  /*****RESET FUNCTIONS*****/

  void ResetParamEstimateState(int correlate, long dim, PatternSearch &PS1);
  /****
  INPUT: Correlation family number (see above), dimension of space, and a PatternSearch.
  EFFECT: Resets all the major state variables.
  ****/
  
  void ResetPatternSearch(PatternSearch &PS1);
  /****
  INPUT: PatternSearch Object.
  EFFECT: Sets PS to be this PatternSearch.
  ****/
  
  bool ResetCorrelationFamily(int correlate);
  /****
  INPUT: Correlation family number (see above).
  OUTPUT: True if number is valid.
  EFFECT: Sets function pointer r to correct function.
  ****/

  bool ResetDim(long dim);
  /****
  INPUT: Dimension of the space.
  OUTPUT: True if dimension is valid.
  EFFECT: P=dim;
  ****/

  /****QUERY FUNCTIONS****/
  
  long GetDim() const;
  /****
  OUTPUT: Returns the dimension of the space.
  ****/

  long GetThetaSize() const;
  /****
  OUTPUT: Returns the size of the theta vector
  ****/     
  
  void debug(char *code) const;
  /****
  EFFECT: Does a complete state dump to stdout.
  ****/ 

  static double GetThetaUpperBound(double delta);
  /****
  INPUT: delta.     
  OUTPUT: Returns a double representing the upper bound for theta = -log(0.8) / curr_delta^2.
  ****/ 
  
  /*****ESTIMATION FUNCTIONS*****/
  
  bool EstimateConstantTrend(pt_collect &pts, mp_vector &fvals, double delta, double* beta,
                             double* sigma2, mp_vector* theta, pt_collect*
                             MLE_Matrix, mp_vector* v);
  /****
  INPUT: Collection of evaluated points, their values, delta , and pointers for the
  parameters to be returned.  theta should come in with a 'guess' value for the
  optimizer to chew on.
  OUTPUT: Returns true if the operation was successful.
  EFFECT: This estimates the parameters, assuming a constant trend, and returns them.
  ****/
  
  bool EstimateQuadraticTrend(pt_collect &pts, mp_vector &fvals, double delta,
                              mp_vector* beta,                              
                              double* sigma2, mp_vector* theta, pt_collect*
                              MLE_Matrix, mp_vector* v);
  /****
  INPUT: Collection of evaluated points, their values, and pointers for the
  parameters to be returned. theta should come in with a 'guess' value for the
  optimizer to chew on. 
  OUTPUT: Returns true if the operation was successful.
  EFFECT: This estimates the parameters, using a quadratic trend constructed
  using generalized least squares.  Then it returns the parameters.   
  ****/
 
  bool EstimateCustomTrend(pt_collect &pts, mp_vector &fvals, double delta, void
                           (*fcn2)(const mp_vector &, mp_vector&), long
                           dimf, mp_vector* beta, double* sigma2, mp_vector* theta,
                           pt_collect* MLE_Matrix, mp_vector* v);
  /****
  INPUT: Collection of evaluated points, their values, and pointers for the
  parameters to be returned.  There is also a function pointer for the
  'a' function, and the dimension of the vector it returns.  theta should come
  in with a 'guess' value for the optimizer to chew on.
  OUTPUT: Returns true if the operation was successful.
  EFFECT: This estimates the function.
  ****/ 

  friend void OptimizeMLEConstant(long dimtheta, mp_vector &x, double &f, bool &success,void*PE) 
    {(*(ParamEstimate*)PE).MLEConstant(dimtheta,x,f,success);}
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Indirection to get around C++'s inability to let function pointer
  point to class member functions.  This gets optimized by PatternSearch.
  ****/

  friend void OptimizeMLEQuadratic(long dimtheta, mp_vector &x, double &f, bool &success,void *PE) 
    {(*(ParamEstimate*)PE).MLEQuadratic(dimtheta,x,f,success);}
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Indirection to get around C++'s inability to let function pointer
  point to class member functions.  This gets optimized by PatternSearch.
  ****/

  friend void OptimizeMLECustom(long dimtheta, mp_vector&x, double &f, bool &success, void *PE)
  {(*(ParamEstimate*)PE).MLECustom(dimtheta,x,f,success);}
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Indirection to get around C++'s inability to let function pointer
  point to class member functions.  This gets optimized by PatternSearch.
  ****/

  double EvaluateMyCorrelation(long dim, const mp_vector &x1, const mp_vector &x2,
                               const mp_vector &theta);
  /****
  INPUT: Dimesion of space, two points, and a 'vector' theta.  For the Isotropic functions, this
  theta is a one dimensional vector, for the Product functions, this is a
  P-dimensional vector.
  OUTPUT: Correlation between the two points.
  ****/       

private:
  
  void SetDefault(long dim);
  /****
  INPUT: Dimension of the space.
  EFFECTS: This will use a happy default set of values to set up the parameter
  estimation. 
  ****/ 

  /*For the Optimization*/

  void MLEConstant(long dimtheta, mp_vector &x, double &f, bool &success);
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Evaluates the MLE for a given theta.  Uses the constant trend.  This
  is to be optimized by the Pattern Search, via a little bit of indirection.
  ****/

  void MLEQuadratic(long dimtheta, mp_vector &x, double &f, bool &success);
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Evaluates the MLE for a given theta.  Uses the quadratic trend.  This
  is to be optimized by the Pattern Search, via a little bit of indirection.
  ****/  

  void MLECustom(long dimtheta, mp_vector &x, double &f, bool &success);
  /****
  INPUT: PatternSearch calling scheme.  See PatternSearch.h
  EFFECT: Evaluates the MLE for a given theta.  Uses the custom trend.  This is
  to be optimized by the Pattern Search, via a little bit of indirection.
  ****/
  
  bool pseudoinvert(pt_collect &M, double *log_det);
  /****
  INPUT: A matrix M, and a pointer to a mp_vector.
  OUTPUT: true if svd can be computed, else false.
  EFFECT: Replaces M with the inverse (if possible) or pseudoinverse of M.  If
  d!=NULL, an mp_vector with the one over the singular values, with those singular
  values below TOLERENCE instead to zero, is returned.  The matrix has to be
  *SYMMETRIC AND SQUARE* for this to work.
  ****/

  inline void GenerateCorrelationMatrix(pt_collect &R, const mp_vector &theta);
  /****
  INPUT: A np x np pt_collect, and the correlation family parameter vector, theta.
  EFFECT: Replaces R with the appropriate correlation matrix.
  ****/
  
  bool ParamEstimate::ComputeCholeskyInverse(pt_collect *A,double* log_det);
  /****
  INPUT: Matrix A,  which should be positive definite, of which to find the
  inverse using the Cholesky factorization and a pointer to the matrix in which to
  store the inverse INV.  A * INV == I.  Also a pointer to a double.
  OUTPUT: true if CLAPACK's dpotri_ returns OK.
  EFFECT: Calls CLAPACK's dpotrf_ routine to compute the Cholesky factorization
  of A, and then calls doprti_ to find the inverse from it.  Computes the
  log determinant of A iff log_det!=NULL, and stores it in log_det.
  Based roughly on Anthony Padula's code, modified by Chris Siefert.
  ****/ 

  
  bool ComputeSVDInverse(pt_collect *A, double *log_det);
  /****
  INPUT: Matrix A, to take the SVD of and a pointer for the log_det.
  OUTPUT: true if CLAPACK's dgesvd_ returns ok.
  EFFECT: Calls CLAPACK's dgesvd routine to compute the Singular Value
  Decomposition(SVD) of A.  Recall that a SVD of A yields U, D, and V',
  s.t. A=U*D* V', where U and V are orthogonal matrices, and D is a diagonal
  matrix with the singular values on its main diagonal.  The pseudoinverse is
  computed, (Ai = V*Di*U') using the tolerance from Matlab. The log determinant is
  calculated if log_det!=NULL, and it is stored there.   
  ****/  

/*
  Correlation Functions
  INPUT: Dimesion of space, two points, and a 'vector' theta.  For the Isotropic functions, this
  theta is a one dimensional vector, for the Product functions, this is a
  P-dimensional vector.
  NOTE: You specify these by use of the integer #defined correlation family
  numbers, and ResetCorrelationFamily().  These are static class member functions.
*/

static double Cor_Exp_Isotropic(long dim, const mp_vector &x1, const mp_vector &x2,
                                  const mp_vector &theta);
  /*Exponential Isotropic Correlation Family
    r_theta(s,t) = exp(-theta ||s-t||)
  */
static double Cor_Gauss_Isotropic(long dim, const mp_vector &x1, const mp_vector &x2,
                           const mp_vector &theta);
  /*Gaussian Isotropic Correlation Family
    r_theta(s,t) = exp(-theta ||s-t||)^2
  */  
static double Cor_Exp_Product(long dim, const mp_vector &x1, const mp_vector &x2, const
                       mp_vector &theta);
  /*Exponential Product Correlation Family
    r_theta(s,t) = \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|)
  */
static double Cor_Gauss_Product(long dim, const mp_vector &x1, const mp_vector &x2,
                         const mp_vector &theta);
  /*Gaussian Product Correlation Family
    r_theta(s,t) = \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|^2)
  */
static double Cor_Cubic_Isotropic(long dim, const mp_vector &x1, const mp_vector &x2,
                         const mp_vector &theta);
  /*Cubic Product Correlation Family
    r_theta(s,t) = { 1-1.5*(||s-t||)/theta + 0.5 ((||s-t||)/theta)^3, for ||s-t|| < theta }
                   { 0 otherwise }
  */
  
private:
  
  pt_collect* points;   /*pointer to the evaluated points*/
  mp_vector* values;    /*pointer to function values*/
  long np;              /*number of point evaluated*/

  PatternSearch* PS;    /*some kind of optimizer*/
  long P;               /*dimension of the space*/
  double (*r)(long dim,const mp_vector&, const mp_vector&, const mp_vector& );
  /*for pointing to the correlation function*/
  long fdim;            /*dimension of the vector returned by fcn*/
  void (*fcn)(const mp_vector &, mp_vector&);
  /*for pointing to the 'A' function in custom estimation*/
  double curr_delta;    /*Current delta. Used for bounds on theta*/
  double MACH_EPS;      /*Machine Epsilon*/
};/*end class*/

double GetMachineEpsilon();

#endif
