/****
ParamEstimate Implementation - 6/29/00 cmsief
It slices, it dices, it does singular value decomposition!
****/

#include "maps_general.h"
#include "PatternSearch.h"
#include "CompassSearch.h"
#include "f2c.h"
#include "ParamEstimate.h"
#include <cmath>          /*e^whatever and sqrt*/
#include <malloc.h>
#include <cstdlib>
#include <iostream>

// #define COR_EXP_ISOTROPIC 0
// #define COR_GAUSS_ISOTROPIC 1
// #define COR_EXP_PRODUCT 2
// #define COR_GAUSS_PRODUCT 3
// #define COR_CUBIC_ISOTROPIC 4



/**********************************CONSTRUCTORS + DESTRUCTOR**************************************/

ParamEstimate::ParamEstimate() {
  PS=NULL;
  P=0;
  np=0;
  r=NULL;
  fdim=0;
  fcn=NULL;
  curr_delta=0;
  points=NULL;
  values=NULL;
  MACH_EPS=GetMachineEpsilon();
}/*end default constructor*/

ParamEstimate::ParamEstimate(ParamEstimate &PE){
  P=PE.P;
  r=PE.r;/*I know that I'm doing this*/
  np=PE.np;
  fcn=PE.fcn;
  fdim=PE.fdim;
  mp_vector temp(fdim);temp=0;
  PS=new CompassSearch(fdim,temp);
  curr_delta=PE.curr_delta;
  points=PE.points;/*I know that I'm doing this.*/
  values=PE.values;/*I know that I'm doing this.*/
  MACH_EPS=PE.MACH_EPS;
}/*end copy constructor*/

ParamEstimate::ParamEstimate(long dim) {
/*Prep work and call AutoGenerate*/
  PS=NULL;
  points=NULL;
  values=NULL;
  MACH_EPS=GetMachineEpsilon();  
  SetDefault(dim);
}/*end special constructor*/

ParamEstimate::ParamEstimate(int correlate, long dim, PatternSearch &PS1) {
  P=dim;
  points=NULL;
  values=NULL;
  np=0;
  fcn=NULL;
  fdim=0;
  curr_delta=0;
  PS=NULL;
  r=NULL;
  ResetCorrelationFamily(correlate);  
  ResetPatternSearch(PS1);
  MACH_EPS=GetMachineEpsilon();
}/*end special constructor*/

ParamEstimate::~ParamEstimate(){
  points=NULL;/*I know that I'm doing this*/
  values=NULL;
  if(PS!=NULL) {delete PS;PS=NULL;}//Am I going to lose memory here?
}/*end destructor*/

ParamEstimate& ParamEstimate::operator=(ParamEstimate &PE){
  if (this!=&PE) {
    if(PE.PS!=NULL) ResetPatternSearch(*PE.PS);
    else if(PS!=NULL) {delete PS;PS=NULL;}    
    P=PE.P;
    r=PE.r;
    fcn=PE.fcn;
    fdim=PE.fdim;
    np=PE.np;
    curr_delta=PE.curr_delta;
    points=PE.points;/*I know that I'm doing this.*/
    values=PE.values;/*I know that I'm doing this.*/
    MACH_EPS=PE.MACH_EPS;
  }/*end if*/
  return (*this);
}/*end overloaded=*/


/**************************************RESET FUNCTIONS*****************************************/

void ParamEstimate::ResetParamEstimateState(int correlate, long dim, PatternSearch &PS1){
/****
INPUT: Correlation family number, dimension of space, and a PatternSearch.
EFFECT: Resets all the major state variables.
****/
  ResetCorrelationFamily(correlate);
  ResetDim(dim);
  ResetPatternSearch(PS1);  
}/*end ResetParamEstimateState*/
  
void ParamEstimate::ResetPatternSearch(PatternSearch &PS1) {
/****
INPUT: PatternSearch Object
EFFECT: Sets PS to be this PatternSearch
****/
  mp_vector temp(fdim);temp=0;
  if(PS==NULL) PS=new CompassSearch(fdim,temp);
  else (*PS)=(PS1);  

}/*end ResetPatternSearch*/
  
bool ParamEstimate::ResetCorrelationFamily(int correlate){
/****
INPUT: Correlation family number
OUTPUT: True if number is valid
EFFECT: Sets function pointer r to correct function
****/
  switch(correlate) {
  case COR_EXP_ISOTROPIC:
    r=Cor_Exp_Isotropic;
    return true;   
  case COR_GAUSS_ISOTROPIC:
    r=Cor_Gauss_Isotropic;
    return true;    
  case COR_EXP_PRODUCT:
    r=Cor_Exp_Product;
    return true;       
  case COR_GAUSS_PRODUCT:
    r=Cor_Gauss_Product;
    return true;
  case COR_CUBIC_ISOTROPIC:
    r=Cor_Cubic_Isotropic;
    return true;
  default:
    return false;
  };/*end switch*/
}/*end ResetCorrelationFamily*/
  
bool ParamEstimate::ResetDim(long dim) {
/****
INPUT: Dimension of the space.
OUTPUT: True if dimension is valid.
EFFECT: P=dim;
****/
  if(dim>0) {P=dim;return true;}
  else return false;
}/* ResetDim*/

/***************************************QUERY FUNCTIONS**************************************/

long ParamEstimate::GetDim() const {
/****
OUTPUT: Returns the dimension of the space.
****/
  return P;
}/*end GetDim*/


long ParamEstimate::GetThetaSize() const {
/****
OUTPUT: Returns the size of the theta vector
****/     
  if(r==Cor_Exp_Isotropic || r==Cor_Gauss_Isotropic||r==Cor_Cubic_Isotropic) return 1;
  else if(r==Cor_Exp_Product || r==Cor_Gauss_Product) return P; 
  else return ERROR;
}/*end GetThetaSize*/

double ParamEstimate::GetThetaUpperBound(double delta){
/****
INPUT: delta.       
OUTPUT: Returns a double representing the upper bound for theta = -log(0.8) / delta^2.
****/
  return (-log(.8) / sqr(delta));
}/*end GetThetaUpperBound*/

void ParamEstimate::debug(char *code) const{
/****
EFFECT: Does a complete state dump to stdout.
****/
  cout<<code<<":***ParamEstimate::debug()***\n";
  cout<<code<<":P="<<P<<" curr_delta="<<curr_delta<<" fdim="<<fdim<<" np="<<np<<endl;   
  if(r==Cor_Exp_Isotropic) cout<<code<<":Exponential Isotropic Correlation\n";
  else if(r==Cor_Gauss_Isotropic) cout<<code<<":Gaussian Isotropic Correlation\n";
  else if(r==Cor_Exp_Product) cout<<code<<":Exponential Product Correlation\n";
  else if(r==Cor_Gauss_Product) cout<<code<<":Gaussian Product Correlation\n";
  else if(r==Cor_Cubic_Isotropic) cout<<code<<":Cubic Isotropic Correlation\n";
  else cout<<code<<":Correlation unknown\n";
  if(points!=NULL) cout<<code<<":points="<<(*points); else cout<<code<<":points=NULL\n";
  if(values!=NULL) cout<<code<<":values="<<(*values); else cout<<code<<":values=NULL\n";
  cout<<code<<":***end PE::debug()***\n";
}/*end debug*/
  

/**************************************ACTION FUNCTIONS**************************************/
  
bool ParamEstimate::EstimateConstantTrend(pt_collect &pts, mp_vector &fvals,
                                          double delta, double* beta,
                                          double* sigma2, mp_vector* theta, pt_collect*
                                          MLE_Matrix, mp_vector* v) {  
/****
INPUT: Collection of evaluated points, their values, and pointers for the
parameters to be returned.  theta should come in with a 'guess' value for the
optimizer to chew on.
OUTPUT: Returns true if the operation was sucessful.
EFFECT: This estimates the parameters, assuming a constant trend, and returns them.
****/
#if DEBUG>=2
  cout<<"PEECT  :Started\n";  
#endif
  bool success;
  curr_delta=delta;
  np=pts.num_rows();
  points=&pts;
  values=&fvals;
  pt_collect R(np,np);
  mp_vector A(np);
  A=1;/*fills A with 1's*/
  long i,j;/*counter*/

  /*Pattern Search w/bounds on OptimizeMLEConstant*/
  PS->CleanSlate(theta->dim(), *theta,START_THETA_STEP, STOP_THETA_STEP,
                 OptimizeMLEConstant, (void *) this);


  PS->BeginSearch();
  PS->GetMinPoint(*theta);
  
  for(i=0;i<theta->dim();i++)
    if((*theta)[i]<TOLERANCE) (*theta)[i]=TOLERANCE;

#if DEBUG>=2
  cout<<"PEECT   :Upperbound="<<GetThetaUpperBound(delta)<<" theta="<<*theta;
#endif  

  /*Generate the Correlation Matrix*/
  GenerateCorrelationMatrix(R,(*theta));

  /*Generate the MLE Matrix*/
  (*MLE_Matrix)(0,0)=0.0; 
  for(i=1;i<np+1;i++) {
    (*MLE_Matrix)(0,i)=A[i-1];
    (*MLE_Matrix)(i,0)=A[i-1];
    for(j=1;j<np+1;j++) 
      (*MLE_Matrix)(i,j)=R(i-1,j-1);
  }/*end for*/
#if DEBUG >=3  
  cout<<"PEECT  :MLE Matrix Constructed:"<<(*MLE_Matrix);
#endif

  success=pseudoinvert(R,NULL);
  if(!success) return (success); /*break out!*/
  success=pseudoinvert(*MLE_Matrix,NULL);
  if(!success) return (success); /*break out!*/  
#if DEBUG >=3  
  cout<<"PEECT  :Pseudoinversions complete\n";
#endif
  (*v)=A*R;/*this multiply assumes that A is actually AT */
  (*beta)=((*v) * (*values)) / ((*v)*A); /*a hack for the constant trend case*/
  (*v)=((*values) - A*(*beta));
  (*sigma2)=((*v) * R * (*v))/np; /*again, the special multiply assumes that the first v
                                    is vT*/
  (*v)=(*v)*R;/*special multiply which assumes v=vT*/

#if DEBUG >=2
  cout<<"PEECT  :terminating success="<<success<<"\n";
  cout<<"PEECT  :theta ="<<(*theta)[0]<<endl; 
  cout<<"PEECT  :beta ="<<(*beta)<<endl;
#endif

  return (success);
}/*end EstimateConstantTrend*/


bool ParamEstimate::EstimateQuadraticTrend(pt_collect &pts, mp_vector &fvals,
                                           double delta, mp_vector* beta,
                                           double* sigma2, mp_vector* theta,
                                           pt_collect* MLE_Matrix, mp_vector*
                                           v){
/****
INPUT: Collection of evaluated points, their values, and pointers for the
parameters to be returned. theta should come in with a 'guess' value for the
optimizer to chew on. 
OUTPUT: Returns true if the operation was successful.
EFFECT: This estimates the parameters, using a quadratic trend constructed
using generalized least squares.  Then it returns the parameters.   
****/

#if DEBUG>=2
  cout<<"PEEQT  :Started\n";  
#endif
  bool success;
  curr_delta=delta;
  np=pts.num_rows();
  points=&pts;
  values=&fvals;
  long i, j, g, h,z;/*counter*/

#if DEBUG>=3
  cout<<"PEMQ   :MLEQuadratic Called\n";
  cout<<"PEMQ   :trial_theta=(";
  for(i=0;i<dimtheta;i++) cout<<x[i]<<" ";
  cout<<")\n";    
#endif  

  double sum=0.0;
  pt_collect R(np,np);
  long k = 1+P+(P*(P+1))/2;
  mp_vector tempv(k);
  pt_collect X(np,k), tempmat(k,k); 

  /*Pattern Search w/bounds on OptimizeMLEQuadratic*/
  PS->CleanSlate(theta->dim(), *theta,START_THETA_STEP, STOP_THETA_STEP,
                 OptimizeMLEQuadratic, this);

  PS->BeginSearch();
  
  PS->GetMinPoint(*theta);

  for(i=0;i<(*theta).dim();i++)
    if((*theta)[i]<TOLERANCE) (*theta)[i]=TOLERANCE; 
  
  /* Steps:
     1) Generate the correlation matrix
     2) Compute the trend
     3) Subtract off the trend from the function values.
     4) Krig them and compute sigma2, sum
  */  
  
  /*1) Generate the correlation matrix*/
  GenerateCorrelationMatrix(R,*theta);
  
  /*2) Compute the trend - code from A. Padula's krigifier*/
  for( i = 0; i < np; i++) {
    X[i][0] = 1;
    g = P+1;
    for( j = 1; j <=P; j++) {
      X[i][j] = (*points)[i][j-1];
      
      for(h = j-1;h < P; h++) {
	X[i][g] = (*points)[i][j-1]*(*points)[i][h];
	g++;
      }/*end for*/
    }/*end for*/
  }/*end for*/


  /*2a) Compute MLE_Matrix for MSE calculations*/

  /* Layout:
     [ 0 XT ]
     [ X R  ]
  */
  for(i=0;i<k+np;i++)
    for(j=0;j<k+np;j++) {
      if(i<k&&j<k) (*MLE_Matrix)[i][j]=0.0;
      else if (i<k&&j>=k) (*MLE_Matrix)[i][j]= X[j-k][i];//X transpose
      else if (i>=k&&j<k) (*MLE_Matrix)[i][j] = X[i-k][j];
      else /*i,j>=k*/
        (*MLE_Matrix)[i][j]=R[i-k][j-k];
    }/*end for*/

  /*2b) Invert the correlation, MSE matrices*/
  success=pseudoinvert(R,&sum);
  success=pseudoinvert(*MLE_Matrix,NULL);
  if(!success) return (success); /*break out!*/  



  /*2c) Finish Trend Calculations*/  
  (*beta) = transpose(X) * (R * (*values));
  tempmat = (transpose(X) * R) *  X;
  success=pseudoinvert(tempmat,NULL);
  (*beta) = tempmat * (*beta); //This is actually what I want!
  
  /*3) Subtract off the trend from the function values.*/
  for(z = 0; z < np; z++) {
    tempv[0] = 1;
    g = P+1;
    for( j = 1; j <=P; j++) {
      tempv[j] = (*points)[z][j-1];
      
      for(h = j-1;h < P; h++) {
        tempv[g] = (*points)[z][j-1]*(*points)[z][h];
        g++;
      }/*end for*/
    }/*end for*/
    (*v)[z] = (*values)[z] - tempv * (*beta);
  }/*end for*/

  
  /*4) Krig them and compute sigma2, sum*/
  (*sigma2)=((*v)*R*(*v))/np;
  (*v)=(*v)*R;
  return (success);
  
}/*end EstimateQuadraticTrend*/

    

bool ParamEstimate::EstimateCustomTrend(pt_collect &pts, mp_vector &fvals,
                                        double delta, void
                                        (*fcn2)(const mp_vector &, mp_vector&), long
                                        dimf, mp_vector* beta, double* sigma2, mp_vector* theta,
                                        pt_collect* MLE_Matrix, mp_vector* v){
/****
INPUT: Collection of evaluated points, their values, and pointers for the
parameters to be returned.  There is also a function pointer for the
'a' function, and the dimension of the vector it returns.  theta should come
in with a 'guess' value for the optimizer to chew on.
OUTPUT: Returns true if the operation was successful.
EFFECT: This estimates the function.
****/
  
#if DEBUG >=2
  cout<<"PEECUT :Started\n";  
#endif
  bool success;
  curr_delta=delta;
  np=pts.num_rows();
  fcn=fcn2;
  points=&pts;
  values=&fvals;
  fdim=dimf;
  
  pt_collect R(np,np);
  pt_collect A(np,fdim);
  pt_collect AT(fdim,np); /*A transpose*/
  mp_vector ai(fdim);
  pt_collect temp(fdim,fdim);
  long i,j;/*counters*/
  
  /*Pattern Search w/bounds on MLECustom*/
  PS->CleanSlate(theta->dim(), *theta,START_THETA_STEP, STOP_THETA_STEP,
                   OptimizeMLECustom, this);

  PS->BeginSearch();
  PS->GetMinPoint(*theta);

  for(i=0;i<theta->dim();i++)
    if((*theta)[i]<TOLERANCE) (*theta)[i]=TOLERANCE; 
  
  /*Generate the Correlation Matrix*/
  GenerateCorrelationMatrix(R,(*theta));
  
  /*Generate A and AT*/
  for(i=0;i<np;i++) {
    (*v)=(*points).row(i);
    (*fcn)(*v,ai);
    for(j=0;j<fdim;j++) {
      A[i][j]=ai[j];
      AT[j][i]=ai[j];
    }/*end for*/
  }/*end for*/
  
  /*Generate MLE_Matrix*/

  for(i=0;i<fdim+np;i++)
    for(j=0;j<fdim+np;j++)
      if(i<fdim && j<fdim) (*MLE_Matrix)[i][j]=0;
      else if(i<fdim) (*MLE_Matrix)[i][j]=AT[i][j-fdim]; /*and j>fdim*/
      else if(j<fdim) (*MLE_Matrix)[i][j]=A[i-fdim][j]; /*and i>fdim*/
      else (*MLE_Matrix)[i][j]=R[i-fdim][j-fdim];
#if DEBUG >=3  
  cout<<"PEECUT :A, AT and MLE Matrix Constructed:"<<(*MLE_Matrix);
#endif
  success=pseudoinvert(R,NULL);
  if(!success) return (success); /*break out!*/
  success=pseudoinvert(*MLE_Matrix,NULL);
  if(!success) return (success); /*break out!*/
#if DEBUG >=3  
  cout<<"PEECUT :Pseudoinversions complete\n";
#endif
  /*Get the Params*/
  temp=AT*R*A;
  success=(int)pseudoinvert(temp,NULL);
  if(!success) return (success); /*break out!*/  
  (*beta)=temp*AT*R*fvals;
  (*v)=fvals- A*(*beta);
  (*sigma2)=((*v)*R*(*v))/np; /*the hacked multiply assumes that the first v is vT*/
  (*v)=(*v) * R;
#if DEBUG >=2
  cout<<"PEECUT :terminating success="<<success<<"\n";
#endif   
  return(success);
}/*end EstimateCustomTrend*/
  
  
void ParamEstimate::SetDefault(long dim){
/****
INPUT: Dimension of the space.
EFFECTS: This will use a happy default set of values to set up the parameter estimation.
Specifically, it chooses Gaussian Isotropic Correlation and uses a Compass Search.
****/ 
  P=dim;
  r=Cor_Gauss_Isotropic;
  np=0;
  fdim=1;
  mp_vector temp(fdim);temp=0;
  PS=new CompassSearch(fdim,temp);
  curr_delta=0;
  fcn=NULL;
}/*end SetDefault*/

void ParamEstimate::MLEConstant(long dimtheta, mp_vector &x, double &f, bool &success){
/****
INPUT: PatternSearch calling sequence. See PatternSearch.h.
EFFECT: Evaluates the MLE for a given theta.  Uses the constant trend.   
NOTES: 5/31/99 - tested and functions well.
****/
  long i;/*counter*/
#if DEBUG>=3
  cout<<"PEMC   :MLEConstant Called\n";
  cout<<"PEMC   :trial_theta=(";
  for(i=0;i<dimtheta;i++) cout<<x[i]<<" ";
  cout<<")\n";    
#endif  
  double beta, sigma2=GetThetaUpperBound(curr_delta), sum=0.0;

  for(i=0;i<dimtheta;i++)
    if((x[i]>sigma2)||(x[i]<TOLERANCE)){
      success=false;return;
    }/*end if*/

  mp_vector v(np);
  mp_vector theta=x;
  pt_collect R(np,np);
  mp_vector A(np);
  A=1;/*fills A with 1's*/

  GenerateCorrelationMatrix(R,theta);
  success=pseudoinvert(R,&sum);
#if DEBUG>=4
  cout<<"PEMC   :Pseudoinverse Complete\n";
#endif
  if(success) {
 

    v=A*R; /*this multiply assumes that A is actually AT */
    beta=(v * (*values)) / (v*A); /*a hack for the constant trend case*/
    v=((*values) - A*beta);
    sigma2=(v * R * v)/np; /*again, the hacked multiply assumes that the first v
                               is vT*/ 
#if DEBUG>=4
    cout<<"PEMC   :beta="<<beta<<endl;
    cout<<"PEMC   :sigma2="<<sigma2<<endl;
#endif
    if(sigma2<TOLERANCE) {
#if DEBUG>=2
      cout<<"PEMC   :ERROR(non-fatal) sigma2="<<sigma2<<" resetting to TOLERANCE..\n";
#endif
      sigma2=TOLERANCE;
      }/*a tempory hack. do I want to do this?*/    
    f=np*log(sigma2) +sum;/*returned value*/
#if DEBUG>=3
    cout<<"PEMC   :Call Complete\n";
#endif
  }/*end if*/  
}/*end MLEConstant*/


void ParamEstimate::MLEQuadratic(long dimtheta, mp_vector &x, double &f, bool &success){
/****
INPUT: PatternSearch calling scheme.  See PatternSearch.h
EFFECT: Evaluates the MLE for a given theta.  Uses the quadratic trend.  This
is to be optimized by the Pattern Search, via a little bit of indirection.
****/  
  /* Steps:
     1) Generate the inverse correlation matrix
     2) Compute the trend
     3) Subtract off the trend from the function values.
     4) Krig them and compute sigma2, sum
     5) Return f
  */  
  long i, j, g, h,z;/*counter*/
#if DEBUG>=3
  cout<<"PEMQ   :MLEQuadratic Called\n";
  cout<<"PEMQ   :trial_theta=(";
  for(i=0;i<dimtheta;i++) cout<<x[i]<<" ";
  cout<<")\n";    
#endif  

  double sigma2=GetThetaUpperBound(curr_delta), sum=0.0;

  for(i=0;i<dimtheta;i++)
    if((x[i]>sigma2)||(x[i]<TOLERANCE)){
      success=false;return;
    }/*end if*/

  mp_vector v(np);
  mp_vector theta=x;
  pt_collect R(np,np);
  long k = 1+P+(P*(P+1))/2;
  mp_vector beta(k), tempv(k);
  pt_collect X(np,k), tempmat(k,k);

  
  /*1) Generate the inverse correlation matrix*/
  GenerateCorrelationMatrix(R,theta);
  success=pseudoinvert(R,&sum);

  /*2) Compute the trend - code from A. Padula's krigifier*/
  for( i = 0; i < np; i++) {
    X[i][0] = 1;
    g = P+1;
    for( j = 1; j <=P; j++) {
      X[i][j] = (*points)[i][j-1];
      
      for(h = j-1;h < P; h++) {
	X[i][g] = (*points)[i][j-1]*(*points)[i][h];
	g++;
      }/*end for*/
    }/*end for*/
  }/*end for*/

  beta = transpose(X) * (R * (*values));
  tempmat = (transpose(X) * R) *  X;
  pseudoinvert(tempmat,NULL);
  beta = tempmat * beta; //This is actually what I want!
  
  /*3) Subtract off the trend from the function values.*/
  for(z = 0; z < np; z++) {
    tempv[0] = 1;
    g = P+1;
    for( j = 1; j <=P; j++) {
      tempv[j] = (*points)[z][j-1];
      
      for(h = j-1;h < P; h++) {
        tempv[g] = (*points)[z][j-1]*(*points)[z][h];
        g++;
      }/*end for*/
    }/*end for*/
    v[z] = (*values)[z] - tempv * beta;
  }/*end for*/

  
  /*4) Krig them and compute sigma2, sum*/
  sigma2=(v*R*v)/np;


  if(sigma2<TOLERANCE) {
#if DEBUG>=2
    cerr<<"PEMQ   :ERROR(non-fatal) sigma2="<<sigma2<<" resetting to TOLERANCE..\n";
#endif
    sigma2=TOLERANCE;
  }/*end if*/

  /*5) Return f*/

  f=np*log(sigma2)+sum;
  
}/*MLEQuadratic*/


void ParamEstimate::MLECustom(long dimtheta, mp_vector &x, double &f, bool &success){
/****
INPUT: PatternSearch calling sequence.  See PatternSearch.h
EFFECT: Evaluates the MLE for a given theta.  Uses the custom trend.   
****/
  long i;/*counter*/
#if DEBUG>=2
  cout<<"PEMCU  :MLECustom Called\n";
  cout<<"PEMCU  :trial_theta=(";
  for(i=0;i<dimtheta;i++) cout<<x[i]<<" ";
  cout<<")\n";    
#endif
  double sigma2=GetThetaUpperBound(curr_delta), sum=0.0;
  for(i=0;i<dimtheta;i++)
    if(x[i]>sigma2||x[i]<TOLERANCE){
      success=false;return;
    }/*end if*/
  
  pt_collect R(np,np);
  mp_vector v(np);
  mp_vector theta=x;
  mp_vector beta(fdim);
  pt_collect temp(fdim,fdim);
  mp_vector ai(fdim);
  pt_collect A(np,fdim);
  pt_collect AT(fdim,np); /*A transpose*/


  GenerateCorrelationMatrix(R,theta);
  success=pseudoinvert(R,&sum);
#if DEBUG>=3
  cout<<"PEMCU  :Pseudoinverse Complete\n";
#endif

  if(success) {

    /*Generate A and AT*/
    for(i=0;i<np;i++) {
      v=(*points).row(i);
      (*fcn)(v,ai);
      for(long j=0;j<fdim;j++) {
        A[i][j]=ai[j];
        AT[j][i]=ai[j];
      }/*end for*/
    }/*end for*/
#if DEBUG>=3
  cout<<"PEMCU  :A,AT computed\n";
#endif    
    temp=AT*R*A;
    success=pseudoinvert(temp,NULL);
    if(success) {
      beta=temp*AT*R*(*values);
      v=(*values)- A*beta;
      sigma2=(v*R*v)/np; /*multiply assumes that the first v is vT*/
#if DEBUG>=3
    cout<<"PEMCU  :beta="<<beta<<endl;
    cout<<"PEMCU  :sigma2="<<sigma2<<endl;
#endif     
    if(sigma2<TOLERANCE) {
#if DEBUG>=2
      cerr<<"PEMCU  :ERROR(non-fatal) sigma2="<<sigma2<<" resetting to TOLERANCE..\n";
#endif
      sigma2=TOLERANCE;
    }/*end if*/

    f=np*log(sigma2)+sum;
    }/*end if*/      
  }/*end if*/
#if DEBUG>=2
    cout<<"PEMCU   :Call Complete\n";
#endif 
}/*end MLECustom*/

inline void ParamEstimate::GenerateCorrelationMatrix(pt_collect &R, const mp_vector &theta) {
/****
INPUT: A np x np pt_collect, and the correlation family parameter vector, theta.
EFFECT: Replaces R with the appropriate correlation matrix.
NOTES: 5/31/99 TESTED OK
****/
  mp_vector temp(P);
  for(long i=0;i<np;i++) {
    R[i][i]=1.0;
    temp=(*points).row(i);
    for(long j=i+1;j<np;j++) {
      R[i][j]=(*r)(P,temp, (*points).row(j), theta);
      R[j][i]=R[i][j];      
    }/*end for*/
  }/*end for*/
}/*end GenerateCorrelationMatrix*/


bool ParamEstimate::pseudoinvert(pt_collect &M, double * log_det) {
/****
INPUT: A matrix M, and a pointer to a mp_vector.
OUTPUT: true if svd can be computed, else false.
EFFECT: Replaces M with the inverse (if possible) or pseudoinverse of M.  If
d!=NULL, an mp_vector with the one over the singular values, with those singular
values below TOLERENCE instead to zero, is returned.  The matrix has to be
*SYMMETRIC AND SQUARE* for this to work.
****/
  if(M.num_rows()!=M.num_cols()) return false;
  bool ret_val;
  long md=M.num_rows();
#if DEBUG >=3
  cout<<"PEPI   :Pseudoinverse Commencing md="<<md<<"\n";      
#endif
  /*We're only going to attempt to do the Cholesky on non MLE_Matrix matricies.
    This is because the Cholesky will almost always fail on them (first column starts
    with a 0).
  */
  if(md==np && ComputeCholeskyInverse(&M,log_det)) ret_val=true;
  else {
    /*Try Iterative Cholesky*/
    pt_collect Mb(M);
    double TOADD = md*(md*MACH_EPS)/(1-md*MACH_EPS);
    for( long i = 0; i < md; i++)
      M[i][i] += TOADD;
      
    /*Calculate the Cholesky Decomposition*/
    ret_val = ComputeCholeskyInverse(&M, log_det);
    if(ret_val==false) {
      M=Mb;
      ret_val=ComputeSVDInverse(&M,log_det);
    }/*end if*/
  }/*end else*/

  return ret_val;
}/*end pseudoinvert*/

bool ParamEstimate::ComputeCholeskyInverse(pt_collect *A,double* log_det) {
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
      
  if(A->num_rows()!=A->num_cols()) return false;
  long i,j,md=A->num_rows(), idx;
  double* MA=(double*)malloc(sizeof(double)*md*md);
  char ch='U';
  if(log_det!=NULL) (*log_det)=0.0;
 
  for(i=0;i<md;i++)
    for(j=0;j<md;j++)
      MA[md*j+i]=(*A)[i][j]; /*This should flip the Matrix into column-major order
                              ala Fortran.*/
      
  dpotrf_(&ch, &md, MA, &md, &idx); 
  if(idx==0) {   
    /*Compute the log determinent*/
    if(log_det!=NULL) {
      for(j=0;j<md;j++) {
        if(MA[md*j+j] < TOLERANCE) (*log_det)+=log(TOLERANCE);
        else (*log_det)+=log(MA[md*j+j]);
      }/*end for*/      
      (*log_det)*=2.0;/*Double to account for Cholesky*/

    }/*end if*/
    dpotri_(&ch, &md, MA, &md, &idx);     
  }/*end if*/
    
  if(idx<0) {free(MA);return false;}
  else if (idx>0) {
#if DEBUG>=2
    if(log_det ==NULL) {
    cout<<"PECC   : WARNING - A is not positive definite at row "
        <<idx<<"... switching to SVD"<<endl;
    }
#endif
    free(MA);
    return false;
  }  
  /*Copy data to outgoing matrix*/
  for(j=0;j<md;j++) {/*cols*/
    idx=md*j;/*Partial computation of index... for speed*/

    for(i=0;i<=j;i++) {/*rows*/
      (*A)[i][j]=MA[idx+i];
      (*A)[j][i]=MA[idx+i];
    }/*end for*/
  }/*end for*/ 

  /*Free memory!*/
  free(MA);
  return true;  
}/*end ComputeCholesky */

bool ParamEstimate::ComputeSVDInverse(pt_collect* A, double * log_det) {
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

  if(A->num_rows()!=A->num_cols()) return false;
  
  long md=A->num_rows();
  long idx, LWORK=md*md+5*md;
  double* MA=(double*)malloc(sizeof(double)*md*md);
  double* MU=(double*)malloc(sizeof(double)*md*md);  
  double* MVT=(double*)malloc(sizeof(double)*md*md);
  double* MS=(double*)malloc(sizeof(double)*md);
  double* WORK=(double*)malloc(sizeof(double)*LWORK);
  pt_collect U(md,md);
  pt_collect VT(md,md);
  pt_collect D(md,md);

  char ch='A';
  long i,j;/*counters*/

  for(i=0;i<md;i++)
    for(j=0;j<md;j++)
      MA[md*j+i]=(*A)[i][j]; /*This should flip the Matrix into column-major order
                               ala Fortran.*/
      
  dgesvd_(&ch,&ch,&md,&md,MA,&md,MS,MU,&md,MVT,&md,WORK,&LWORK,&idx);

  if(idx<0) {free(MA);free(MU);free(MVT);free(MS);free(WORK);return false;}
#if DEBUG>0
  else if (idx>0) cerr<<"PECSVD : WARNING - SVD has "<<idx<<" unconverged superdiagonal values\n";
#endif
  
  /*Copy data to all the outgoing matrices*/
  for(j=0;j<md;j++) {/*cols*/
    idx=md*j;/*Partial computation of index... for speed*/
    for(i=0;i<md;i++) {/*rows*/
      U[i][j]=MU[idx+i];
      VT[i][j]=MVT[idx+i];
      if(i==j) D[i][j]=MS[i];
      else D[i][j]=0.0;
    }/*end for*/
  }/*end for*/

  if(log_det!=NULL) {
    (*log_det)=0.0;
    /*Compute log determinant of matrix via singular values = eigenvalues*/
    for(i=0;i<md; i++)
      if (fabs(D[i][i])<TOLERANCE) (*log_det)+=log(TOLERANCE);
      else (*log_det)+=log(D[i][i]);
  }/*end if*/
  
  double t_tolerance= MACH_EPS * A->num_rows() * D(0,0);
    
  for(i=0;i<md;i++) {
    if (fabs(D[i][i])<t_tolerance) D[i][i]=0;
    else D[i][i]=1/D[i][i];
  }/*end for*/
      
  (*A)=transpose(VT) * D * transpose(U);
  
#if DEBUG>=4
  cout<<"PECSVD :A="<<(A);
  cout<<"PECSVD :U="<<(*U);
  cout<<"PECSVD :D="<<(*D);
  cout<<"PECSVD :VT="<<(*VT);
  cout<<"PECSVD :U*D*VT="<<((*U)*(*D)*(*VT));
#endif  
  /*Free memory!*/
  free(MA);free(MU);free(MVT);free(MS);free(WORK);
  return true;  
}/*end ComputeSVD*/
                  

/**************************CORRELATION FUNCTIONS*****************************************/

double ParamEstimate::Cor_Exp_Isotropic(long dim, const mp_vector &x1, const mp_vector
                         &x2, const mp_vector &theta) {
  /*Exponential Isotropic Correlation Family
    r_theta(s,t) = exp(-theta ||s-t||)
  */
  /*This assumes that the theta of value is theta[0]*/
  mp_vector temp(dim);
  temp=x2-x1;
  return exp(-theta[0] * temp.l2norm());  
}/*end Cor_Exp_Isotropic*/

double ParamEstimate::Cor_Gauss_Isotropic(long dim, const mp_vector &x1, const mp_vector
                           &x2, const mp_vector &theta) {
  /*Gaussian Isotropic Correlation Family
    r_theta(s,t) = exp(-theta ||s-t||)^2
  */ 
  /*This assumes that the theta of value is theta[0]*/
  mp_vector temp(dim);
  temp=x2-x1;
  return exp(-theta[0] * temp.l2norm_sqr());    
}/*end Cor_Gauss_Isotropic*/

double ParamEstimate::Cor_Exp_Product(long dim, const mp_vector &x1, const mp_vector &x2,
                       const mp_vector &theta) {
  /*Exponential Product Correlation Family
    r_theta(s,t) = \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|)
  */
  double sum=1.0;
  for(long i=0;i<dim;i++) sum*=exp(-theta[i]*fabs(x2[i]-x1[i]));
  return sum;
}/*end Cor_Exp_Product*/

double ParamEstimate::Cor_Gauss_Product(long dim,const mp_vector &x1, const mp_vector
                         &x2, const mp_vector &theta) {
  /*Gaussian Product Correlation Family
    r_theta(s,t) = \prod_{j=1}^{P} exp(-theta_j * |s_j - t_j|^2)
  */
  double sum=1.0;
  for(long i=0;i<dim;i++) sum*=exp(-theta[i]*sqr(x2[i]-x1[i]));
  return sum;
}/*end Cor_Gauss_Product*/

double ParamEstimate::Cor_Cubic_Isotropic(long dim, const mp_vector &x1,
                                const mp_vector &x2,const mp_vector &theta) {
  /*Cubic Product Correlation Family
    r_theta(s,t) = { 1-1.5*(||s-t||)/theta + 0.5 ((||s-t||)/theta)^3, for ||s-t|| < theta }
                   { 0 otherwise }
  */
 /*This assumes that the theta of value is theta[0]*/  
  mp_vector temp(dim);
  temp=x2-x1;
  double norm=temp.l2norm();
  double norm3theta= (norm*norm*norm)/(theta[0]*theta[0]*theta[0]);
  if(norm < theta[0]) 
    return 1.0 - 1.5*(norm / theta[0]) + 0.5 * norm3theta;
  else
    return 0.0;
  
}/*end Cor_Cubic_Isotropic*/

double ParamEstimate::EvaluateMyCorrelation(long dim, const mp_vector &x1, const mp_vector &x2,
                             const mp_vector &theta) {
  /****
  INPUT: Dimesion of space, two points, and a 'vector' theta.  For the Isotropic functions, this
  theta is a one dimensional vector, for the Product functions, this is a
  P-dimensional vector.
  OUTPUT: Correlation between the two points.
  ****/     
  return (*r)(dim,x1,x2,theta);
}/*end EvaluateMyCorrelation*/

double GetMachineEpsilon() {
  /*From Heath's Scientific Computing*/

  return fabs( 3.0 * (4.0/3.0 - 1.0) - 1.0 );
}/*end GetMachineEpsilon*/
