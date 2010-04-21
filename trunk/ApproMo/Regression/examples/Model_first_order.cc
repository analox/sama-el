#include "Model_first_order.h"
#include "iostream.h"

int Model_first_order(Array<double> InputData,
		      Array<double> TargetData,
		      int start,
		      int ende,
		      Array<double> &coefficient) {
  
  //  int N = InputData.dim(0);

  int D = InputData.dim(1);
  int p = D+1;
  
  int N = ende-start;

  /* 
     printf("Number of data points: %d\n  Dimension %d\n",N,p);
     
     FILE  *fout1,*fout2;
     fout2 = fopen("t2.dat","wt");
     fout1 = fopen("t3.dat","wt");
  */

  double **a, *b, d;
  int *indx;
  double **X;
  double *Y;
  
  // index for resorting
  indx = ivector(1,p);
  
  // input data 
  a = dmatrix(1,N,1,p);
  // a'*a
  X = dmatrix(1,p,1,p);
  
  // results
  b = dvector(1,N);
  
  // X'*b
  Y = dvector(1,p);
  
  // initialise a
  for (int i = 1; i <= N; i++) { 
    a[i][1] = 1.0;
    for (int pp = 1; pp <= D; pp++) {
      a[i][pp+1] = InputData(start + i-1, pp-1);
    }
    b[i] = TargetData(start + i - 1,0);
  }

  // calculate X = a'*a
  for (int i = 1; i <= p; i++) {
    for (int j = 1; j <= p; j++) {
      X[i][j] = 0.0;
    }
  }
  
  for (int ii = 1; ii <= p; ii++) {
    for (int jj = 1; jj <= p; jj++) {
      for (int i = 1; i <= N; i++) {
	X[ii][jj] += a[i][ii] * a[i][jj];
      }
    }
  }

  for (int ii = 1; ii <= p; ii++) {
    for (int jj = 1; jj <= p; jj++) {
      cout << ii << jj << ": " <<X[ii][jj] << " ";
    }
    cout << endl;
  }
   
  // calculating Y = X'*b
  for (int i = 1; i <= p; i++) {
    Y[i] = 0.0;
  }
  
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= p; j++) {
      Y[j] += a[i][j] * b[i];
    }
  }
  
  for (int j = 1; j <= p; j++) {
    cout << j << ": " << Y[j] << endl;
  }

  /*
  printf("the equation system\n");
  for (int ii = 1; ii <= p; ii++) {
    for (int jj = 1; jj <= p; jj++) {
      printf("%lf  " ,X[ii][jj]);
    }
    printf("   %lf\n",Y[ii]);
  }
  */

  ludcmp(X,p,indx,&d);
  lubksb(X,p,indx,Y);
   
  //printf("coefficients \n");
  
  coefficient.resize((unsigned int)p,(unsigned int)1);
  for (int i = 1; i <= p; i++) {
    coefficient(i-1,0) = Y[i];
    cout << Y[i] << endl;
    //printf("%lf   ",Y[i]);
  }

  /*  for (int i = 1; i <= p; i++) {
 
    printf("%lf   ",Y[i]);
  }
  printf("\n");
  */

  /*
  for (int ii = 1; ii <= N; ii++){ 
     for (int jj = 1; jj <= N; jj++){ 
       int i = (ii-1) *N + jj;
       double x = 0.+((double)ii/N);
       double y = 0.+((double)jj/N);
       fprintf(fout2,"%lf %lf %lf\n",x,y,Y[1]+Y[2]*x+Y[3]*y+ Y[4]*x*x + Y[5] *y*y);
     }
  }

  fclose(fout1);
  fclose(fout2);
  */
  free_ivector(indx,1,p);
  free_dmatrix(a,1,N,1,p);
  free_dmatrix(X,1,p,1,p);
  free_dvector(b,1,N);
  free_dvector(Y,1,p);

  return 0;
}

// ************************************************************

