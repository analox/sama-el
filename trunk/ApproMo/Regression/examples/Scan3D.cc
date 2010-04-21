#include "Array.h"
#include "RegrApp.h"


using namespace std;

int main( int argc, char **argv )
{

  Database myData;
  RegrApp myNet(10,false);
 
  int m = 75;    
  int p = 2;           

  Array < double > InputData(m,p);
  Array < double > TargetData(m,1);  

  FILE* fin; 
  char Inputfilename[100];
  strcpy(Inputfilename,((string)"../data/input1" + ".inp").c_str());
  fin = fopen(Inputfilename,"rt");
  double tmpx, tmpy, tmpz;
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpx);
    //    X[a][0L] = tmpx;
    InputData(a,0) = tmpx;
  }
  fclose(fin);

  strcpy(Inputfilename,((string)"../data/input2" + ".inp").c_str());
  fin = fopen(Inputfilename,"rt");
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpy);
    //    X[a][1L] = tmpy;
    InputData(a,1) = tmpy;
  }
  fclose(fin);
  strcpy(Inputfilename,((string)"../data/target" + ".inp").c_str());
  fin = fopen(Inputfilename,"rt");
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpz);
    //    y[a] = tmpz;
    TargetData(a,0) = tmpz;
    cout << InputData(a,0) << " " << InputData(a,1) << " " << TargetData(a,0) << endl;
  }
  fclose(fin);

  myData.AddTrainingData(InputData,TargetData);
  myNet.Train(myData);

  FILE* fout1;
  char Outputfilename[100];
  strcpy(Outputfilename,((string)"results" + ".out").c_str());
  double result;
  fout1 = fopen(Outputfilename,"wt");

  // for evaluating an array of InputDatas use:
  // myApprox.Evaluate(InputData,TargetData);

  int jj = 0;
  int kk = 0;
  // for scanning the approximation of an area use
  myNet.ScanSquare3D(0,100,TargetData);

  for ( int ii = 0; ii < TargetData.dim(0); ii++) {
    jj = (int) TargetData(ii,0);
    kk = (int) TargetData(ii,1);
    result = TargetData(ii,2);
    fprintf(fout1, "%f %f %f \n" , (double)jj, (double)kk, result);
    fflush(fout1);
  }
  fclose(fout1);

  cout << "Error: " << myNet.MSE(myData) << endl;

}
