#include "Array/ArrayOp.h"
#include "Database.h"
#include "EvoNet.h"


using namespace std;

int main( int argc, char **argv )
{

  Database myData;

  int in = 2;
  int out = 1;
  int dimension = 1;

  Array < int > cmat((unsigned int)(3*dimension),(unsigned int)(3*dimension+1));
  for ( unsigned i = 0; i < 3*dimension; i++) {
    for ( int j = 0; j < 3*dimension+1; j++) {
      cmat(i,j) = 1;
    }
  }
  
  Array < double > wmat((unsigned int)(3*dimension),(unsigned int)(3*dimension+1));
  for ( unsigned i = 0; i < 3*dimension; i++) {
    for ( int j = 0; j < 3*dimension+1; j++) {
      wmat(i,j) = Rng::uni(-0.01,0.01);
    }
  }

  EvoNet myNet(2,1, cmat, wmat);
 
  int m = 75;    
  int p = 2;           

  Array < double > InputData(m,p);
  Array < double > TargetData(m,1);  

  FILE* fin; 
  char Inputfilename[100];
  strcpy(Inputfilename,((string)"input1" + ".inp").c_str());
  fin = fopen(Inputfilename,"rt");
  double tmpx, tmpy, tmpz;
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpx);
    //    X[a][0L] = tmpx;
    InputData(a,0) = tmpx;
  }
  fclose(fin);

  strcpy(Inputfilename,((string)"input2" + ".inp").c_str());
  fin = fopen(Inputfilename,"rt");
  for ( long a = 0; a < 75; a++) {
    fscanf(fin,"%lf", &tmpy);
    //    X[a][1L] = tmpy;
    InputData(a,1) = tmpy;
  }
  fclose(fin);
  strcpy(Inputfilename,((string)"target" + ".inp").c_str());
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
  myNet.StructureOptimization(myData);
  myNet.Train(myData);
  myNet.StructureOptimization(myData);
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

  cout << "Saving the net..." << endl;
  myNet.Save( "Net.sav");

}
