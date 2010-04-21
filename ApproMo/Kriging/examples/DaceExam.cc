/* DaceExam.cc
   This is the demo main program for the Krigifier and the Kriging approximation classes
   This program will approximate the data of the matlab dace algorithms to verify the functions of the C++ package. The data are loaded from input1.inp, input2.inp and target.inp. The results are output to the file "results.out" */

#include "KrigApprox.h" 
#include "krig.h"
#include "krigify.h"
#include "rngs.h"


using namespace std;

int main( void ) {
  
  FILE* fstuetz;
  char fstuetzname[500];
  sprintf ( fstuetzname, "stuetzstellen_red.dat" );
  fstuetz = fopen ( fstuetzname, "rt" );
  Array < double > Dat ( 87*4, 3 );
  int counterr = 0;
  fscanf ( fstuetz, "x\n\n" );
  for ( int i = 0; i < 87*4; i++ ) {
    int ii1;
    fscanf ( fstuetz, "%i", &ii1 );
    Dat(counterr,0) = (double ) ii1;
    cout << Dat(counterr,0) << endl;
    counterr++ ;
  }
  fscanf ( fstuetz, "\n\ny\n\n" );
  counterr = 0;
  for ( int i = 0; i < 87*4; i++ ) {
    int ii1;
    fscanf ( fstuetz, "%i", &ii1 );
    Dat(counterr,1) = (double ) ii1;
    cout << Dat(counterr,1) << endl;
    counterr++ ;
  }
  fscanf ( fstuetz, "\n\nz\n\n" );
  counterr = 0;
  for ( int i = 0; i < 87*4; i++ ) {
    double ii1;
    fscanf ( fstuetz, "%lf", &ii1 );
    Dat(counterr,2) = ii1;
    cout << Dat(counterr,2) << endl;
    counterr++ ;
  }
  for ( int i = 0; i < 87*4; i++ ) {
    cout << Dat(i,0) << " " << Dat(i,1) << " " << Dat(i,2) << " " << endl;
  }
  fclose ( fstuetz );
  FILE* fout11;
  char fout11name[500];
  sprintf ( fout11name, "current_mesh.dat" );
  fout11 = fopen ( fout11name, "wt" );
  for ( int i = 0; i < 87*4; i++ ) {
    fprintf ( fout11, "%lf %lf %lf\n",  Dat(i,0),  Dat(i,1),  Dat(i,2) );
  }
  fclose ( fout11 );
 
  Database myData1;
  KrigApprox myApprox1;
  
  int m1 = 87*4;       // The number of initial design sites used to create the approximation
  int p1 = 2;                       // The dimension of the space ( Given in "krig.dat" )

  Matrix<double> X1(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y1(m1);         // Vector to store the function values at the initial design sites 
  X1.newsize(m1,p1);

  Array < double > InputData1(m1,p1);
  Array < double > TargetData1(m1,1);

  for ( long a = 0; a < m1; a++) {
    InputData1(a,0) = Dat(a,0);
    if ( Dat (a,0) == 20 ) InputData1(a,0) = 0;
    if ( Dat (a,0) == 52 ) InputData1(a,0) = 0.2;
    if ( Dat (a,0) == 77 ) InputData1(a,0) = 0.4;
    if ( Dat (a,0) == 109 ) InputData1(a,0) = 0.6;
    InputData1(a,1) = Dat(a,1)/100;
    TargetData1(a,0) = Dat(a,2)*100;
  }

  fout11 = fopen ( fout11name, "wt" );
  for ( int i = 0; i < 87*4; i++ ) {
    fprintf ( fout11, "%lf %lf %lf\n",  InputData1(i,0),   InputData1(i,1),   TargetData1(i,0) );
  }
  fclose ( fout11 );
 

  for ( long a = 0; a < m1; a++) {
    cout << InputData1(a,0) << " " << InputData1(a,1) << " " << TargetData1(a,0) << endl;
  }
 
  cout << "creating approximation..." << endl;
  myData1.AddTrainingData(InputData1,TargetData1);
  myApprox1.Train(myData1);

  FILE* fout12;
  char Outputfilename2[100];
  strcpy(Outputfilename2,((string)"initial_mesh" + ".dat").c_str());
  double result2;
  fout12 = fopen(Outputfilename2,"wt");

  // for evaluating an array of InputDatas use:
  // myApprox.Evaluate(InputData,TargetData);
  Array < double > Inpu ( 2, 2 );
  Inpu(0,0) = 21;
  Inpu(0,1) = 286;
  Inpu(1,0) = 19;
  Inpu(1,1) = 286; //0.0171
  myApprox1.Evaluate(Inpu,TargetData1);
  cout << Inpu(0,0) << " " << Inpu(0,1) << "  " << TargetData1(0,0) << endl;
  cout << Inpu(1,0) << " " << Inpu(1,1) << "  " << TargetData1(1,0) << endl;
  
  /*
  for ( long a = 0; a < m1; a++) {
    cout << InputData1(a,0) << " " << InputData1(a,1) << " " << TargetData1(a,0) << endl;
  }
  */

  Array < double > ResultDat ( 400*400,1);
  Array < double > TestDat ( 400*400,2);
  int ccc = 0;
  for ( int j = 0; j < 400; j++ ) {
    for ( int i = 0; i < 400; i++ ) {
      TestDat (ccc, 0 ) = (double)i/100;
      TestDat ( ccc,1 ) = (double)j/100;
      ccc++;
    }
  }
  myApprox1.Evaluate(TestDat,ResultDat);

  double jj22 = 0;
  double kk22 = 0;
  double resz;
  for ( int ii = 0; ii < (int)ResultDat.dim(0); ii++) {
    jj22 =  TestDat(ii,0);
    kk22 =  TestDat(ii,1);
    resz = ResultDat(ii,0);
    fprintf(fout12, "%f %f %f \n" , (double)jj22, (double)kk22, resz);
    fflush(fout12);
  }
  fclose(fout12);

  exit(0);
  
  /*
  int ParametersToKeepConstant = m1;
  Array < double > Inputs ( 1, ParametersToKeepConstant+2 );
  Vector < double > InputVector( 2 );
  Vector < double > DataToEvaluate( 2 + ParametersToKeepConstant );
  int resolution = 100;
  int upperBorder = 100;
  int lowerBorder = 0;
  int diff = upperBorder - lowerBorder;
  TargetData1.resize( (unsigned int)sqr(resolution+1), (unsigned int)3);
  int outputEntry = 0;
  for ( int i = 0; i <= resolution; i++ ) {
    for ( int j = 0; j <= resolution; j++ ) {
      InputVector[(long)0] = (double)i*diff/resolution+lowerBorder;
      InputVector[(long)1] = (double)j*diff/resolution+lowerBorder;
      OutputData(outputEntry,0) = (double)i*diff/resolution+lowerBorder;
      OutputData(outputEntry,1) = (double)j*diff/resolution+lowerBorder;

      if ( Inputs.dim(1) > 2 ) {
	int counter = 0;
	for ( int k = 0; k < (int)Inputs.dim(1); k++ ) {
	  if ( k == ind1 ) {
	    DataToEvaluate [ k ] = InputVector [ 0 ];
	  }
	  else {
	    if ( k == ind2 ) {
	      DataToEvaluate [ k ] = InputVector [ 1 ];
	    }
	    else {
	      DataToEvaluate [ k ] = ParametersToKeepConstant ( counter );
	      counter++;
	    }
	  }
	}
      }
      else {
	DataToEvaluate [ 0 ] = InputVector [ 0 ];
	DataToEvaluate [ 1 ] = InputVector [ 1 ];
      }

      OutputData(outputEntry,2) = KrigApp->evalf(DataToEvaluate);
      outputEntry++;
    }
  }
  */

  
  int jj2 = 0;
  int kk2 = 0;
  // for scanning the approximation of an area use
  myApprox1.ScanSquare3D(0,100,TargetData1);

  for ( int ii = 0; ii < (int)TargetData1.dim(0); ii++) {
    jj2 = (int) TargetData1(ii,0);
    kk2 = (int) TargetData1(ii,1);
    result2 = TargetData1(ii,2);
    fprintf(fout12, "%f %f %f \n" , (double)jj2, (double)kk2, result2);
    fflush(fout12);
  }
  fclose(fout12);
	

  exit(0);
  

  Database myData;
  KrigApprox myApprox;
  
  int m = 75;       // The number of initial design sites used to create the approximation
  int p = 2;                       // The dimension of the space ( Given in "krig.dat" )

  Matrix<double> X(1,1);       // The matrix to store the initial design sites
 
  Vector<double> y(m);         // Vector to store the function values at the initial design sites 
  X.newsize(m,p);

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

  cout << "creating approximation..." << endl;
  myData.AddTrainingData(InputData,TargetData);
  myApprox.Train(myData);

  FILE* fout1;
  char Outputfilename[100];
  strcpy(Outputfilename,((string)"initial_mesh" + ".dat").c_str());
  double result;
  fout1 = fopen(Outputfilename,"wt");

  // for evaluating an array of InputDatas use:
  // myApprox.Evaluate(InputData,TargetData);

  int jj = 0;
  int kk = 0;
  // for scanning the approximation of an area use
  myApprox.ScanSquare3D(0,100,TargetData);

  for ( int ii = 0; ii < (int)TargetData.dim(0); ii++) {
    jj = (int) TargetData(ii,0);
    kk = (int) TargetData(ii,1);
    result = TargetData(ii,2);
    fprintf(fout1, "%f %f %f \n" , (double)jj, (double)kk, result);
    fflush(fout1);
  }
  fclose(fout1);
	
  //cout << "Error: " << myApprox.Error_MSE(myData) << endl;

  return 0;
} 
