#include "RegrApp.h"
#include <iostream.h>
#include "Model_first_order.h"
#include "Model_second_order.h"
#include <string.h>

using namespace std;

int main ( int argc, char **argv ) {

  int a = 7; int b = 2; int c = 1;
  Array < double > InputData(a,b);
  Array < double > OutputData(a,c);
  Array < double > TargetData(a,c);

  Array < double > Input1(1,b);
  Array < double > Output1(1,1);
  for (int i = 0; i < b; i++) {
    Input1(0,i) = i*i*5*(i+1)+15-3*i;
  }

  cout << "InputData:" << endl;
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      InputData(i,j) = i*i+j*j;
      InputData(0,0) = 25.3;
      InputData(1,1) = 17.5;
      InputData(2,1) = 2.0;
      InputData(3,0) = 38.2;
      InputData(4,0) = 39;
      InputData(4,1) = 39;
      InputData(5,0) = 0;
      InputData(5,1) = 39;
      InputData(5,0) = 20;
      InputData(5,1) = 39;
      cout << "Point: " << i << " " << j << ": " << InputData(i,j) << endl;
    }
  }

  cout << "TargetData" << endl;
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < c; j++) {
      TargetData(i,j) = i*5*(j+1)+j/3;
       TargetData(0,0) = 22.1;
       TargetData(1,0) = -12.3;
       TargetData(2,0) = 32;
       TargetData(3,0) = -200;
       TargetData(4,0) = 20;
       TargetData(5,0) = 20;
       TargetData(5,0) = -20;
      cout << "Point: " << i << " " << j << ": " << TargetData(i,j) << endl;
    }
  }

  cout << "Model first order!" << endl;

  int start = 0;
  int end = a;
  Array<double> coefficient(1,1);
  Array < double > InputDatainvers(a,b);
  Array < double > TargetDatainvers(a,c);
  for ( int i = 0; i < InputData.dim(0); i++) {
    for (int j = 0; j < InputData.dim(1);j++) {
      InputDatainvers( (InputData.dim(0)-i-1), j) = InputData (i,j);
    }
    TargetDatainvers( (InputData.dim(0)-i-1), 0) = TargetData (i,0);
  }

  /*
    Model_first_order(InputDatainvers,
		    TargetDatainvers,
		    start,
		    end,
		    coefficient);
  
  cout << "Model end!" << endl;

  cout << "Model second order!" << endl;
   
  Model_second_order(InputData,
		     TargetData,
		     start,
		     end,
		     coefficient); */
  Model_second_order_interaction(InputDatainvers,
				 TargetDatainvers,
				 start,
				 end,
				 coefficient);

  double y1;

  double x1 = Input1(0,0);
  double x2 = Input1(0,1);
  y1 = 1.0*coefficient(0,0) + x1 * coefficient(1,0) + x2 * coefficient(2,0) + x1*x1*coefficient(3,0) + x2*x2 *coefficient(4,0) + x1*x2*coefficient(5,0);
  cout << "Ergebnis: " << y1 << endl;
  //  x1 = 25.3;
  //  x2 = 1;

  FILE* fo;
  char Filename[100];
  strcpy(Filename,((string)"ModelFO" + ".out").c_str());
  fo = fopen(Filename,"wt");

  for (int i = 0; i < InputData.dim(0); i++ ) {
    fprintf(fo,"%f %f %f \n", InputData(i,0), InputData(i,1), TargetData(i,0) );
  }
   
  for ( int i = 0; i <= 40; i++ ) {
    for ( int j = 0; j <= 40; j++ ) {
      x1 = i;
      x2 = j;
      y1 = 1.0*coefficient(0,0) + x1 * coefficient(1,0) + x2 * coefficient(2,0) + x1*x1*coefficient(3,0) + x2*x2 *coefficient(4,0) + x1*x2*coefficient(5,0);
      // cout << "Ergebnis: " << y1 << endl;
      fprintf(fo,"%f %f %f \n", x1,x2,y1);
      fflush(fo);
    }
  }
  fclose(fo);
  
  cout << "Model end!" << endl;

  Database myData;

  myData.AddTrainingData(InputData,TargetData);
  
  cout << "Konstruktor!" << endl;

    RegrApp myRegression(1,false);

  cout << "Regression 1 start!" << endl;

  myRegression.Train(myData);

  /*
  RegrApp myRegression2(2,false);

  myRegression2.Train(myData);*/

  /*
  RegrApp myRegression3(2,true);

  myRegression3.Train(myData);

  myRegression3.Evaluate(Input1,Output1);

  cout << "Error: " << myRegression3.Error_MSE(myData) << endl;

  for ( int i = 0; i < Output1.dim(0); i++) {
    cout << "Result " << i << ": " << Output1(i,0)<< endl;
  }

  FILE* fo2;
  char Filename1[100];
  strcpy(Filename1,((string)"ModelRegrApp" + ".out").c_str());
  fo2 = fopen(Filename1,"wt");
   
  for (int i = 0; i < InputData.dim(0); i++ ) {
    fprintf(fo2,"%f %f %f \n", InputData(i,0), InputData(i,1), TargetData(i,0) );
  }


  int counter =0;
  for ( int i = 0; i <= 40; i++ ) {
    for ( int j = 0; j <= 40; j++ ) {
      Input1(0,0) = (double) i;
      Input1(0,1) = (double) j;
      myRegression3.Evaluate(Input1,Output1);
      fprintf(fo2,"%f %f %f \n", (double)i,(double)j,Output1(0,0));
      counter++;
      fflush(fo2);
    }
  }
  fclose(fo2);

  Array < double > Scan(1,1);
  int z = 0;
  int zz = 40;
  myRegression3.ScanSquare3D( z, zz, Scan);
  //  myRegression3.ScanSquare3D( z, zz, Scan);
  */
  cout << "Everything is done!" << endl;

  // Regression 1. order with 3 parameters and less data
  cout << endl << endl;
  Array < double > InputData1(2,3);
  Array < double > OutputData1(2,1);
  Array < double > TargetData1(2,1);

  InputData1(0,0) = 32;
  InputData1(1,0) = 45;
  //  InputData1(2,0) = 52;
  //InputData1(3,0) = 35;
  //  InputData1(4,0) = 67;
  //  InputData1(5,0) = 2;
  InputData1(0,1) = 89;
  InputData1(1,1) = 3;
  //  InputData1(2,1) = 43;
  //  InputData1(3,1) = 2;
  //  InputData1(4,1) = 34;
  //  InputData1(5,1) = 22;
  InputData1(0,2) = 34;
  InputData1(1,2) = 98;
  //  InputData1(2,2) = 87;
  //  InputData1(3,2) = 43;
  //  InputData1(4,2) = 23;
  //  InputData1(5,2) = 23;


  TargetData1(0,0) = 21;
  TargetData1(1,0) = 114;
  //  TargetData1(2,0) = 235;
  //  TargetData1(3,0) = 47;
  //  TargetData1(4,0) = 231;
  //  TargetData1(5,0) = 341;


  InputData1.resize((unsigned int)14,(unsigned int)2);
  TargetData1.resize((unsigned int)14,(unsigned int)1);

  InputData1(0,0) = 195.;
  InputData1(1,0) = 255.;
  InputData1(2,0) = 195.;
  InputData1(3,0) = 255.;
  InputData1(4,0) = 225;
  InputData1(5,0) = 225;
  InputData1(6,0) = 225.;
  InputData1(7,0) = 195.;
  InputData1(8,0) = 255.;
  InputData1(9,0) = 225.;
  InputData1(10,0) = 225.;
  InputData1(11,0) = 225.;
  InputData1(12,0) = 225.;
  InputData1(13,0) = 230.;

  InputData1(0,1) = 4;
  InputData1(1,1) = 4;
  InputData1(2,1) = 4.6;
  InputData1(3,1) = 4.6;
  InputData1(4,1) = 4.2;
  InputData1(5,1) = 4.1;
  InputData1(6,1) = 4.6;
  InputData1(7,1) = 4.3;
  InputData1(8,1) = 4.3;
  InputData1(9,1) = 4;
  InputData1(10,1) = 4.7;
  InputData1(11,1) = 4.3;
  InputData1(12,1) = 4.72;
  InputData1(13,1) = 4.3;

  //InputData1.resize((unsigned int)1,(unsigned int)2);
  //TargetData1.resize((unsigned int)1,(unsigned int)1);

  // InputData1(0,0) = 195;
  /*  InputData1(1,0) = 255;
  InputData1(2,0) = 195;
  InputData1(3,0) = 255;
  InputData1(4,0) = 225;
  InputData1(5,0) = 225;
  InputData1(6,0) = 225;
  InputData1(7,0) = 195;
  InputData1(8,0) = 255;
  InputData1(9,0) = 225;
  InputData1(10,0) = 225;
  InputData1(11,0) = 225;
  InputData1(12,0) = 225;
  InputData1(13,0) = 230;*/

  //InputData1(0,1) = 4;
  /*InputData1(1,1) = 4;
  InputData1(2,1) = 4.6;
  InputData1(3,1) = 4.6;
  InputData1(4,1) = 4.2;
  InputData1(5,1) = 4.1;
  InputData1(6,1) = 4.6;
  InputData1(7,1) = 4.3;
  InputData1(8,1) = 4.3;
  InputData1(9,1) = 4;
  InputData1(10,1) = 4.7;
  InputData1(11,1) = 4.3;
  InputData1(12,1) = 4.72;
  InputData1(13,1) = 4.3;*/

  TargetData1(0,0) = 1004;
  TargetData1(1,0) = 1636;
  TargetData1(2,0) = 852;
  TargetData1(3,0) = 1506;
  TargetData1(4,0) = 1272;
  TargetData1(5,0) = 1270;
  TargetData1(6,0) = 1269;
  TargetData1(7,0) = 903;
  TargetData1(8,0) = 1555;
  TargetData1(9,0) = 1260;
  TargetData1(10,0) = 1146;
  TargetData1(11,0) = 1276;
  TargetData1(12,0) = 1225;
  TargetData1(13,0) = 1321;

  //InputData1.resize((unsigned int)6,(unsigned int));
  //TargetData1.resize((unsigned int)6,(unsigned int)1);



  Database myData2;
  myData2.AddTrainingData(InputData1,TargetData1);
  //myData2.LoadTrainingData("TrainingDataCFD.dat",0);
  myData2.PrintTrainingData();
  RegrApp myRegression4(1,false);
  cout << "Regression 1 start!" << endl;
  myRegression4.Train(myData2);

  for ( int i = 0; i < myRegression4.columnDimension; i++ ) {
    cout << myRegression4.ModelCoefficients[i] << endl;
  }
  cout << myRegression4.columnDimension << endl;

  myRegression4.CalcNaturalCoefficients();
  for ( int i = 0; i < myRegression4.columnDimension; i++ ) {
    cout << myRegression4.NaturalModelCoefficients[i] << endl;
  }
  cout << myRegression4.columnDimension << endl;

  /*
  int start1 = 0;
  int end1 = InputData1.dim(0);
  Array<double> coefficient1(1,1);
  Array < double > InputDatainvers1(InputData1.dim(0),InputData1.dim(1));
  Array < double > TargetDatainvers1(TargetData1.dim(0),1);
  for ( int i = 0; i < InputData1.dim(0); i++) {
    for (int j = 0; j < InputData1.dim(1);j++) {
      InputDatainvers1( (InputData1.dim(0)-i-1), j) = InputData1 (i,j);
    }
    TargetDatainvers1( (InputData1.dim(0)-i-1), 0) = TargetData1 (i,0);
  }

  Model_first_order(InputDatainvers1,
		    TargetDatainvers1,
		    start1,
		    end1,
		    coefficient1);
  */



}
   
