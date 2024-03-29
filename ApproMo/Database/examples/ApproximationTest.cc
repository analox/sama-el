#include <stdio.h>
#include <Database.h>

using namespace std;

int main( int argc, char **argv)
{

  int a = 5; int b = 6; int c = 1;
  Array < double > InputData(a,b);
  Array < double > OutputData(a,c);
  Array < double > TargetData(a,c);

  for (int i = 0; i < a; i++) {
    for (int j = 0; j < b; j++) {
      InputData(i,j) = i+j;
      cout << "Point: " << i << " " << j << ": " << InputData(i,j) << endl;
    }
  }
  
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < c; j++) {
      TargetData(i,j) = i+j;
      cout << "Point: " << i << " " << j << ": " << TargetData(i,j) << endl;
    }
  }

  Database myData;

  myData.AddTrainingData(InputData,TargetData);

  myData.SaveTrainingData("TestData.dat");
  
  cout << "Data saved!" << endl;

  myData.LoadTrainingData("TestData.dat",1);

  cout << "Data loaded!" << endl;

  myData.SaveTrainingData("TestData.dat");

  myData.PrintTrainingData();
  cout << endl;

  myData.LoadTrainingData("TestData.dat",1);

  myData.PrintTrainingData();
  cout << endl;

  myData.DeleteLastData();

  myData.DeleteDataAtPosition(5);

  myData.PrintTrainingData();

  cout << myData.GetMaxArchiveLength() << endl;

  cout << myData.GetNumTrainingData() << endl;
  
  cout << myData.GetMaxTrainingData() << endl;

  cout << myData.GetUsedTrainingData() << endl;

  cout << myData.GetNewestEntry() << endl;

  myData.SetMaxArchiveLength(10);

  myData.PrintTrainingData();

  cout << myData.GetMaxArchiveLength() << endl;

  cout << myData.GetNumTrainingData() << endl;
  
  cout << myData.GetMaxTrainingData() << endl;

  cout << myData.GetUsedTrainingData() << endl;

  cout << myData.GetNewestEntry() << endl;

  myData.SetMaxTrainingData(400);

  myData.PrintTrainingData();

  cout << myData.GetMaxArchiveLength() << endl;

  cout << myData.GetNumTrainingData() << endl;
  
  cout << myData.GetMaxTrainingData() << endl;

  cout << myData.GetUsedTrainingData() << endl;

  cout << myData.GetNewestEntry() << endl;

  cout << "Everything is done!" << endl;

  Array < double > A1 (1,1);
  Array < double > A2 (1,1);
  myData.GetLastData (A1,A2);
  for ( int i = 0; i < (int)A1.dim(0); i++ ) {
    for ( int j = 0; j < (int)A1.dim(1); j++ ) {
      cout << A1(i,j);
    }
    cout << A2(i,0) << endl;
  }
}
