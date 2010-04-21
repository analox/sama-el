#include "Database.h"
#include <iostream>

using namespace std;

// ==========================================================================

Database::Database() {

  NumTrainingData = 0;
  UsedTrainingData = 0;
  MaxTrainingData = 400;
  MaxArchiveLength = 1000;
  NewestEntry = -1;
  Array < double > I (0,1);
  Array < int > J (0,1);
  Input.resize ( I );
  Target.resize( I );
  Information.resize ( J );

}

// ==========================================================================

int Database::AddTrainingData ( Array < double > NewInput,
			        Array < double > NewTarget ) {

  Array < int > Generation ( NewInput.dim(0), 1 );
  for ( unsigned int i = 0; i < NewInput.dim(0); i++ ) {
    Generation ( i, 0 ) = -1;
  }
  AddTrainingData ( NewInput, 
		    NewTarget,
		    Generation );
  return 0;
}

// ==========================================================================

int Database::AddTrainingData ( Array < double > NewInput,
			        Array < double > NewTarget,
				int numberOfGeneration ) {

  Array < int > Generation ( NewInput.dim(0), 1 );
  for ( unsigned int i = 0; i < NewInput.dim(0); i++ ) {
    Generation ( i, 0 ) = numberOfGeneration;
  }
  AddTrainingData ( NewInput, 
		    NewTarget,
		    Generation );
  return 0;
}
// ==========================================================================

int Database::AddTrainingData ( Array < double > NewInput,
			        Array < double > NewTarget,
				Array < int > Generation ) {

  // if dimensions do not fit the dimension of the old training data return with an error
  if  ( Input.nelem() !=0 ) {
    if ((NewInput.dim(1) != Input.dim(1)) || (NewTarget.dim(1) != Target.dim(1))) {
      return 1;
    }
  }

  Array <double> d1( NewInput.dim(1) );
  Array <double> d2( NewTarget.dim(1) );
  Array <int> d3( 1 );
  
  for (unsigned int i = 0; i < NewInput.dim(0); i++) {
    for (unsigned int j = 0; j < NewInput.dim(1); j++) {
      d1(j) = NewInput(i,j);
    }
    for (unsigned int j = 0; j < NewTarget.dim(1); j++) {
      d2(j) = NewTarget(i,j);
    }
    for (unsigned int j = 0; j < 1; j++) {
      d3(j) = Generation(i,j);
    }

    AddElement (d1,d2,d3);
  }
  return 0;
}

// ==========================================================================

int Database::AddTrainingData ( vector < double > NewInput, 
				vector < double > NewTarget ) {

  int numberOfGeneration = -1;
  AddTrainingData ( NewInput,
		    NewTarget,
		    numberOfGeneration );
  return 0;
}

// ==========================================================================

int Database::AddTrainingData ( vector < double > NewInput, 
				vector < double > NewTarget,
				int numberOfGeneration ) {

  // if dimensions do not fit the dimension of the old training data return with an error
  if  ( Input.nelem() !=0 ) {
    if ((NewInput.size() != Input.dim(1)) || (NewTarget.size() != Target.dim(1))) {
      return 1;
    }
  }

  unsigned InputDimension = NewInput.size();
  unsigned TargetDimension = NewTarget.size();
  Array < double > InputArray(InputDimension);
  Array < double > TargetArray(TargetDimension);
  Array < int > Generation(1);
  for (unsigned int i = 0; i < InputDimension; i++) {
    InputArray(i) = NewInput[i];
  }
  for (unsigned int i = 0; i < TargetDimension; i++) {
    TargetArray(i) = NewTarget[i];
  }
  Generation(0) = numberOfGeneration;

  AddElement (InputArray, TargetArray, Generation);   
  return 0;
}

// ==========================================================================

int Database::AddElement ( Array < double > NewInput, 
			   Array < double > NewTarget,
			   Array < int > Generation ) {
    
  // --- if the data base is empty add the first set  
  if  (Input.dim(0)==0) {
    Input.resize ( (unsigned int)0, (unsigned int)NewInput.dim(0) );
    Target.resize ( (unsigned int)0, (unsigned int)NewTarget.dim(0) );
    Information.resize ( (unsigned int)0, (unsigned int)1 );
 
    Input.append_rows( NewInput );
    Target.append_rows( NewTarget );
    Information.append_rows ( Generation );

  }
  else {
    // --- if there are already data sets stored, then check whether there is still space left
    // --- if not, remove the oldest element from the archive
    if ((Input.dim(0) >= (unsigned int)MaxArchiveLength) ){
      Input.remove_row(0);
      Target.remove_row(0);
      Information.remove_row(0);
    } 
    // --- add the new data to the archive

    Input.append_rows( NewInput );
    Target.append_rows ( NewTarget );
    Information.append_rows ( Generation );

  }

  // --- store the actual number of available data points
  NumTrainingData = Input.dim(0);
  NewestEntry = Input.dim(0)-1;
  
  if (NumTrainingData > MaxTrainingData) 
    UsedTrainingData = MaxTrainingData;
  else
    UsedTrainingData = NumTrainingData;
  
  return 0;
}

// ==========================================================================

int Database::DeleteFirstData() {

  if ( Input.nelem() == 0 ) 
    return 1;
  DeleteDataAtPosition (0);
  return 0;
}

// ==========================================================================

int Database::DeleteLastData() {

  if ( Input.nelem() == 0 ) 
    return 1;
  
  Input.remove_row(Input.dim(0)-1);
  Target.remove_row(Target.dim(0)-1);
  Information.remove_row(Information.dim(0)-1);

  NumTrainingData = Input.dim(0);
  NewestEntry = Input.dim(0)-1;
  
  if (NumTrainingData > MaxTrainingData) 
    UsedTrainingData = MaxTrainingData;
  else
    UsedTrainingData = NumTrainingData;

  return 0;
}

// ==========================================================================

int Database::DeleteDataAtPosition( int Index ) {

  if ( (unsigned int)Index > (Input.dim(0)-1) ) 
    return 1;
  
  Input.remove_row(Index);
  Target.remove_row(Index);
  Information.remove_row(Index);

  NumTrainingData = Input.dim(0);
  NewestEntry = Input.dim(0)-1;
  
  if (NumTrainingData > MaxTrainingData) 
    UsedTrainingData = MaxTrainingData;
  else
    UsedTrainingData = NumTrainingData;

  return 0;
}

// ==========================================================================

int Database::DeleteDataFromIndex1ToIndex2 ( int Index1,
					     int Index2 ) {

  if ( Index1 < 0 ) return 1;
  if ( (unsigned int) Index2 > Input.dim(0)-1 ) return 1;

  for ( int i = Index2; i >= Index1; i-- ) {
    DeleteDataAtPosition ( i );
  }

  return 0;
}

// ==========================================================================

int Database::DeleteGeneration ( int numberOfGeneration ) {

  int upper = Information.dim(0)-1;
  for ( int i = upper; i >= 0; i-- ) {
    if ( Information(i,0) == numberOfGeneration )
      DeleteDataAtPosition(i);
  }

  return 0;
}

// ==========================================================================

int Database::GetLastData ( Array < double > &InputData, 
			    Array < double > &TargetData ) {
  
  Array < int > Generation (1,1);
  GetLastData ( InputData, TargetData, Generation );
  return 0;
}

// ==========================================================================

int Database::GetLastData ( Array < double > &InputData, 
			    Array < double > &TargetData,
			    Array < int > &numberOfGeneration ) {
  
  if (Input.ndim() == 0) {
    cout << "No training data available!" << endl;
    return 1;
  }
  InputData.resize ((unsigned int)(UsedTrainingData), (unsigned int)Input.dim(1), true);
  TargetData.resize ((unsigned int)(UsedTrainingData), (unsigned int)Target.dim(1), true);
  numberOfGeneration.resize ((unsigned int)(UsedTrainingData), (unsigned int)Information.dim(1), true);
  
  for (int i = (int)UsedTrainingData-1; i >= 0; i--) {
    for (unsigned int j = 0; j < Input.dim(1); j++)
      InputData(i,j) = Input(Input.dim(0) - i - 1,j);
    for (unsigned int j = 0; j < Target.dim(1); j++)
      TargetData(i,j) = Target(Target.dim(0) - i - 1, j);
    for ( unsigned int j = 0; j < Information.dim(1); j++ ) 
      numberOfGeneration(i,j) = Information ( Information.dim(0) - i - 1, j );
  }
  return 0;
}

// ==========================================================================

int Database::GetDataAtPosition( Array < double > &InputData, 
				 Array < double > &TargetData,
				 int Index ) {

  int number;
  GetDataAtPosition ( InputData, TargetData, number, Index );
  return 0;
}

// ==========================================================================

int Database::GetDataAtPosition( Array < double > &InputData, 
				 Array < double > &TargetData,
				 int &numberOfGeneration,
				 int Index ) {

  if ( (unsigned int) Index > (Input.dim(0)-1) ) 
    return 1;

  InputData.resize ( 1, (unsigned int)Input.dim(1), true);
  TargetData.resize ( 1, (unsigned int)Target.dim(1), true);

  for (unsigned int j = 0; j < Input.dim(1); j++)
    InputData(0,j) = Input(Index,j);
  for (unsigned int j = 0; j < Target.dim(1); j++)
    TargetData(0,j) = Target(Index,j);
  numberOfGeneration = Information ( Index, 0 );
  return 0;
}

// ==========================================================================

int Database::GetDataOfGeneration ( int numberOfGeneration,
				    Array < double >& InputData,
				    Array < double >& TargetData ) {

  int total = 0;
  for ( unsigned int i = 0; i < Information.dim(0); i++ ) {
    if ( Information(i,0) == numberOfGeneration )
      total++;
  }
  if ( total == 0 ) 
    return 1;
  
  InputData.resize( total, Input.dim(1) );
  TargetData.resize( total, Target.dim(1) );
  int numm = 0;
  for ( unsigned int i = 0; i < Information.dim(0); i++ ) {
    if ( Information(i,0) == numberOfGeneration ) {
      for ( unsigned int j = 0; j < Input.dim(1); j++ ) {
	InputData(numm,j) = Input(i,j);
      }
      for ( unsigned int j = 0; j < Target.dim(1); j++ ) {
	TargetData(numm,j) = Target(i,j);
      }
      numm++;
    }
  }
  
  return 0;
}

// ==========================================================================

int Database::LoadTrainingData ( string Filename, 
				 int add) {
  
  FILE *fin;
  
  if ((fin = fopen(Filename.c_str(), "rt")) == NULL) {
    return 1;
  }
    
  int NumberOfDataToRead;
  int NewInputDimension;
  int NewOutputDimension;
  
  fscanf(fin,"%i %i %i", &NumberOfDataToRead, &NewInputDimension, &NewOutputDimension);
  
  if (add) {
    // if dimensions do not fit the dimension of the old training data
    // return with an error
    if (((unsigned int)NewInputDimension != Input.dim(1)) || ((unsigned int)NewOutputDimension != Target.dim(1))) {
      return 2;
    }
  }
  // --- remove old training data if add == FALSE
  if ( add == 0 ) {
    int inp = Input.dim(0);
    for ( int i = 0; i < inp; i++) {
      Input.remove_row(0);
      Target.remove_row(0);
      Information.remove_row(0);
    }
    NumTrainingData = 0;
    NewestEntry = -1;
    UsedTrainingData = 0;
  }

  Array < double > NewInput(NewInputDimension);
  Array < double > NewTarget(NewOutputDimension);
  Array < int > Generation(1);

  // --- read new data
  double ftmp;
  for(unsigned int i = 0; i < (unsigned int)NumberOfDataToRead; i++){
    // --- if there are already data sets stored, then check whether there is still space left
    // --- if not, remove the oldest element from the archive
    if (NumTrainingData > (unsigned int)MaxArchiveLength) {
      Input.remove_row(0);
      Target.remove_row(0);
      Information.remove_row(0);
    }

    for (unsigned int j = 0; j < (unsigned int)NewInputDimension; j++) {
      fscanf(fin, "%lf", &ftmp);
      NewInput(j) = ftmp;
    }
    for (unsigned int j = 0; j < (unsigned int)NewOutputDimension; j++) {
      fscanf(fin, "%lf", &ftmp);
      NewTarget(j) = ftmp;
    } 
    int iftmp;
    fscanf(fin, "%i", &iftmp );
    Generation(0) = iftmp;

    AddElement(NewInput,NewTarget,Generation);
  }    

  if (NumTrainingData > MaxTrainingData) 
    UsedTrainingData = MaxTrainingData;
  else
    UsedTrainingData = NumTrainingData;

  fclose(fin);
  return 0;
}

// ==========================================================================

int Database::SaveTrainingData( string Filename ) {

  FILE* fout;

  if ((fout = fopen(Filename.c_str(), "wt")) == NULL) {
    return 1;
  }
  
  if (Input.ndim() == 0) {
    cout << "No data in database!" << endl;
    return 2;
  }
  
  int InputDimension = Input.dim(1);
  int OutputDimension = Target.dim(1);

  fprintf(fout, "%i %i %i \n", Input.dim(0), InputDimension, OutputDimension);

  for (unsigned int i = 0; i < Input.dim(0); i++) {
    for (unsigned int j = 0; j < Input.dim(1); j++) {
      fprintf(fout, "%f ", Input(i,j));
    }
    fprintf(fout, "\n");
    for (unsigned int j = 0; j < Target.dim(1); j++) {
      fprintf(fout, "%f ", Target(i,j));
    }
    fprintf(fout, "\n");
    fprintf(fout, "%i", Information(i,0));
    fprintf(fout, "\n");
    fflush(fout);
  }

  fclose(fout);
  return 0;
}

// ==========================================================================

int Database::PrintTrainingData() {

  if ( Input.dim(0) == 0 )
    cout << "Database is empty!" << endl;
  for (unsigned i = 0; i < Input.dim(0); i++) { 
    cout << "Number " << i << ": ";
    for(unsigned j = 0; j < Input.dim(1); j++){
      cout << Input(i,j) << " ";
    }
    cout << "    ";
    for(unsigned j = 0; j < Target.dim(1); j++){
      cout << Target(i,j) << " ";
    }
    cout << "    ";
    cout << Information(i,0);
    cout << endl;
  }
  return 0;
}

// ==========================================================================

int Database::GetNumTrainingData() {
  
  return NumTrainingData;

}

// ==========================================================================

int Database::GetMaxTrainingData() {
  
  return MaxTrainingData;

}

// ==========================================================================

int Database::GetUsedTrainingData() {
  
  return UsedTrainingData;

}

// ==========================================================================

int Database::GetMaxArchiveLength() {
  
  return MaxArchiveLength;

}

// ==========================================================================

int Database::GetNewestEntry() {
  
  return NewestEntry;

}

// ==========================================================================

int Database::SetMaxArchiveLength(int NewLength) {
  
  MaxArchiveLength = NewLength;

  // check if the maximal length is greater than the actual number of data
  if ( Input.dim(0) > (unsigned int)MaxArchiveLength ) {
    while ( Input.dim(0) > (unsigned int)MaxArchiveLength ) {
      Input.remove_row(0);
      Target.remove_row(0);
      Information.remove_row(0);
    }
    NumTrainingData = Input.dim(0);
    NewestEntry = Input.dim(0)-1;
  
    if (NumTrainingData > MaxTrainingData) 
      UsedTrainingData = MaxTrainingData;
    else
      UsedTrainingData = NumTrainingData;
    return 1;
  }

  return 0;
}

// ==========================================================================

int Database::SetMaxTrainingData(int NewMax) {
  
  MaxTrainingData = NewMax;

  // check if the new maximal number is greater than the actual number

  if (NumTrainingData > MaxTrainingData) 
    UsedTrainingData = MaxTrainingData;
  else
    UsedTrainingData = NumTrainingData;

  return 0;
}


