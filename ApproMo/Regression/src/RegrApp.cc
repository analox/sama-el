#include "RegrApp.h"
#include "sqr.h"

using namespace std;

RegrApp::RegrApp() {
  
}

//=======================================================================

RegrApp::RegrApp( int degree ) {

  degOfModel = degree;
  interactionTerms = false;
  InputDataIsNormalized = false;
  ModelCoefficients = new double[0];
  maxInColumn = new double[0];
  minInColumn = new double[0];
  ModelIsTrained = false;
}

//=======================================================================

RegrApp::RegrApp( int degree,
		  bool interact) {

  degOfModel = degree;
  interactionTerms = interact;
  InputDataIsNormalized = false;
  ModelCoefficients = new double[0];
  maxInColumn = new double[0];
  minInColumn = new double[0];
  ModelIsTrained = false;
}

//=======================================================================

RegrApp::~RegrApp() {

  delete[] ModelCoefficients;
  delete[] maxInColumn;
  delete[] minInColumn;
}

//=======================================================================

double RegrApp::MSE( Database Data ) {

  Array < double > InputData(1,1);
  Array < double > TargetData(1,1);
  Data.GetLastData(InputData,TargetData);
  
  Array < double > OutputData(1,1);
  Evaluate ( InputData, OutputData );

  double rms_LB = 0;
  for ( int i = 0; i < (int)InputData.dim(0); i++) {
    rms_LB += sqr( OutputData(i,0) - TargetData(i,0) );
  }

  return sqrt(rms_LB/(double)InputData.dim(0));
}

//=======================================================================

double RegrApp::MSE( Array < double > Input,
		     Array < double > Target ) {

  Array < double > InputData(1,1);
  Array < double > TargetData(1,1);
  InputData.resize ( Input.dim (0), Input.dim(1) );
  TargetData.resize ( Target.dim(0), Target.dim(1) );
  for ( int i = 0; i < (int)Input.dim(0); i++ ) {
    for ( int j = 0; j < (int)Input.dim(1); j++ ) {
      InputData (i,j) = Input (i,j);
      //TargetData (i,j) = Target (i,j);
    }
	TargetData(i,0) = Target( i, 0 );
  }
  
  Array < double > OutputData(InputData.dim(0),1);
  Evaluate ( InputData, OutputData );

  double rms_LB = 0;
  for ( int i = 0; i < (int)InputData.dim(0); i++) {
    rms_LB += sqr( OutputData(i,0) - TargetData(i,0) );
  }

  return sqrt(rms_LB/(double)InputData.dim(0));
}

//=======================================================================

double RegrApp::MSE( Population offsprings ) {

  // the parameters have to be in the first Chromosome!!!
  unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > (offsprings[0][0]).size();
  unsigned numberOfOffsprings = offsprings.size();

  Array < double > InputData (numberOfOffsprings,InputDimension);
  Array < double > TargetData (numberOfOffsprings, 1);

  for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
    for ( unsigned j = 0; j < InputDimension; j++ ) {
      InputData(i,j) = dynamic_cast < ChromosomeT < double >& > (offsprings[i][0] )[j];
    }
    TargetData (i,0) = offsprings[i].fitnessValue();
  }  

  double error = MSE ( InputData, TargetData );

  return error;
}

//======================================================================

int RegrApp::Evaluate( Population& offsprings ) {

  if ( !ModelIsTrained ) {
    cout << "Regression model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
    return 1;
  }
  unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > (offsprings[0][0]).size();
  unsigned numberOfOffsprings = offsprings.size();

  Array < double > NNInput (numberOfOffsprings,InputDimension);
  Array < double > NNOutput (numberOfOffsprings, 1);
  
  for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
    for ( unsigned j = 0; j < InputDimension; j++ ) {
      NNInput(i,j) = dynamic_cast < ChromosomeT < double >& > (offsprings[i][0] )[j];
    }
  }
  Evaluate(NNInput, NNOutput);

  for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
    offsprings[i].setFitness(NNOutput(i,0));
  }
  
  return 0;
}

//=======================================================================

int RegrApp::Evaluate( Individual& offspring ) {

  if ( !ModelIsTrained ) {
    cout << "Regression model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
    return 1;
  }
  unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > (offspring[0]).size();
  Array < double > NNInput (1,InputDimension);
  Array < double > NNOutput (1,1);

  for ( unsigned j = 0; j < InputDimension; j++ ) {
    NNInput(0,j) = dynamic_cast < ChromosomeT < double >& > (offspring[0])[j];
  }
  Evaluate(NNInput, NNOutput);

  offspring.setFitness(NNOutput(0,0));
  
  return 0;
}

//=======================================================================

int RegrApp::Evaluate(  Array < double > InputData,
			Array < double >& OutputData ) {

  if ( !ModelIsTrained ) {
    cout << "Regression model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
    return 1;
  }
  Array < double > ExtInputData (InputData.dim(0), (InputData.dim(1)+1) );
  Array < double > CompleteInputData ( InputData.dim(0) , columnDimension );
  OutputData.resize ( (unsigned int)InputData.dim(0), (unsigned int)1 );
  bool interact = interactionTerms;
  int degree = degOfModel;
  
  if ( !InputDataIsNormalized ) {
    for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
      ExtInputData(i,0) = 1.0;
      for ( int j = 0; j < (int)InputData.dim(1); j++ ) {
	ExtInputData(i,j+1) = InputData(i,j);
      }
    }
  }
  else {
    for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
      ExtInputData(i,0) = 1.0;
      for ( int j = 0; j < (int)InputData.dim(1); j++ ) {
	ExtInputData(i,j+1) = (InputData(i,j)-(maxInColumn[j] + minInColumn[j])/2)/((maxInColumn[j]-minInColumn[j])/2);
	// cout << "ExtData: " << ExtInputData(i,j+1) << endl;
      }
    }
  }

  if ( !interact ) {
    // format in matrix: 1, x0, x0p2, x0p3, ..., x1, x1p2, x1p3, ..., xn, xnp2, xnp3, ...
    double tempValue;
    int arrayPosition;
    for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
      CompleteInputData(i,0) = ExtInputData(i,0);
      arrayPosition = 1;
      for ( int j = 1; j < (int)ExtInputData.dim(1); j++ ) {
	tempValue = ExtInputData(i,j);
	CompleteInputData(i,arrayPosition) = tempValue;
	arrayPosition++;
	for ( int k = 1; k < degree; k++ ) {
	  tempValue *= tempValue;
	  CompleteInputData(i,arrayPosition) = tempValue;
	  arrayPosition++;
	}
      }
    }
  }
  else {
    int arrayPosition;
    vector < int > counterOnInputData(degree);
    double tempNumber;
    for ( int h = 0; h < (int)InputData.dim(0); h++ ) {
      for ( int i = 0; i < degree; i++ ) 
	counterOnInputData[i] = 0;
      bool DimensionOfColumnsDetermined = false;
      arrayPosition = 0;
      while ( !DimensionOfColumnsDetermined ) {
	tempNumber = 1;
	for ( int i = 0; i < degree; i++ ) {
	  tempNumber *= ExtInputData( h, counterOnInputData[i] );
	}
	CompleteInputData( h, arrayPosition ) = tempNumber;
	arrayPosition++;
	counterOnInputData[degree-1]++;
	for ( int i = (degree-1); i >= 1; i-- ) {
	  if ( counterOnInputData[i] > (int)InputData.dim(1) ) {
	    counterOnInputData[i-1]++;
	  }
	}
	for ( int i = 0; i < degree; i++ ) {
	  if ( counterOnInputData[i] > (int)InputData.dim(1) ) {
	    counterOnInputData[i] = counterOnInputData[i-1];
	  }
	}
	if ( counterOnInputData[0] > (int)InputData.dim(1) )
	  DimensionOfColumnsDetermined = true;
      }
    }
  }

  //for ( int i = 0; i < CompleteInputData.dim(1); i++ )
  // cout << i << " " << CompleteInputData(0,i) << endl;

  double result;
  for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
    result = 0;
    for ( int j = 0; j < (int)CompleteInputData.dim(1); j++ ) { 
      result += CompleteInputData(i,j) * ModelCoefficients[j];
    }
    OutputData(i,0) = result;
    // cout << OutputData(i,0) << endl;
  }

  return 0;
}

//=======================================================================

int RegrApp::Train( Array < double > InputData,
		    Array < double > TargetData ) {

  // default 
  Train_LSM ( InputData, TargetData );
  return 0;
}

//=======================================================================

int RegrApp::Train_LSM ( Array < double > InputData,
			 Array < double > TargetData ) {

  Database tempForTraining;
  tempForTraining.AddTrainingData ( InputData,
				    TargetData );
  Train_LSM ( tempForTraining );
  return 0;
}

//=======================================================================

int RegrApp::Train( Database Data ) {

  // default 
  Train_LSM ( Data );
  return 0;
}

//=======================================================================

int RegrApp::Train_LSM( Database Data ) {

  if ( degOfModel < 1 ) {
    cout << "The degree of the regression model must be greater 0!";
    return 1;
  }
  int degree = degOfModel;
  bool interact = interactionTerms;
  //int UsedData = Data.GetUsedTrainingData();
  Array < double > InputData(1,1);
  Array < double > TargetData(1,1);

  Data.GetLastData(InputData,TargetData);

  Array < double > NormInputData(InputData.dim(0),InputData.dim(1));

  NormalizeData(InputData, NormInputData);
  /* for ( int j = 0; j < (int)InputData.dim(1); j++ ) {
    for ( int i = 0; i < (int) InputData.dim(0); i++ ) {
      NormInputData(i,j) = InputData(i,j);
      cout << NormInputData(i,j) << endl;
    }
  }
  */
  // Matrix with a first column with the values "1" for convenience
  Array < double > ExtNormInputData( InputData.dim(0), (InputData.dim(1)+1) );
  for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
    ExtNormInputData(i,0) = 1.0;
    for ( int j = 1; j < ( (int)InputData.dim(1)+1); j++ ) {
      ExtNormInputData(i,j) = NormInputData(i,j-1);
      //cout << i << ", " << j << ": " << ExtNormInputData(i,j) << endl;
    }
  }

  // determine matrix dimensions
  int rowDimension = 0;
  columnDimension = 0;
  rowDimension = InputData.dim(0);

  vector <int> counterOnInputData(degree);
  double tempNumber;
  for ( int i = 0; i < degree; i++ ) 
    counterOnInputData[i] = 0;
  bool DimensionOfColumnsDetermined = false;

  if ( !interact ) {
    columnDimension = InputData.dim(1) * degree + 1;
  }
  else {
    while ( !DimensionOfColumnsDetermined ) {
      tempNumber = 1;
      for ( int i = 0; i < degree; i++ ) {
	tempNumber *= ExtNormInputData( 0, counterOnInputData[i] );
      }
      columnDimension++;
      counterOnInputData[degree-1]++;
      for ( int i = (degree-1); i >= 1; i-- ) {
	if ( counterOnInputData[i] >  (int)InputData.dim(1) ) {
	  counterOnInputData[i-1]++;
	}
      }
      for ( int i = 0; i < degree; i++ ) {
	if ( counterOnInputData[i] >  (int)InputData.dim(1) ) {
	  counterOnInputData[i] = counterOnInputData[i-1];
	}
      }
      if ( counterOnInputData[0] > (int)InputData.dim(1) )
	DimensionOfColumnsDetermined = true;
    }
  }
      
  // generating input matrix X, pay attention to size: 1..r, 1..c
  double **X;
  X = dmatrix((long)1, (long)InputData.dim(0),(long)1, (long)(columnDimension));
  cout << "Matrix dimension of columns: " << columnDimension << endl;
  
  if ( !interact ) {
    // format in matrix: 1, x0, x0p2, x0p3, ..., x1, x1p2, x1p3, ..., xn, xnp2, xnp3, ...
    double tempValue;
    int matrixPosition;
    for ( int i = 0; i < (int)ExtNormInputData.dim(0); i++ ) {
      X[i+1][1] = ExtNormInputData(i,0);
      matrixPosition = 2;
      for ( int j = 1; j < (int)ExtNormInputData.dim(1); j++ ) {
	tempValue = ExtNormInputData(i,j);
	X[i+1][matrixPosition] = tempValue;
	matrixPosition++;
	for ( int k = 1; k < degree; k++ ) {
	  tempValue *= tempValue;
	  X[i+1][matrixPosition] = tempValue;
	  matrixPosition++;
	}
      }
    }
  }
  else {
    // example for the format in matrix: x0x0x0, x0x0x1, x0x0x2, x0x1x1, x0x1x2, x0x2x2, x1x1x1, x1x1x2, x1x2x2, x2x2x2
    int matrixPosition;
    for ( int h = 0; h < (int)InputData.dim(0); h++ ) {
      for ( int i = 0; i < degree; i++ ) 
	counterOnInputData[i] = 0;
      DimensionOfColumnsDetermined = false;
      matrixPosition = 1;
      while ( !DimensionOfColumnsDetermined ) {
	tempNumber = 1;
	for ( int i = 0; i < degree; i++ ) {
	  tempNumber *= ExtNormInputData( h, counterOnInputData[i] );
	}
	X[h+1][matrixPosition] = tempNumber;
	matrixPosition++;
	counterOnInputData[degree-1]++;
	for ( int i = (degree-1); i >= 1; i-- ) {
	  if ( counterOnInputData[i] >  (int)InputData.dim(1) ) {
	    counterOnInputData[i-1]++;
	  }
	}
	for ( int i = 0; i < degree; i++ ) {
	  if ( counterOnInputData[i] >  (int)InputData.dim(1) ) {
	    counterOnInputData[i] = counterOnInputData[i-1];
	  }
	}
	if ( counterOnInputData[0] >  (int)InputData.dim(1) )
	  DimensionOfColumnsDetermined = true;
      }
    }
  }

  for ( int i = 1; i <= (int)InputData.dim(0); i++) {
    for ( int j = 1; j <= columnDimension; j++ ) {
      //cout << "Matrix: " << i << " " << j << ": " << X[i][j] << endl;
    }
  }

  // generating target vector, pay attention to size: 1..r
  double *y;
  y = dvector( 1, InputData.dim(0) );
  for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
    y[i+1] = TargetData(i,0);
    //cout << "Target " << i << ": " << TargetData(i,0)<< endl;
  }

  // calculate A = X'X
  double **A;
  A = dmatrix (1, columnDimension , 1, columnDimension );
  for ( int i = 1; i <= columnDimension; i++ ) {
    for ( int j = 1; j <= columnDimension; j++ ) {
      A[i][j] = 0.0;
    }
  }
  for ( int i = 1; i <= columnDimension; i++ ) {
    for ( int j = 1; j <= columnDimension; j++ ) {
      for ( int k = 1; k <= (int)InputData.dim(0); k++ ) {
	A[i][j] += X[k][i]*X[k][j];
      }
    }
  }

  /*
  for ( int i = 1; i <= columnDimension; i++ ) {
    for ( int j = 1; j <= columnDimension; j++ ) {
      cout << i << j << ": " << A[i][j] << " ";
    }
    cout << endl;
  }
  */

  // calculate b = X'y
  double *b;
  b = dvector (1, columnDimension);
  for ( int i = 1; i <= columnDimension; i++ ) {
    b[i] = 0.0;
  }
  for ( int i = 1; i <= (int)InputData.dim(0); i++ ) {
    for ( int j = 1; j <= columnDimension; j++ ) {
      // cout << " Matrixwert " << i << j << ": " << X[i][j] << endl;
      // cout << " Target: " << y[i] << endl;
      b[j] += X[i][j] * y[i];
    }
  }
  
  for ( int j = 1; j <= columnDimension; j++ ) {
    //cout << j << ": " << b[j] << endl;
  }

  int *indx;
  indx = ivector(1,columnDimension);
  double d;
  delete[] ModelCoefficients;
  ModelCoefficients = new double[columnDimension];

  ludcmp(A, columnDimension, indx, &d);
  lubksb(A, columnDimension, indx, b);

  for ( int i = 1; i <= columnDimension; i++ ) {
    ModelCoefficients[i-1] = b[i];
    ModelIsTrained = true;
    // cout << ModelCoefficients[i-1] << endl;
  }

  free_ivector(indx,1,columnDimension);
  free_dmatrix(X,1,InputData.dim(0),1,columnDimension);
  free_dmatrix(A,1,columnDimension,1,columnDimension);
  free_dvector(y,1,InputData.dim(0));
  free_dvector(b,1,columnDimension);

  return 0;
}

//=======================================================================

int RegrApp::NormalizeData( Array < double > InputData,
			    Array < double >& NormInputData ) {

  int numberOfParameters = InputData.dim(1);
  delete[] maxInColumn;
  delete[] minInColumn;
  maxInColumn = new double[numberOfParameters];
  minInColumn = new double[numberOfParameters];
  for ( int i = 0; i < (int)InputData.dim(1); i++ ) {
    maxInColumn[i] = -1e10;
    minInColumn[i] = 1e10;
    for ( int j = 0; j < (int)InputData.dim(0); j++ ) {
      if ( InputData(j,i) > maxInColumn[i] ) maxInColumn[i] = InputData(j,i);
      if ( InputData(j,i) < minInColumn[i] ) minInColumn[i] = InputData(j,i);
    }
    //cout << "Maximum in " << i << ": " << maxInColumn[i] << endl;
    //cout << "Minimum in " << i << ": " << minInColumn[i] << endl;
  }
  
  InputDataIsNormalized = true;
  for ( int i = 0; i < (int)InputData.dim(1); i++ ) {
    if ( InputDataIsNormalized )
      if ( maxInColumn[i] <= minInColumn[i] )
	InputDataIsNormalized = false;
  }
  
  // InputDataIsNormalized = false;
  if ( !InputDataIsNormalized ) {
    for (int i = 0; i < (int)InputData.dim(0); i++ ) {
      for ( int j = 0; j < (int)InputData.dim(1); j++ ) {
	NormInputData(i,j) = InputData(i,j);
      }
    }
  }
  else {
    for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
      for ( int j = 0; j < (int)InputData.dim(1); j++ ) {
	NormInputData(i,j) = (InputData(i,j)-(maxInColumn[j] + minInColumn[j])/2)/((maxInColumn[j]-minInColumn[j])/2);
	//if ( (NormInputData(i,j) < 1e-10) && (NormInputData(i,j) > -1e-10) ) NormInputData(i,j) = 0.;
	//cout << InputData(i,j) << " " << maxInColumn[j] << " " << minInColumn[j] << " " << NormInputData(i,j) << " : " <<  (InputData(i,j)-(maxInColumn[j] + minInColumn[j])/2.0) <<  endl;
      }
    }
  }

  return 0;
}

//=======================================================================

int RegrApp::ScanSquare3D( int lowerBorder,
			   int upperBorder,
			   Array < double >& OutputData ) {

  int i = 0;
  int j = 1;
  Array < double > values (0);
  ScanSquare3D ( lowerBorder,
		 upperBorder,
		 OutputData,
		 i,
		 j,
		 values );
  return 0;
}

//=======================================================================

int RegrApp::CalcNaturalCoefficients() {

  NaturalModelCoefficients = new double[columnDimension];
  for ( int i = 1; i < columnDimension; i++ ) {
    //cout << maxInColumn[i-1] << " " << minInColumn[i-1] << endl;
    NaturalModelCoefficients[i] = ModelCoefficients[i]*1./((maxInColumn[i-1]-minInColumn[i-1])/2);
  }
  NaturalModelCoefficients[0] = ModelCoefficients[0];
  for ( int i = 1; i < columnDimension; i++ ) { 
    NaturalModelCoefficients[0] = NaturalModelCoefficients[0]-ModelCoefficients[i]* ((maxInColumn[i-1] + minInColumn[i-1])/2)/((maxInColumn[i-1]-minInColumn[i-1])/2);
    cout << ModelCoefficients[i] << endl;
  }
  return 0;
}

//=======================================================================

int RegrApp::ScanSquare3D( int lowerBorder,
			   int upperBorder,
			   Array < double >& OutputData,
			   int ind1,
			   int ind2,
			   Array < double > ParametersToKeepConstant ) {

  Array < double > InputData( 1,2 );
  Array < double > Output( 1,1 );
  Array < double > DataToEvaluate ( 1, ParametersToKeepConstant.dim(0)+2 );

  int resolution = 100;
  int diff = upperBorder - lowerBorder;
  OutputData.resize( (unsigned int)sqr(resolution+1), (unsigned int)3);
  int outputEntry = 0;
  for ( int i = 0; i <= resolution; i++ ) {
    for ( int j = 0; j <= resolution; j++ ) {
      InputData(0,0) = (double)i*diff/resolution+lowerBorder;
      InputData(0,1) = (double)j*diff/resolution+lowerBorder;
      OutputData(outputEntry,0) = (double)i*diff/resolution+lowerBorder;
      OutputData(outputEntry,1) = (double)j*diff/resolution+lowerBorder;
      
      if ( DataToEvaluate.dim(1) > 2 ) {
	int counter = 0;
	for ( int k = 0; k < (int)DataToEvaluate.dim(1); k++ ) {
	  if ( k == ind1 ) {
	    DataToEvaluate ( 0, k ) = InputData (0,0);
	  }
	  else {
	    if ( k == ind2 ) {
	      DataToEvaluate ( 0, k ) = InputData ( 0, 1 );
	    }
	    else {
	      DataToEvaluate ( 0, k ) = ParametersToKeepConstant ( counter );
	      counter++;
	    }
	  }
	}
      }
      else {
	DataToEvaluate ( 0, 0 ) = InputData ( 0, 0 );
	DataToEvaluate ( 0, 1 ) = InputData ( 0, 1 );
      }

      Evaluate(DataToEvaluate,Output);
      OutputData(outputEntry,2) = Output(0,0);
      outputEntry++;
    }
  }

  return 0;
}

//=======================================================================

int RegrApp::GetDegreeOfRegrModel() {

  return degOfModel;
}

//=======================================================================

int RegrApp::SetDegreeOfRegrModel( int newDegree ) {

  degOfModel = newDegree;
  return 0;
}

//=======================================================================

bool RegrApp::interactionTermsConsidered() {

  return interactionTerms;
}

//=======================================================================

int RegrApp::SetInteractionTerms ( bool choice ) {

  interactionTerms = choice;
  return 0;
}

//=======================================================================
