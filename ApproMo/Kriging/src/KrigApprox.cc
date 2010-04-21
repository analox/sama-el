#include "KrigApprox.h"

using namespace std;

KrigApprox::KrigApprox() {

  KrigApp = new krigapprox();
  ModelIsTrained = false;
}

//=======================================================================

KrigApprox::~KrigApprox() {
  
  delete KrigApp;
}

//=======================================================================

int KrigApprox::Train ( Array < double > InputData,
			Array < double > TargetData ) {

  // default
  Train_MLM ( InputData, TargetData );
  return 0;
}


//=======================================================================

int KrigApprox::Train_MLM ( Array < double > InputData,
			    Array < double > TargetData ) {

  Database tempForTraining;
  tempForTraining.AddTrainingData ( InputData,
				    TargetData );
  Train_MLM ( tempForTraining );
  return 0;
}

//=======================================================================

int KrigApprox::Train ( Database Data ) {

  // default
  Train_MLM ( Data );
  return 0;
}

//=======================================================================

int KrigApprox::Train_MLM ( Database Data ) {

  Array < double > InputData(1,1);
  Array < double > TargetData(1,1);
  Data.GetLastData( InputData,
		    TargetData );

  Matrix < double > InputMatrix( InputData.dim(0), InputData.dim(1) );
  Vector < double > TargetVector( TargetData.dim(0) );
  for ( int i = 0; i < (int)InputData.dim(0); i++) {
    for ( int j = 0; j < (int)InputData.dim(1); j++) {
      InputMatrix[(long)i][(long)j] = InputData(i,j);
    }
    TargetVector[(long)i] = TargetData(i,0);
  }
  
  KrigApp->setvalues( InputMatrix, TargetVector, 1 );
  ModelIsTrained = true;

  return 0;
}

//=======================================================================

double KrigApprox::MSE( Database Data ) {

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

double KrigApprox::MSE( Array < double > input, Array < double > target ) {
	Array < double > inputData( 1, 1 );
	Array < double > targetData( 1, 1 );
	inputData.resize ( input.dim ( 0 ), input.dim( 1 ) );
	targetData.resize ( target.dim( 0 ), target.dim( 1 ) );
	for ( int i = 0; i < ( int ) input.dim( 0 ); i++ ) {
		for ( int j = 0; j < ( int ) input.dim( 1 ); j++ ) {
			inputData ( i, j ) = input ( i, j );
		}
		targetData ( i, 0 ) = target ( i, 0 );
	}
		

	cout << "in Kriging MSE" << endl;

	Array < double > outputData ( inputData.dim( 0 ), 1 );
	Evaluate ( inputData, outputData );
			
	cout << " in Kriging MSE , after Evaluate " << endl;

	double rms_LB = 0;
	for ( int i = 0; i < ( int ) inputData.dim( 0 ); i++ ) {
		rms_LB += sqr( outputData( i, 0 ) - targetData( i, 0 ) );
	}
	

	cout << "in Kriging MSE, before returning " << endl;

	return ( rms_LB / ( double ) inputData.dim( 0 ) );
}


//=======================================================================

double KrigApprox::MSE( Population offsprings ) {

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

//=======================================================================

int KrigApprox::Evaluate( Population& offsprings ) {

  if ( !ModelIsTrained ) {
    cout << "Kriging model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
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

int KrigApprox::Evaluate( Individual& offspring ) {

  if ( !ModelIsTrained ) {
    cout << "Kriging model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
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

int KrigApprox::Evaluate( Array < double > InputData,
			  Array < double >& OutputData ) {

  if ( !ModelIsTrained ) {
    cout << "Kriging model has to be trained first before it can be used for evaluation! Please train it first!" << endl;
    return 1;
  }
  Vector < double > InputVector( InputData.dim(1) );
  OutputData.resize( InputData.dim(0) , (unsigned int)1);
  for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
    for ( int j = 0; j < (int)InputData.dim(1); j++) {
      InputVector[(long)j] = InputData(i,j);
    }
    OutputData(i,0) = KrigApp->evalf(InputVector);
  }

  return 0;
}

//=======================================================================

int KrigApprox::ScanSquare3D( int lowerBorder,
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

int KrigApprox::ScanSquare3D( int lowerBorder,
			      int upperBorder,
			      Array < double >& OutputData,
			      int ind1,
			      int ind2,
			      Array < double > ParametersToKeepConstant ) {

  Array < double > Inputs ( 1, ParametersToKeepConstant.dim(0)+2 );
  Vector < double > InputVector( 2 );
  Vector < double > DataToEvaluate( 2 + ParametersToKeepConstant.dim(0) );
  int resolution = 400;
  int diff = upperBorder - lowerBorder;
  OutputData.resize( (unsigned int)sqr(resolution+1), (unsigned int)3);
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

  return 0;
}

//=======================================================================
