#include "InterpolateRBF.h"
#include <iostream.h>
#include "EALib/sqr.h"
#include "LinAlg/linalg.h"
#include <iomanip.h>


//=======================================================================
InterpolateRBF::InterpolateRBF() { }

//=======================================================================
InterpolateRBF::~InterpolateRBF() { }

//=======================================================================
InterpolateRBF::InterpolateRBF( unsigned inDim, kernelType aKernel, Array < double > kernParams ) {   //interpolation
	setParams( inDim, aKernel, kernParams );
}

//=======================================================================
void InterpolateRBF::setParams( unsigned inDim, kernelType aKernel, Array < double > kernParams ) {
	inputDim = inDim;
	kernelFunction = aKernel;
	kernelParams = kernParams;

	switch( kernelFunction ) {
		case gaussian: //kernelFuncPtr = &gaussFunc;
			break;
		case linear_spline: //kernelFuncPtr = &linearFunc;
			break;		
		case cubic_spline: //kernelFuncPtr = &cubicFunc;
			break;		
		case multiquadric: //kernelFuncPtr = &multiquadricFunc;
			break;
		case inverse_multiquadric: //kernelFuncPtr = &inverseMultiquadricFunc;
			break;
		case cauchy: //kernelFuncPtr = &cauchyFunc;
			break;
		default: //kernelFuncPtr = &linearFunc;
			break;
	}
}

//=======================================================================
double InterpolateRBF::mse( Array < double > input, Array < double > target ) {
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
		
	Array < double > outputData ( inputData.dim( 0 ), 1 );
	evaluate ( inputData, outputData );
			
	double rms_LB = 0;
	for ( int i = 0; i < ( int ) inputData.dim( 0 ); i++ ) {
		rms_LB += sqr( outputData( i, 0 ) - targetData( i, 0 ) );
	}
	
	return ( rms_LB / ( double ) inputData.dim( 0 ) );
}

//=======================================================================
double InterpolateRBF::mse( Population offsprings ) {
	// the parameters have to be in the first Chromosome!!!
	unsigned inputDimension = dynamic_cast < ChromosomeT < double >& > (offsprings[0][0]).size();
	unsigned numberOfOffsprings = offsprings.size();
	
	Array < double > inputData (numberOfOffsprings,inputDimension);
	Array < double > targetData (numberOfOffsprings, 1);
	
	for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
		for ( unsigned j = 0; j < inputDimension; j++ ) {
			inputData(i,j) = dynamic_cast < ChromosomeT < double >& > (offsprings[i][0] )[j];
		}
		targetData (i,0) = offsprings[i].fitnessValue();
	}  
	
	double error = mse ( inputData, targetData );
	
	return error;
}


//=======================================================================
void InterpolateRBF::train( Array < double > inputData, Array < double > targetData ) {
	unsigned i,j,k;
	
	nBasis = inputData.dim( 0 );
	weightMatrix.resize( nBasis, 1, true );
	designMatrix.resize( inputData.dim( 0 ), nBasis, true );
	baseCentres.resize( nBasis, inputDim, true );
	initArrayToZero( weightMatrix );
	initArrayToZero( designMatrix );
	initArrayToZero( baseCentres );
	
	for( i = 0; i < nBasis; i++ ) {
		for( j = 0; j < inputDim; j++ ) {
			baseCentres( i,j ) = inputData( i,j );		
		}
	}
		
	cout << "Setting design matrix..." << endl;	
	setDesignMatrix( inputData );

	cout<< "Setting weight matrix..." << endl;
	setWeightMatrix( targetData );

	cout<<"Training done..."<<endl;	
}

//=======================================================================
void InterpolateRBF::setWeightMatrix( Array < double >& targetData ) {
	//Array < double > A = innerProduct( transpose( designMatrix ), designMatrix );

	unsigned i, j;

	//A = invert( A );
	//weightMatrix = innerProduct( innerProduct( A, transpose( designMatrix ) ), targetData ); 

        weightMatrix = innerProduct( invert(designMatrix), targetData );

}


//=======================================================================
void InterpolateRBF::setDesignMatrix( Array < double >& inputData ) {
	unsigned i, j;

	double distance;
	
	for ( i = 0; i < designMatrix.dim( 0 ); i++ ) {
		for ( j = 0; j < designMatrix.dim( 1 ); j++ ) {
			//designMatrix( i,j ) = kernelFuncPtr( inputData, i, baseCentres, j, kernelParams );
		
			distance = distance = getDistance( inputData, i, baseCentres, j );

			switch( kernelFunction ) {
				case gaussian:
                                        designMatrix(i,j) = exp( -1.0 * ( distance / kernelParams(0,0) ) * ( distance /
 kernelParams(0,0) ) );
                                        break;
                                case linear_spline:
                                        designMatrix(i,j) = distance;
                                        break;
                                case cubic_spline:
                                        designMatrix(i,j) = sqr( distance ) * distance;
                                        break;
                                case multiquadric:
                                        designMatrix(i,j) = sqrt( sqr( distance )+sqr( kernelParams( 0, 0 ) ) )/kernelParams( 0,0 );
                                        break;
                                case inverse_multiquadric:
                                        designMatrix(i,j) = kernelParams( 0, 0 ) / sqrt( sqr( distance ) + sqr( kernelParams( 0, 0 ) ) );
                                        break;
                                case cauchy:
                                        designMatrix(i,j) =  sqr( kernelParams( 0, 0 ) ) / ( sqr( distance ) + sqr( kernelParams( 0, 0 ) ) );

			}
		}
	}
}


//=======================================================================
void InterpolateRBF::setDesignVector( Array < double > input, unsigned index, Array < double > centres, Array < double >& outputVect ) {
	unsigned i, j;
	
	for ( i = 0; i < input.dim( 0 ); i++ ) {
		outputVect( i, 0 ) = kernelFuncPtr( input, i, centres, index, kernelParams );
	}
}

//=======================================================================
double InterpolateRBF::gaussFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params ) {
	double sigmaParam = params( 0, 0 );
	double distance = getDistance( input, i, centre, j );
	
	return exp( -1.0 * ( distance / sigmaParam ) * ( distance / sigmaParam ) );
}

//=======================================================================
double InterpolateRBF::cubicFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params ) {
	double distance = getDistance( input, i, centre, j );
	return sqr( distance ) * distance;
}

//=======================================================================
double InterpolateRBF::linearFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params) {
	return  getDistance( input, i, centre, j);
}

//=======================================================================
double InterpolateRBF::multiquadricFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params ) {
	double distance = getDistance( input, i, centre, j );
	return sqrt( sqr( distance )+sqr( params( 0, 0 ) ) )/params( 0,0 );
}

//=======================================================================
double InterpolateRBF::inverseMultiquadricFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params ) {
	double distance = getDistance( input, i, centre, j );
	return params( 0, 0 ) / sqrt( sqr( distance ) + sqr( params( 0, 0 ) ) );
}

//=======================================================================
double InterpolateRBF::cauchyFunc( Array < double > input, unsigned i, Array < double > centre, unsigned j, Array < double > params ) {
	double distance = getDistance( input, i, centre, j );
	return sqr( params( 0, 0 ) ) / ( sqr( distance ) + sqr( params( 0, 0 ) ) );
}

//=======================================================================
Array < double > InterpolateRBF::getDesignMatrix( ) {
	return designMatrix;
}

//=======================================================================
Array < double > InterpolateRBF::getWeightMatrix( ) {
	return weightMatrix;
}

//=======================================================================
Array < double > InterpolateRBF::getBaseCentres( ) {
	return baseCentres;
}

//=======================================================================
double InterpolateRBF::getDistance( Array < double >& inputData, unsigned index1, Array <double>& baseCentres, unsigned index2 ) {
	int i = inputData.dim( 1 );
	int j = baseCentres.dim( 1 );
	if ( i != j ) { 
		cout << "error in get distance, not same dimensionality!" << endl; 
		exit(1);
	}

	int k;
	double sum=0;
	for ( k = 0; k < i; k++ ) {
		double diff = inputData( index1, k ) - baseCentres( index2, k );
		sum += diff * diff;
	}

	return sqrt( sum );
}


//======================================================================
void InterpolateRBF::evaluate( Population& offsprings ) {
	unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > ( offsprings[ 0 ][ 0 ] ).size( );
	unsigned numberOfOffsprings = offsprings.size( );
	
	Array < double > NNInput ( numberOfOffsprings, InputDimension );
	Array < double > NNOutput ( numberOfOffsprings, 1 );
	
	for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
		for ( unsigned j = 0; j < InputDimension; j++ ) {
			NNInput( i, j ) = dynamic_cast < ChromosomeT < double >& > ( offsprings[ i ][ 0 ] )[ j ];
		}
	}
	evaluate( NNInput, NNOutput );
	
	for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
		offsprings[ i ].setFitness( NNOutput( i, 0 ) );
	}
}

//=======================================================================
void InterpolateRBF::evaluate( Individual& offspring ) {
	unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > ( offspring[ 0 ] ).size( );
	Array < double > NNInput ( 1, InputDimension );
	Array < double > NNOutput ( 1, 1 );
	
	for ( unsigned j = 0; j < InputDimension; j++ ) {
		NNInput( 0, j ) = dynamic_cast < ChromosomeT < double >& > ( offspring[ 0 ] )[ j ];
	}
	evaluate( NNInput, NNOutput );
	offspring.setFitness( NNOutput( 0, 0 ) );
}


void InterpolateRBF::evaluate(  Array < double > inputData, Array < double >& outputData ) {
        unsigned nToEvaluate = inputData.dim( 0 );
        Array < double > localDesignMatrix( nToEvaluate, nBasis );
        unsigned i, j;


	double distance;
        for ( i = 0; i < nToEvaluate; i++ ) {
                for ( j = 0; j < nBasis; j++ ) {
		         //localDesignMatrix( i, j ) = getDistance( inputData, i, baseCentres, j);		// for linear spline
		         //localDesignMatrix( i, j ) = linearFunc( inputData, i, baseCentres, j, kernParams);
                
			distance = getDistance( inputData, i, baseCentres, j );

			switch( kernelFunction ) {
				case gaussian: 
					localDesignMatrix(i,j) = exp( -1.0 * ( distance / kernelParams(0,0) ) * ( distance / kernelParams(0,0) ) );
					break;
				case linear_spline:
					localDesignMatrix(i,j) = distance;
					break;
				case cubic_spline:
					localDesignMatrix(i,j) = sqr( distance ) * distance;
					break;
				case multiquadric:
					localDesignMatrix(i,j) = sqrt( sqr( distance )+sqr( kernelParams( 0, 0 ) ) )/kernelParams( 0,0 );
					break;
				case inverse_multiquadric:
					localDesignMatrix(i,j) = kernelParams( 0, 0 ) / sqrt( sqr( distance ) + sqr( kernelParams( 0, 0 ) ) );
					break;
				case cauchy:
					localDesignMatrix(i,j) =  sqr( kernelParams( 0, 0 ) ) / ( sqr( distance ) + sqr( kernelParams( 0, 0 ) ) );
					break;
			}

		}
        }

        outputData = innerProduct( localDesignMatrix, weightMatrix );
}



//=======================================================================
void InterpolateRBF::KMean( unsigned nCluster, unsigned maxIter, Array < double > input, Array < double >& centresFound, Array< unsigned >& clustMap ) {
	cout << endl << "Start K-Mean clustering to find " << nCluster << " centres..." << endl;
	unsigned i, j, k;
	Rng::seed( time( NULL ) );
	unsigned nInput = input.dim( 0 );
	unsigned probDimension = input.dim( 1 );
	Array < double > clusterCentres( nCluster, probDimension );
	
	//randomly select one input for each cluster centre
	for( i = 0; i < nCluster; i++ ) {
		Rng::uni.low( i * nInput / nCluster );
		Rng::uni.high( ( i + 1 ) * nInput / nCluster );
		
		unsigned randNumb = ( int )( floor( Rng::uni( ) ) );
		for( j = 0; j < probDimension; j++ ) {
			clusterCentres( i, j ) = input( randNumb, j ); 
		}
	}

	bool converged = false;

	Array < unsigned > clusterMap( nInput, 1 );
	Array < unsigned > clusterMapCopy( nInput, 1 );
	initArrayToZero( clusterMap );
	initArrayToZero( clusterMapCopy );

	for( i = 0; i < maxIter && !converged ; i++ ) {
		clusterMapCopy = clusterMap;
		for( j = 0; j < nInput; j++ ) {
			double minimDistance = getDistance( input, j, clusterCentres, clusterMap( j, 0 ) );
			for ( k = 0; k < nCluster ; k++ ) {
				if( k!=clusterMap( j, 0 ) ) {
					double distance = getDistance( input, j, clusterCentres, k);
					if( distance < minimDistance ) {
						minimDistance = distance;
						clusterMap(j,0) = k;
					}
				}
			}
		}

		//compare map
		converged = true;
		for( j = 0; j < nInput; j++ ) {
			if( clusterMap( j, 0 ) != clusterMapCopy( j, 0 )) {
				converged = false;
				break;
			}
		}
		
		if( !converged ) {
			//find new centres
			initArrayToZero( clusterCentres );
			Array < unsigned > sumArray( nCluster, 1 );
			initArrayToZero( sumArray );
			for ( j = 0; j < nInput; j++ ) {
				unsigned content = clusterMap( j, 0 );
				sumArray( content, 0)++;
				for ( k = 0; k < probDimension; k++ ) {
					clusterCentres( content, k ) += input( j, k );
				}
			}
			for ( j = 0; j < nCluster; j++ ) {
				if( sumArray( j, 0) > 0 ) {
					for ( k = 0; k < probDimension; k++ ) {
						clusterCentres( j, k ) /= sumArray( j, 0 );
					}
				}
			}
		}
		cout << "Running K-Mean iteration " << i << "..." << endl;
	}

	centresFound = clusterCentres;
	clustMap = clusterMap;

	cout << "K-Mean algorithm found the following centres : " << endl;
	cout << "-----------------------------------------------" << endl;
	
	for ( i = 0; i < clusterCentres.dim( 0 ); i++ ) {
		cout << "Centre " << i << " ==> ";
		for ( j = 0; j < clusterCentres.dim( 1 ); j++ ) {
			cout << clusterCentres( i, j ) << "  ";
		}
		cout << endl;
	}
}

