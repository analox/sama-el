#include <iostream.h>
#include <time.h>
#include <fstream.h>
#include <vector.h>
#include "InterpolateRBF.h"
#include "LinAlg/linalg.h"
#include "ArrayOp.h"


/*

To run:

Just do interpolation, one step of matrix inverse.
./InterpolateRBFTest [testFunction] [numbOfInputSamples] [kernelFunction] [kernelWidth]
e.g.: ./InterpolateRBFTest ackley 50 gaussian 1.0

testFunction : { ackley, sphere, griewank, rastrigin }

kernelFunction : {  gaussian, linear_spline, cubic_spline,  multiquadric, inverse_multiquadric, cauchy }


*/

double sphere( vector < double > v ) {
	double sum = 0;
	int i;
	int inDim = v.size( );
	for ( i = 0; i < inDim; i++ ) {
		sum += v[ i ] * v[ i ];
	}
	return sum;
}

double rastrigin( vector < double > l_array ) {
	int varsize = l_array.size( );
	int w, i, j, a = 2; 
	float c;
	double l_value, unNormVar;
  	
	l_value = 0;
	for ( w = 0; w < varsize; w++ ) {
		l_value = l_value + ( ( ( l_array[ w ] ) * ( l_array[ w ] ) ) - 10.0 * cos( 2 * M_PI * ( l_array[ w ] ) ) );
	}
	l_value = ( l_value + 10.0 * varsize );
	return l_value;
}

double griewank(  vector< double > l_array ) {
        double l_value, l_Sumobj, l_Productobj;
        int w;
        int varsize = l_array.size( );

        l_value = 0;
        l_Sumobj = 0;
        l_Productobj = 1;

        for ( w = 0; w < varsize; w++ ) {
                l_Sumobj = l_Sumobj + ( ( l_array[ w ] * l_array[ w ] ) / 4000 );
                l_Productobj = l_Productobj * cos( l_array[ w ] / ( 4 * sqrt( w + 1 ) ) );
        }
        l_value = ( l_Sumobj + 1 - l_Productobj );
	return l_value;
}

double ackley( vector < double > l_array ) {
	int w;
	int varsize = l_array.size( );
	double l_value, sqterm, costerm, a, b;
	l_value = 0;
	sqterm = 0;
	costerm = 0;
	
	for ( w = 0; w < varsize; w++ ) {
		sqterm = sqterm + ( l_array[ w ] * l_array[ w ] );
		costerm = costerm + cos( 2 * M_PI * l_array[ w ] );
	}

	l_value = ( ( -20 ) * exp( ( -0.2 ) * sqrt( sqterm / varsize ) ) ) - ( exp( costerm / varsize ) ) + 20 + exp( 1 );
	return l_value;
}

int main( int argc, char** argv ) {
	int i, j, k;

	if( argc<5 ) {
		cout << "Usage: ./InterpolateRBFTest [testFunction] [numbOfInputSamples] [kernelFunction] [kernelWidth] e.g.: ./InterpolateRBFTest ackley 50 gaussian 1.0" << endl;
		exit(1);
	}


	char lm[ 50 ]; 
	char tf[ 20 ];
	unsigned inDim = 2,  nInput;
	char kernFunc[ 40 ]; kernelType kernFuncU;
	double loBound = -4, wholeRange = 8, kernWidth;

	double ( *evaluate ) ( vector < double > );

	sprintf( tf, "%s", argv[ 1 ] );
	if( !strcmp( tf, "ackley" ) ) { evaluate = &ackley; }
	else if( !strcmp( tf, "griewank" ) ) { evaluate = &griewank; }
	else if( !strcmp( tf, "sphere" ) ) { evaluate = &sphere; }
	else if( !strcmp( tf, "rastrigin" ) ) { evaluate = &rastrigin; }

	nInput = atoi( argv[ 2 ] );
	
	sprintf( kernFunc, "%s", argv[ 3 ] );
	if( !strcmp( kernFunc, "gaussian" ) ) { kernFuncU = gaussian; }
	else if( !strcmp( kernFunc, "linear_spline" ) ) { kernFuncU = linear_spline; }
	else if( !strcmp( kernFunc, "cubic_spline" ) ) { kernFuncU = cubic_spline; }
	else if( !strcmp( kernFunc, "multiquadric" ) ) { kernFuncU = multiquadric; }
	else if( !strcmp( kernFunc, "inverse_multiquadric" ) ) { kernFuncU = inverse_multiquadric; }
	else if( !strcmp( kernFunc, "cauchy" ) ) { kernFuncU = cauchy; }

	kernWidth = atof( argv[ 4 ] );
	
	Array < double > inArray( nInput, inDim );
	Array < double > targArray( nInput, 1);
	
	Rng::seed(time(NULL));
	Rng::uni.low(0.0);
	Rng::uni.high(1.0);
	for ( i = 0; i < nInput; i++ ) {
		for ( j = 0; j < inDim; j++ ) 
			inArray( i, j ) = loBound + wholeRange * Rng::uni();
	}
		
	for ( i = 0; i < nInput; i++ ) {
		vector < double > arr( inDim );
		for ( j = 0; j < inDim; j++ ) {
			arr[ j ] = inArray( i, j );
		}
		targArray( i, 0 ) = evaluate( arr );
	}

	InterpolateRBF aRBF;
	
	Array < double > kernParams( 1, 1 );
	kernParams( 0, 0 ) = kernWidth;
		
	aRBF = InterpolateRBF( inDim, kernFuncU, kernParams );

	aRBF.train( inArray, targArray );

	cout << endl << "MSE of training input is " << aRBF.mse( inArray, targArray ) << endl;;

	return 0;
}
