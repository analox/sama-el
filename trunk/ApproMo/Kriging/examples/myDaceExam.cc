#include "Array/ArrayOp.h"
#include "EALib/Population.h"
#include <vector.h>
#include "KrigApprox.h" 
#include "krig.h"
#include "krigify.h"
#include "rngs.h"


using namespace std;


double sphere( vector < double > v ) {
	double sum = 0;
	int i;
	int inDim = v.size( );
	for ( i = 0; i < inDim; i++ ) {
		sum += v[ i ] * v[ i ];
	}
	return sum;
}

double rosenbrock( vector < double > l_array ) {
	int varsize = l_array.size( );
	int w;
	double l_value;
	
	l_value = 0;

	for ( w = 1; w < varsize; w++ ) {
		l_value = l_value + ( 100 * ( ( l_array[ w ] - ( l_array[ w - 1 ] * l_array[ w - 1 ] ) ) * ( l_array[ w ] - ( l_array[ w - 1 ] * l_array[ w - 1 ] ) ) ) + ( ( 1 - l_array[ w - 1 ] ) * ( 1 - l_array[ w - 1 ] ) ) );
	}
	return l_value;
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
                l_Productobj = l_Productobj * cos( l_array[ w ] / ( 4 * sqrt( (double)w + 1.0 ) ) );
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

	l_value = ( ( -20 ) * exp( ( -0.2 ) * sqrt( sqterm / varsize ) ) ) - ( exp( costerm / varsize ) ) + 20 + exp( 1.0 );
	return l_value;
}


int main( void ) {
	unsigned i, j, k;

	char tf[ 20 ];
	unsigned inDim = 2,  nInput = 50;

	double loBound = -4, wholeRange = 8;

	double ( *evaluate ) ( vector < double > );
	
	sprintf( tf, "ackley" );
	if( !strcmp( tf, "ackley" ) ) { evaluate = &ackley; }
	else if( !strcmp( tf, "griewank" ) ) { evaluate = &griewank; }
	else if( !strcmp( tf, "rosenbrock" ) ) { evaluate = &rosenbrock; }	
	else if( !strcmp( tf, "sphere" ) ) { evaluate = &sphere; }
	else if( !strcmp( tf, "rastrigin" ) ) { evaluate = &rastrigin; }  


	Array < double > inArray( nInput, inDim );
	Array < double > targArray( nInput, 1);
	
	//might need to change to lhs later
	srand( time( NULL ) );
	for ( i = 0; i < nInput; i++ ) {
		for ( j = 0; j < inDim; j++ ) 
			inArray( i, j ) = loBound + wholeRange * ( rand( ) % 12134 ) / 12133.0;
	}
		
	for ( i = 0; i < nInput; i++ ) {
		vector < double > arr( inDim );
		for ( j = 0; j < inDim; j++ ) {
			arr[ j ] = inArray( i, j );
		}
		targArray( i, 0 ) = evaluate( arr );
	}

	KrigApprox myApprox1;
	
	cout<<"before train"<<endl;
	

	myApprox1.Train( inArray, targArray );
	
	//myApprox1.Evaluate(TestDat,ResultDat);
	
	cout<<"after train"<<endl;

	cout << endl << "MSE of training input is " << myApprox1.MSE( inArray, targArray ) << endl;


	Array<double> test(1,2); test(0,0) = 0.0; test(0,1) = 0.0;
	Array<double> testOut(1,1);	
	myApprox1.Evaluate(test,testOut);
	cout<<"Evaluate 0,0  = " << testOut(0,0)<<endl;

	return 0;
} 
