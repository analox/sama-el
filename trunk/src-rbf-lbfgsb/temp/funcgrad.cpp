#include<ZNminmax.h>
#include<vector.h>
#include<Array/Array.h>
#include<ap.h>
#include<lbfgsb.h>

extern void lbfgsbminimize(const int& n,
     const int& m,
     ap::real_1d_array& x,
     const double& epsg,
     const double& epsf,
     const double& epsx,
     const int& maxits,
     const ap::integer_1d_array& nbd,
     const ap::real_1d_array& l,
     const ap::real_1d_array& u,
     int& info);
extern unsigned dimension;

void fool( double x, double y, double *z)
{
	*z = x*y + y;
   	return;
}

double small()
{
   	double one,two,z,tsmall;

   	one=1.e0;
   
   	two=2.e0;
   	tsmall=one;
   	do {
      		tsmall=tsmall/two;
      		fool(tsmall,one,&z);
   	} 
	while (z>1.e0);
   
   	return tsmall*two*two;
}

double max(double a, double b){
	return a>b?a:b;
}

//
// objective and gradient function (finite differencing is used here)
//
void funcgrad(const ap::real_1d_array& x, double& f, ap::real_1d_array& g) {

	double rteps = sqrt( small() );

	double delta;
	unsigned i, j;
	
	vector<double> xx(dimension);
	for(i=1; i<=dimension; i++) {
		xx[i-1] = x(i);
	}
	LSEvaluate( xx, f );
	
	vector<double> temp( dimension );
	for( i=1; i<=dimension; i++ ) {
		
      		delta = rteps * max( 1.e0, fabs( x(i) ) );
      		if ( x(i) < 0.e0 ) delta = -delta;
      
      		for( j=1; j<=dimension; j++ ) temp[ j-1 ] = x( j );
		temp[ i-1 ] = x( i ) + delta;
		
	  	LSEvaluate( temp, g( i ) );
	  
      		g( i )=( g( i ) - f )/delta;
		
	}		
}

void findLocalOpt( vector<double> xx, vector<double>& newX, Array<double> boundTR ) {


	unsigned i, j;
	int info;

	ap::real_1d_array startX;
	startX.setbounds(1,dimension);
	ap::integer_1d_array nbd;
	nbd.setbounds(1,dimension);
	ap::real_1d_array lb;
	lb.setbounds(1,dimension);
	ap::real_1d_array ub;
	ub.setbounds(1,dimension);	

	for(i=1; i<=dimension; i++) {
		startX( i ) = xx[ i-1 ];
		nbd( i ) = 2;
		lb( i ) = boundTR( i-1, 0 );
		ub( i ) = boundTR( i-1, 1 );
	}

	lbfgsbminimize( dimension, 5, startX, 0.0, 0.0, 0.0, 100, nbd, lb, ub, info);

	for(i=1; i<=dimension; i++) {
		newX[ i-1 ] = startX( i );
	}

}

