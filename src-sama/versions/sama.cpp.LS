/*
*	EXAMPLE OF SURROGATE-ASSISTED MEMETIC ALGORITHM (using GA+L-BFGS-B+RBF MODEL)
*	
*/

//
// ansi c/c++ header files
//
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
using namespace std;
//
// shark header files
//
#include <EALib/Population.h>
#include <Rng/RNG.h>
#include <Array/Array.h>

//
// Distribution file
//
#include <limits>
#include <Rng/Uniform.h>
#include "tools/Matrix.cpp"
#include "tools/distribution/DistributionFunc.h"
#include "tools/database/ILDatabase.h"
#include "tools/distribution/WeightCal.h"

//
// Appromo header files
//
#include <InterpolateRBF.h>
#include <fifo.h>

// For multiple models
#include <RegrApp.h>
#include <KrigApprox.h>
#include <krig.h>
#include <krigify.h>

//
// L-BFGS-B header file
//
#include <lbfgsb.h>
//
// some definitions
//
#define PRECISION 1.E-5

//
// fitness functions prototypes
//
double sphere( vector<double>& x );
double ackley( vector<double>& x );
double rastrigin( vector<double>& x );
double griewank( vector<double>& x );

//
// functions related to local search
//
double localSearchApproximate( vector<double>& xStart, double& locFit, unsigned trIteration);
void getInitialBoundTR( Array<double>& boundTR, Array<double> input, Array<double> target );
void getInitialDeltaTR( Array<double>& boundTR, vector<double> xStart, vector<double>& delta );
bool atBoundary( vector<double>& xOptLoc, Array<double>& boundTR );
double decidePsi( double rho, Array<double>& boundTR, vector<double>& xOptLoc );
void decideNewX( double rho, vector<double>& xStart, vector<double> xOptLoc, double& fitStart, double fitOptLoc );
void getNewDeltaAndBoundTR( double psi, vector<double>& delta, Array<double>& boundTR, vector<double> xStart );
//void LSEvaluate( vector<double> xStart, double& predicted);

//
// functions related to database
//
void appendDatabase( vector<double>& a, double fitness );
double getDistance( vector<double>& a, vector<double>& b );

//
// function related to learning of evolvability
//
_Matrix<double> metricMeasure(DistributionFunc & func, ILDatabase & db);
void printVec(ostream & output, vector<double>& vec);
//
// global variables
//
double(*evaluate)( vector<double>& );
double objectiveFunction( vector<double>& x );
unsigned dimension=30;
unsigned mcount = 0;
double currBest = HUGE;
vector<double> currBestVec( dimension );
vector<double> lowestBound(dimension);
vector<double> highestBound(dimension);
//double MAX_VALUE = numeric_limits<double>::max();
double MAX_VALUE = 1.7e+308;



FIFO myData;

ILDatabase database;
	

//
// model specific
//
InterpolateRBF aRBF;
//Nghia: add more models here (Neural net?)
RegrApp aRegression(2);
KrigApprox aKrig;
unsigned iModel = 0;

/**
Nghia: Local search on surrogate model function. Number of iterations = number of exact function
calls in trusted-region framework. Need ID for which model it should run on.
*/
int main( int argc, char** argv ) 
{
//double localSearchApproximate( vector<double>& xStart, double& locFit, unsigned trIteration) {
	
	unsigned i, j;
	double locFit = 0;

	cout << "Read input " << endl;
	ifstream fin1("logs/input.dat");
	int w, h;
	fin1 >> h >> w;
	cout << h << ", " << w << endl;
	// set array size
	Array<double> input(h, w);
	for(int i=0; i < h; i++)
	{
		for(int j=0; j < w; j++)
		{
			double tmp;
			fin1 >> tmp;
			// set array i,j
			input[i][j] = tmp;
		}
	}
	fin1.close();
	

	cout << "Read target " << endl;
	ifstream fin2("logs/target.dat");
	fin2 >> h >> w;
	cout << h << ", " << w << endl;
	// set array size
	Array<double> target(h, w);
	for(int i=0; i<h; i++)
	{
		for(int j=0; j<w; j++)
		{
			double tmp;
			fin2 >> tmp;
			// set array i,j
			target[i][j] = tmp;
		}
	}
	fin2.close();

	cout << "Read startPoint " << endl;
	ifstream fin3("logs/startPoint.dat");
	fin3 >> dimension;
	// set array size
	vector<double> xStart(dimension);
	for(int j=0; j<dimension; j++)
	{
		double tmp;
		fin3 >> tmp;
		// set array i,j
		xStart[j] = tmp;
	}

	fin3.close();


	int trIteration = 2;
	
	evaluate = &griewank;
	for(i=0; i<dimension; i++) {
		lowestBound[i] = -600;
		highestBound[i] = 600;
	}
	
	//======== build optimal model for xStart  =======================//	
	// Nghia: build more model here //
	int temp = iModel;
	iModel = 1;
	if (iModel == 0) {		
		cout << endl << "--> BUILDING RBF..." << endl;
		Array<double> kernParams(1,1);				// for completeness sake, no use for rbf linear spline
		kernParams(0,0) = 1;
		aRBF= InterpolateRBF( dimension, linear_spline, kernParams );	// see rbf.h for more options
		aRBF.train(input, target);
		
	} 
	else if (iModel == 1) {
		cout << endl << "--> BUILDING PR..." << endl;
		//unsigned polyOrder = 2;
		//RegrApp tempPoly( polyOrder );
		//tempPoly.Train( input, target );
		//aRegression = tempPoly;
		aRegression.Train( input, target );
	}
	else {
		cout << endl << "--> BUILDING KRIG..." << endl;
		aKrig.Train(input, target);
	}
	iModel = temp;
	cout << "Done building model" << endl;
	//==================== end model specific code ====================//

	double fitStart;
	double fitStartApprox, fitOptLocApprox;
	
	//
	// determine initial tr
	//
	Array<double> boundTR(dimension, 2);
	vector<double> delta(dimension);
	getInitialBoundTR( boundTR, input, target );
	getInitialDeltaTR( boundTR, xStart, delta );
	
	//
	// do trust region here
	//	
	fitStart = objectiveFunction( xStart ); 		
	double predicted;
	double rho;
	vector<double> xOptLoc( dimension );

	// Database: create start point
	ChromosomeT<double> chrom1(xStart);
	Individual startPoint(chrom1);
	startPoint.setFitness(fitStart);
	
	// debug code
	cout << "Start point" << endl; printVec(cout, xStart);	
	// Start optimizing in TRF
	for(i=0; i<trIteration; i++) {	
		cout << "Iteration: " << i << endl;		
		//
		// perform local search
		//
		findLocalOpt( xStart, xOptLoc, boundTR ); 
		// Nghia: exact function value	
		double fitOptLoc = objectiveFunction( xOptLoc );

		cout << "End point" << endl; printVec(cout, xOptLoc);	
		// debug code 
		if (fitOptLoc != fitOptLoc) {
			cout << "DIE here!" << endl;	
			exit(0);
		}
		// Nghia: approximate function at start point
		LSEvaluate( xStart, predicted);
		fitStartApprox = predicted;
	        // Nghia: approximate value at end point	
		LSEvaluate( xOptLoc, predicted);
		fitOptLocApprox = predicted;
	
		//
		// calculate rho,  small number added to prevent divide by zero
		//
		rho = ( fitStart - fitOptLoc ) / ( fitStartApprox - fitOptLocApprox + PRECISION );

		cout << "fitStart = " << fitStart << endl;
                cout << "fitOptLoc = " << fitOptLoc << endl;
                cout << "fitStartApprox = " << fitStartApprox << endl;
                cout << "fitOptLocApprox = " << fitOptLocApprox << endl;
                cout << "RHO = " << rho << endl;

		//
		// determine next tr
		//
		double psi = decidePsi( rho, boundTR, xOptLoc );
	
		decideNewX( rho, xStart, xOptLoc, fitStart, fitOptLoc );
		getNewDeltaAndBoundTR( psi, delta, boundTR, xStart );
	       // debug code
		cout << "New start point" << endl; printVec(cout, xStart);
		for(int k =0;  k < dimension; k++) 
		{
			cout << boundTR(k,0) << " , " << boundTR(k,1) << endl;;
		}	
	}


	// Return fitness to reference variable	
	locFit = fitStart;
	cout << "End here" << endl;
}


//
// model specific
//
void LSEvaluate( vector<double> xStart, double& predicted) {
	unsigned i;
	Array<double> xx( 1, dimension );
	Array<double> result(1,1);
	for(i=0; i<dimension; i++) {
		xx(0,i) = xStart[ i ];
	}	
	int temp = iModel;
	iModel = 1;
	//cout << "Approx Eval - " << iModel << endl;
	if (iModel == 0)
	{
		aRBF.evaluate(xx, result);	
		predicted = result(0,0);

	}
	else if (iModel == 1) {
		aRegression.Evaluate(xx, result);	
		predicted = result(0,0);
	}
	else {
		aKrig.Evaluate(xx, result);	
		predicted = result(0,0);
	}
	iModel = temp;
}

void getNewDeltaAndBoundTR( double psi, vector<double>& delta, Array<double>& boundTR, vector<double> xStart ) {
	unsigned i;
	
	for(i=0; i<dimension; i++) {
		delta[i] *= psi;		
		boundTR(i,0) = max( lowestBound[i] , xStart[i] - delta[i] );
		boundTR(i,1) = min( highestBound[i] , xStart[i] + delta[i] );
	}
}
// Nghia: basically copy xOptLoc to xStart if rho > 0
void decideNewX( double rho, vector<double>& xStart, vector<double> xOptLoc, double& fitStart, double fitOptLoc ) {
	if(rho > 0) {
		unsigned i;
		for( i=0; i<dimension; i++ ) {
			xStart[i] = xOptLoc[i];
		}
		fitStart = fitOptLoc;
	}
}

double decidePsi( double rho, Array<double>& boundTR, vector<double>& xOptLoc ) {

	if(rho<=0.25) return 0.25;
	if(rho<=0.75) return 1;
	else {
		//
		// simplified to only return 2
		//
		return 2;
	
		//
		// original version
		//
		// if( atBoundary( xOptLoc, boundTR ) ) return 2;
		// else return 1
	}
}

bool atBoundary( vector<double>& xOptLoc, Array<double>& boundTR ) {

	unsigned i;

	for( i=0; i<dimension; i++ ) {
		if(  
			( xOptLoc[ i ] > boundTR(i,0) - PRECISION  &&  xOptLoc[ i ] < boundTR( i, 0 )+PRECISION )  ||  
			( xOptLoc[ i ] > boundTR(i,1) - PRECISION  &&  xOptLoc[ i ] < boundTR( i, 1 )+PRECISION )  
		) {
			return true;
		}
	}

	return false;
}


void getInitialBoundTR( Array<double>& boundTR, Array<double> input, Array<double> target ) {
	
	unsigned i, j, size = input.dim(0);
	
	for(i=0; i<dimension; i++) {
		double min = input(0,i), max = input(0,i);
		
		for(j=1; j<size; j++) {
			if(input(j,i)<min) min = input(j,i);
			if(input(j,i)>max) max = input(j,i);
		}	

		boundTR(i,0) = min;
		boundTR(i,1) = max;
	}
}

void getInitialDeltaTR( Array<double>& boundTR, vector<double> xStart, vector<double>& delta ) {
	
	unsigned i;
	double upper, lower;
	for(i=0; i<dimension; i++) {
		upper = boundTR(i,1) - xStart[i];
		lower = xStart[i] - boundTR(i,0);
		
		delta[i] = min(upper,lower);
	}
	for(i=0; i<dimension; i++) {
		boundTR(i,0) = xStart[i] - delta[i];
		boundTR(i,1) = xStart[i] + delta[i];
	}
	
}

double objectiveFunction( vector<double>& x ) {
	unsigned i;
	double tempFit = evaluate( x );

	mcount++;
        if(tempFit<currBest) {
        	currBest = tempFit;
                for(i=0 ;i<dimension; i++) {
                	currBestVec[ i ] =  x[ i ];
                        
                }
	}

	ofstream bestFile("logs/bestFile.dat", ios::app);
	ofstream bestVecFile("logs/bestVecFile.dat", ios::app);

        bestFile << mcount << "\t" << currBest << endl;
        bestVecFile << mcount << "\t";
        for( i=0; i<dimension; i++ ) {
        	bestVecFile << currBestVec[ i ] << "\t";
        }
        bestVecFile << endl;
	
	bestFile.close();
	bestVecFile.close();

	//appendDatabase( x, tempFit );

	return tempFit;
}

void appendDatabase( vector<double>& a, double fitness ) {

	vector<double> target(1);

	target[0] = fitness;

	bool equalDatabase = false;

	unsigned i, j,
		currDataSize = myData.getCurrentSize();
	

	Array<double> tempInp = myData.getInput();
	for(i=0; i<currDataSize; i++) {
		vector<double> b(dimension);
		for(j=0; j<dimension; j++) {
			b[j] = tempInp(i,j);
		}
		if( getDistance( a, b ) < PRECISION ) {
			equalDatabase = true;
			break;
		}
	}

	if( !equalDatabase ) myData.append( a, target );

}

double getDistance( vector<double>& a, vector<double>& b ) {

	//assert( a.size() == b.size() );
	unsigned i, j;

	double sum = 0;
	for(i=0; i<a.size(); i++) {
		sum += ( a[i]-b[i] )*( a[i]-b[i] );

	}

	return sqrt(sum);

}

double sphere( vector<double>& x ) {

	unsigned size = x.size();
	unsigned i;
	double result=0;
	
	for(i=0; i<size; i++) {
		result += x[i]*x[i];
	
	}
	
	return result;
	
}


double ackley( vector< double >& x )
{
    const double A = 20.;
    const double B = 0.2;
    const double C = znPi2C;

    unsigned i, n;
    double   a, b;

    for( a = b = 0., i = 0, n = x.size( ); i < n; ++i ) {
        a += x[ i ] * x[ i ];
	b += cos( C * x[ i ] );
    }

    return -A * exp( -B * sqrt( a / n ) ) - exp( b / n ) + A + znEC;
}



double rastrigin( vector < double >& x ) {
        int varsize = x.size( );
        int w, i, j, a = 2;
        float c;
        double result, unNormVar;

        result = 0;
        for ( w = 0; w < varsize; w++ ) {
                result = result + ( ( ( x[ w ] ) * ( x[ w ] ) ) - 10.0 * cos( 2 * M_PI * ( x[ w ] ) ) );
        }
        result = ( result + 10.0 * varsize );
        return result;
}


double griewank(vector<double>& l_array)
{
        int varsize = l_array.size();

        double l_value, l_Sumobj, l_Productobj;
        int w;

        l_value = 0;
        l_Sumobj = 0;
        l_Productobj = 1;


        for (w = 0; w < varsize; w++) {

                l_Sumobj = l_Sumobj + ((l_array[w] * l_array[w])/4000);
                l_Productobj = l_Productobj * cos(l_array[w]/sqrt(w+1.0));
        }
        l_value = (l_Sumobj + 1 - l_Productobj) ;

        return l_value;
}

/**
Evolvability Learning on stochastic variation/distribution & DB of individual learning
[IL Index] [0] -> expected resulted fitness
[IL Index] [1] -> expected cost required
*/
_Matrix<double> metricMeasure(DistributionFunc & func, ILDatabase & db)
{
	_Matrix<double> val(db.nMethods, 2);
	
	// for each individidual learning method
	for (int i = 0; i < db.nMethods; i++)
	{
		list<DBEntry> listDB = db.getListForMethod(i);
		// check if list empty
		if (listDB.empty())
		{
			cout << "List DB is empty " << endl;
			val.at(i,0) = MAX_VALUE; // some big value
			val.at(i,1) = 0;	
		}
		else
		{
			_Matrix<double> array = db.getListArrayForMethod(i);
			vector<double> w = WeightCal::weightCal(&func, array);

			double sum = 0;
			//cout << "Weights: "<< w.size() << endl << array << endl;
			cout << "Weights: "<< w.size() << endl;
			for (int k = 0; k < w.size(); k++) 
			{
				cout << w[k] <<  ", ";
				sum += w[k];
			}
			cout << " | " << sum << endl;
				
			// measure expect fitness obtained and cost
			// if no relevant data obtained
			if (sum <= 1E-9) 
			{
//			 	System.out.println("No relevance samples");
				val.at(i, 0) = MAX_VALUE;
				val.at(i, 1) = 0;
			}
			else 
			{ // measure expected fitness
//				System.out.println("Have relevance samples");
				val.at(i, 0) = 0;
				val.at(i, 1) = 0;

				list<DBEntry>::iterator lIter = listDB.begin();	
				
				double j = 0;
				for (lIter=listDB.begin(); lIter != listDB.end(); ++lIter)
				{
					double fitness = ((*lIter).getEndPoint()).fitnessValue();
					// fitness
					val.at(i, 0) += w[j] * fitness;
					// cost
					//cout << "Cost " << j << ": " << ((*lIter).getCost()) << endl;
					val.at(i, 1) += w[j] * ((*lIter).getCost());

					j++;
				}
				cout << "FI : " << val.at(i, 0) << ", C : " << val.at(i, 1) << endl;
				
			}
		}
	}
	return val;
}

void printVec(ostream & output, vector<double>& vec)
{

	for(int i=0; i<vec.size(); i++)
	 	output<<vec[i] << " ";
	output << endl;

}
// ==================================================   eof  ============================================ //
int _main()
{
// UNIT TEST
	// Create database
	ILDatabase db(2);
	cout << "Double max: " << MAX_VALUE << endl;
	double _a[] = {2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
	vector<double> a(_a, _a + 12);
	ChromosomeT<double> chrom1(a);
	Individual startPoint(chrom1);
	startPoint.setFitness(11.0);

	double _b[] = {	1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3};
	vector<double> b(_b, _b + 12);
	ChromosomeT<double> chrom2(b);
	Individual endPoint(chrom2);
	endPoint.setFitness(10);

	DBEntry entry1(0, startPoint, endPoint, 0.123, 10);
	
	double _a2[] = {2, 4, 2, 3, 2, 1, 2, 3, 2, 3, 2, 3};
        vector<double> a2(_a2, _a2 + 12);
        ChromosomeT<double> chrom12(a2);
        Individual startPoint2(chrom12);
	startPoint2.setFitness(20);

        double _b2[] = { 1, 3, 1, 3, 1, 3, 0, 3, 0, 3, 1, 3};
        vector<double> b2(_b2, _b2 + 12);
        ChromosomeT<double> chrom22(b2);
        Individual endPoint2(chrom22);
	endPoint2.setFitness(15);

	DBEntry entry2(1, startPoint2, endPoint2, 0.5, 15);

	db.add(entry1, 0); db.add(entry2,0);
	db.prepareListArray();
	cout << db.toString();
	
	// Create distribution
	double _mean[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	vector<double> mean(_mean, _mean+12);
	double _sigma[] = {10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0,
     			   0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0,
			   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10
							};
	_Matrix<double> sigMat(_sigma, 12, 12);
	//cout << "Create Gaussian distribution " << endl;
	GaussianDistFunc dist(mean, sigMat);
	fflush(stdout);
	
	_Matrix<double> mat = metricMeasure( dist, db);
	cout << "Weighted sample: "<< endl << mat << endl;
	return 0;	
}
