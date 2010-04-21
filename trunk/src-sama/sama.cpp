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
#include <string>
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
#include "objFunctions.cpp"
#include "tools/utilities.h"

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
//unsigned dimension=30;
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
unsigned iModel;
//
// EL specified
//
int nRep = 1; // number of stochastic variations
int nIL = 3;  // number of surrogate models considerred

//
// Extra variable
//
long long startT = 0;

int main( int argc, char** argv ) {

	srand( time( NULL ) );
	Rng::seed(time(NULL));
	
	unsigned i, j;
	
	//
	// default settings are listed below for minimizing 30D griewank(range -600 to 600), 8000 evaluations, 1 elitist, popsize of 100, mutateProb 0.01, crossProb 0.90
	// Local search performed on all offsprings, with intensity of 3 exact evaluations (reflected from the number of trust region iteration);
	//
	evaluate = &griewank;

	bool rotated, shifted;
	unsigned popSize = 100,
		 nElitist = 1,
		 maxCount = 8000;
		 
	double crossProb = 0.90;
	double mutateProb = 0.01;
	double memeticProb = 1.0;
	unsigned trIteration = 3;	
	unsigned databaseBuildingPhase = 1000;		// in terms of function evaluation mcount
		 
	//Population parent(popSize, ChromosomeT<double>(dimension));

	for(i=0; i<dimension; i++) {
		lowestBound[i] = -600;
		highestBound[i] = 600;
	}
	

	//
	// read parameter file, by default the filename is param.dat , 
	// otherwise put it as the first argument of the main executable file, 
	// e.g.  ./ga paramFileName
	//
	char paramFileName[30];
	if( argc>1 ) {
		sprintf( paramFileName, "%s", argv[1] ); 
	}
	else sprintf( paramFileName, "param.dat" );
	bool useParamFile;

	//
	// read param file
	//
	char dummy[20];
	ifstream paramFile( paramFileName );
		

	paramFile >> dummy >> useParamFile;
	


	if( useParamFile ) {
		char funcName[20];
		bool sameBoundFlag;
		paramFile >> dummy >> dimension
			  >> dummy >> popSize
			  >> dummy >> maxCount
			  >> dummy >> mutateProb
			  >> dummy >> crossProb
			  >> dummy >> memeticProb
			  >> dummy >> trIteration
			  >> dummy >> databaseBuildingPhase
			  >> dummy >> funcName
			  >> dummy >> sameBoundFlag;

		if     ( !strcmp(funcName, "griewank") ) evaluate = &griewank;
		else if( !strcmp(funcName, "rastrigin") ) evaluate = &rastrigin;
		else if( !strcmp(funcName, "ackley") ) evaluate = &ackley;
		else if( !strcmp(funcName, "sphere") ) evaluate = &sphere;
		else if( !strcmp(funcName, "bump") ) evaluate = &bump;
		else if( !strcmp(funcName, "elliptic") ) evaluate = &elliptic;
		else if( !strcmp(funcName, "griewankrosenbrock") ) evaluate = &griewank_rosenbrock;
		else if( !strcmp(funcName, "rosenbrock") ) evaluate = &rosenbrock;
		else if( !strcmp(funcName, "scaffer") ) evaluate = &scaffer;
		else if( !strcmp(funcName, "schwefel102") ) evaluate = &schwefel102;
		else if( !strcmp(funcName, "step") ) evaluate = &step;
		else if( !strcmp(funcName, "weierstrass") ) evaluate = &weierstrass;
		else if( !strcmp(funcName, "schwefel206") ) evaluate = &schwefel206;
		else if( !strcmp(funcName, "schwefel213") ) evaluate = &schwefel213;
		else {
			cout << "Fitness function undefined !" << endl;
			exit( 1 );
		}

		lowestBound.resize( dimension );
		highestBound.resize( dimension );
		if( sameBoundFlag ) {
			paramFile >> dummy >> lowestBound[ 0 ]
				  >> dummy >> highestBound[ 0 ];

			for( i=1; i<dimension; i++ ) {
				lowestBound[ i ] = lowestBound[ 0 ];
				highestBound[ i ] = highestBound[ 0 ]; 
			}
		}	
		else {
			paramFile >> dummy;
			for( i=0; i<dimension; i++ ) {
				paramFile >> lowestBound[ i ];
			}
			paramFile >> dummy;
			for( i=0; i<dimension; i++ ) {
				paramFile >> highestBound[ i ];
			}
		}
		
		initTransform();
 		paramFile >> dummy >> rotated;
		paramFile >> dummy >> shifted;
		cout << "Shifted: " << shifted << " Rotated: " << rotated << endl;
		char rotFile[100], shiftFile[100];
		if (shifted) {
			sprintf(shiftFile, "input_data/%s_func_data.txt", funcName);
			loadTranslation(shiftFile, dimension);
		}
		if (rotated) {
			sprintf(rotFile, "input_data/%s_M_D%d.txt", funcName, dimension);
			loadRotation(rotFile, dimension);
		}

	}
	paramFile.close();
	

	//
	// bestFile.dat logs the best fitness found so far, bestVecFile.dat logs the best design vector found so far
	//
	ofstream bestFile("logs/bestFile.dat",ios::trunc);
	ofstream bestVecFile("logs/bestVecFile.dat",ios::trunc);
	bestFile.close();
	bestVecFile.close();
	//
	// store frequency of profile usage
	ofstream freqFile("logs/freqFile.dat",ios::trunc);

	myData = FIFO( dimension*100, dimension, 1 );

	ChromosomeT<double> tempLoBound(lowestBound);
        ChromosomeT<double> tempHiBound(highestBound);
	
	currBestVec.resize(dimension); 
        
	Population parent(popSize, ChromosomeT<double>(dimension));

	for(i=0; i<parent.size(); i++) {
                dynamic_cast<ChromosomeT<double>&>(parent[i][0]).initialize(tempLoBound,tempHiBound);
		
		double tempFit =  objectiveFunction( dynamic_cast<ChromosomeT<double>&>(parent[i][0]) );
                parent[i].setFitness( tempFit );

        }

	Population offspring(popSize, ChromosomeT<double>(dimension));
	
	parent.setMinimize();
	offspring.setMinimize();

	// NOTE: need to create a copy of offspring for crossover i & (i+1)
	//Population offspringBuffer(popSize, ChromosomeT<double>(dimension));
	//offspringBuffer.setMinimize();
	
	unsigned t = 1;
	double bestFitness =  parent.best().fitnessValue();
        cout << "Best fitness at generation " << t << " = " << bestFitness << endl;

	/* NGHIA: Global variable */
	database = ILDatabase(nIL);
	int nIndivs = popSize;
	// prepare for Uniform distribution
	vector<double> uBound(dimension);
	vector<double> lBound(dimension);
	UniformDistFunc distUniform(lBound, uBound);
	// prepare for Gaussian distribution
	// Create distribution to pre-calculate inverse matrix 
	double mutVar = 1.0;
	vector<double> mean(dimension);
	for (i = 0; i < dimension; i++)
		mean[i] = 0;
	_Matrix<double> sigMat(dimension, dimension);
	for (i = 0; i < dimension; i++)
		for (j =0; j < dimension; j++)
			sigMat.at(i,j) = ((i==j) ? mutVar : 0);
	GaussianDistFunc distGaussian(mean, sigMat);
	// Matrix to store evolvability of (variation, model)
	vector< _Matrix<double> > evol_metric(nRep);
	// As calculation is only required once for crossover
	bool isCalculated;	
	_Matrix<double> v1(nIL,2);
	evol_metric.push_back(v1);
	_Matrix<double> v2(nIL,2);
	evol_metric.push_back(v2);
	// _Matrix to store frequency of usage
	_Matrix<double> freqF(nRep, nIL);
	// Probability for building database phase
	double crossProbInit = 0.9;
	double mutateProbInit = 0.01;	
	//
	// optimization starts here !!!
	//
	while( mcount < maxCount )  {
	
		/* ****************** database building phase ********************/
		if (mcount < databaseBuildingPhase) 
		{
			cout << "Database building phase..." << endl;
			//
			// decide which allele to be crossed
			//
			vector<bool> pos(dimension);

			//
			// do uniform crossover
			//
			offspring = parent;
			for( i = 0; i < offspring.size( )-1; i += 2 ) {
	 			for(j=0; j<dimension; j++) {
		               		pos[j] = (rand()%11)>5?0:1;
		        	}
			 	if( Rng::coinToss( crossProbInit ) ) {
		        		offspring[ i ][ 0 ].crossoverUniform( offspring[ i+1 ][ 0 ],pos);
				}
			}		

			//
			// do uniform mutation
			//
			for(i=0; i<offspring.size(); i++) {
		               dynamic_cast< ChromosomeT< double >& >( offspring[ i ][ 0 ] ).mutateUniform(tempLoBound, tempHiBound, mutateProbInit);
			}

			for(i=0; i<offspring.size(); i++) {
				double tempFit;
				tempFit = objectiveFunction( dynamic_cast<vector<double>&>( offspring[ i ][ 0 ])  );
				offspring[i].setFitness( tempFit );
			}
		}
		else
		{
			/* ============ Learning of Evolvability ==========================*/
			// Check if performance is good to increase factor
	
			int ilIndex, vIndex;
		
			// Reset frequency counter
			for (vIndex = 0; vIndex < nRep; vIndex++)
				for (ilIndex = 0; ilIndex < nIL; ilIndex++)
					freqF.at(vIndex, ilIndex) = 0;
			// Estimate density function for crossover
			for (i = 0; i < dimension; i++)	{
				uBound[i] = -MAX_VALUE;
				lBound[i] = MAX_VALUE;
			}
				
			for (i = 0; i < parent.size(); i++)
			{
				vector <double> geneVals = dynamic_cast<vector<double>&>( parent[ i ][ 0 ]);
				for (j = 0; j < dimension; j++)
				{
					if (uBound[j] <= geneVals[j]) uBound[j] = geneVals[j];
					if (lBound[j] >= geneVals[j]) lBound[j] = geneVals[j];
				}
			}
			distUniform.setBounds(lBound, uBound);

			ChromosomeT<double> tempUBound(uBound);
		        ChromosomeT<double> tempLBound(lBound);	

			// Set calculation for crossover to FALSE
			isCalculated = false;
			// For each individual in offspring pool
			for (i = 0; i < parent.size(); i++)
			{
				// Prepare database
				cout << "Prepare database" << endl;
				database.prepareListArray();
				cout << "Database: " << endl << database.toString() << endl;
				// For each stochastic distribution: 0 for crossover, 1 for mutation
				for (vIndex = 0; vIndex < nRep; vIndex++)
				{
					// Estimate mutation density distribution of xStart
					if (vIndex == 1) 
					{
						// Calculate density function once
						vector<double> meanDist = 
								dynamic_cast< vector<double> &>( parent[i][0]);
						distGaussian.setMean(meanDist);
						// Calculate the fitness/cost matrix for all models
						evol_metric[vIndex] = metricMeasure(distGaussian, database);
						cout << "Finish building mutation distribution" << endl;
					}
					else
					{
						// check if FI/C for crossover has been calculated in the generation
						//if (!isCalculated) 
						//{
							evol_metric[vIndex] = metricMeasure(distUniform, database);
						//	isCalculated = true;
							cout << "Finish building uniform distribution" << endl;
						//}
					}
				
				}
				// Display evol_metric
				cout << "Computed evolvability metric" << endl;
				for (vIndex = 0; vIndex < nRep; vIndex++)
					cout << evol_metric[vIndex] ;
				cout << endl;
				// Compute Evolvability Metric & Select best variation + model to use
				double best = -MAX_VALUE;
				bool randomMode = true;
				int bestV, bestIL; // best variation + model to use
				for (vIndex = 0; vIndex < nRep; vIndex++)
					for (ilIndex = 0; ilIndex < nIL; ilIndex++)
					{
						double parentFitness = parent[i].getFitness();
						double fi = parentFitness - evol_metric[vIndex].at(ilIndex, 0);
						double cost = evol_metric[vIndex].at(ilIndex, 1);
						cout << "EFI : " << fi << ", EC : " << cost << endl;
						if ((fi >=0) && cost > 0)
						{
							randomMode = false;
							double test = fi/ cost;	
							if (test >= best) 
							{
								bestV = vIndex;
								bestIL = ilIndex;
								best = test;
							}
						}
					}
				cout << "Finish choosing best profile" << endl;
				// check if random mode is required
				if (randomMode) {
					// Choose nRep
					double threshold = 1.0/ nRep;
					double rand = Rng::uni() ;
					for (int v = 0;  v < nRep; v++)
						if ((v+1) * threshold > rand)  
						{
							bestV = v;
							break;
						}
		    				
					// Choose model
					threshold = 1.0/ nIL;
					rand = Rng::uni() ;
					for (int v = 0;  v < nIL; v++)
						if ((v+1) * threshold > rand) 
						{
							bestIL = v;
							break;
						}
					cout << i << ".Random mode" << "( " << bestV << ", " << bestIL << ")" << endl;
						
				}
				else 
					cout << i << ".Best mode" << "( " << bestV << ", " << bestIL << ")" << endl;
				// increase counter
				freqF.at(bestV, bestIL) = freqF.at(bestV, bestIL) + 1;
				// If (crossover + IL)
				if (bestV == 0)
				{ 
					iModel = bestIL;
					cout << "Starting crossover + IL-" << iModel << endl << "--------------------------" << endl;
					//
					// decide which allele to be crossed
					//
					vector<bool> pos(dimension);

					//
					// do uniform crossover
					//
					// NOTE: need to create a copy of offspring for crossover i & (i+1) 
					// aka "offspringBuffer"
					offspring[i] = parent[i];
					for(j=0; j<dimension; j++) 
					{
		       	       			pos[j] = (rand()%11)>5?0:1;
		        		}
				 	if( Rng::coinToss( crossProb ) ) 
					{
	       		        		offspring[ i ][ 0 ].crossoverUniform( parent[ (i+1) % popSize ][ 0 ],pos);
					}
						
					// do extra mutation of only crossover is considerred
					if (nRep == 1)
					{
						
						cout << "Add uniform mutation" << endl;
						dynamic_cast< ChromosomeT< double >& >( offspring[ i ][ 0 ] ).mutateUniform(tempLBound, tempUBound, mutateProbInit);
						//dynamic_cast< ChromosomeT< double >& >( offspring[ i ][ 0 ] ).mutateUniform(tempLoBound, tempHiBound, mutateProbInit);
					}
					// do TRF learning
				
					double tempFit; // exact function evaluation

					if( mcount > databaseBuildingPhase && Rng::coinToss( memeticProb ) ) {
						localSearchApproximate( dynamic_cast<vector<double>&>( offspring[ i ][ 0 ] ), tempFit, trIteration);
					}
					else 
						tempFit = objectiveFunction( dynamic_cast<vector<double>&>( offspring[ i ][ 0 ])  );
			
					offspring[i].setFitness( tempFit );
					cout << "Finish crossover + IL-" << iModel << endl << "--------------------------" << endl;
				}
				// If (mutation + IL)
				else
				{
					iModel = bestIL;
					cout << "Starting mutation  + IL-" << iModel << endl << "--------------------------" << endl;
					//
					// do uniform mutation
					//
					offspring[i] = parent[i];
				        //dynamic_cast< ChromosomeT< double >& >( offspring[ i ][ 0 ] ).mutateUniform(tempLoBound, tempHiBound, mutateProb);
				        if ( Rng::coinToss(mutateProb) )
						dynamic_cast< ChromosomeT< double >& >( offspring[ i ][ 0 ] ).mutateNormal(mutVar);
					// do TRF learning
	
					double tempFit;
				
					if( mcount > databaseBuildingPhase && Rng::coinToss( memeticProb ) ) {
						localSearchApproximate( dynamic_cast<vector<double>&>( offspring[ i ][ 0 ] ), tempFit, trIteration);
					}
					else 
						tempFit = objectiveFunction( dynamic_cast<vector<double>&>( offspring[ i ][ 0 ])  );
			
					offspring[i].setFitness( tempFit );
					cout << "Finish mutation  + IL-" << iModel << endl << "--------------------------" << endl;
				}
			}
			// scale counter to frequency
			for (vIndex = 0; vIndex < nRep; vIndex++)
				for (ilIndex = 0; ilIndex < nIL; ilIndex++)
				{
					freqF.at(vIndex, ilIndex) = freqF.at(vIndex, ilIndex)/(1.0 * popSize);
					freqFile << freqF.at(vIndex, ilIndex) << " ";
				}
			freqFile << endl;
		}
		// Merge back to offspring for selection stage
		/*==================== END - NGHIA ============================= */

		//
	        //scaling and selection   
		//
		const unsigned Omega       = 5;
		vector< double > window( Omega );
		offspring.linearDynamicScaling( window, t );
	        parent.selectProportional( offspring, nElitist );
		cout << "Finish selection & merging" << endl;	
	
		bestFitness = parent.best().fitnessValue();
		cout << "Best fitness at generation " << t+1 << " = " << bestFitness << endl;
		
		t++;
	}	

	bestFile.close();
	bestVecFile.close();
	freqFile.close();
	return 0;
}
/**
Nghia: Local search on surrogate model function. Number of iterations = number of exact function
calls in trusted-region framework. Need ID for which model it should run on.
*/
double localSearchApproximate( vector<double>& xStart, double& locFit, unsigned trIteration) {
	
	unsigned i, j;
	
	Array<double> input = myData.getInput();
	Array<double> target = myData.getTarget();
		
	
	//
	// use k-nearest
	//
	//
	unsigned currSize = myData.getCurrentSize();
	unsigned kNearest = (unsigned) ( (dimension+1)*(dimension+2)/2 );
        // Nghia: Creat new population of Individual of 2 chromosomes. Use Population structure to sort distance.	
	Population pop( currSize, ChromosomeT<double>(dimension), ChromosomeT<double>(1) );
	pop.setMinimize();
	
		
	//
	// utilize the population structure to find nearest points
	//
	vector<double> tempVec(dimension);
	vector<double> tempVec2(dimension);
        for( j=0; j<dimension; j++ ) {
               tempVec2[j] = xStart[j];
        }
	for( i=0; i<pop.size(); i++ ) {
		for( j=0; j<dimension; j++ ) {
			tempVec[j] = input(i,j);
		}
		
		ChromosomeT<double> tempChrom(tempVec);
		dynamic_cast< ChromosomeT<double>& > (pop[i][0]).initialize( tempChrom, tempChrom );
		dynamic_cast< ChromosomeT<double>& > (pop[i][1]).initialize( target(i,0), target(i,0) );
		
		vector<double> tempVec2(dimension);
		for( j=0; j<dimension; j++ ) {
			tempVec2[j] = xStart[j];
		}
		pop[i].setFitness( getDistance( tempVec, tempVec2  ) );
	}
	// Nghia: sort population by distance. Only choose k-nearest
	pop.sort();
	input.resize( kNearest, dimension, false );
	target.resize( kNearest, 1, false );
	for( i=0; i<kNearest; i++ ) {
		for( j=0; j<dimension; j++ ) {
			input( i, j ) = dynamic_cast<vector<double>&>( pop[i][0] )[ j ];
		}
		target( i, 0 ) = dynamic_cast<vector<double>&>( pop[i][1] )[ 0 ];
	}
	
	
	//======== build optimal model for xStart  =======================//	
	// Nghia: build more model here //
	int temp = iModel;
	if (nIL == 1) iModel = 1;
	if (iModel == 0) {		
		cout << endl << "--> BUILDING RBF..." << endl;
		Array<double> kernParams(1,1);				// for completeness sake, no use for rbf linear spline
		kernParams(0,0) = 1;
		aRBF= InterpolateRBF( dimension, linear_spline, kernParams );	// see rbf.h for more options

		startT = getTime(); // Get start time
		aRBF.train(input, target);
		cout << "RBF Training Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;
	} 
	else if (iModel == 1) {
		cout << endl << "--> BUILDING PR..." << endl;
		//unsigned polyOrder = 2;
		//RegrApp tempPoly( polyOrder );
		//tempPoly.Train( input, target );
		//aRegression = tempPoly;

		startT = getTime(); // Get start time
		aRegression.Train( input, target );
		cout << "PR Training Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;	
	}
	else {
		cout << endl << "--> BUILDING KRIG..." << endl;
		startT = getTime(); // Get start time		
		aKrig.Train(input, target);
		cout << "GP Training Time: " << (getTime() - startT)/1000.0 << "ms"  << endl;
	}
	iModel = temp;
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
	
	// Start optimizing in TRF
	for(i=0; i<trIteration; i++) {
			
		// debug code
		cout << "Start point" << endl; printVec(cout, xStart);
		cout << "End point" << endl; printVec(cout, xOptLoc);			
		//
		// perform local search
		//
		findLocalOpt( xStart, xOptLoc, boundTR ); 
		// Nghia: exact function value	
		double fitOptLoc = objectiveFunction( xOptLoc );

		// debug code 
		if (fitOptLoc != fitOptLoc) {
			cout << "DIE here!" << endl;
			ofstream inputFile("logs/input.dat",ios::trunc); inputFile << input; inputFile.close();
			ofstream targetFile("logs/target.dat",ios::trunc); targetFile << target ; targetFile.close();
			ofstream startFile("logs/startPoint.dat",ios::trunc); printVec(startFile, xStart); startFile.close();			
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

		
	}
	
	// Database: create end point
	ChromosomeT<double> chrom2(xStart);
	Individual endPoint(chrom2);
	endPoint.setFitness(fitStart);
	//double fi = startPoint.getFitness() - endPoint.getFitness();
	DBEntry entry(0, startPoint, endPoint, 1, 1);
	cout << "Add point to List - " << iModel << endl;
	database.add(entry, iModel); 

	// Return fitness to reference variable	
	locFit = fitStart;
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
	if (nIL == 1) iModel = 1;
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

	appendDatabase( x, tempFit );

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
