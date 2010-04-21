
#include "sqr.h"
#include "ZNconfig.h"
#include "Population.h"
#include "Array.h"
#include "EvoNet.h"


using namespace std;

//=======================================================================
//
// fitness function: Ackley's function
//
double ackley( const std::vector< double >& x )
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

//=======================================================================
//
// fitness function: sphere model
//
double sphere( const vector< double >& x )
{
    unsigned i;
    double   sum;
    for( sum = 0., i = 0; i < x.size( ); i++ )
        sum += sqr( i * ( x[ i ] + i ) );
    return sum;
}
double sphere1( const vector< double >& x )
{
    unsigned i;
    double   sum;
    for( sum = 0., i = 0; i < x.size( ); i++ )
        sum += sqr( ( x[ i ] ) );
    return sum;
}

//=======================================================================
//
// main program
//

int main( int argc, char **argv )
{
  cout << "Starting" << endl;

  Database myData;

  FILE *fout, *fout2, *fout3;
  fout = fopen("result1.dat" ,"wt");
  fout2 = fopen("path1.dat" ,"wt");
  fout3 = fopen("path2.dat" ,"wt");


  FILE *fNet1, *fNet2, *fNet3;
  fNet1 = fopen("NetOut1.dat","wt");
  fNet2 = fopen("NetOut2.dat","wt");
  fNet3 = fopen("NetOut3.dat","wt");
  
  // test();
    
  //cout << "fertig" << endl;

  //exit(0);

 //
    // constants
    //
    const unsigned Mu           = 2;
    const unsigned Lambda       = 10;
    const unsigned Dimension    = 30;
    const unsigned Iterations   = 500;
    const unsigned Interval     = 1;
    const unsigned NSigma       = 1;

    const double   MinInit      = -3;
    const double   MaxInit      = +15;
    const double   SigmaInit    = 3;

    const bool     PlusStrategy = false;

    unsigned       i, t;

    //
    // initialize random number generator
    //
    Rng::seed( argc > 1 ? atoi( argv[ 1 ] ) : 1234 );

    //
    // define populations
    //
    Population parents   ( Mu,     ChromosomeT< double >( Dimension ),
			           ChromosomeT< double >( NSigma    ) );
    Population offsprings( Lambda, ChromosomeT< double >( Dimension ),
			           ChromosomeT< double >( NSigma    ) );

    //
    // minimization task
    //
    parents   .setMinimize( );
    offsprings.setMinimize( );

    //
    // initialize parent population
    //
    for( i = 0; i < parents.size( ); ++i ) {
        dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize( MinInit,   MaxInit   );
	dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 1 ] ).initialize( SigmaInit, SigmaInit );
    }

    //
    // selection parameters (number of elitists)
    //
    unsigned numElitists = PlusStrategy ? Mu : 0;

    //
    // standard deviations for mutation of sigma
    //
    double     tau0 = 1. / sqrt( 2. * Dimension );
    double     tau1 = 1. / sqrt( 2. * sqrt( ( double )Dimension ) );

    //
    // evaluate parents (only needed for elitist strategy)
    //
    if( PlusStrategy )
        for( i = 0; i < parents.size( ); ++i )
	    parents[ i ].setFitness( sphere1( dynamic_cast< std::vector< double >& >( parents[ i ][ 0 ] ) ) );

    unsigned inDimension = Dimension;
    unsigned outDimension = 1;
    unsigned numberOfHiddenNeurons = 5;
    unsigned numberOfNeurons = numberOfHiddenNeurons+inDimension+outDimension;

    Array< int > cmat( numberOfNeurons, numberOfNeurons+1 );
    Array<double> wmat( numberOfNeurons, numberOfNeurons+1 );
      //Initialize a feed-forward network

   for ( unsigned ii = inDimension; ii < inDimension + numberOfHiddenNeurons; ii++ ){
     for ( unsigned jj = 0; jj < inDimension ; jj++ ){
       cmat( ii, jj ) = 1;
       wmat( ii, jj) = Rng::uni(-0.01,0.01);
     }
   }
   // Set bias values:
    for ( unsigned ii = inDimension; ii < inDimension + numberOfHiddenNeurons+outDimension; ii++ ){
      cmat( ii, cmat.dim( 1 )-1 ) = 1;
      wmat( ii, cmat.dim( 1 )-1) = Rng::uni(-0.01,0.01);;
    }
  
    //Set connections from hidden to output
    
    for ( unsigned ii = inDimension+numberOfHiddenNeurons; ii < inDimension + numberOfHiddenNeurons+outDimension; ii++ ){
      for (unsigned jj =inDimension; jj < inDimension+numberOfHiddenNeurons; jj++ ){
	cmat( ii, jj ) = 1;
	wmat( ii, jj) = Rng::uni(-0.01,0.01);
      }
    }

    
    EvoNet Network((int)Dimension, 1, cmat, wmat);
    myData.SetMaxTrainingData(300);
    // EvoNet Network("../data/EvoTest0.net");
    
    //
    // iterate
    //
    for( t = 0; t < Iterations; ++t ) {
        //
        // generate new offsprings
        //
        for( i = 0; i < offsprings.size( ); ++i ) {
	    //
	    // select two random parents
	    //
	    Individual& mom = parents.random( );
	    Individual& dad = parents.random( );

	    //
	    // define temporary references for convenience
	    //
	    ChromosomeT< double >& objvar = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] );
	    ChromosomeT< double >& sigma  = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 1 ] );

	    //
	    // recombine object variables discrete, step sizes intermediate
	    //
	    objvar.recombineDiscrete       ( mom[ 0 ], dad[ 0 ] );
	    sigma .recombineGenIntermediate( mom[ 1 ], dad[ 1 ] );

	    //
	    // mutate object variables normal distributed,
	    // step sizes log normal distributed
	    //
	    sigma .mutateLogNormal( tau0,  tau1 );
	    objvar.mutateNormal   ( sigma, true );
	}

	//
	// evaluate objective function (parameters in chromosome #0)
	//
	

	// test = 1 // test with Shark interfaces
	// test = 2 // test with Array interfaces
	int test = 2;
	int fraction = 2000; // fraction of network calls
	//int fraction2 = 20; // fraction of network optimisations
	if ((t % (unsigned int)fraction) != (unsigned int)fraction - 1) {
	  // real fitness function  
	  if (test == 2) {
	    // allocate memory for training data
	    Array <double> Input(offsprings.size( ), Dimension);
	    Array <double> Target(offsprings.size( ), 1);
	    Array <double> Output(offsprings.size( ), 1);
	    
	    // store training data in temporary variables
	    for( i = 0; i < offsprings.size( ); ++i ) {
	      
	      for (unsigned int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
		Input(i,j) = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ])[j];
		//cout << Input(i,j) << " ";
	      }
	      
	      
	      Target(i,0) = sphere1( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] ));
	      //cout << " - " << Target(i,0) << endl;
	      
	    }
	

	    // ----------------------------------------
	    // update the network

	    // evaluate before training

	    if (t != 0)	    Network.Evaluate(Input, Output);
	    double NNerror = 0.0;
	    
	    for( i = 0; i < offsprings.size( ); ++i ){
	      fprintf(fNet2,"%lf ",Output(i,0));
	      fprintf(fNet3,"%lf ",Target(i,0));
	      NNerror += sqr( (Output(i, 0)-Target(i,0)));
	    }

	    NNerror /= (double)offsprings.size();	    
	    NNerror = sqrt(NNerror);
	    
	    cout << "Before Training: " << endl;
	    cout << "  on population        : " << NNerror << endl;
	    if (t != 0) 
	      cout << "  on training data set : " << Network.MSE(myData) << endl;
	    if (t != 0) 
	      //	    cout << "  on training data set " <<  Network.meanSquaredError(Input, Target) << endl;
	    
  

	    fprintf(fNet2,"\n");
	    fprintf(fNet3,"\n");
	    

	    cout << "Before Training: " << NNerror << endl;
	    fprintf(fNet1,"%lf ", NNerror);
	    
	    // add training data to database
	    if (t != 5) myData.AddTrainingData(Input,Target);
	    else {
	      cout << "This time we do not add the new training data" << endl;
	    }
	    cout << "Train starts!" << endl;
	    Network.Train_Rprop(myData);
	    cout << "Train finished!" << endl;
	    
	    // evaluate after training 
	    Network.Evaluate(Input, Output);
	    NNerror = 0.0;
	    
	    for( i = 0; i < offsprings.size( ); ++i ){
	      NNerror += sqr( (Output(i, 0)-Target(i,0) ));
	    }
	    NNerror /= (double)(offsprings.size( ));
	    NNerror = sqrt(NNerror);
	    
	    cout << "After Training: " << endl;
	    cout << "  on population        : " << NNerror << endl;
	    cout << "  on training data set : " << Network.MSE(myData) << endl;
	    fprintf(fNet1,"%lf \n", NNerror);
	  }
	  fflush(fNet1);
	  fflush(fNet2);
	  fflush(fNet3);
	  

	  // assign the evaluated fitness values to the offpsrings
	  
	  for( i = 0; i < offsprings.size( ); ++i )
	    offsprings[ i ].setFitness( sphere1( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] )));
	  
	  // the Shark interface
	  if (test == 1) {
	    // train network
	    //	    Network.Train(offsprings);	  
	  }


	} else {
	  // network approximation

	  if (test == 1) {
	    
	    //	    Network.Evaluate(offsprings);
	  }
	
	  if (test == 2) { 
	      Array <double> Input(offsprings.size( ), Dimension);
	      Array <double> Output(offsprings.size( ), 1);
	      for( i = 0; i < offsprings.size( ); ++i )
		for (unsigned int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
		  Input(i,j) = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ])[j];
		}
	      
	      Network.Evaluate(Input, Output);
	      for( i = 0; i < offsprings.size( ); ++i ){
		offsprings[ i ].setFitness(Output(i,0)); 
		
	      }   
	  }
	}

	for( i = 0; i < offsprings.size( ); ++i ) {
	  
	  for (unsigned int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
	    fprintf(fout2,"%e ", dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ])[j]);
	  }
	  fprintf(fout2,"\n");
	  
	}

 

	//
	// select (mu,lambda) or (mu+lambda)
	//
	parents.selectMuLambda( offsprings, numElitists );

	//
	// print out best value found so far
	//
	if( t % Interval == 0 )
	    std::cout << t << "\tbest value = "
		      << parents.best( ).fitnessValue( ) 
		      << " " << dynamic_cast< ChromosomeT< double >& >( parents.best()[ 1 ])[ 0 ] 
		      << std::endl;
	fprintf(fout, "%d %e\n", t, log(parents.best().fitnessValue()));
	fflush (fout);

	cout << "----------------------------" << endl;
	
    }

    fclose(fout2);
    fclose(fout);
    return 0;
}





