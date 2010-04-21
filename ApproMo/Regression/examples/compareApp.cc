
#include "sqr.h"
#include "ZNconfig.h"
#include "Population.h"
#include "Array.h"
#include "RegrApp.h"

// #include "EvoNet.h"


using namespace std;

//=======================================================================
//
// fitness function: sphere model
//
double sphere( const vector< double >& x )
{
    unsigned i;
    double   sum;
    for( sum = 0., i = 0; i < x.size( ); i++ )
      sum += sqr( (x[i]) );
    return sum;
}

//=======================================================================

//
// main program
//

int main( int argc, char **argv )
{

  // Database
  Database myData;

  // Regression models
  RegrApp Network(1);

  // LogFile
  FILE *fout;
  fout = fopen("LogFile.dat" ,"wt");

  const unsigned Mu           = 1;
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
  Population parents   ( Mu,    
			 ChromosomeT< double >( Dimension ),
			 ChromosomeT< double >( NSigma    ) );
  Population offsprings( Lambda, 
			 ChromosomeT< double >( Dimension ),
			 ChromosomeT< double >( NSigma    ) );
  
  //
  // minimization task
  //
  parents   .setMinimize( );
  offsprings.setMinimize( );
  
  cout << "Starting..." << endl;

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
      parents[ i ].setFitness( sphere( dynamic_cast< std::vector< double >& >( parents[ i ][ 0 ] ) ) );
  
  Array < int > cmat(Dimension*3,Dimension*3+1);
  for ( unsigned i = 0; i < Dimension; i++) {
    for ( int j = 0; j < Dimension; j++) {
      cmat(i,j) = 1;
    }
  }
  Array < double > wmat(Dimension*3,Dimension*3+1);
  for ( unsigned i = 0; i < Dimension; i++) {
    for ( int j = 0; j < Dimension; j++) {
      wmat(i,j) = Rng::uni(-0.1,0.1);
    }
  }

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
    
    int fraction = 2000; // fraction of network calls
    int fraction2 = 20; // fraction of network optimisations
    
    // real fitness function  
    // allocate memory for training data
    Array <double> Input(offsprings.size( ), Dimension);
    Array <double> Target(offsprings.size( ), 1);
    Array <double> Output(offsprings.size( ), 1);
    
    // store training data in temporary variables
    for( i = 0; i < offsprings.size( ); ++i ) {
      for (int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
	Input(i,j) = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ])[j];
	Target(i,0) = sphere( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] ));
      }
    }

    // Evaluation
    for( i = 0; i < offsprings.size( ); ++i ) {
      offsprings[ i ].setFitness( sphere( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] )));
    }

    // Selection
    parents.selectMuLambda( offsprings, numElitists );

    // Training
    myData.AddTrainingData( Input,Target);
    Network.Train(myData);
    cout << "Network Error: " << Network.Error_MSE(myData) << endl;

    cout  << t << " Best value = "
	  << parents.best( ).fitnessValue( ) 
	  << " " << dynamic_cast< ChromosomeT< double >& >( parents.best()[ 1 ])[ 0 ] 
	  << endl;
    fprintf( fout, "%i %f \n", t, parents.best().fitnessValue() );
    
    double compi = parents.best().fitnessValue();
    if ( compi < 0.00001 ) {
      cout << "Evolution done!" << endl;
      fclose (fout);
      exit (0);
    }
  }
  fclose(fout);
}

	/*	Target(i,0) = sphere( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] ));
      }
    
    

    
    // ----------------------------------------
    // update the model
    
    // evaluate before training
    
      Network.Evaluate(Input, Output);
      Network.Train(myData);


	cout << "done!" << endl;	double NNerror = 0.0;
	
	for( i = 0; i < offsprings.size( ); ++i ){
	  //cout << Output(i, 0) << " " << Target(i,0) << "  -> " << Output(i, 0)-Target(i,0)<<endl;
	  
	  fprintf(fNet2,"%lf ",Output(i,0));
	  fprintf(fNet3,"%lf ",Target(i,0));
	  
	  NNerror += sqr( Output(i, 0)-Target(i,0));
	  
	  //Network.Train(Input,Output);
	}
	fprintf(fNet2,"\n");
	fprintf(fNet3,"\n");
	
	NNerror /= (offsprings.size( ));
	cout << "Before Training: " << NNerror << endl;
	fprintf(fNet1,"%f ", NNerror);
	
	// add training data to database
	if (t != 5) myData.AddTrainingData(Input,Target);
	else {
	  cout << "This time we do not add the new training data" << endl;
	}
	// do the structure optimisation
	if ((t % fraction2) == fraction2 - 1) {
	  //Network.StructureOptimization(myData);//, Input,Target);
	}
	else {
	  cout << "train only" << endl;
	  Network.Train(myData);
	  
	}
	
	// evaluate after training 
	Network.Evaluate(Input, Output);
	NNerror = 0.0;
	
	for( i = 0; i < offsprings.size( ); ++i ){
	  //cout << Output(i, 0) << " " << Target(i,0) << "  -> " << Output(i, 0)-Target(i,0)<<endl;
	  NNerror += sqr( Output(i, 0)-Target(i,0));
	  
	  //Network.Train(Input,Output);
	}
	NNerror /= (offsprings.size( ));
	cout << "After Training: " << NNerror << endl;
	fprintf(fNet1,"%f \n", NNerror);
      }
      fflush(fNet1);
      fflush(fNet2);
      fflush(fNet3);
      
      
      // assign the evaluated fitness values to the offpsrings
      
      for( i = 0; i < offsprings.size( ); ++i )
	//offsprings[ i ].setFitness( ackley( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] )));
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
	  for (int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
	    Input(i,j) = dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ])[j];
	  }
	
	Network.Evaluate(Input, Output);
	for( i = 0; i < offsprings.size( ); ++i ){
	  offsprings[ i ].setFitness(Output(i,0)); 
	  
	}   
      }
    }
    
    for( i = 0; i < offsprings.size( ); ++i ) {
      
      for (int j = 0; j < dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ]).size(); j++) {
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

      */



