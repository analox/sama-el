//===========================================================================
/*!
 *  \file simpleFFNet.cpp
 *
 *
 *  \par Copyright (c) 1999, 2003:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *      <BR>
 *
 *  \par Project:
 *      ReClaM
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: simpleFFNet.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: simpleFFNet.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:13  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:05  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of ReClaM. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 */
//===========================================================================


#include "Array/Array.h"
#include "ReClaM/FFNet.h"
#include "ReClaM/Rprop.h"
#include "ReClaM/MeanSquaredError.h"
#include "ReClaM/ErrorMeasures.h"
#include "Rng/GlobalRng.h"

using namespace std;


// Define own network class to combine a feed forward
// neural network with a mean squared error measure for training
// and the RProp algorithm for optimization by resilent backpropagation.
// The mean squared error is only used for training, but the
// error percentage will be used to evaluate the fitness of the network.
// Then the best network depending on the error percentage will
// be used as final result:
//
class MyNet :	public FFNet,
		public IRpropPlus,
		public MeanSquaredError,
		public ErrorMeasures
{	

    public:
  
        // Constructor for a network based on a given connection
        // matrix:
        MyNet( unsigned nIn, unsigned nOut, const Array< int >& cmat ) : 
            FFNet( nIn, nOut, cmat ) { }

        // Constructor for reading the network structure from
        // a file:
        MyNet( const string &filename ) : FFNet( filename ) { }
  

        // Define alternative sigmoid activation function for 
        // the hidden units:
        double g( double a ) 
        { 
            return a / ( 1 + fabs( a ) ); 
        }
  
        // A way to calculate the derivative of the
        // activation function, when given the result
        // of the activation:
        double dg( double ga ) 
        { 
            return  ( 1 - sgn( ga ) * ga ) * ( 1 - sgn( ga ) * ga );
        }
  
  
        // Use linear output neurons:
        double gOutput( double a ) 
        { 
            return a; 
        }
  
        // Because of linearity, the derivative
        // the output neurons will always be a constant:
        double dgOutput(double ga) 
        { 
            return 1;
        }

};


int main()
{
    // Just counter variables:
    unsigned i,j;


    unsigned numberOfHiddenNeurons = 5;

    // We use a network here with two input neurons, one
    // output neuron and five hidden neurons:
    unsigned numberOfNeurons = numberOfHiddenNeurons+2+1;

    unsigned numberOfLearningCycles = 100;


    // Create connection matrix for the neurons:
    Array< int > con( numberOfNeurons, numberOfNeurons+1 );

    // Establish connections from the input neurons to each 
    // hidden neuron and from each hidden neuron to the output neuron:
    for ( i = 2; i < con.dim( 0 ); i++ )
    {
        for ( j = 0; j < i; j++ )
	{
            con( i, j ) = 1;
        }
    }

    // Set bias values:
    for ( i = 2; i < con.dim( 0 ); i++ )
    {
        con( i, con.dim( 1 )-1 ) = 1;
    }
  
    // Construct two new network objects on the base of the
    // connection matrix. The first one will be used for
    // training, the second one will store the network
    // with the smallest validation error:
    MyNet net( 2, 1, con ), 
          netMin( 2, 1, con );


    // Initialize the weights of the neuron connections
    // by uniformally distributed random numbers between
    // "-0.1" and "0.1":
    net.initWeights( -0.1, 0.1 );

    // The untrained network is the best one until now.
    netMin = net;
  
    // Use two distinct data sets. One will only be used for
    // training, the other will be used to examine
    // the fitness of the trained network:
    Array< double > inputTrain( 121, 2u ), targetTrain( 121, 1u );
    Array< double > inputValidate( 121, 2u ), targetValidate( 121,1u );

    // Construct data:
    for ( i = 0; i < 6; i++ )
    {
        for ( j = 0; j < 6; j++ )
        {
            inputTrain( i * 6 + j, 0 ) = i * 0.4 - 1;
            inputTrain( i * 6 + j, 1 ) = j * 0.4 - 1;

            inputValidate( i * 6 + j, 0 ) = i * 0.4 - 1;
            inputValidate( i * 6 + j, 1 ) = j * 0.4 - 1;

            targetTrain( i * 6 + j, 0 ) = 
                ( i * 0.4 - 1 ) * ( j * 0.4 - 1 ) + Rng::gauss( 0, 0.2 );

            targetValidate( i * 6 + j, 0 ) = 
                ( i * 0.4 - 1 ) * ( j * 0.4 - 1 ) + Rng::gauss( 0, 0.2 );
        }
    }


    // Current error percentage of the training set:
    double t; 
      
    // Current error percentage of the validation set:
    double v;

    // The smallest error percentage for the validation set:
    double vMin = net.errorPercentage( inputValidate, targetValidate );

    // The training epoch, where the network with the
    // smallest error percentage for the validation set 
    // occured:
    unsigned epochMin = 0;

    // Initialize the RProp optimization algorithm:
    net.initRprop( 0.01 );


    cout << "Train network:\n" << endl 
         << "epoch:\ttraining error:\tvalidation error:" << endl
         << "-----------------------------------------" << endl;

    // Training of the network:
    for ( unsigned epoch = 1; epoch <= numberOfLearningCycles; epoch++ ) 
    {
        // Train the net with RProp ...
        net.rprop( inputTrain, targetTrain, 1.2, .5, 50, 0 );

        //  ... and calculate the (monitoring) errors:
        t = net.errorPercentage( inputTrain, targetTrain );
        v = net.errorPercentage( inputValidate, targetValidate );

        // Monitor the results:
        std::cout << epoch << "\t" << t << "\t\t" << v << endl;

        // Memorize the network with smallest validation error:
        if ( v < vMin ) 
        {
            vMin = v;
            netMin = net;
            epochMin = epoch;
        }
    }
    
    //  Output of the performance values of the best network:
    t = netMin.errorPercentage( inputTrain, targetTrain );
    v = netMin.errorPercentage( inputValidate, targetValidate );

    cout << endl << "\n\nError of network with best validation error:\n" 
         << endl
         << "epoch:\ttraining error:\tvalidation error:" << endl
         << "-----------------------------------------" << endl;
    cout << epochMin << "\t" << t << "\t" << v << endl;
    
    //  Output of the structure of the best network:
    cout << "\n\n\nStructure of this network\n(no. of input neurons, "
         << "no. of output neurons, connection matrix, weight matrix):\n" 
         << endl;
    cout.precision( 6 );
    cout.setf( ios::fixed | ios::showpos );
    cout << netMin << endl;
    
    // Show network behaviour for one input pattern:
    cout << "Processing a single input pattern:\n" << endl;
    Array< double > in( 2 ), out( 1 );
    in( 0 ) = 0.3;
    in( 1 ) = -0.1;
    netMin.model( in, out );
    cout << "Input:\t( " << in( 0 ) << ", " << in( 1 ) << " ) " << endl;
    cout << "Output:\t" << out( 0 ) << endl;
    
    exit( EXIT_SUCCESS );
}





