//===========================================================================
/*!
 *  \file simpleRNNet.cpp
 *
 *
 *  \par Copyright (c) 1999,2003:
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
 *      $RCSfile: simpleRNNet.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: simpleRNNet.cpp,v $
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


#include "ReClaM/MSERNNet.h"
#include "ReClaM/Rprop.h"
#include "Array/Array.h"

using namespace std;

// Define own network class to combine a recurrent
// neural network with a mean squared error measure
// and the RProp algorithm for optimization:
//
class MyNet : public IRpropPlus, public MSERNNet 
{
    public:

        // Standard network constructor:
        MyNet( ) : MSERNNet( ) { };

        // Constructor for a network based on a given connection
        // matrix:
        MyNet( Array< int > con ) : MSERNNet( con ) { };
    
        // This method will return the derivatives of the error
        // of the network.
        Array< double> getDedw( ){ return dedw; }

  
    private:

        // Define own activation function for the
        // hidden neurons, corresponding to the
        // prediction task (see below):
        double g( double a )  { return tanh( a ); }

        // A way to calculate the derivative of the
        // activation function, when given the result
        // of the activation:
        double dg( double ga ) { return 1 - ga * ga; }

};


int main()
{
    // Just counter variables:
    unsigned i, j, k; 

    unsigned numberOfHiddenNeurons( 5 );

    // We use a network here with one input, one
    // output and five hidden neurons:
    unsigned numberOfNeurons = numberOfHiddenNeurons + 1 + 1;

    // We will use the network for a mini time series,
    // where we will predict the value for one step in
    // the future:
    unsigned numberOfTimeSteps = 2;

    // We will use 100 input patterns...
    unsigned numberOfDataPoints = 100;

    unsigned numberOfLearningCycles = 100;

    // ..., where 90 will be used for training
    // and the rest for initialization of the
    // network:
    unsigned warmUpLength = 10;


    // Create connection matrix for two time steps and
    // 7 neurons:
    Array< int > con( numberOfTimeSteps, numberOfNeurons, numberOfNeurons );

    // For the initial time step establish connections
    // from the input neuron to each hidden neuron and
    // from each hidden neuron to the output neuron:
    for ( i = 1; i < con.dim( 1 ); i++ )
    {
        for ( j = 0; j <= i; j++ )
	{
            con( 0,i,j ) = 1;
        }
    }

    // For the second time step establish connections
    // between all neurons:
    for ( k = 1; k < con.dim( 0 ); k++ )
    {
        for ( i = 0; i < con.dim( 1 ); i++ )
	{
            for ( j = 0; j < con.dim( 2 ); j++ )
	    {
	        con( k, i, j ) = 1;
            }
        } 
    } 

    // Construct new network object on the base of the
    // connection matrix:
    MyNet net( con );

    // Initialize the weights of the neuron connections
    // by uniformally distributed random numbers between
    // "-0.1" and "0.1":
    net.initWeights( -0.1, 0.1 );

    // Use the first ten input patterns for initializing
    // the internal state of the network:
    net.includeWarmUp( warmUpLength );

    // Create the data set.  
    Array< double > input(  numberOfDataPoints, 1 ), 
                    target( numberOfDataPoints, 1 );

    // The task is to predict the next amplitude of an sin-like 
    // oszillation:
    for ( i = 0; i < input.dim( 0 ); i++ )
    {
        input( i, 0 )  = sin( 0.1 * i );
        target( i, 0 ) = sin( 0.1 * ( i + 1 ) );
    }


    // Initialize the RProp optimization algorithm:
    net.initRprop( 0.01 );

   
    cout << "Train network:\n" << endl 
         << "No. of cycles\tIteration Error:" << endl
         << "--------------------------------" << endl;

    // Train the network by using the RProp algorithm,
    // monitor error each 10 training cycles:
    for ( i = 1; i <= numberOfLearningCycles; i++ )
    {
        net.rprop( input, target );
        if ( i % 10 == 0 ) 
	{
            cout << i << "\t\t" 
                 << net.error( input, target ) << endl;
        }
    }
   
    // Get the prediction results of the trained network:
    Array< double > output( target.dim( 0 ), target.dim( 1 ) );
    net.model( input, output );
  
    cout << endl << "Evaluate trained network:\n" << endl 
         << "x( t ):\t\tx( t + 1 ):\tprediction:" << endl
         << "-------------------------------------------" << endl; 

    // Output of each data input value, its value
    // in the next time step and the prediction of the
    // trained network:
    for ( i = warmUpLength; i < target.dim( 0 ); i++ ) 
    {
        cout << input( i, 0 ) << "\t" << target( i, 0 ) << "\t" 
             << output( i, 0 ) << endl;
    }

    exit( EXIT_SUCCESS );
}
