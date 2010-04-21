//===========================================================================
/*!
 *  \file simpleRBFNet.cpp
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
 *      $RCSfile: simpleRBFNet.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: simpleRBFNet.cpp,v $
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


#include <ReClaM/Rprop.h>
#include <ReClaM/RBFNet.h>
#include <ReClaM/ErrorMeasures.h>
#include <ReClaM/MeanSquaredError.h>
#include <Array/ArrayOp.h>
#include <Rng/GlobalRng.h>

using namespace std;


// Define own network class to combine a radial basis function
// neural network with a mean squared error measure for training
// and the RProp algorithm for optimization:
// For monitoring purposes, there is also included the
// binary criterion error measure.
//
class MyRBFN : public RBFNet,
               public IRpropPlus,
	       public MeanSquaredError,
	       public ErrorMeasures
{
    public:

        // Constructor for a network with basis structure but no content:
        MyRBFN( unsigned numInput, unsigned numOutput, unsigned numHidden )
            : RBFNet( numInput, numOutput, numHidden ) { }

};



int main( )
{ 
    // Just counter variables:
    unsigned i, t;

    // No. of neurons (input, output, hidden):
    const unsigned inputDimension = 1;
    const unsigned outputDimension = 1; 
    const unsigned nHidden = 2;

    // Use 100 input patterns:
    const unsigned numSamples = 100;
  
    // Define input and target patterns:
    Array< double > inTrain ( numSamples, inputDimension ),
                    outTrain( numSamples, outputDimension );


    // Create data set using gauss distributed random numbers:
    for ( i = 0; i < numSamples / 2; i++ )
    {
        inTrain( i, 0 ) = Rng::gauss( 2.0, 0.3 );
        outTrain( i, 0 ) = 0;
        inTrain( i + numSamples / 2, 0 ) = Rng::gauss( 5.0, 0.5 );
        outTrain( i + numSamples / 2, 0 ) = 1;
    }
  
    // Create an empty network with two input and hidden neurons and
    // one output neuron:
    MyRBFN rbfn( inputDimension, outputDimension, nHidden );     

    // Use the input and target patterns for initialization of
    // the network:
    rbfn.initRBFNet( inTrain, outTrain );

    // Show performance of the network after initialization:
    cout << "After initialization of RBFN:\n" << endl;
    cout << "Mean squared error:\t" 
         << rbfn.meanSquaredError( inTrain, outTrain ) << endl;
    cout << "Classification error:\t" 
         << rbfn.binaryCriterion( inTrain, outTrain ) << endl;

    // Extract center, variance and weights of connection to
    // output neurons for the hidden neurons...
    Array< double > c = rbfn.getCenter( );
    Array< double > v = rbfn.getVariance( );
    Array< double > w = rbfn.getWeights( );

    // ...and monitor the data:
    cout << "Data of the " << c.dim( 0 ) << " hidden neurons:" << endl;
    for ( i = 0; i < c.dim( 0 ); i++ )
    {
        cout << i + 1 << ". neuron:\n\tcenter = " << c( i, 0 ) 
             << ",\n\tvariance = " 
	     << v( i, 0 ) << ",\n\tweight to output neuron = " << w( 0 , i ) 
             << "." << endl;
    }

    cout << "\n\nTraining of the network:\n" << endl;

    // Initialize the RProp optimization algorithm:
    rbfn.initRprop( 0.0125 );

    cout << "epoch:\ttraining error:\tclassification error:" << endl;
    cout << "---------------------------------------------" << endl;
    // Train the network for ten epochs:
    for ( t = 0; t < 10; ++t ) 
    {
        // Use RProp for training with all patterns:
        rbfn.rprop( inTrain, outTrain );

        // Output of training and monitoring error:
        cout << t << '\t'
      	     << rbfn.error( inTrain, outTrain ) << "\t"
	     << rbfn.binaryCriterion( inTrain, outTrain ) << "\t"
	     << endl;
    }

    cout << "\n\nAfter Training of RBFN:\n" << endl;

    // Show performance of the network after training:
    cout << "Mean squared error:\t" 
         << rbfn.meanSquaredError( inTrain, outTrain ) << endl;
    cout << "Classification error:\t" 
         << rbfn.binaryCriterion( inTrain, outTrain ) << endl;

    // Extract center, variance and weights of connection to
    // output neurons for the hidden neurons...
    c = rbfn.getCenter( );
    v = rbfn.getVariance( );
    w = rbfn.getWeights( );

    // ...and monitor the data:
    cout << "Data of the " << c.dim( 0 ) << " hidden neurons:" << endl;
    for ( i = 0; i < c.dim( 0 ); i++ )
    {
        cout << i + 1 << ". neuron:\n\tcenter = " << c( i, 0 ) 
             << ",\n\tvariance = " << v( i, 0 ) 
             << ",\n\tweight to output neuron = " 
             << w( 0, i ) << "." << endl;
    }

    exit( EXIT_SUCCESS );
}

