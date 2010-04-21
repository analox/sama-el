/*! saveState.cpp
* ======================================================================
*
*  File    :  stateState.cpp
*  Created :  12.08.1995
*
*  Copyright (c) 1995,1998 Martin Kreutz
*
*  Institut fuer Neuroinformatik
*  Ruhr-Universitaet Bochum
*  44780 Bochum, Germany<BR>
*  Phone: +49-234-32-25558<BR>
*  Fax:   +49-234-32-14209<BR>
*  eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
*  www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
*  
*    \par Project:
*        EALib
*  
*    \par Language and Compiler:
*        C++, egcs (Linux)
*  
*    \par File and Revision:
*        $RCSfile: saveState.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: saveState.cpp,v $
*        Revision 2.1  2004/06/17 23:05:26  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.0  2003/11/28 16:23:08  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:03  shark-admin
*        INI Administration
*  <BR>
*  
*  
*  
*  Last update: 17.08.1998
*
* ----------------------------------------------------------------------
*
*  This file is part of the EALib. This library is free software;
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
* ======================================================================
*/

#include <cerrno>
#include <cstdio>
#include <iomanip>
#include <fstream>

#include "sqr.h"
#include "MathConst.h"
#include "Population.h"

#ifdef __UNIX__
#include <unistd.h>
#else
extern "C" int unlink( const char * );
#endif

using namespace std;

//=======================================================================
//
// fitness function: Ackley's function
//
double ackley( const vector< double >& x )
{
    const double A = 20.;
    const double B = 0.2;
    const double C = Pi2;

    unsigned i, n;
    double   a, b;

    for( a = b = 0., i = 0, n = x.size( ); i < n; ++i ) {
        a += x[ i ] * x[ i ];
	b += cos( C * x[ i ] );
    }

    return -A * exp( -B * sqrt( a / n ) ) - exp( b / n ) + A + E;
}


//=======================================================================
//
// main program
//
int main( int argc, char **argv )
{
    //
    // constants
    //
    const unsigned Mu           = 15;
    const unsigned Lambda       = 100;
    const unsigned Dimension    = 30;
    const unsigned Iterations   = 500;
    const unsigned Interval     = 10;
    const unsigned NSigma       = 1;

    const double   MinInit      = -3;
    const double   MaxInit      = +15;
    const double   SigmaInit    = 3;

    const bool     PlusStrategy = false;

    const char*    LockFile     = "saveState.lock";
    const char*    TmpFile      = "saveState.tmp";

    const unsigned Precision    = 16;

    unsigned       i, t;

    ifstream       loadFile;
    ofstream       saveFile;
    ofstream       logFile;

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
    double   tau0 = 1. / sqrt( 2. * Dimension );
    double   tau1 = 1. / sqrt( 2. * sqrt( ( double )Dimension ) );

    //
    // evaluate parents (only needed for elitist strategy)
    //
    if( PlusStrategy )
        for( i = 0; i < parents.size( ); ++i )
	    parents[ i ].setFitness( ackley( dynamic_cast< vector< double >& >( parents[ i ][ 0 ] ) ) );

    //
    // check whether a process was/is running (lock file exists)
    //
    loadFile.open( LockFile );

    if( loadFile ) {
	//
	// load current state
	//
	loadFile >> t;
	Rng::readFrom( loadFile );
	loadFile >> parents >> offsprings;

	if( ! loadFile ) {
	    cerr << "error: read failed on file " << LockFile << endl;
	    exit( 1 );
	}

	loadFile.close( );
    } else
        t = 0; // reset number of generations

    //
    // iterate
    //
    for( ; t < Iterations; ++t ) {
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
	for( i = 0; i < offsprings.size( ); ++i )
	    offsprings[ i ].setFitness( ackley( dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] ) ) );

	//
	// select (mu,lambda) or (mu+lambda)
	//
	parents.selectMuLambda( offsprings, numElitists );

	if( t % Interval == 0 ) {
	    //
	    // save current state,
	    // use temporary file in order to avoid inconsistent states
	    //
	    saveFile.open( TmpFile );
	    if( ! saveFile ) {
	        cerr << "error: can't open file " << TmpFile << endl;
		exit( 1 );
	    }
	    saveFile << setprecision( Precision );
	    saveFile << t+1 << endl;
	    Rng::writeTo( saveFile );
	    saveFile << parents << offsprings;
	    if( ! saveFile ) {
	        cerr << "error: can't write to file " << TmpFile << endl;
		exit( 1 );
	    }
	    saveFile.close( );
#ifdef _WIN32
            //
            // windows doesn't overwrite existing files
            // so remove old file before renaming new file
            //
            if( t > 0 )   // lock file exists
                if( unlink( LockFile ) != 0 ) {
                    perror( argv[ 0 ] );
                    cerr << "error: can't unlink file " << LockFile << endl;
	            exit( 1 );
                }
#endif
	    if( rename( TmpFile, LockFile ) != 0 ) {
	        perror( argv[ 0 ] );
		cerr << "error: can't rename temporary file to " << LockFile << endl;
		exit( 1 );
	    }

	    //
	    // print out best value found so far
	    //
	    cout << t << "\tbest value = "
		 << parents.best( ).fitnessValue( ) << endl;
	}
    }

    //
    // clean up and unlink an existing lock file
    //
    if( unlink( LockFile ) != 0 ) {
        perror( argv[ 0 ] );
        cerr << "error: can't unlink file " << LockFile << endl;
	exit( 1 );
    }

    return 0;
}
