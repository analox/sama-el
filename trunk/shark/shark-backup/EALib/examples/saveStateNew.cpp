/*! saveState.cpp
* ======================================================================
*
*  File    :  stateState.cpp
*  Created :  12.08.1995
*
*  Copyright (c) 1995,1999 Martin Kreutz
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
*        $RCSfile: saveStateNew.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: saveStateNew.cpp,v $
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
*  Last update: 13.03.1999
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


#include "sqr.h"
#include "MathConst.h"
#include "EvoAlg.h"

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
// 
//
class ES : public EvoAlg
{
  public:
    ES( unsigned mu,
	unsigned lambda,
	unsigned dimension,
	unsigned nsigma,
	double   min,
	double   max,
	double   sigma,
	bool     plus,
	const string& saveFile = "" )
      : EvoAlg( mu,     Individual( ChromosomeT< double >( dimension ),
				    ChromosomeT< double >( nsigma    ) ),
		lambda, Individual( ChromosomeT< double >( dimension ),
				    ChromosomeT< double >( nsigma    ) ),
		saveFile )
    {
        minInit      = min;
	maxInit      = max;
	sigmaInit    = sigma;
	plusStrategy = plus;

	//
	// selection parameters (number of elitists)
	//
	numElitists = plusStrategy ? parents.size( ) : 0;

	//
	// standard deviations for mutation of sigma
	//
	tau0 = 1. / sqrt( 2. * parents[ 0 ][ 0 ].size( ) );
	tau1 = 1. / sqrt( 2. * sqrt( ( double )parents[ 0 ][ 0 ].size( ) ) );

	if( generation == 0 )   // no lock file loaded (, probably)
	    initialize( );
    }

    void initialize( )
    {
	//
	// initialize parent population
	//
	for( unsigned i = 0; i < parents.size( ); ++i ) {
	    dynamic_cast< ChromosomeT< double >& >
	        ( parents[ i ][ 0 ] ).initialize( minInit,   maxInit   );
	    dynamic_cast< ChromosomeT< double >& >
	        ( parents[ i ][ 1 ] ).initialize( sigmaInit, sigmaInit );
	}

	//
	// evaluate parents (only needed for elitist strategy)
	//
	if( plusStrategy ) {
	    for( unsigned i = 0; i < parents.size( ); ++i )
	        parents[ i ].setFitness( ackley(
                    dynamic_cast< vector< double >& >( parents[ i ][ 0 ] ) ) );
	}
    }

    void generate( )
    {
        //
        // generate new offsprings
        //
        for( unsigned i = 0; i < offsprings.size( ); ++i ) {
	    //
	    // select two random parents
	    //
	    Individual& mom = parents.random( );
	    Individual& dad = parents.random( );

	    //
	    // define temporary references for convenience
	    //
	    ChromosomeT< double >& objvar =
	        dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] );
	    ChromosomeT< double >& sigma  =
	        dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 1 ] );

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
    }

    void evaluate( )
    {
	//
	// evaluate objective function (parameters in chromosome #0)
	//
	for( unsigned i = 0; i < offsprings.size( ); ++i )
	    offsprings[ i ].setFitness( ackley(
                dynamic_cast< ChromosomeT< double >& >
		    ( offsprings[ i ][ 0 ] ) ) );
    }

    void select( )
    {
	//
	// select (mu,lambda) or (mu+lambda)
	//
	parents.selectMuLambda( offsprings, numElitists );
    }

    void display( )
    {
        //
        // print out best value found so far
        //
        cout << generation << "\tbest value = "
	     << parents.best( ).fitnessValue( ) << endl;
    }

  private:
    double   minInit;
    double   maxInit;
    double   sigmaInit;
    bool     plusStrategy;

    unsigned numElitists;
    double   tau0;
    double   tau1;
};


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
    const unsigned NSigma       = 1;
    const unsigned Iterations   = 500;
    const unsigned Interval     = 10;

    const double   MinInit      = -3;
    const double   MaxInit      = +15;
    const double   SigmaInit    = 3;

    const bool     PlusStrategy = false;

    const string   SaveFile     = "save.state";

    ES es( Mu, Lambda, Dimension, NSigma,
	   MinInit, MaxInit, SigmaInit, PlusStrategy, SaveFile );

    es.setMaxGeneration( Iterations );
    es.setDisplayInterval ( Interval );
    es.setMinimize     ( );

    //
    // iterate
    //
    while( ! es.terminated( ) )
        es.iterate( );

    return 0;
}
