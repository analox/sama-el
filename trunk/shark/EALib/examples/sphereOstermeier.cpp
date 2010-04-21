/*! sphereOstermeier.cpp
* ======================================================================
*
*  File    :  sphereOstermeier.cpp
*  Created :  01.01.1995
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
*        $RCSfile: sphereOstermeier.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: sphereOstermeier.cpp,v $
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
*  Last update: 10.08.1998
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

#include "EALib/sqr.h"
#include "ZNconfig.h"
#include "EALib/Population.h"

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
        sum += sqr( i * ( x[ i ] + i ) );
    return sum;
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
    const unsigned Mu           = 4;
    const unsigned Lambda       = 20;
    const unsigned Dimension    = 30;
    const unsigned Iterations   = 2000;
    const unsigned Interval     = 10;

    const double   MinInit      = -3;
    const double   MaxInit      = +5;
    const double   SigmaInit    = 1;

    unsigned       i, t;

    //
    // initialize random number generator
    //
    Rng::seed( argc > 1 ? atoi( argv[ 1 ] ) : 1234 );

    //
    // parameters for the derandomized step size adaptation
    //
    ChromosomeT< double >::DerandomConst derandom( Dimension );

    //
    // define populations
    //
    Population parents   ( Mu,     ChromosomeT< double >( derandom.nobj( ) ),
			           ChromosomeT< double >( derandom.npar( ) ) );
    Population offsprings( Lambda, ChromosomeT< double >( derandom.nobj( ) ),
			           ChromosomeT< double >( derandom.npar( ) ) );

    //
    // minimization task
    //
    parents   .setMinimize( );
    offsprings.setMinimize( );

    //
    // initialize parent population
    //
    for( i = 0; i < parents.size( ); ++i ) {
        dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 0 ] ).initialize( MinInit, MaxInit );
	dynamic_cast< ChromosomeT< double >& >( parents[ i ][ 1 ] ).initializeDerandom( SigmaInit, SigmaInit );
    }

    //
    // iterate
    //
    for( t = 0; t < Iterations; ++t ) {
        //
        // generate new offsprings: no recombination, only reproduction
        //
        offsprings.reproduce( parents );

	//
	// mutate object variables and adapt strategy parameter
	//
	for( i = 0; i < offsprings.size( ); ++i )
	    dynamic_cast< ChromosomeT< double >& >( offsprings[ i ][ 0 ] ).mutateDerandom( offsprings[ i ][ 1 ], derandom );

	//
	// evaluate objective function (parameters in chromosome #0)
	//
	for( i = 0; i < offsprings.size( ); ++i )
	    offsprings[ i ].setFitness( sphere( dynamic_cast< vector< double >& >( offsprings[ i ][ 0 ] ) ) );

	//
	// select (mu,lambda)
	//
	parents.selectMuLambda( offsprings );

	//
	// print out best value found so far
	//
	if( t % Interval == 0 )
	    cout << t << "\tbest value = "
		 << parents.best( ).fitnessValue( ) << endl;
    }

    return 0;
}
