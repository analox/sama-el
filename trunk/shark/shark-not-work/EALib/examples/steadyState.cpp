/*! steadyState.cpp
* ======================================================================
*
*  File    :  steadyState.cpp
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
*        $RCSfile: steadyState.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: steadyState.cpp,v $
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
#include "EALib/Population.h"

//=======================================================================
//
// fitness function: sphere model
//
double sphere( const std::vector< double >& x )
{
    unsigned i;
    double   sum;
    for( sum = 0., i = 0; i < x.size( ); i++ )
        sum += sqr( x[ i ] );
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
    const unsigned PopSize     = 50;
    const unsigned Dimension   = 20;
    const unsigned NumOfBits   = 10;
    const unsigned Iterations  = 10000;
    const unsigned DspInterval = 100;
    const unsigned Omega       = 5;
    const unsigned CrossPoints = 2;
    const double   CrossProb   = 0.6;
    const double   FlipProb    = 1. / ( Dimension * NumOfBits );
    const bool     UseGrayCode = true;
    const Interval RangeOfValues( -3, +5 );

    unsigned i, t;

    //
    // initialize random number generator
    //
    Rng::seed( argc > 1 ? atoi( argv[ 1 ] ) : 1234 );

    //
    // define populations
    //
    Individual kid( ChromosomeT< bool >( Dimension * NumOfBits ) );
    Population pop( PopSize, kid );

    //
    // scaling window
    //
    std::vector< double > window( Omega );

    //
    // temporary chromosome for decoding
    //
    ChromosomeT< double > dblchrom;

    //
    // minimization task
    //
    pop.setMinimize( );

    //
    // initialize all chromosomes of the population
    //
    for( i = 0; i < pop.size( ); ++i ) {
        dynamic_cast< ChromosomeT< bool >& >( pop[ i ][ 0 ] ).initialize( );
	dblchrom.decodeBinary( pop[ i ][ 0 ], RangeOfValues, NumOfBits, UseGrayCode );
	pop[ i ].setFitness( sphere( dblchrom ) );
    }

    //
    // iterate
    //
    for( t = 0; t < Iterations; ++t ) {
        //
        // scale fitness values and use proportional selection
        //
        pop.linearDynamicScaling( window, t );

        //
        // recombine by crossing over two parents
        //
        if( Rng::coinToss( CrossProb ) )
	    kid[ 0 ].crossover( pop.selectOneIndividual( )[ 0 ],
				pop.selectOneIndividual( )[ 0 ],
				CrossPoints );
	else
	    kid = pop.selectOneIndividual( );

        //
        // mutate by flipping bits
        //
	dynamic_cast< ChromosomeT< bool >& >( kid[ 0 ] ).flip( FlipProb );

        //
        // evaluate objective function
        //
	dblchrom.decodeBinary( kid[ 0 ], RangeOfValues, NumOfBits, UseGrayCode );
	kid.setFitness( sphere( dblchrom ) );

	//
	// replace the worst individual in the population
	//
	pop.worst( ) = kid;

        //
        // print out best value found so far
        //
        if( t % DspInterval == 0 )
            std::cout << t << "\tbest value = "
                 << pop.best( ).fitnessValue( ) << "\n";
    }

    return 0;
}
