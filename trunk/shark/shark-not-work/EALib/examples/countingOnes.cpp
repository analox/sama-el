/*! countingOnes.cpp
* ======================================================================
*
*  File    :  coutingOnes.cpp
*  Created :  01.01.1995
*
*  Copyright (c) 1995,1998 Martin Kreutz
*  EMail - Martin.Kreutz@neuroinformatik.ruhr-uni-bochum.de
*
*  Institut fuer Neuroinformatik
*  Ruhr-Universitaet Bochum
*  44780 Bochum, Germany
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
*        $RCSfile: countingOnes.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: countingOnes.cpp,v $
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
*  Last update: 12.08.1998
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

#include "EALib/Population.h"

//=======================================================================
//
// fitness function: counting ones problem
//
double ones( const std::vector< bool >& x )
{
    unsigned i;
    double   sum;
    for( sum = 0., i = 0; i < x.size( ); sum += x[ i++ ] );
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
    const unsigned PopSize     = 20;
    const unsigned Dimension   = 200;
    const unsigned Iterations  = 1000;
    const unsigned Interval    = 10;
    const unsigned NElitists   = 1;
    const unsigned Omega       = 5;
    const unsigned CrossPoints = 2;
    const double   CrossProb   = 0.6;
    const double   FlipProb    = 1. / Dimension;

    unsigned i, t;

    //
    // initialize random number generator
    //
    Rng::seed( argc > 1 ? atoi( argv[ 1 ] ) : 1234 );

    //
    // define populations
    //
    Population parents   ( PopSize, ChromosomeT< bool >( Dimension ) );
    Population offsprings( PopSize, ChromosomeT< bool >( Dimension ) );

    //
    // scaling window
    //
    std::vector< double > window( Omega );

    //
    // maximization task
    //
    parents   .setMaximize( );
    offsprings.setMaximize( );

    //
    // initialize all chromosomes of parent population
    //
    for( i = 0; i < parents.size( ); ++i )
        dynamic_cast< ChromosomeT< bool >& >( parents[ i ][ 0 ] ).initialize( );

    //
    // evaluate parents (only needed for elitist strategy)
    //
    if( NElitists > 0 )
        for( i = 0; i < parents.size( ); ++i )
            parents[ i ].setFitness( ones( dynamic_cast< std::vector< bool >& >( parents[ i ][ 0 ] ) ) );

    //
    // iterate
    //
    for( t = 0; t < Iterations; ++t ) {
        //
        // recombine by crossing over two parents
        //
        offsprings = parents;
        for( i = 0; i < offsprings.size( )-1; i += 2 )
            if( Rng::coinToss( CrossProb ) )
                offsprings[ i ][ 0 ].crossover( offsprings[ i+1 ][ 0 ],
						CrossPoints );

        //
        // mutate by flipping bits
        //
        for( i = 0; i < offsprings.size( ); ++i )
            dynamic_cast< ChromosomeT< bool >& >( offsprings[ i ][ 0 ] ).flip( FlipProb );

        //
        // evaluate objective function
        //
        for( i = 0; i < offsprings.size( ); ++i )
            offsprings[ i ].setFitness( ones( dynamic_cast< std::vector< bool >& >( offsprings[ i ][ 0 ] ) ) );

        //
        // scale fitness values and use proportional selection
        //
        offsprings.linearDynamicScaling( window, t );
	parents.selectProportional( offsprings, NElitists );

        //
        // print out best value found so far
        //
        if( t % Interval == 0 )
            std::cout << t << "\tbest value = "
                 << parents.best( ).fitnessValue( ) << "\n";
    }

    return 0;
}
