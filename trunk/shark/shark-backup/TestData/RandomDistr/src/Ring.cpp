//===========================================================================
/*!
 *  \file Ring.cpp
 *
 *
 *  \author  Martin Kreutz
 *  \date    20.08.1995
 *
 *  \par Copyright (c) 1999-2003:
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
 *      TestData
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: Ring.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: Ring.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:20  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of TestData. This library is free software;
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



#ifdef __GNUC__
#pragma implementation
#endif

#include <algorithm>
#include "EALib/sqr.h"
#include "ZNconfig.h"
#include "TestData/RandomDistr/Ring.h"

//===========================================================================

Ring::Ring( double xvar,
	    double yvar,
	    double xrad,
	    double yrad )
  : uni    ( -znPiC, znPiC, RNG::globalRng ),
    gauss  ( 0, 1, RNG::globalRng ),
    xStdDev( sqrt( xvar ) ),
    yStdDev( sqrt( yvar ) ),
    xRad   ( xrad ),
    yRad   ( yrad )
{
}

Ring::Ring( double xvar,
	    double yvar,
	    double xrad,
	    double yrad,
	    RNG& rng )
  : RandomVector< double >( rng ),
    uni    ( -znPiC, znPiC, rng ),
    gauss  ( 0, 1, rng ),
    xStdDev( sqrt( xvar ) ),
    yStdDev( sqrt( yvar ) ),
    xRad   ( xrad ),
    yRad   ( yrad )
{
}

//===========================================================================

Array< double > Ring::operator( )( )
{
    Array< double > x( 2 );
    double alpha = uni( );
    x( 0 ) = gauss( xRad * cos( alpha ), xStdDev );
    x( 1 ) = gauss( yRad * sin( alpha ), yStdDev );
    return x;
}

//===========================================================================

double Ring::p( const Array< double >& x ) const
{
    unsigned i;
    double   n, sum;

    n   = ceil( znPi2C * std::max( xRad, yRad ) / std::min( xStdDev, yStdDev ) );
    sum = 0;

    for( i = unsigned( ceil( n ) ); i--; )
        sum += exp( -( sqr( ( x( 0 ) - xRad * cos( i*znPi2C/n ) ) / xStdDev ) +
		       sqr( ( x( 1 ) - yRad * cos( i*znPi2C/n ) ) / yStdDev ) ) / 2 );

    return sum / ( znPi2C * xStdDev * yStdDev );
}

//===========================================================================
