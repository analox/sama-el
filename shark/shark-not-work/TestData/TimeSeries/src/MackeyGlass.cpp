//===========================================================================
/*!
 *  \file MackeyGlass.cpp
 *
 *
 *  \author  Martin Kreutz
 *  \date    16.09.1998
 *
 *  \par Copyright (c) 1998-2003:
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
 *      $RCSfile: MackeyGlass.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: MackeyGlass.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
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

#include <cmath>
#include "Array/ArrayOp.h"
#include "TestData/TimeSeries/MackeyGlass.h"

MackeyGlass::MackeyGlass( double   a,
			  double   b,
			  double   c,
			  unsigned tau,
			  double   stepsize,
			  unsigned subsmpl,
			  double   init )
  : RK4_1D( stepsize ),
    A( a ), B( b ), C( c ), T( tau ), subsample( subsmpl ),
    initial( init ), x( tau )
{
    reset( );
}

void MackeyGlass::reset( )
{
    t = T;
    x = initial;
    y = initial;
}

double MackeyGlass::operator ( ) ( )
{
//
// der KCC optimiert die Variable index weg, dies fuehrt dazu, dass
// der alte Wert von t vergessen wird...
//
#ifdef __CRAY_T3E__
    volatile
#endif
    unsigned index;
    unsigned i = 0;
    do{
        index = t % T;
        x( index ) = RK4_1D::operator( )( );
    } while( ++i < subsample );
    return x( index );
}

double MackeyGlass::derivative( double y ) const
{
    double xT = x( t % T );
    return A * xT / ( 1 + pow( xT, B ) ) - C * y;
}

Generator< double >* MackeyGlass::clone( ) const
{
    return new MackeyGlass( *this );
}
