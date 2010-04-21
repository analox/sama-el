//===========================================================================
/*!
 *  \file Lorenz63.cpp
 *
 *  \author  Martin Kreutz
 *  \date    16.09.1998
 *
 *  \par Copyright (c) 1998, 2003:
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
 *      $RCSfile: Lorenz63.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: Lorenz63.cpp,v $
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

#include "TestData/TimeSeries/Lorenz63.h"

Lorenz63::Lorenz63( double sigpar,
		    double bpar,
		    double rpar,
		    double steppar )
  : RK4  ( 3, steppar ),
    sigma( sigpar ),
    b    ( bpar   ),
    r    ( rpar   )
{
    reset( );
}

void Lorenz63::derivative( const Array< double >& y,
			   Array< double >& dy ) const
{
    dy( 0 ) = -sigma * ( y( 0 ) - y( 1 ) );
    dy( 1 ) = -y( 0 ) * y( 2 ) + r * y( 0 ) - y( 1 );
    dy( 2 ) = y( 0 ) * y( 1 ) - b * y( 2 );
}

Generator< Array< double > >* Lorenz63::clone( ) const
{
    return new Lorenz63( *this );
}
