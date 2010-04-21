//===========================================================================
/*!
 *  \file SelectComponent.h
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
 *      $RCSfile: SelectComponent.h,v $<BR>
 *
 *  \par Changes:
 *      $Log: SelectComponent.h,v $
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:07  shark-admin
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


#ifndef __cplusplus
#error Must use C++.
#endif

#ifndef __SELECTCOMPONENT_H
#define __SELECTCOMPONENT_H

#include "Array/Array.h"
#include "TestData/TimeSeries/Generator.h"

//===========================================================================

template < class T >
class SelectComponent : public Generator< T >
{
  public:
    SelectComponent( const Generator< Array< T > >& gen, unsigned n = 0 )
      : generator( gen.clone( ) ), num( n )
    {
    }

    SelectComponent( const SelectComponent< T >& s )
      : Generator< T >( s ),
        generator( s.generator->clone( ) ),
        num      ( s.num )
    {
    }

    ~SelectComponent( )
    {
        delete generator;
    }

    void reset( )
    {
        generator->reset( );
    }

    T operator ( ) ( )
    {
        x = ( *generator )( );
	SIZE_CHECK( x.ndim( ) == 1 )
	RANGE_CHECK( num < x.dim( 0 ) )
	return x( num );
    }

  private:
    Generator< Array< T > >* generator;
    unsigned num;
    Array< T > x;

    Generator< T >* clone( ) const
    {
        return new SelectComponent< T >( *this );
    }
};

//===========================================================================

#endif /* !__SELECTCOMPONENT_H */
