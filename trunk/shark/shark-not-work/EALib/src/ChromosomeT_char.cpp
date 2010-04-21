/*! ChromosomeT_char.cpp
* ======================================================================
*
*  File    :  ChromosomeT_char.cpp
*  Created :  17.08.1998
*
*  Copyright (c) 1995-2000 Martin Kreutz
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
*        $RCSfile: ChromosomeT_char.cpp,v $<BR>
*  
*    \par Changes:
*        $Log: ChromosomeT_char.cpp,v $
*        Revision 2.2  2004/06/17 23:08:59  saviapbe
*        An error in the INT's URL was corrected.
*
*        Revision 2.1  2004/06/17 22:53:19  saviapbe
*        The standard file header was for doxygen adapted.
*
*        Revision 2.0  2003/11/28 16:23:09  shark-admin
*        Revision tag reset to revision tag 2.x
*  
*        Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
*        INI Administration
*  <BR>
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

#include "EALib/ChromosomeT.h"

//===========================================================================

/*
#ifndef __NO_GENERIC_IOSTREAM

void ChromosomeT< char >::writeTo( ostream& os ) const
{
    os << "ChromosomeT<" << typeid( char ).name( ) << ">(" << size() << ")" << endl;
    for( unsigned i = 0; i < size( ); i++ )
        os.put( ( *this )[ i ] );
    os << endl;
}

void ChromosomeT< char >::readFrom( istream& is )
{
    string s;
  //is.getline( s );
    is >> s;
    is.get( );   // skip end of line

    if( is.good( ) &&
	s.substr( 0, 12 ) == "ChromosomeT<" &&
	s.find( '>' ) != string::npos &&
	s.substr( 12, s.find( '>' ) - 12 ) == typeid( char ).name( ) ) {

        resize( atoi( s.substr( s.find( '>' ) + 2 ).c_str( ) ) );
	for( unsigned i = 0; i < size( ); i++ )
	    is.get( ( *this )[ i ] );
    } else
        is.setf( ios::failbit );
}

#endif // !__NO_GENERIC_IOSTREAM
*/

//===========================================================================
