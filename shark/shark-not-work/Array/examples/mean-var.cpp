//===========================================================================
/*!
 *  \file mean-var.cpp
 *
 *
 *  \par Copyright (c) 1999, 2003:
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
 *      Array
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: mean-var.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: mean-var.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:07  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:02  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of Array. This library is free software;
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

 
#include "Array/ArrayIo.h"
#include "Array/ArrayOp.h"

using namespace std;

int main( )
{
    unsigned i, num;
    Array< double > data, m, v;

    //
    // read the data file
    // items must be separated by white space characters
    // and records by newline characters
    //
    readArray( data, cin );

    m.resize( data.dim( 1 ) );
    v.resize( data.dim( 1 ) );
    m = 0.;
    v = 0.;
    for( num = i = 0; i < data.dim( 0 ); ++i ) {
        if( data( i, 0 ) < 10 ) {
	    m += data[ i ];
	    v += data[ i ] * data[ i ];
	    num++;
	}
    }

    m /= double( num );
    v /= double( num );
    v -= m * m;

    cout << "mean = " << m << endl;
    cout << "var  = " << v << endl;

    for( i = 0; i < m.dim( 0 ); ++i )
        cout << "$" << m( i ) << " \\pm " << sqrt( v( i ) ) << "$" << endl;

    cout << "num = " << num << endl;

    return 0;
}
