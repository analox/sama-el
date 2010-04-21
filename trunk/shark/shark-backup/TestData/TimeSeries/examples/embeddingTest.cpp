//===========================================================================
/*!
 *  \file embeddingTest.cpp
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
 *      $RCSfile: embeddingTest.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: embeddingTest.cpp,v $
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


#include <iostream>
#include "TestData/TimeSeries/Counter.h"
#include "TestData/TimeSeries/Embedding.h"
#include "TestData/TimeSeries/IOSamples.h"
#include "Array/ArrayIo.h"
#include "Array/ArrayOp.h"

using namespace std;

int main( )
{
    unsigned i;
    Embedding< int > embed( Counter< int >( 0, 1 ), 4, 3, 10 );
    IOSamples< int > iosam( Counter< int >( 0, 1 ), 3, 3, 2, 2 );

    for( i = 0; i < 10; ++i )
        cout << embed( ) << endl;


    for( i = 0; i < 10; ++i ) {
        Array< int > in, out;
	iosam( in, out );
        cout << in << out << endl;
    }

    return 0;
}
