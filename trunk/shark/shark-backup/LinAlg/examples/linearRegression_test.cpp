//===========================================================================
/*!
 *  \file linearRegression_test.cpp
 *
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
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: linearRegression_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: linearRegression_test.cpp,v $
 *      Revision 2.1  2005/10/11 13:56:57  christian_igel
 *      added ()
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *
 *
 *  This file is part of LinAlg. This library is free software;
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
#include "LinAlg/LinearRegression.h"

using namespace std;

// Data for the x-vectors.
double data_x[8][2] = {
  {0.5, 1.0},
  {1.5, 1.0},
  {2.0, 1.0},
  {2.5, 1.0},
  {3.0, 1.0},
  {3.0, 1.0},
  {3.5, 1.0},
  {4.0, 1.0}
};

// Data for the corresponding y-vectors.
double data_y[8][2] = {
  {0.5, 1.0},
  {1.0, 1.0},
  {0.5, 1.0},
  {1.5, 1.0},
  {0.5, 1.0},
  {1.0, 1.0},
  {2.0, 1.0},
  {1.5, 1.0}
};

// Data for x-vectors for the model test.
double data_test[3][2] = {
  {2.5, 1.0},
  {1.0, 1.0},
  {1.75, 1.0}
};

 
int main()
{
    Array< double >   x( 8, 2 ),      // Matrix of all x-vectors.
                      y( 8, 2 ),      // Matrix of all y-vectors.
                      x_test( 3, 2 ), // Matrix of all x-vectors used
                                      // for testing the model.
                      yout,           // Current y-vector-approximation
                                      // of the model.
                      A,              // Regression coefficient A.
                      b;              // Regression coefficient b.
    unsigned          i, j;           // Counter variables.


    // Copy data of all x/y-vectors to matrices:
    for ( i = 0; i < 8; i++ ) 
    {
        for (j = 0; j < 2; j++ ) {
            x( i, j ) = data_x[ i ][ j ];
            y( i, j ) = data_y[ i ][ j ];
            if ( i < 3 ) x_test( i, j ) = data_test[ i ][ j ];
        }
    }

    // Create empty regression model:
    LinearRegression model;

    // Add all x/y-vector-pairs to the model:
    model.train( x, y );

    // Extract regression coefficients of the final model:
    A = model.A( );
    b = model.b( );
    cout << "A = ";
    writeArray( A, cout );
    cout << "b = ";
    writeArray( b, cout );
    cout << "\n\n";

    // Test y-vector-approximations for all x-vectors
    // given by the model without these x/y-pairs:
    for ( i = 0; i < 8; i++ )
    {
        model.leaveOneOut( x[ i ], y[ i ], yout );
        cout << "data vector x:                         ";
        writeArray( x[ i ], cout );
        cout << "corresponding data vector y:           ";
        writeArray( y[ i ], cout );
        cout << "model approximation for data vector y: ";
        writeArray( yout, cout );
        cout << endl;
    }

    // Testing final regression model for some
    // vector pairs:
    cout << "\n\n";
    for ( i = 0; i < 3; i++ )
    {
        model.recall( x_test[ i ], yout );
        cout << "Given x-vector:                      ";
        writeArray( x_test[ i ], cout );
        cout << "y-vector-approximation of the model: ";
        writeArray( yout, cout );
        cout << endl;
    }

    return 0;
}





