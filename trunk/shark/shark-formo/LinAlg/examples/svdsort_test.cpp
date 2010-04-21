//===========================================================================
/*!
 *  \file svdsort_test.cpp
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
 *      $RCSfile: svdsort_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: svdsort_test.cpp,v $
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

#include <iostream>
#include "Array/ArrayIo.h"
#include "LinAlg/linalg.h"

using namespace std;

int main() 
{
    Array2D< double > A(2,2);    // input matrix 
    Array2D< double > U(2,2);    // the two orthogonal
    Array2D< double > V(2,2);    // matrices
    Array  < double > w(2);      // vector for singular values
    unsigned          curr_row,  // currently considered matrix
                      curr_col;  // row/column

    // Initialization values for input matrix:
    double input_matrix[4] = 
    { 
         1.,  3.,  

        -4.,  3.
    };

    // Initializing matrices and vector:
    for (curr_row = 0; curr_row < 2; curr_row++)
    {
        for (curr_col = 0; curr_col < 2; curr_col++)
        {
	        A(curr_row, curr_col) = input_matrix[curr_row*2+curr_col];
	        U(curr_row, curr_col) = 0.;
	        V(curr_row, curr_col) = 0.;
        }
	  w(curr_row) = 0.;
    }

    // Singular value decomposition:
    svd(A, U, V, w);


    // Output before sorting:
    cout << "matrices and singular values before sorting:\n" << endl;
    cout << "input matrix:" << endl;
    writeArray(A, cout);
    cout << "column-orthogonal matrix U:" << endl;
    writeArray(U, cout);
    cout << "orthogonal matrix V:" << endl;
    writeArray(V, cout);
    cout << "vector of singular values:" << endl;
    writeArray(w, cout);

    // Sorting singular values, adjusting matrices:
    svdsort(U, V, w);

    // Output after sorting:
    cout << "\n\nmatrices and singular values after sorting:\n" << endl;
    cout << "column-orthogonal matrix U:" << endl;
    writeArray(U, cout);
    cout << "orthogonal matrix V:" << endl;
    writeArray(V, cout);
    cout << "vector of singular values:" << endl;
    writeArray(w, cout);

    return 0;
}
