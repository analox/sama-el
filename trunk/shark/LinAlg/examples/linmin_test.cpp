//===========================================================================
/*!
 *  \file linmin_test.cpp
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
 *      $RCSfile: linmin_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: linmin_test.cpp,v $
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
#include "LinAlg/arrayoptimize.h"


using namespace std;


// The testfunction: 
double testfunc(const Array< double >& x)
{
    return x(0) * x(0) + 2;
}


int main()
{
    double          fret(0.); // function value at the point that
                              // is found by the function
    Array< double > p(1),     // initial search starting point
                    xi(1);    // direction of search

    // Initialize search starting point and direction:
    p(0) = -3.;
    xi(0) = 3.;

    cout << "Performing a line minimization for function f(x) = x*x + 2 "
         << "with starting point (-3, 11) and direction (3, 11):" 
         << endl << endl;

    // Minimize function:
    linmin(p, xi, fret, testfunc);

    // Output of results:
    cout << "x-value of final point:           ";
    writeArray(p, cout);
    cout << "y-value of final point = minimum: ";
    cout << fret << endl;	
    cout << "starting point has moved to the   " << endl
         << "predefined direction by:          ";
    writeArray(xi, cout);

    return 0;
}
