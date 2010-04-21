//===========================================================================
/*!
 *  \file bfgs2_test.cpp
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
 *      $RCSfile: bfgs2_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: bfgs2_test.cpp,v $
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


// The derivation of the testfunction: 
void dtestfunc(const Array< double >& x, Array< double >& y)
{
    y = 2*x(0);
    return;
}


int main()
{
    double          fret(0.),      // minimum of the function
                    gtol(0.);      // convergence requirement
    unsigned        iter(0),       // no. of iterations performed
                                   // by the bfgs-function
                    iterMax(200),  // max. number of allowed iterations
                    i;             // counter value
    Array< double > p(1),          // initial starting point
                    xi(1);         // search direction


    cout << "Performing minimization for function f(x) = x*x + 2 "
         << "with starting point (-3,11) and direction (3,11):" << endl 
         << endl;

    // Perform minimization with all four line minimization
    // algorithms as subroutines:
    for (i = 0; i < 4; i++)
    {
        // Initialize starting point and search direction:
        p(0) = -3.;
        xi(0) = 3.;

        // Output of subroutine types:
        switch( i ) {
	    case 0 : cout << "Using linesearch as subroutine:" << endl;
	             break;
	    case 1 : cout << "Using cubelinesearch as subroutine:" << endl;
	             break;
            case 2 : cout << "Using linear minimization as subroutine:" 
                          << endl;
	             break;
            case 3 : cout << "Using derived linear minimization as "
                          << "subroutine:" << endl;
	             break;
        }
      
        // Minimization of test function:
        bfgs2(p, gtol, iter, fret, testfunc, dtestfunc, iterMax, 
             (LinMinTypes) i);

        // Output of results:
        cout << "current starting point = point where function " 
             << "takes it minimum value:" << endl;
        writeArray(p, cout);
        cout << "minimum of function:" << endl;
        cout << fret << endl;	
        cout << "number of iterations performed:" << endl;
        cout << iter << endl;
        cout << endl << endl;

    }

    return 0;
}
