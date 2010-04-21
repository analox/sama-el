//===========================================================================
/*!
 *  \file lnrsrch_test.cpp
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
 *      $RCSfile: lnsrch_test.cpp,v $<BR>
 *
 *  \par Changes:
 *      $Log: lnsrch_test.cpp,v $
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

#include<iostream>
#include "Array/Array.h"
#include "Array/ArrayIo.h"
#include "LinAlg/arrayoptimize.h"

using namespace std;


// Function to be decreased, implements f(x) = x*x
//
double my_func(const Array< double > &x)
{
    return x(0) * x(0);
}


// Derivative function f'(x) = 2x of f(x) = x*x
//
double gradient(double x)
{
    return 2*x;
}


int main()
{

    Array< double > x_old(1),       // old point
                    x_new(1),       // new point
                    p(1),           // Newton direction 
                    grad(1);        // gradient of function at old point
    double          f_old(0.),      // function-value of old point 
                    f_new(0.),      // function-value of new point
                    no_iter(1.);    // no. of iterations, that
                                    // "lnsrch" will do internally
    bool            status(false);  // status of "lnsrch"-function
    unsigned        i;              // number of "lnsrch" iterations
 

    // input points for "my_func":
    //
    double x_values[19] = 
    {
          -9., -8., -7., -6., 
          -5., -4., -3., -2., -1.,  
	   0.,  1.,  2.,  3.,  4.,
	   5.,  6.,  7.,  8.,  9.
    };

    // Initialize values for first call of "lnsrch":
    //
    x_old(0) = x_values[0];       // starting point is -10.0 
    f_old = my_func(x_old);       // Calculate function value of starting point
    grad(0) = gradient(x_old(0)); // Calculate gradient of starting point
    p(0) = 10.;                   // decrease function by going to
                                  // the "right" on the X-axis.
    x_new(0) = x_old(0);          // No new values calculated yet 
    f_new = f_old;

    cout << "Performing line search for function f(x) = x*x "
         << "with starting point (-10, 100) and direction (10, 100):" 
         << endl << endl;

    // Head for output of results:
    cout << "status:\t\tnew point:\tfunction-value:\n";
    cout << "-------------------------------------------------------------\n";
 
    // Presettings for formatted output:
    cout.setf(ios::fixed | ios::right | ios::internal | ios::showpos);
    cout.precision(2);
    cout.width(5);

    // Output of starting values:
    cout << "false" << "\t\t" << x_new(0) << "\t\t" << f_new << endl;     

    // Take 15 "lnsrch" iterations:
    for (i = 0; i < 14; i++) 
    {
        // Calculate next new point:
        lnsrch(x_old, f_old, grad, p, x_new, f_new, no_iter, status,my_func);

        // Output of calculated results:
        if (status == 0) cout << "false"; else cout << "true";
	    cout << "\t\t" << x_new(0) << "\t\t" << f_new << endl;     

        // Setting values for next iteration:
        x_old(0) = x_new(0);
        f_old = f_new;
        grad(0) = gradient(x_old(0));
    }

    return 0;
}
