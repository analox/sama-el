//===========================================================================
/*!
 *  \file bfgs2.cpp
 *
 *  \brief Minimization of a function by using a modified 
 *         Broyden-Fletcher-Goldfarb-Shanno-algorithm.
 *         
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Copyright (c) 1998-2000:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 
 
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: bfgs2.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: bfgs2.cpp,v $
 *      Revision 2.1  2004/06/16 15:22:20  saviapbe
 *
 *      Some bugs in the doxygen documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.6  2002/11/19 15:46:54  kreutz
 *      definition of min and max templates removed, include ZNminmax included instead
 *
 *      Revision 1.5  2002/05/16 13:55:04  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2001/11/30 13:25:45  rudi
 *      Doxygen comments added.
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
 *  The algorithm was taken from the library:
 *
 *  ============================================================<BR>
 *  COOOL           version 1.1           ---     Nov,  1995    <BR>
 *  Center for Wave Phenomena, Colorado School of Mines         <BR>
 *  ============================================================<BR>
 *
 *  This code is part of a preliminary release of COOOL (CWP
 *  Object-Oriented Optimization Library) and associated class 
 *  libraries. 
 *
 *  The COOOL library is a free software. You can do anything you want
 *  with it including make a fortune.  However, neither the authors,
 *  the Center for Wave Phenomena, nor anyone else you can think of
 *  makes any guarantees about anything in this package or any aspect
 *  of its functionality.
 *
 *  Since you've got the source code you can also modify the
 *  library to suit your own purposes. We would appreciate it 
 *  if the headers that identify the authors are kept in the 
 *  source code.
 *
 *  ================================<BR>
 *  author:  H. Lydia Deng, 06/17/96<BR>
 *  ================================<BR>
 */
//===========================================================================//


#include "ZNconfig.h"
#include "ZNminmax.h"

#include <cmath>
#include "LinAlg/arrayoptimize.h"

#define DLINMIN
//#define LINMIN

/*
#ifdef _WIN32
#ifndef __MIN_MAX__
#define __MIN_MAX__
namespace std {
//
// undefine macros min and max to avoid conflicts with template names
//
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

//
// template functions 'min' and 'max' are not defined for _WIN32
// due to name conflicts
//
template < class T > inline T min( T a, T b ) { return a < b ? a : b; }
template < class T > inline T max( T a, T b ) { return a > b ? a : b; }
}
#endif
#endif
*/

using namespace std;

//===========================================================================
/*!
 *  \brief Minimizes function "func" by using a modified
 *         Broyden-Fletcher-Goldfarb-Shanno-algorithm.
 *
 *  Given a function \em func of \em N variables, from which
 *  the gradient or first partial derivatives at arbitrary
 *  points can be computed and which can be
 *  locally approximated by the quadratic form of equation
 *
 *  \f$
 *      f(x) \approx  c - b x + \frac{1}{2} x A x
 *  \f$
 *
 *  then the number of unknown parameters in \f$ f \f$ is equal
 *  to the number of free parameters in \f$ A \f$ and \f$ b \f$,
 *  which is \f$ \frac{1}{2} N (N + 1) \f$.
 *  To gather information about these parameters, function evaluations
 *  and line minimizations are used.
 *  \em BFGS as a variable metric method builds up, iteratively,
 *  a good approximation to the inverse Hessian matrix \f$ A^{-1} \f$,
 *  that is, to construct a sequence of matrices \f$ H_i \f$ with
 *  the property,
 *
 *  \f$
 *      \lim_{i\to\infty} H_i = A^{-1}
 *  \f$
 *
 *  The number of iterations can be restricted.
 * 
 *      \param p          N-dimensional initial starting point.
 *      \param gtol       Convergence requirement on zeoring the gradient.
 *      \param iter       Number of iterations that were
 *                        performed by the function.
 *      \param fret       Minimum value of the function.
 *      \param func       The function that will be minimized.
 *      \param dfunc      The derivation of \em func.
 *      \param iterMax    Max. number of allowed iterations.
 *      \param linmintype Type of algorithm used for line minimization,
 *                        possible values are \em LNSRCH
 *                        for #lnsrch, \em CBLNSRCH
 *                        for #cblnsrch,
 *                        \em LINMIN for #linmin and
 *                        \em DERIVLINMIN for #dlinmin .
 *      \return           none.
 *
 *  \example bfgs2_test.cpp
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
 *
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *  
 *  \sa lnsrch.cpp, cblnsrch.cpp, linmin.cpp, dlinmin.cpp, bfgs.cpp
 */
void bfgs2( Array< double >& p,
	    double gtol,
	    unsigned& iter,
	    double&   fret,
	    double  (*func)( const Array< double >& ),
	    void   (*dfunc)( const Array< double >&, Array< double >& ),
	    unsigned iterMax,
            LinMinTypes linmintype )
{

    const double   Lambda  = 0.025;
    const double   STPMX = 100.;

    bool check(false);
    unsigned i, j, n = p.nelem( );
    double   d, err, scale, sum, stpmax;
    double   test, temp;

    Array< double > p1( n );
    Array< double > g0( n );
    Array< double > g1( n );
    Array< double > s ( n );
    Array< double > g ( n );

    Array< double > H ( n, n );

    Array< double > gamma( n );
    Array< double > delta( n );

    iter = 0;

    dfunc( p, g0 );

    for( sum = 0., i = 0; i < n; i++ ) {
        for( H( i, i ) = 1., j = 0; j < i; j++ )
	    H( i, j ) = H( j, i ) = 0.;
	s( i ) = -g( i );
	sum += p( i ) * p( i );
    }

    stpmax = STPMX * std::max( sqrt( sum ), ( double ) n );

    for( err = 0., i = 0; i < n; ++i )
        err += g0( i ) * g0( i );
    if( err < gtol*gtol ) {
        fret = func( p );
	return;
    }

    for( i = 0; i < n; ++i ) {
        for( j = 0; j < n; ++j )
	    H( i, j ) = 0.;
	H( i, i ) = 1.;
    }

    for( ; ; ) {
        for( i = 0; i < n; ++i )
	    for( s( i ) = 0., j = 0; j < n; ++j )
	        s( i ) -= H( i, j ) * g0( j );

	//
	// line search
	//
        switch (linmintype) 
	{
	    case LNSRCH :
	        lnsrch( p, func( p ), g0, s, p1, fret, stpmax, check, func );
	        cout << "check = " << check << endl;
                break;
	    case CBLNSRCH :
	        cblnsrch( p, func( p ), g0, s, p1, fret, func, Lambda );
                break;
	    case LINMIN :
	        p1 = p;
	        linmin( p1, s, fret, func );
                break;
	    case DERIVLINMIN :
	        p1 = p;
	        dlinmin( p1, s, fret, func, dfunc );
                break;
	    default :
                p1 = p;
	        dlinmin( p1, s, fret, func, dfunc );
                break;
        }

	//cblnsrch( p, func( p ), g0, s, p1, fret, func, Lambda );
	dfunc( p1, g1 );
	for( err = 0., i = 0; i < n; ++i )
	    err += g1( i ) * g1( i );
	if( ++iter >= iterMax || err <= gtol*gtol )
	    break;

	//
	// test for convergence
	//
	test = 0.;
	for( i = 0; i < n; i++ ) {
	    temp = fabs( p1( i ) - p( i ) ) / std::max( fabs( p( i ) ), 1.);
	    if( temp > test )
	        test = temp;
	}
	if( test < gtol )
	    return;

        for( d = 0., i = 0; i < n; ++i )
	    d += ( gamma(i) = g1(i)-g0(i) ) * ( delta(i) = p1(i)-p(i) );

	for( i = 0; i < n; ++i )
	    for( s( i ) = 0., j = 0; j < n; ++j )
	        s( i ) += H( i, j ) * gamma( j );

	if( d < 1e-8 ) {
	    for( i = 0; i < n; ++i ) {
	        for( j = 0; j < n; ++j )
		    H( i, j ) = 0.;
		H( i, i ) = 1.;
	    }
	} else {
	    for( scale = 0., i = 0; i < n; ++i )
	        scale += gamma( i ) * s( i );
	    scale = ( scale / d + 1 ) / d;

	    for( i = 0; i < n; ++i ) {
	        g0( i ) = g1( i );
		p ( i ) = p1( i );
	        for( j = 0; j < n; ++j )
		    H( i, j ) += scale * delta(i) * delta(j)
		      - ( s(i) * delta(j) + s(j) * delta(i) ) / d;
	    }
	}
    }

    for( i = 0; i < n; ++i )
        p( i ) = p1( i );
}





