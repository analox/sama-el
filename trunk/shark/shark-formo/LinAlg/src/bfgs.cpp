//===========================================================================
/*!
 *  \file bfgs.cpp
 *
 *  \brief Minimization of a function by using the 
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
 *      <BR>
 *
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: bfgs.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: bfgs.cpp,v $
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
 *      Revision 1.5  2002/05/16 14:15:00  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2002/02/15 13:22:37  rick
 *      adapt DataAnalysis 3.5 for Linux
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
//===========================================================================//

#include "ZNconfig.h"
#include "ZNminmax.h"

#include <cmath>
#include "LinAlg/arrayoptimize.h"

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
 *  \brief Minimizes function "func" by using the 
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
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em p is not
 *             one-dimensional 
 *
 *  \example bfgs_test.cpp  
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
 *  \sa lnsrch.cpp, cblnsrch.cpp, linmin.cpp, dlinmin.cpp, bfgs2.cpp
 */
void bfgs
(
    Array< double >& p,
    double    gtol,
    unsigned& iter,
    double&   fret,
    double  (*func)( const Array< double >& ),
    void   (*dfunc)( const Array< double >&, Array< double >& ),
    unsigned  iterMax,
    LinMinTypes linmintype
)
{
    SIZE_CHECK( p.ndim( ) == 1 )

    const double   EPS   = 1e-12; // 3.0e-8;
    const double   TOLX  = 4 * EPS;
    const double   STPMX = 100.;

    unsigned i,its,j;
    double   den,fac,fad,fae,fp,stpmax,sum,sumdg,sums,temp,test;

    unsigned        n = p.nelem( );
    Array< double > dg  ( n );
    Array< double > g   ( n );
    Array< double > Hdg ( n );
    Array< double > H   ( n, n );
    Array< double > pnew( n );
    Array< double > s   ( n );

    fp = func( p );
    dfunc( p, g );
    SIZE_CHECK( g.samedim( p ) )

    for( sum = 0., i = 0; i < n; i++ ) {
        for( H( i, i ) = 1., j = 0; j < i; j++ )
	    H( i, j ) = H( j, i ) = 0.;
	s( i ) = -g( i );
	sum += p( i ) * p( i );
    }

    stpmax = STPMX * std::max( sqrt( sum ), ( double ) n );

    for( its = 1; its <= iterMax; its++ ) {
        iter = its;

	//
	// line search
	//
	switch (linmintype) {
	  case LNSRCH : 
              bool check;
	      lnsrch( p, fp, g, s, pnew, fret, stpmax, check, func );
	      cout << "check = " << check << endl;
              break;
	  case CBLNSRCH :
	      cblnsrch( p, fp, g, s, pnew, fret, func );
              break;
	  case LINMIN :
	      pnew = p;
	      linmin( pnew, s, fret, func );
              break;
          case DERIVLINMIN :
	      pnew = p;
              dlinmin( pnew, s, fret, func, dfunc );
              break;
	  default : // DERIVLINMIN
	      pnew = p;
              dlinmin( pnew, s, fret, func, dfunc );
              break;
        }

	fp = fret;
	for( i = 0; i < n; i++ ) {
	    s( i ) = pnew( i ) - p( i );
	    p ( i ) = pnew( i );
	}

	test = 0.;
	for( i = 0; i < n; i++ ) {
	    temp = fabs( s( i ) ) / std::max( fabs( p( i ) ), 1.);
	    if( temp > test )
	        test = temp;
	}

	if( test < TOLX )
	    return;

	for( i = 0; i < n; i++ )
	    dg( i ) = g( i );

	dfunc( p, g );
	test = 0.;
	den  = std::max( fret, 1. );
	for( i = 0; i < n; i++ ) {
	    temp = fabs( g( i ) ) * std::max( fabs( p( i ) ), 1. ) / den;
	    if( temp > test )
	        test = temp;
	}

	if( test < gtol )
	    return;

	for( i = 0; i < n; i++ )
	    dg( i ) = g( i ) - dg( i );

	for( i = 0; i < n; i++ )
	    for( Hdg( i ) = 0., j = 0; j < n; j++ )
	        Hdg( i ) += H( i, j ) * dg( j );

	fac = fae = sumdg = sums = 0.;
	for( i = 0; i < n; i++ ) {
	    fac   += dg( i ) * s  ( i );
	    fae   += dg( i ) * Hdg( i );
	    sumdg += dg( i ) * dg ( i );
	    sums  += s ( i ) * s  ( i );
	}

	if( fac * fac > EPS * sumdg * sums ) {
	    fac = 1. / fac;
	    fad = 1. / fae;

	    for( i = 0; i < n; i++ )
	        dg( i ) = fac * s( i ) - fad * Hdg( i );

	    for( i = 0; i < n; i++ )
	        for( j = 1; j < n; j++ )
		    H( i, j ) += fac * s( i ) * s( j )
		      - fad * Hdg( i ) * Hdg( j ) + fae * dg( i ) * dg( j );
	}

	for( i = 0; i < n; i++ )
	    for( s( i ) = 0., j = 0; j < n; j++ )
	        s( i ) -= H( i, j ) * g( j );
    }
}



