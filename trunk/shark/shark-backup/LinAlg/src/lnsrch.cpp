//===========================================================================
/*!
 *  \file lnsrch.cpp
 *
 *  \brief Line search algorithm for finding a
 *         minimum value of a function.
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
 *      $RCSfile: lnsrch.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: lnsrch.cpp,v $
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
 *      Revision 1.5  2002/11/19 15:46:54  kreutz
 *      definition of min and max templates removed, include ZNminmax included instead
 *
 *      Revision 1.4  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 13:26:06  rudi
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
 *
 */
//===========================================================================

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

//===========================================================================
/*!
 *  \brief Does a line search, i.e. given a nonlinear function, a starting 
 *         point and a direction,
 *         a new point is calculated where the function has
 *         decreased "sufficiently". 
 *
 *  The line search algorithms are based on the Newton method of
 *  approximating root values of monotone decreasing
 *  functions. When the derivative \f$ f' \f$ of the function \f$ f \f$ 
 *  at a starting 
 *  point \f$ x \f$ on the X-axis can be calculated, the intersection
 *  \f$ x' \f$ of the tangent at point \f$ f(x) \f$ (with gradient 
 *  \f$ f'(x) \f$) with the X-axis can be used to get a better approximation
 *  of the minima at \f$ x_{min} \f$.
 *  
 *  This function is based on this idea: Given a nonlinear function 
 *  \f$ f \f$, a n-dimensional starting point
 *  \f$ x_{old} \f$ and a direction \f$ p \f$ (known as
 *  Newton direction), a new point \f$ x_{new} \f$ is calculated
 *  as
 *
 *  \f$
 *      x_{new} = x_{old} + \lambda p, \hspace{2em} 0 < \lambda \leq 1
 *  \f$
 *
 *  in a way that \f$ f(x_{new})\f$ has decreased sufficiently.
 *  Sufficiently means that
 *
 *  \f$
 *      f(x_{new}) \leq f(x_{old}) + \alpha \nabla f \cdot (x_{new} - x_{old})
 *  \f$
 *
 *  where \f$ 0 < \alpha < 1 \f$ is a fraction of the initial
 *  rate of decrease \f$ \nabla f \cdot p \f$.
 * 
 *  This function can be used for minimization or solving a set of
 *  nonlinear equations of the form \f$ F(x) = 0 \f$.
 *  Finding the root value (the x-value at which the related function
 *  will intersect the X-axis) will then solve the equations.
 *  In contrast to cubic line search (cblnsrch.cpp), here the
 *  root value can be calculated by only three gradient
 *  of function evaluations. 
 *
 *  \param xold   n-dimensional starting point.
 *  \param fold   The value of function \em func at point \em xold.
 *  \param g      The n-dimensional gradient of function \em func at point 
 *                \em xold.
 *  \param p      n-dimensional point to specify the search direction (the 
 *                Newton direction).
 *  \param x      New n-dimensional point along the direction \em p
 *                from \em xold where the function \em func decreases
 *                "sufficiently".
 *  \param f      The new function value for point \em x.
 *  \param stpmax No. of the steps \em lnsrch.
 *                will do internally, to limit the number
 *                avoids evaluating the function in regions where
 *                it is undefined or subject to overflow.
 *  \param check  Set to \em false on a normal
 *                exit and to \em true when \em x is to close to
 *                \em xold (signals convergence for minimization
 *                algorithms).
 *  \param func   Function to decrease.
 *  \return       none.
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em xold is not
 *         one-dimensional or has not the same size than
 *         \em g or \em p or \em x
 *
 *  \example lnsrch_test.cpp
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
 *
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
void lnsrch
(
    Array< double >& xold,
    double fold,
    Array< double >& g,
    Array< double >& p,
    Array< double >& x,
    double& f,
    double stpmax,
    bool& check,
    double (*func)( const Array< double >& )
)
{
    SIZE_CHECK( xold.ndim( ) == 1 &&
		xold.samedim( g ) &&
		xold.samedim( p ) &&
		xold.samedim( x ) )

    const double ALF  = 1.0e-4;
    const double TOLX = 1.0e-7;

    unsigned n, i;
    double a,alam,alam2(0.),alamin,b,disc,f2(0.),rhs1,rhs2,slope,sum,temp,
      test,tmplam;

    n     = xold.nelem( );
    check = false;

    for( sum = 0., i = 0; i < n; i++ )
        sum += p( i ) * p( i );
    sum = sqrt( sum );

    if( sum > stpmax )
        for( i = 0; i < n; i++ )
	    p( i ) *= stpmax / sum;

    /*
     * dot product of search direction and gradient
     */
    for( slope = 0., i = 0; i < n; i++ )
        slope += g( i ) * p( i );

    test = 0.;
    for( i = 0; i < n; i++ ) {
        temp = fabs( p( i ) ) / std::max( fabs( xold( i ) ), 1. );
	if( temp > test )
	    test = temp;
    }

    alamin = TOLX / test;
    alam   = 1.;
    for( ; ; ) {
        for( i = 0; i < n; i++ )
	    x( i ) = xold( i ) + alam * p( i );
	f = func( x );

	if( alam < alamin ) {
	    for( i = 0; i < n; i++ )
	        x( i ) = xold( i );
	    check = true;
	    return;
	} else if( f <= fold + ALF * alam * slope )
	    return;
	else {
	    if( alam == 1. )
	        tmplam = -slope / ( 2. * ( f - fold - slope ) );
	    else {
	        rhs1 = ( f - fold - alam * slope ) / ( alam * alam );
		rhs2 = ( f2 - fold - alam2 * slope ) / ( alam2 * alam2 );
		a = 3. * ( rhs1 - rhs2 ) / ( alam - alam2 );
		b = ( -alam2 * rhs1 + alam * rhs2 ) / ( alam - alam2 );
		if( a == 0. )
		    tmplam = -slope / ( 2. * b );
		else {
		    disc = b * b - a * slope;
		    if( disc < 0. ) {
                        std::cerr << "Roundoff problem in lnsrch." << std::endl;
			//ASSERT( false );
		    } else
		        tmplam = ( -b + sqrt( disc ) ) / a;
		}
		if( tmplam > 0.5 * alam )
		    tmplam = 0.5 * alam;
	    }
	}

	alam2 = alam;
	f2    = f;
	alam  = std::max( tmplam, 0.1 * alam );
    }
}
