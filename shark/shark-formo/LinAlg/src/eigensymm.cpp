//===========================================================================
/*!
 *  \file eigensymm.cpp
 *
 *  \brief Used to calculate the eigenvalues and the eigenvectors of a
 *         symmetric matrix.
 *
 *  Here the eigenvectors and eigenvalues are calculated by using
 *  Givens and Householder reduction of the matrix to tridiagonal form.
 *
 *  \author  M. Kreutz
 *  \date    1998
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
 *      $RCSfile: eigensymm.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: eigensymm.cpp,v $
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
 *
 */
//===========================================================================

#include "ZNconfig.h"
#include <cmath>
#include "LinAlg/linalg.h"
//===========================================================================
/*!
 *  \brief Calculates the eigenvalues and the normalized
 *   eigenvectors of a symmetric matrix "amatA" using the Givens
 *   and Householder reduction.  
 *
 *  Given a symmetric \f$ n \times n \f$ matrix \em A, this function 
 *  calculates the eigenvalues \f$ \lambda \f$ and the eigenvectors \em x, 
 *  defined as
 *
 *  \f$
 *      Ax = \lambda x
 *  \f$
 *
 *  where \em x is a one-column matrix and the matrix multiplication
 *  is used for \em A and \em x.
 *  Here, the Givens reduction as a modification of the Jacobi method
 *  (s.a. #eigensymmJacobi) is used. Instead of trying to reduce the
 *  matrix all the way to diagonal form, we are content to stop
 *  when the matrix is tridiagonal. This allows the function
 *  to be carried out in a finite number of steps, unlike the
 *  Jacobi method, which requires iteration to convergence.
 *  So in comparison to the Jacobi method, this function is
 *  faster for matrices with an order greater than 10.
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which must be symmetric, so
 *                    only the bottom triangular matrix must contain values.
 *                    Values below the diagonal will be destroyed. The method 
 *                    uses this matrix as a buffer for intermediate results.
 *      \param	vmatA \f$ n \times n \f$ matrix with the calculated normalized 
 *                    eigenvectors, each column will contain an
 *                    eigenvector.
 *	\param  dvecA n-dimensional vector with the calculated 
 *                    eigenvalues in descending order.
 *      \return       none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em amatA is not a square matrix
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *     previously 'eigensymm', renamed by S. Wiegand 2003/10/01 
 *
 *  \par Status
 *      stable
 *
 */
void eigensymm_obsolete
(
	const Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK( amatA.dim( 0 ) == amatA.dim( 1 ) )
	vmatA.resize( amatA );
	dvecA.resize( amatA.dim( 0 ) );

	const unsigned maxIterC = 50;

	double **const x = amatA.ptrArr( );
	double **a = vmatA.ptrArr( );
	double *d = dvecA.elemvec( );
	unsigned n = dvecA.nelem( );

    unsigned i, j, k, l, m;
    double   b, c, f, g, h, hh, p, r, s, scale, *e, tmp;

    //
    // special case n = 1
    //
    if( n == 1 )
	{
        a[ 0 ][ 0 ] = 1;
		d[ 0 ] = x[ 0 ][ 0 ];
		return;
    }

    //
    // copy matrix
    //
    for( i = 0; i < n; i++ )
	{
        for( j = 0; j <= i; j++ )
		{
			a [i] [j] = x [i] [j];
		}
	}
    tmp = x[ n-1 ][ n-1 ];
    e   = x[ n-1 ];

    //
    // reduction to tridiagonal form
    //
    for( i = n; i-- > 1; )
    {
        h = scale = 0.0;

        if( i > 1 )
		{
            // scale row
            for( k = 0; k < i; k++ )
			{
                scale += fabs( a [i] [k] );
			}
		}

        if( scale == 0.0 )
		{
            e [i] = a [i] [i-1];
		}
        else
        {
            for( k = 0; k < i; k++ )
            {
                a [i] [k] /= scale;
                h += a [i] [k] * a [i] [k];
            }

            f           = a [i] [i-1];
            g           = f > 0 ? -sqrt( h ) : sqrt( h );
            e [i]       = scale * g;
            h          -= f * g;
            a [i] [i-1] = f - g;
            f           = 0.0;

            for( j = 0; j < i; j++ )
            {
                a [j] [i] = a [i] [j] / (scale * h );
                g = 0.0;

                // form element of a*u
                for( k = 0; k <= j; k++ )
				{
                    g += a [j] [k] * a [i] [k];
				}

                for( k = j + 1; k < i; k++ )
				{
                    g += a [k] [j] * a [i] [k];
				}

                // form element of p
                f += ( e [j] = g / h ) * a [i] [j];
            }

            hh = f / (h + h );

            // form reduced a
            for( j = 0; j < i; j++ )
            {
                f     = a [i] [j];
                g     = e [j] - hh * f;
                e [j] = g;

                for( k = 0; k <= j; k++ )
				{
                    a [j] [k] -= f * e [k] + g * a [i] [k];
				}
            }

            for( k = i; k--; )
			{
                a [i] [k] *= scale;
			}
        }

        d [i] = h;
    }

    d [0] = e [0] = 0.0;

    // accumulation of transformation matrices
    for( i = 0; i < n; i++ )
    {
        if( d [i] )
        {
            for( j = 0; j < i; j++ )
            {
                g = 0.0;

                for( k = 0; k < i; k++ )
				{
                    g += a [i] [k] * a [k] [j];
				}

                for( k = 0; k < i; k++ )
				{
                    a [k] [j] -= g * a [k] [i];
				}
            }
        }

        d [i]     = a [i] [i];
        a [i] [i] = 1.0;

        for( j = 0; j < i; j++ )
		{
            a [i] [j] = a [j] [i] = 0.0;
		}
    }

	//
	// eigenvalues from tridiagonal form
	//
    if( n <= 1 )
	{
		return;
	}

    for( i = 1; i < n; i++ )
	{
        e [i-1] = e [i];
	}

    e [n-1] = 0.0;

    for( l = 0; l < n; l++ )
    {
        j = 0;

        do
        {
            // look for small sub-diagonal element
            for( m = l; m < n-1; m++ )
            {
                s = fabs( d [m] ) + fabs( d [m+1] );
                if( fabs( e [m] ) + s == s )
				{
					break;
				}
            }

            p = d [l];

            if( m != l )
            {
                if( j++ == maxIterC )
				{
					throw 1;
				}

                // form shift
                g = (d [l+1] - p ) / (2.0 * e [l] );
                r = sqrt( g * g + 1.0 );
                g = d [m] - p + e [l] / (g + ((g ) > 0 ? fabs(r ) : -fabs(r )));
                s = c = 1.0;
                p = 0.0;

                for( i = m; i-- > l; )
                {
                    f = s * e [i];
                    b = c * e [i];

                    if( fabs( f ) >= fabs( g ) )
                    {
                        c       = g / f;
                        r       = sqrt( c * c + 1.0 );
                        e [i+1] = f * r;
                        s       = 1.0 / r;
                        c      *= s;
                    }
                    else
                    {
                        s       = f / g;
                        r       = sqrt( s * s + 1.0 );
                        e [i+1] = g * r;
                        c       = 1.0 / r;
                        s      *= c;
                    }

                    g       = d [i+1] - p;
                    r       = (d [i] - g ) * s + 2.0 * c * b;
                    p       = s * r;
                    d [i+1] = g + p;
                    g       = c * r - b;

                    // form vector
                    for( k = 0; k < n; k++ )
                    {
                        f           = a [k] [i+1];
                        a [k] [i+1] = s * a [k] [i] + c * f;
                        a [k] [i  ] = c * a [k] [i] - s * f;
                    }
                }

                d [l] -= p;
                e [l]  = g;
                e [m]  = 0.0;
            }
        }
        while( m != l );
    }

	//
	// sorting eigenvalues
	//
    eigensort( vmatA, dvecA );

	//
	// normalizing eigenvectors
	//
    for( j = n; j--; )
    {
        s = 0.0;
        for( i = n; i--; )
		{
			s += a [i] [j] * a [i] [j];
		}
        s = sqrt( s );

        for( i = n; i--; )
		{
			a [i] [j] /= s;
		}
    }

    x [n-1] [n-1] = tmp;
}

//===========================================================================
/*!
 *  \brief Used as frontend for 
 *  #eigensymm_obsolete(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA) 
 *  for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" 
 *  using the Givens and Householder reduction.
 *
 *  Frontend for function 
 *  #eigensymm_obsolete(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA), 
 *  when using type \em Array instead of type \em Array2D.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom triangular matrix must contain values.
 *                Values below the diagonal will be destroyed.The method 
 *                uses this matrix as a buffer for intermediate results.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column will contain one eigenvector.
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order.
 *      \return   none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em A is not a square matrix
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *     previously 'eigensymm', renamed by S. Wiegand 2003/10/01 
 *
 *  \par Status
 *      stable
 *
 */
void eigensymm_obsolete
(
    const Array< double >& A,
          Array< double >& G,
	  Array< double >& l
)
{
    SIZE_CHECK( A.ndim( ) == 2 && A.dim( 0 ) == A.dim( 1 ) )

	G.resize( A );
	l.resize( A.dim( 0 ) );

	Array2DReference< double > amatL( const_cast< Array< double >& >( A ) );
	Array2DReference< double > vmatL( G );

    eigensymm_obsolete( amatL, vmatL, l );
}


//===========================================================================
/*!
 *  \brief Used as frontend for 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA) 
 *  for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix 
 *  'amatA' using the Givens and Householder reduction. Each time this frontend is called additional 
 *  memory is allocated for intermediate results.
 *
 *  Frontend for function 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA), 
 *  when memory for 'odvecA' should be allocated by the call of this frontend.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom triangular matrix must contain values.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column will contain one eigenvector.
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order.
 *      \return   none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em A is not a square matrix
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */

void eigensymm
(
    const Array2D< double >& A,
          Array2D< double >& G,
	  Array< double >& l
)
{
    SIZE_CHECK( A.dim( 0 ) == A.dim( 1 ) )

    G.resize( A );
    l.resize( A.dim( 0 ) );

    Array<double> od;
    od.resize( A.dim( 0 ) );

    Array2DReference< double > amatL( const_cast< Array2D< double >& >( A ) );
    Array2DReference< double > vmatL( G );

    eigensymm( amatL, vmatL, l, od);
}


//===========================================================================
/*!
 *  \brief Used as frontend for 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA) 
 *  for calculating the eigenvalues and the normalized eigenvectors of a symmetric 
 *  matrix "amatA" using the Givens and Householder reduction. Each time this frontend is 
 *  called additional memory is allocated for intermediate results.
 *
 *  Frontend for function
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA), 
 *  when using type \em Array instead of type \em Array2D. Each time this frontend is 
 *  called additional memory is allocated for intermediate results.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom
 *                triangular matrix must contain values.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column will contain one
 *                eigenvector.
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order.
 *
 *      \return   none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em A is not a square matrix
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
void eigensymm
(
    const Array< double >& A,
	Array< double >& G,
	Array< double >& l
)
{
    SIZE_CHECK( A.ndim( ) == 2 && A.dim( 0 ) == A.dim( 1 ) )

	G.resize( A );
	l.resize( A.dim( 0 ) );

	Array2DReference< double > amatL( const_cast< Array< double >& >( A ) );
	Array2DReference< double > vmatL( G );

    eigensymm( amatL, vmatL, l );
}


//===========================================================================
/*!
 *  \brief Calculates the eigenvalues and the normalized
 *  eigenvectors of a symmetric matrix "amatA" using the Givens
 *  and Householder reduction without corrupting  "amatA" during application.  
 *
 *  Given a symmetric \f$ n \times n \f$ matrix \em A, this function 
 *  calculates the eigenvalues \f$ \lambda \f$ and the eigenvectors \em x, 
 *  defined as
 *
 *  \f$
 *      Ax = \lambda x
 *  \f$
 *
 *  where \em x is a one-column matrix and the matrix multiplication
 *  is used for \em A and \em x.
 *  Here, the Givens reduction as a modification of the Jacobi method
 *  (s.a. #eigensymmJacobi) is used. Instead of trying to reduce the
 *  matrix all the way to diagonal form, we are content to stop
 *  when the matrix is tridiagonal. This allows the function
 *  to be carried out in a finite number of steps, unlike the
 *  Jacobi method, which requires iteration to convergence.
 *  So in comparison to the Jacobi method, this function is
 *  faster for matrices with an order greater than 10.
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which must be symmetric, so
 *                    only the bottom triangular matrix must contain values.
 *      \param	vmatA \f$ n \times n \f$ matrix with the calculated normalized 
 *                    eigenvectors, each column will contain an
 *                    eigenvector.
 *	\param  dvecA n-dimensional vector with the calculated 
 *                    eigenvalues in descending order.
 *	\param  odvecA n-dimensional vector with the calculated 
 *                     offdiagonal of the Householder transformation.
 *
 *      \return       none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em amatA is not a square matrix
 *
 *  \example eigensymm_test.cpp  
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *      2003/09/19 18:55:22  S. Wiegand
 *      lower triangular part of 'amatA' will not corrupted anymore
 *
 *  \par Status
 *      stable
 *
 */
void eigensymm
(
	const Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA,
	Array  < double >& odvecA
)
{
    SIZE_CHECK( amatA.dim( 0 ) == amatA.dim( 1 ) )

    vmatA.resize( amatA );
    dvecA.resize( amatA.dim( 0 ) );
    odvecA.resize( amatA.dim( 0 ) );

    const unsigned maxIterC = 50;

    double **const x = amatA.ptrArr( );
    double **a = vmatA.ptrArr( );
    double *e = odvecA.elemvec( );
    double *d = dvecA.elemvec( );
    unsigned n = dvecA.nelem( );

    unsigned i, j, k, l, m;
    double   b, c, f, g, h, hh, p, r, s, scale;

    //
    // special case n = 1
    //
    if( n == 1 )
	{
        a[ 0 ][ 0 ] = 1;
		d[ 0 ] = x[ 0 ][ 0 ];
		return;
    }

    //
    // copy matrix
    //
    for( i = 0; i < n; i++ )
	{
        for( j = 0; j <= i; j++ )
		{
			a [i] [j] = x [i] [j];
		}
	}
 
    //
    // reduction to tridiagonal form
    //
    for( i = n; i-- > 1; )
    {
        h = scale = 0.0;

        if( i > 1 )
		{
            // scale row
            for( k = 0; k < i; k++ )
			{
                scale += fabs( a [i] [k] );
			}
		}

        if( scale == 0.0 )
		{
            e [i] = a [i] [i-1];
		}
        else
        {
            for( k = 0; k < i; k++ )
            {
                a [i] [k] /= scale;
                h += a [i] [k] * a [i] [k];
            }

            f           = a [i] [i-1];
            g           = f > 0 ? -sqrt( h ) : sqrt( h );
            e [i]       = scale * g;
            h          -= f * g;
            a [i] [i-1] = f - g;
            f           = 0.0;

            for( j = 0; j < i; j++ )
            {
                a [j] [i] = a [i] [j] / (scale * h );
                g = 0.0;

                // form element of a*u
                for( k = 0; k <= j; k++ )
				{
                    g += a [j] [k] * a [i] [k];
				}

                for( k = j + 1; k < i; k++ )
				{
                    g += a [k] [j] * a [i] [k];
				}

                // form element of p
                f += ( e [j] = g / h ) * a [i] [j];
            }

            hh = f / (h + h );

            // form reduced a
            for( j = 0; j < i; j++ )
            {
                f     = a [i] [j];
                g     = e [j] - hh * f;
                e [j] = g;

                for( k = 0; k <= j; k++ )
				{
                    a [j] [k] -= f * e [k] + g * a [i] [k];
				}
            }

            for( k = i; k--; )
			{
                a [i] [k] *= scale;
			}
        }

        d [i] = h;
    }

    d [0] = e [0] = 0.0;

    // accumulation of transformation matrices
    for( i = 0; i < n; i++ )
    {
        if( d [i] )
        {
            for( j = 0; j < i; j++ )
            {
                g = 0.0;

                for( k = 0; k < i; k++ )
				{
                    g += a [i] [k] * a [k] [j];
				}

                for( k = 0; k < i; k++ )
				{
                    a [k] [j] -= g * a [k] [i];
				}
            }
        }

        d [i]     = a [i] [i];
        a [i] [i] = 1.0;

        for( j = 0; j < i; j++ )
		{
            a [i] [j] = a [j] [i] = 0.0;
		}
    }

	//
	// eigenvalues from tridiagonal form
	//
    if( n <= 1 )
	{
		return;
	}

    for( i = 1; i < n; i++ )
	{
        e [i-1] = e [i];
	}

    e [n-1] = 0.0;

    for( l = 0; l < n; l++ )
    {
        j = 0;

        do
        {
            // look for small sub-diagonal element
            for( m = l; m < n-1; m++ )
            {
                s = fabs( d [m] ) + fabs( d [m+1] );
                if( fabs( e [m] ) + s == s )
				{
					break;
				}
            }

            p = d [l];

            if( m != l )
            {
                if( j++ == maxIterC )
				{
					throw 1;
				}

                // form shift
                g = (d [l+1] - p ) / (2.0 * e [l] );
                r = sqrt( g * g + 1.0 );
                g = d [m] - p + e [l] / (g + ((g ) > 0 ? fabs(r ) : -fabs(r )));
                s = c = 1.0;
                p = 0.0;

                for( i = m; i-- > l; )
                {
                    f = s * e [i];
                    b = c * e [i];

                    if( fabs( f ) >= fabs( g ) )
                    {
                        c       = g / f;
                        r       = sqrt( c * c + 1.0 );
                        e [i+1] = f * r;
                        s       = 1.0 / r;
                        c      *= s;
                    }
                    else
                    {
                        s       = f / g;
                        r       = sqrt( s * s + 1.0 );
                        e [i+1] = g * r;
                        c       = 1.0 / r;
                        s      *= c;
                    }

                    g       = d [i+1] - p;
                    r       = (d [i] - g ) * s + 2.0 * c * b;
                    p       = s * r;
                    d [i+1] = g + p;
                    g       = c * r - b;

                    // form vector
                    for( k = 0; k < n; k++ )
                    {
                        f           = a [k] [i+1];
                        a [k] [i+1] = s * a [k] [i] + c * f;
                        a [k] [i  ] = c * a [k] [i] - s * f;
                    }
                }

                d [l] -= p;
                e [l]  = g;
                e [m]  = 0.0;
            }
        }
        while( m != l );
    }

	//
	// sorting eigenvalues
	//
    eigensort( vmatA, dvecA );

	//
	// normalizing eigenvectors
	//
    for( j = n; j--; )
    {
        s = 0.0;
        for( i = n; i--; )
		{
			s += a [i] [j] * a [i] [j];
		}
        s = sqrt( s );

        for( i = n; i--; )
		{
			a [i] [j] /= s;
		}
    }
}


//===========================================================================
/*!
 *  \brief Used as frontend for 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA) 
 *  for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the 
 *  Givens and Householder reduction without corrupting the input matrix "amatA" during application.
 *
 *  Frontend for function 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA), 
 *  when using type \em Array instead of type \em Array2D.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom triangular matrix must contain values.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column will contain one
 *                eigenvector.
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order.
 *	\param  od n-dimensional vector with the calculated 
 *              offdiagonal of the Householder transformation.
 *
 *      \return   none.
 *
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em A is not a square matrix
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes
 *         none
 *
 *  \par Status
 *      stable
 *
 */
void eigensymm
(
    const Array< double >& A,
          Array< double >& G,
          Array< double >& l,
          Array< double >& od
)
{
    SIZE_CHECK( A.ndim( ) == 2 && A.dim( 0 ) == A.dim( 1 ) )

	G.resize( A );
	l.resize( A.dim( 0 ) );
	od.resize( A.dim( 0 ) );

	Array2DReference< double > amatL( const_cast< Array< double >& >( A ) );
	Array2DReference< double > vmatL( G );

    eigensymm( amatL, vmatL, l, od);
}















