//===========================================================================
/*!
 *  \file svd.cpp
 *
 *  \brief Used for singular value decomposition of rectangular and
 *         square matrices.
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
 *      $RCSfile: svd.cpp,v $<BR>
 *      $Revision: 2.3 $<BR>
 *      $Date: 2005/10/06 13:23:24 $
 *
 *  \par Changes:
 *      $Log: svd.cpp,v $
 *      Revision 2.3  2005/10/06 13:23:24  christian_igel
 *      bug found by Thomas removed
 *
 *      Revision 2.2  2004/06/16 15:22:20  saviapbe
 *
 *      Some bugs in the doxygen documentation were removed.
 *
 *      Revision 2.1  2004/06/01 15:41:44  saviapbe
 *
 *      The bugs in the doxygen's documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.8  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.7  2001/11/30 13:26:06  rudi
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
#include <cstdlib>
#include "LinAlg/linalg.h"

#ifndef min
#define min( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )
#endif

#ifndef max
#define max( a, b ) ( ( a ) > ( b ) ? ( a ) : ( b ) )
#endif

static double SIGN (double a, double b)
{
    return b > 0 ? fabs (a) : -fabs (a);
}

//===========================================================================
/*!
 *  \brief Determines the singular value decomposition of a rectangular
 *  matrix "amatA".
 *
 *  Given a \f$ m \times n \f$ matrix \em amatA, this routine computes its 
 *  singular value decomposition, defined as
 *
 *  \f$
 *      A = UWV^T
 *  \f$
 *
 *  where W is an \f$ n \times n \f$ diagonal matrix with positive
 *  or zero elements, the so-called \em singular \em values.
 *  The matrices \em U and \em V are each orthogonal in the sense
 *  that their columns are orthonormal, i.e.
 *
 *  \f$
 *      UU^T = VV^T = V^TV = 1
 *  \f$   
 *
 *      \param  amatA The input matrix \em A, with size \f$ m \times n \f$ and
 *                    \f$ m \geq n \f$.
 *      \param  umatA The \f$ m \times n \f$ column-orthogonal matrix \em U 
 *                    determined by the function.
 *      \param  vmatA The \f$ n \times n \f$ orthogonal matrix \em V
 *                    determined by the function.
 *      \param  w     n-dimensional vector with the calculated singular values.
 *      \return       none
 *
 *  \example svd_test.cpp
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
 */
void svd
(
	const Array2D< double >& amatA,
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& w // singular values
)
{
	unsigned m = amatA.rows( ); /* rows */
	unsigned n = amatA.cols( ); /* cols */
	double** const a = amatA.ptrArr( ); /* input matrix */
	double **u = umatA.ptrArr( ); /* left vectors */
	double **v = vmatA.ptrArr( ); /* right vectors */

    int flag;
    unsigned i, its, j, jj, k, l, nm(0);
    double anorm, c, f, g, h, s, scale, x, y, z;

	Array< double > rv1( n );

    /* copy A to U */
    for (i = 0; i < m; i++)
	{
        for (j = 0; j < n; j++)
		{
            u [i] [j] = a [i] [j];
		}
	}

    /* householder reduction to bidiagonal form */
    g = scale = anorm = 0.0;

    for (i = 0; i < n; i++)
	{
        l = i + 1;
        rv1 (i) = scale * g;
        g = s = scale = 0.0;

		if (i < m)
		{
			for (k = i; k < m; k++)
			{
				scale += fabs (u [k] [i]);
			}

			if (scale != 0.0)
			{
				for (k = i; k < m; k++)
				{
					u [k] [i] /= scale;
					s += u [k] [i] * u [k] [i];
				}

				f = u [i] [i];
				g = -SIGN (sqrt (s), f);
				h = f * g - s;
				u [i] [i] = f - g;

				for (j = l; j < n; j++)
				{
					s = 0.0;
					for (k = i; k < m; k++)
					{
						s += u [k] [i] * u [k] [j];
					}

					f = s / h;
					for (k = i; k < m; k++)
					{
						u [k] [j] += f * u [k] [i];
					}
				}

				for (k = i; k < m; k++)
				{
					u [k] [i] *= scale;
				}
			}
		}

        w( i ) = scale * g;
        g = s = scale = 0.0;

        if (i < m && i != n-1)
		{
            for (k = l; k < n; k++)
			{
                scale += fabs (u [i] [k]);
			}

			if (scale != 0.0)
			{
				for (k = l; k < n; k++)
				{
					u [i] [k] /= scale;
					s += u [i] [k] * u [i] [k];
				}

				f = u [i] [l];
				g = -SIGN (sqrt(s), f);
				h = f * g - s;
				u [i] [l] = f - g;

				for (k = l; k < n; k++)
				{
                    rv1( k ) = u [i] [k] / h;
				}

				for (j = l; j < m; j++)
				{
					s = 0.0;
					for (k = l; k < n; k++)
					{
						s += u [j] [k] * u [i] [k];
					}

					for (k = l; k < n; k++)
					{
						u [j] [k] += s * rv1( k );
					}
				}

				for (k = l; k < n; k++)
				{
					u [i] [k] *= scale;
				}
			}
		}

        anorm = max (anorm, fabs (w( i )) + fabs (rv1( i )));
    }

    /* accumulation of right-hand transformations */
    for (l = i = n; i--; l--)
	{
		if (l < n)
		{
			if (g != 0.0)
			{
				for (j = l; j < n; j++)
				{
					/* double division avoids possible underflow */
					v [j] [i] = (u [i] [j] / u [i] [l]) / g;
				}

				for (j = l; j < n; j++)
				{
					s = 0.0;
					for (k = l; k < n; k++)
					{
						s += u [i] [k] * v [k] [j];
					}

					for (k = l; k < n; k++)
					{
						v [k] [j] += s * v [k] [i];
					}
				}
			}

			for (j = l; j < n; j++)
			{
				v [i] [j] = v [j] [i] = 0.0;
			}
		}

		v [i] [i] = 1.0;
		g = rv1( i );
    }

    /* accumulation of left-hand transformations */
    for (l = i = min (m, n); i--; l--)
	{
		g = w( i );

		for (j = l; j < n; j++)
		{
			u [i] [j] = 0.0;
		}

		if (g != 0.0)
		{
			g = 1.0 / g;

			for (j = l; j < n; j++)
			{
				s = 0.0;
				for (k = l; k < m; k++)
				{
					s += u [k] [i] * u [k] [j];
				}

				/* double division avoids possible underflow */
				f = (s / u [i] [i]) * g;

				for (k = i; k < m; k++)
				{
					u [k] [j] += f * u [k] [i];
				}
			}

			for (j = i; j < m; j++)
			{
				u [j] [i] *= g;
			}
		}
		else
		{
			for (j = i; j < m; j++)
			{
				u [j] [i] = 0.0;
			}
		}

		u [i] [i]++;
    }

    /* diagonalization of the bidiagonal form */
    for (k = n; k--; )
	{
		for (its = 1; its <= 30; its++)
		{
			flag = 1;

			/* test for splitting */
//			for (l = k + 1; l--; )
/* Thomas Buecher: 
        Änderung, die fehlerhaft sein kann, aber Absturz verhindert:
	l kann 0 werden --> nm (zwei Zeilen tiefer) wird riesig, da 0-1 auf 
	unsigned typ berechnet wird
*/
			for (l = k; l > 0; l--)
			{
				/* rv1 [0] is always zero, so there is no exit */
				nm = l - 1;

				if (fabs (rv1( l )) + anorm == anorm)
				{
					flag = 0;
					break;
				}

				if (fabs (w( nm )) + anorm == anorm)
				{
					break;
				}
			}

			if (flag)
			{
				/* cancellation of rv1 [l] if l greater than 0 */
				c = 0.0;
				s = 1.0;

				for (i = l; i <= k; i++)
				{
					f = s * rv1( i );
					rv1( i ) *= c;

					if (fabs (f) + anorm == anorm)
					{
						break;
					}

					g = w( i );
					h = hypot(f, g);
					w( i ) = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;

					for (j = 0; j < m; j++)
					{
						y = u [j] [nm];
						z = u [j] [i];
						u [j] [nm] = y * c + z * s;
						u [j] [i ] = z * c - y * s;
					}
				}
			}

			/* test for convergence */
			z = w( k );

			if (l == k)
			{
				if (z < 0.0)
				{
					w( k ) = -z;
					for (j = 0; j < n; j++)
					{
						v [j] [k] = -v [j] [k];
					}
				}
				break;
			}

			if (its == 30)
			{
				throw k;
			}

			/* shift from bottom 2 by 2 minor */
			x = w( l );
			nm = k - 1;
			y = w( nm );
			g = rv1( nm );
			h = rv1( k );
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = hypot(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN (g, f))) - h)) / x;

			/* next qr transformation */
			c = s = 1.0;

			for (j = l; j < k; j++)
			{
				i = j + 1;
				g = rv1( i );
				y = w( i );
				h = s * g;
				g *= c;
				z = hypot(f, h);
				rv1( j ) = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;

				for (jj = 0; jj < n; jj++)
				{
					x = v [jj] [j];
					z = v [jj] [i];
					v [jj] [j] = x * c + z * s;
					v [jj] [i] = z * c - x * s;
				}

				z = hypot(f, h);
				w( j ) = z;

				/* rotation can be arbitrary if z is zero */
				if (z != 0.0)
				{
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}

				f = c * g + s * y;
				x = c * y - s * g;

				for (jj = 0; jj < m; jj++)
				{
					y = u [jj] [j];
					z = u [jj] [i];
					u [jj] [j] = y * c + z * s;
					u [jj] [i] = z * c - y * s;
				}
			}

			rv1( l ) = 0.0;
			rv1( k ) = f;
			w( k ) = x;
		}
    }
}


//===========================================================================
/*!
 *  \brief Used as frontend for function #svd(const Array2D< double >& amatA, Array2D< double >& umatA, Array2D< double >& vmatA, Array  < double >& w),
 *  calculates singular value decomposition for a square matrix "A".
 *
 *  Frontend for function #svd(const Array2D< double >& amatA, Array2D< double >& umatA, Array2D< double >& vmatA, Array  < double >& w), when using
 *  type \em Array instead of type \em Array2D. Additionally, the
 *  function can only be used for square matrices \em A.
 *
 *      \param  A The input matrix \em A, with size \f$ n \times n \f$.
 *      \param  U The \f$ n \times n \f$ orthogonal matrix determined 
 *                by the function.
 *      \param  V The \f$ n \times n \f$ orthogonal matrix determined by the 
 *                function.
 *      \param  W n-dimensional vector with the calculated singular values,
 *                sorted by ascending order.
 *      \return   none.
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em A
 *             is not a square matrix
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
 */
void svd
(
    const Array< double >& A,
	  Array< double >& U,
	  Array< double >& V,
	  Array< double >& W
)
{
    SIZE_CHECK( A.ndim( ) == 2 && A.dim( 0 ) == A.dim( 1 ) )

	Array2DReference< double > amatL( const_cast< Array< double >& >( A ) );
	Array2DReference< double > umatL( U );
	Array2DReference< double > vmatL( V );

    U.resize( A.dim( 0 ), A.dim( 1 ) );
    V.resize( A.dim( 1 ), A.dim( 1 ) );
    W.resize( A.dim( 1 ) );

    svd( amatL, umatL, vmatL, W );
    svdsort( umatL, vmatL, W );
}










