//===========================================================================
/*!
 *  \file eigensymmJacobi2.cpp
 *
 *  \brief Used to calculate the eigenvectors and eigenvalues of
 *         the symmetric matrix "amatA" using a modified Jacobi method.
 *
 *  Here the so-called Jacobi rotation is used to calculate
 *  the eigenvectors and values, but with a modification to
 *  avoid convergence problems.
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
 *      $RCSfile: eigensymmJacobi2.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: eigensymmJacobi2.cpp,v $
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
 *      Revision 1.3  2002/05/16 13:55:03  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.2  2001/11/30 13:25:45  rudi
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
 *         eigenvectors of the symmetric matrix "amatA" using a modified
 *         Jacobi method.
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
 *  This function uses the Jacobi method as in #eigensymmJacobi,
 *  but the method is modificated after J. von Neumann to avoid
 *  convergence problems when dealing with low machine precision.
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which must be symmetric, so 
 *                    only the bottom
 *                    triangular matrix must contain values.
 *                    Values below the diagonal will be destroyed.
 *      \param	vmatA \f$ n \times n \f$ matrix with the calculated 
 *                    normalized 
 *                    eigenvectors, each column will contain one
 *                    eigenvector.
 *	\param  dvecA n-dimensional vector with the calculated 
 *                    eigenvalues in descending order.
 *      \return       none.
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that
 *             	\em amatA is not a square matrix
 *
 *  \example eigensymmJacobi2_test.cpp
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
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
void eigensymmJacobi2
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK( amatA.ndim( ) == 2 && amatA.dim( 0 ) == amatA.dim( 1 ) )

	vmatA.resize( amatA );
	dvecA.resize( amatA.dim( 0 ) );

	double **a   = amatA.ptrArr( );
	double **vec = vmatA.ptrArr( );
	double  *val = dvecA.elemvec( );
	unsigned n   = dvecA.nelem( );

    unsigned ind, i, j, l, m;
    double anorm, thr, aml, all, amm, x, y;
    double *aim, *ail;
    double sinx, sinx2, cosx, cosx2, sincs;

    /*
     * Diagonalelemente von 'a' nach 'val' kopieren,
     * 'val' konvergiert gegen die Eigenwerte
     */
    for( j = 0; j < n; j++ )
	{
		val [j] = a [j] [j];
	}

    /*
     * Einheitsmatrix in 'vec' initialisieren
     * Quadratnorm der Nebendiagonalelemente von 'a' in 'anorm'
     */
    anorm = 0.0;
    for( i = 0; i < n; i++ )
	{
        for( j = 0; j < i; j++ )
		{
            anorm += a [i] [j] * a [i] [j];
            vec [i] [j] = vec [j] [i] = 0.0;
        }
        vec [i] [i] = 1.0;
    }

    if( anorm <= 0.0 )
	{
        goto done; /* Matrix hat keine Nebendiagonalelemente */
	}

    anorm  = sqrt( 4 * anorm );
    thr    = anorm;
    anorm /= n;

    while( anorm + thr > anorm )  /* relative Genauigkeit */
/*  while( 1 + thr > 1 ) */       /* absolute Genauigkeit */
    {
        thr /= n;

        do
		{
            ind = 0;

            for( l = 0; l < n-1; l++ )
			{
                for( m = l+1; m < n; m++ )
				{
					if( fabs( aml = a [m] [l] ) < thr )
					{
						continue;
					}

					ind   = 1;
					all   = val [l];
					amm   = val [m];
					x     = (all - amm ) / 2;
                    y     = -aml / hypot( aml, x );
					if( x < 0.0 )
					{
						y = -y;
					}
					sinx  = y / sqrt( 2 * (1 + sqrt( 1 - y * y )));
					sinx2 = sinx * sinx;
					cosx  = sqrt( 1 - sinx2 );
					cosx2 = cosx * cosx;
					sincs = sinx * cosx;

                    /* Spalten l und m rotieren */
                    for( i = 0; i < n; i++ )
					{
						if( (i != m ) &&( i != l ) )
						{
                            aim  = i > m ? &a [i] [m] : &a [m] [i];
                            ail  = i > l ? &a [i] [l] : &a [l] [i];
                            x    = *ail * cosx - *aim * sinx;
                            *aim = *ail * sinx + *aim * cosx;
                            *ail = x;
						}
						x = vec [i] [l];
						y = vec [i] [m];
						vec [i] [l] = x * cosx - y * sinx;
						vec [i] [m] = x * sinx + y * cosx;
					}

					x = 2.0 * aml * sincs;
					val   [l] = all * cosx2 + amm * sinx2 - x;
					val   [m] = all * sinx2 + amm * cosx2 + x;
					a [m] [l] = ( all - amm ) * sincs + aml * (cosx2 - sinx2 );
				}
			}
		}
		while( ind );
    }

done:	;

    /*
     * Eigenwerte sortieren
     */
    eigensort( vmatA, dvecA );
}

