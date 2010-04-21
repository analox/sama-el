//===========================================================================
/*!
 *  \file rankDecomp.cpp
 *
 *  \brief Used to calculate the rank, the eigenvalues and eigenvectors
 *         of the symmetric matrix "amatA".
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Copyright (c) 1998-2001:
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
 *      $RCSfile: rankDecomp.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: rankDecomp.cpp,v $
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
#include <cmath>
#include "LinAlg/linalg.h"


//===========================================================================
/*!
 *  \brief Calculates the rank of the symmetric matrix "amatA",
 *         its eigenvalues and eigenvectors.
 *
 *  Determines the rank of matrix \em amatA and additionally
 *  calculates the eigenvalues and eigenvectors of the matrix.
 *  Empty eigenvalues (i.e. eigenvalues equal to zero) are
 *  set to the greatest calculated eigenvalue and each
 *  eigenvector \f$ x_j \f$ is scaled by multiplying
 *  it with the scalar value \f$ \frac{1}{\sqrt{\lambda_j}} \f$.
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which is symmetric, so
 *                    only the bottom triangular matrix must contain
 *                    values. At the end of the function \em amatA
 *                    always contains the full matrix.
 *      \param	vmatA \f$ n \times n \f$ matrix, which will
 *                    contain the scaled eigenvectors.
 *	\param  dvecA n-dimensional vector with the calculated 
 *                    eigenvalues in descending order.
 *      \return       The rank of matrix \em amatA.
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em amatA
 *             is not a square matrix
 *
 *  \example  rank_decomp_test.cpp
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
 * 
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *      2003/10/02 by S. Wiegand
 *      due to name change of 'eigensymm';
 *
 *  \par Status
 *      stable
 *
 */
unsigned rankDecomp
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK( amatA.ndim( ) == 2 && amatA.dim( 0 ) == amatA.dim( 1 ) )

	vmatA.resize( amatA );
	dvecA.resize( amatA.dim( 0 ) );

	double **a = amatA.ptrArr( );
	double **t = vmatA.ptrArr( );
	double *d = dvecA.elemvec( );
	unsigned n = dvecA.nelem( );

    unsigned i, j, r;
    double   dt;

    /*
     * Eigenwerte berechnen (Nebendiagonalen sichern)
     */
    for( i = 0; i + 1 < n; i++ )
        for( j = i + 1; j < n; j++ )
            a [i] [j] = a [j] [i];

    //eigensymm( amatA, vmatA, dvecA );
    eigensymm_obsolete( amatA, vmatA, dvecA );

    for( i = 0; i + 1 < n; i++ )
        for( j = i + 1; j < n; j++ )
            a [j] [i] = a [i] [j];

    /*
     * numerischen Rang ermitteln
     */
    r = rank( amatA, vmatA, dvecA );

    /*
     * Eigenwerte 0 auf "beliebige" Konstante setzen (= groesster Eigenwert)
     */
    for( i = r; i < n; i++ ) d [i] = d [0];

    /*
     * Eigenvektoren skalieren
     */
    if( r )
	for( j = n; j--; ) {
	    dt = 1. / sqrt( d[ j ] );

	    for( i = n; i--; ) t[ i ][ j ] *= dt;
	}

    return r;
}













