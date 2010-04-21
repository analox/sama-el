//===========================================================================
/*!
 *  \file detsymm.cpp
 *
 *  \brief Used to calculate the determinate, the eigenvalues and eigenvectors
 *         of the symmetric matrix "amatA".
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
 
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: detsymm.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: detsymm.cpp,v $
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
 *      Revision 1.3  2002/05/16 13:55:04  rudi
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
 *  \brief Calculates the determinate of the symmetric matrix "amatA".
 *
 *  Calculates the determinate of matrix \em amatA by using its
 *  \em n eigenvalues \f$ x_j \f$ that first will be calculated. 
 *  The determinate is then given as:
 *
 *  \f$
 *  \prod_{j=1}^n x_j
 *  \f$
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which is symmetric, so
 *                    only the bottom triangular matrix must contain
 *                    values. At the end of the function \em amatA
 *                    always contains the full matrix.
 *      \param	vmatA \f$ n \times n \f$ matrix, that will
 *                    contain the scaled eigenvectors at the
 *                    end of the function.
 *	\param  dvecA n-dimensional vector that will contain
 *                    the eigenvalues at the end of the
 *                    function.
 *      \return       The determinate of matrix \em amatA.
 *      \throw check_exception the type of the eception will be
 *             "size mismatch" and indicates that \em amatA is
 *             not a square matrix
 *
 *  \example detsymm_test.cpp
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
 *      2003/10/02 by S. Wiegand
 *      due to name change of 'eigensymm';
 *
 *  \par Status
 *      stable
 *
 */
double detsymm
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK( amatA.dim( 0 ) == amatA.dim( 1 ) )

	unsigned n = amatA.dim( 0 );
    unsigned i;
	unsigned j;
    double   det;

    // Fill upper triangular matrix:
    for( i = 0; i + 1 < n; i++ )
	{
        for( j = i + 1; j < n; j++ )
		{
            amatA( i, j ) = amatA( j, i );
		}
	}

    // Calculate eigenvalues:
    //eigensymm( amatA, vmatA, dvecA );
    eigensymm_obsolete( amatA, vmatA, dvecA );

    for( i = 0; i + 1 < n; i++ )
	{
        for( j = i + 1; j < n; j++ )
		{
            amatA( j, i ) = amatA( i, j );
		}
	}

    // Calculate determinate as product of eigenvalues:
    for( i = 0, det = 1; i < n; det *= dvecA( i++ ) );

    return det;
}


