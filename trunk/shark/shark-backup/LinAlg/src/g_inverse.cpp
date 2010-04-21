//===========================================================================
/*!
 *  \file g_inverse.cpp
 *
 *  \brief Determines the generalized inverse matrix of an input matrix
 *         by using singular value decomposition.
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
 *      $RCSfile: g_inverse.cpp,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: g_inverse.cpp,v $
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
 *      Revision 1.6  2002/08/01 10:14:39  arne
 *      g_inverse() now returns rank
 *
 *      Revision 1.5  2002/07/31 10:05:49  arne
 *      bugfix
 *
 *      Revision 1.4  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 13:25:45  rudi
 *      Doxygen comments added.
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
 *  \brief Calculates the generalized inverse matrix of input matrix "amatA".
 *
 *  Given an input matrix \f$ X \f$ this function uses singular value 
 *  decomposition to determine the generalized inverse matrix \f$ X' \f$,
 *  so that
 *
 *  \f$
 *      XX'X = X
 *  \f$
 *
 *  If \f$ X \f$ is singular, i.e. \f$ det(X) = 0 \f$ or \f$ X \f$ is 
 *  non-square then \f$ X' \f$ is not unique.
 *
 *      \param  amatA \f$ m \times n \f$ input matrix.
 *      \param	bmatA \f$ n \times m \f$ generalized inverse matrix.
 *      \return       none.
 *
 *  \example g_inverse_matrix.cpp
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
 *  \sa svd.cpp
 *
 */
unsigned g_inverse
(
	const Array2D< double >& amatA,
	Array2D< double >& bmatA
)
{
    unsigned i, j, k, m, n, r;

	m = amatA.rows( );
	n = amatA.cols( );

	Array2D< double > umatL( m, n );
	Array2D< double > vmatL( m, n );
	Array  < double > wvecL( n );

	bmatA.resize( n, m );

	if( m == 0 || n == 0 )
	{
		return 0;
	}

    svd    ( amatA, umatL, vmatL, wvecL );
    svdsort( umatL, vmatL, wvecL );

    r = svdrank( amatA, umatL, vmatL, wvecL );

    for( i = 0; i < r; i++ )
	{
		wvecL( i ) = 1. / wvecL( i );
	}
    for(      ; i < n; i++ )
	{
		wvecL( i ) = 0.;
	}

    for( i = 0; i < n; i++ )
	{
		for( j = 0; j < m; j++ )
		{
			double* vi = vmatL[ i ];
			double* uj = umatL[ j ];
			double  t  = 0.;

			for( k = 0; k < n; k++ )
			{
				t += vi[ k ] * wvecL( k ) * uj[ k ];
			}

			bmatA( i, j ) = t;
		}
    }

	return r;
}




