//===========================================================================
/*!
 *  \file eigensort.cpp
 *
 *  \brief Used to sort eigenvalues and their corresponding eigenvectors.
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
 *      $RCSfile: eigensort.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: eigensort.cpp,v $
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
 *  \brief Sorts the eigenvalues in vector "dvecA" and the corresponding
 *         eigenvectors in matrix "vmatA".
 *
 *  Given the matrix \em vmatA of eigenvectors and the vector
 *  \em dvecA of corresponding eigenvalues, the values in \em dvecA will
 *  be sorted by descending order and the eigenvectors in
 *  \em vmatA will change their places in a way, that at the
 *  end of the function an eigenvalue at position \em j of
 *  vector \em dvecA will belong to the eigenvector
 *  at column \em j of matrix \em vmatA.
 *  If we've got for example the following result after calling the function:
 *
 *  \f$
 *      \begin{array}{*{3}{r}}
 *          v_{11} & v_{21} & v_{31}\\
 *          v_{12} & v_{22} & v_{32}\\
 *          v_{13} & v_{23} & v_{33}\\
 *          & & \\
 *          v_1 & v_2 & v_3\\
 *      \end{array}
 *  \f$
 *
 *  then eigenvalue \f$ v_1 \f$ has the corresponding eigenvector
 *  \f$ ( v_{11}\ v_{12}\ v_{13} ) \f$ and \f$ v_1 > v_2 > v_3 \f$.
 *
 *
 *      \param	vmatA \f$ n \times n \f$ matrix with eigenvectors (each column
 *                    contains an eigenvector, corresponding to
 *                    one eigenvalue).
 *	\param  dvecA n-dimensional vector with eigenvalues, will
 *                    contain the eigenvalues in descending order
 *                    when returning from the function.
 *      \return       none.
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em dvecA
 *             is not one-dimensional or that the number of
 *             rows or the number of columns in \em vmatA
 *             is different from the number of values
 *             in \em dvecA 
 *
 *  \example eigensort_test.cpp
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
void eigensort
(
	Array2D< double >& vmatA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK
	(
		dvecA.ndim( ) == 1 &&
		dvecA.dim( 0 ) == vmatA.dim( 0 ) &&
		dvecA.dim( 0 ) == vmatA.dim( 1 )
	)

	double** v = vmatA.ptrArr( );
	double*  d = dvecA.elemvec( );
	unsigned n = dvecA.nelem( );
    unsigned i, j, l;
    double t;

    //
    // sort eigen values
    //
    for( i = 0; i < n - 1; i++ )
	{
		t = d[ l = i ];

		for( j = i + 1; j < n; j++ )
		{
			if( d[ j ] >= t )
			{
				t = d[ l = j ];
			}
		}

		if( l != i )
		{
			d[ l ] = d[ i ];
			d[ i ] = t;

			for( j = 0; j < n; j++ )
			{
				t           = v[ j ][ i ];
				v[ j ][ i ] = v[ j ][ l ];
				v[ j ][ l ] = t;
			}
		}
    }
}
