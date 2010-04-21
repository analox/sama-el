//===========================================================================
/*!
 *  \file svdsort.cpp
 *
 *  \brief Used for sorting singular values and the orthogonal
 *         matrices \em U and \em V after a singular value
 *         decomposition of an input matrix \em A.
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
 *      $RCSfile: svdsort.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: svdsort.cpp,v $
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
 *  \brief Sorts the singular values in vector "wvecA" by descending order.
 *
 *  For a singular value decomposition defined as
 *
 *  \f$
 *      A = UWV^T
 *  \f$
 *
 *  the resulting orthogonal matrices \em U and \em V and the singular
 *  values in \em W can be sorted in a way, that the singular values
 *  are given in descending order, when leaving the function. 
 * 
 *      \param  umatA The \f$ m \times n \f$ matrix \em U.
 *      \param  vmatA The \f$ n \times n \f$ matrix \em V.
 *      \param  wvecA n-dimensional vector containing the singular
 *                    values.
 *      \return       none.
 *      \throw check_exception the type of the exception will
 *             be "size mismatch" and indicates that \em wvecA is not
 *             one-dimensional
 *
 *  \example svdsort_test.cpp
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
 *  \sa svd.cpp
 *
 */
void svdsort
(
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& wvecA
)
{
	SIZE_CHECK
	(
		wvecA.ndim( ) == 1
	)

	unsigned m = umatA.rows( ); /* rows */
	unsigned n = umatA.cols( ); /* cols */
	double **u = umatA.ptrArr( ); /* m x n matrix */
	double **v = vmatA.ptrArr( ); /* n x n matrix */
	double  *w = wvecA.elemvec( ); /* n x 1 vector */

    unsigned i, j, k;
    double   p;

    for (i = 0; i < n - 1; i++)
	{
		p = w [k = i];

		for (j = i + 1; j < n; j++)
		{
			if (w [j] >= p)
			{
				p = w [k = j];
			}
		}

		if (k != i)
		{
			w [k] = w [i];
			w [i] = p;

			for (j = 0; j < n; j++)
			{
				p         = v [j] [i];
				v [j] [i] = v [j] [k];
				v [j] [k] = p;
			}

			for (j = 0; j < m; j++)
			{
				p         = u [j] [i];
				u [j] [i] = u [j] [k];
				u [j] [k] = p;
			}
		}
    }
}




