//===========================================================================
/*!
 *  \file svdrank.cpp
 *
 *  \brief Determines the numerical rank of a rectangular matrix,
 *         when a singular value decomposition for this matrix has taken 
 *         place before.
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
 *
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: svdrank.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: svdrank.cpp,v $
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
 *      Revision 1.5  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2001/11/30 13:26:06  rudi
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

#ifdef _WINDOWS
// disable warning C2055: unreferenced formal parameter
#pragma warning(disable: 4100)
#endif

//===========================================================================
/*!
 *  \brief Determines the numerical rank of a rectangular matrix "amatA",
 *         when a singular value decomposition for "amatA" has taken place 
 *         before.
 *
 *  For a singular value decomposition defined as
 *
 *  \f$
 *      A = UWV^T
 *  \f$
 *
 *  the resulting orthogonal matrices \em U and \em V and the singular
 *  values in \em W sorted by descending order are used to determine
 *  the rank of input matrix \em A. 
 * 
 *      \param  amatA The \f$ m \times n \f$ input matrix \em A, with
 *                    \f$ m \geq n \f$.
 *      \param  umatA The \f$ m \times n \f$ column-orthogonal matrix \em U.
 *      \param  vmatA The \f$ n \times n \f$ orthogonal matrix \em V.
 *      \param  wvecA n-dimensional vector containing the singular
 *                    values in descending order.
 *      \return       none.
 *
 *  \example svdrank_test.cpp
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
unsigned svdrank
(
	const Array2D< double >& amatA,
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& wvecA
)
{
        // unsigned m = amatA.rows( );
	unsigned n = amatA.cols( );
	// double** const a = amatA.ptrArr( );
	// double **u = umatA.ptrArr( );
	// double **v = vmatA.ptrArr( );
	double  *w = wvecA.elemvec( );

    unsigned r;
    double   s, t;

    /*
     * numerischen Rang ermitteln:
     *   wenn die letzten Singulaerwerte < 0 sind, ist der Rechenfehler in den
     *   Werten > 0 mindestens so gross wie der Betrag des letzten Wertes,
     *   d.h. Werte unter dieser Schwelle sind zu verwerfen.
     *   Zusaetzliche Schwellen bilden die relative Maschinengenauigkeit
     *   und der relative Fehler in den Singulaerwerten
     */
    for( r = 0; r < n && w[ r ] > 0.; r++ );

    t = r < n ? fabs( w[ n-1 ] ) : 0.0;
    r = 0;
    s = 0.0;
    while( r < n     &&
	   w[ r ] > t &&
	   w[ r ] + s > s /* &&
           w[ r ] > svderr( a, t, d, n, r ) */ )
        s += w[ r++ ];

    return r;
}





