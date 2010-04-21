//===========================================================================
/*!
 *  \file rank.cpp
 *
 *  \brief Used to determine the rank of the symmetric matrix \em amatA.
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
  *
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: rank.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: rank.cpp,v $
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
 *  \brief Determines the numerical rank of the symmetric matrix "amatA".
 *
 *  Given the \f$ n \times n \f$ matrix \em amatA, this function uses the
 *  eigenvectors \em vmatA and the eigenvalues \em dvecA to
 *  calculate the rank of \em amatA.
 *  For the calculation of the eigenvectors and eigenvalues
 *  see also #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA).
 *
 *      \param  amatA \f$ n \times n \f$ matrix, which is symmetric
 *                    and is not changed by the function.
 *      \param  vmatA \f$ n \times n \f$ matrix with the normalized 
 *                    eigenvectors, each column contains one eigenvector.
 *                    The matrix is not changed by the function.
 *      \param  dvecA n-dimensional vector of the eigenvalues,
 *                    given in descending order.
 *                    The vector is not changed by the function.
 *      \return       The rank of matrix \em amatA.
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em dvecA
 *             is not one-dimensional or that \em amatA or
 *             \em vmatA has not the same number of rows or columns
 *             than \em dvecA contains values
 *
 *  \example  rank_test.cpp
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
unsigned rank
(
	const Array2D< double >& amatA,
	const Array2D< double >& vmatA,
	const Array  < double >& dvecA
)
{
	SIZE_CHECK
	(
		dvecA.ndim( ) == 1 &&
		dvecA.dim( 0 ) == amatA.dim( 0 ) &&
		dvecA.dim( 0 ) == amatA.dim( 1 ) &&
		dvecA.dim( 0 ) == vmatA.dim( 0 ) &&
		dvecA.dim( 0 ) == vmatA.dim( 1 )
	)

	unsigned n = dvecA.nelem( );
    unsigned r;
    double   s;
	double   u;

    /*
     * numerischen Rang ermitteln:
     *   wenn die letzten Eigenwerte < 0 sind, ist der Rechenfehler in den
     *   Werten > 0 mindestens so gross wie der Betrag des letzten Wertes,
     *   d.h. Werte unter dieser Schwelle sind zu verwerfen.
     *   Zusaetzliche Schwellen bilden die relative Maschinengenauigkeit
     *   und der relative Fehler in den Eigenwerten
     */
    for( r = 0; r < n && dvecA( r ) > 0; r++ );

    u = r < n ? fabs( dvecA( n-1 ) ) : 0.0;
    r = 0;
    s = 0.0;
    while( r < n     &&
		   dvecA( r ) > u &&
		   dvecA( r ) + s > s &&
           dvecA( r ) > eigenerr( amatA, vmatA, dvecA, r ) )
	{
        s += dvecA( r++ );
	}

    return r;
}
