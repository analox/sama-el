//===========================================================================
/*!
 *  \file linearRegress.cpp
 *
 *  \brief Given the correlations of two data vectors "x" and "y" and
 *         their mean values, the function in this file 
 *         summarizes the data by finding a linear mapping
 *         that will approximate the data. 
 *
 *  The function in this file is used as subroutine for the class
 *  LinearRegression, where you can directly work with the
 *  data vectors \f$x\f$ and \f$y\f$ itself and the correlation
 *  and mean values are calculated for you. <br>
 *  Please refer to the reference of this class for more
 *  detailed information about linear regression.
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
 *      $RCSfile: linearRegress.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:10 $
 *
 *  \par Changes:
 *      $Log: linearRegress.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
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
 *  \brief Given the correlations of the n-dimensional data vector "x" 
 *         and the m-dimensional data vector "y" and also given
 *         their mean values, this function summarizes the data by 
 *         finding a linear mapping that will approximate the data. 
 *
 *  This function is used as subroutine for the class
 *  LinearRegression, where you can directly work with the
 *  data vectors \f$x\f$ and \f$y\f$ itself and the correlation
 *  and mean values are calculated for you. <br>
 *  Please refer to the reference of this class for more
 *  detailed information about linear regression.
 *  
 *      \param  cxxMatA \f$n \times n\f$ matrix, that contains
 *                      the covariances between the single dimensions
 *                      of the vector \f$x\f$.
 *                      Only the lower triangle matrix must contain values.
 *      \param  cxyMatA \f$n \times m\f$ matrix, that contains
 *                      the covariances between the single values of
 *                      \f$x\f$ and \f$y\f$.
 *      \param	mxVecA  \f$n\f$-dimensional vector that contains the mean
 *                      values of vector \f$x\f$.
 *      \param	myVecA  \f$m\f$-dimensional vector that contains the mean
 *                      values of vector \f$y\f$. 
 *      \param  amatA   \f$m \times n\f$ matrix, that will contain the
 *                      transposed transformation matrix \f$A^T\f$,
 *                      that is used for the linear mapping
 *                      \f$y = A \cdot x + b\f$ of the data.
 *      \param	bvecA   The \f$m\f$-dimensional vector \f$b\f$, 
 *                      that is used for the linear mapping
 *                      \f$y = A \cdot x + b\f$ of the data.
 *      \param	dvecA   \f$n\f$-dimensional temporary vector, that
 *                      will contain the eigenvalues of matrix 
 *                      \em cxxMatA in descending order.
 *      \return none
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em mxVecA or
 *             \em myVecA are not one-dimensional or that the dimensions
 *             of the matrices \em cxxMatA or \em cxyMatA don't
 *             correspond to the sizes of \em mxVecA and \em myVecA
 *             (the size of the first dimension of \em cxxMatA and
 *              \em cxyMatA and the size of the second dimension
 *              of \em cxxMatA must be the same than the size of
 *              \em mxVecA and the size of the second dimension
 *              of \em cxyMatA must be the same than the size
 *              of \em myVecA)
 *             
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
void linearRegress
(
	Array2D< double >& cxxMatA,
	Array2D< double >& cxyMatA,
	Array  < double >& mxVecA,
	Array  < double >& myVecA,
	Array2D< double >& amatA,
	Array  < double >& bvecA,
	Array  < double >& dvecA
)
{
	SIZE_CHECK
	(
		mxVecA.ndim( ) == 1 &&
		myVecA.ndim( ) == 1 &&
		mxVecA.dim( 0 ) == cxxMatA.dim( 0 ) &&
		mxVecA.dim( 0 ) == cxxMatA.dim( 1 ) &&
		mxVecA.dim( 0 ) == cxyMatA.dim( 0 ) &&
		myVecA.dim( 0 ) == cxyMatA.dim( 1 )
	)

    unsigned i, j, k;

	Array2D< double > vmatL( cxxMatA.dim( 0 ), cxxMatA.dim( 0 ) );
    rankDecomp( cxxMatA, vmatL, dvecA );

	amatA.resize( myVecA.dim( 0 ), mxVecA.dim( 0 ) );
	bvecA.resize( myVecA );

    for( i = 0; i < cxxMatA.dim( 0 ); ++i )
	{
        for( j = 0; j < cxxMatA.dim( 1 ); ++j )
        {
            double t = 0.;

            for( k = 0; k < vmatL.dim( 1 ); ++k )
			{
                t += vmatL( i, k ) * vmatL( j, k );
			}

            cxxMatA( i, j ) = t;
        }
	}

    for( i = 0; i < amatA.dim( 0 ); ++i )
	{
        for( j = 0; j < amatA.dim( 1 ); ++j )
        {
            double t = 0.;

            for( k = 0; k < cxxMatA.dim( 0 ); ++k )
			{
                t += cxyMatA( k, i ) * cxxMatA( k, j );
			}

            amatA( i, j ) = t;
        }
	}

    for( i = 0; i < bvecA.dim( 0 ); ++i )
    {
        double t = 0.;

        for( j = 0; j < mxVecA.dim( 0 ); ++j )
		{
            t += amatA( i, j ) * mxVecA( j );
		}

        bvecA( i ) = myVecA( i ) - t;
    }
}
