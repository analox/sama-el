//===========================================================================
/*!
 *  \file discrimAnalysis.cpp
 *
 *  \brief Given "m" classes of data, the covariances between
 *         the values within one class and the covariances
 *         between all the classes, the method in this file calculates
 *         the transformation matrix that will project the
 *         data in a way, that maximum separation of the 
 *         different classes is given.
 *
 *  The method in this file is used as subroutine by class LinearClassifier,
 *  so please refer to the reference of this class for more information.
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Copyright (c) 1998:
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
 *      $RCSfile: discrimAnalysis.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:10 $
 *
 *  \par Changes:
 *      $Log: discrimAnalysis.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:55:04  rudi
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
 *  \brief Given "m" classes of data, the covariances between
 *         the values within one class and the covariances
 *         between all the classes, this method calculates
 *         the transformation matrix that will project the
 *         data in a way, that maximum separation of the 
 *         different classes is given.
 *
 *  The within class and the between class scatter are used to
 *  formulate criteria for class separability. While the within-
 *  covariance matrix contains the expected covariance of each of the
 *  classes, the between-covariance matrix can be seen as covariance
 *  of the data set, whose members are the mean vectors of each
 *  class. <br>
 *  For \em m classes \em m optimizing criteria are used for transforming
 *  the data sets independently. <br>
 *  Based on this the eigenvectors are calculated and those with
 *  eigenvalues equal to zero are neglected, so the resulting
 *  eigenvectors as one-dimensional subspaces of the transformation
 *  are invariant. <br>
 *  The data points are then projected onto the maximally discriminating 
 *  axes, represented by the eigenvectors. <br>
 *  Having this transformation it can be used to classify new data
 *  points by projecting them 
 *  into the discriminant space and measure the euclidian distance of the
 *  data point to each of the \em m classes. The data point
 *  is then associated to the class with the smallest euclidian
 *  distance. <br>
 *  This method is used as subroutine for class LinearClassifier,
 *  so please refer to the reference of this class for more
 *  information about the linear discriminant analysis.
 *
 *  \param betweenCovarA The symmetric \f$m \times m\f$ between-covariance
 *                       matrix. Only the bottom
 *                       triangle matrix must contain values.
 *                       The content of the matrix will be destroyed.
 *  \param withinCovarA  The symmetric \f$m \times m\f$ within-covariance
 *                       matrix. Only the bottom
 *                       triangle matrix must contain values.
 *                       The content of the matrix will be destroyed.
 *  \param transMatA     The calculated \f$m \times m\f$ transformation
 *                       matrix.
 *  \param dvecA         Temporary \f$m\f$-dimensional vector that is used
 *                       to store the calculated eigenvalues.
 *  \param m             Dimension of the discriminant space
 *                       \f$ = \#(classes) - 1\f$.
 *  \return none
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
void discrimAnalysis
(
	Array2D< double >& betweenCovarA,
	Array2D< double >& withinCovarA,
	Array2D< double >& transMatA,
	Array  < double >& dvecA,
	unsigned& m
)
{
    unsigned i, j, k;

    // Calculate transformation matrix for standard normal distribution:
    rankDecomp( withinCovarA, transMatA, dvecA );

    // Covariance matrix is symmetric, so only the lower triangle
    // matrix must contain values. Here the whole matrix is filled
    // with values by mirroring the given values:
    for( i = 0; i + 1 < betweenCovarA.dim( 0 ); ++i )
	{
        for( j = i + 1; j < betweenCovarA.dim( 0 ); ++j )
		{
            betweenCovarA( i, j ) = betweenCovarA( j, i );
		}
	}

    for( i = 0; i < betweenCovarA.dim( 0 ); ++i )
	{
        for( j = 0; j < betweenCovarA.dim( 0 ); ++j )
        {
            double t = 0.;

            for( k = 0; k < betweenCovarA.dim( 0 ); ++k )
			{
                t += transMatA( k, i ) * betweenCovarA( k, j );
			}

            withinCovarA( i, j ) = t;
        }
	}

    for( i = 0; i < betweenCovarA.dim( 0 ); ++i )
	{
        for( j = 0; j <= i; ++j )
        {
            double t = 0.;

            for( k = 0; k < betweenCovarA.dim( 0 ); ++k )
			{
                t += withinCovarA( i, k ) * transMatA( k, j );
			}

            betweenCovarA( i, j ) = t;
        }
	}

    // Dimension reduction by KLT on the mean values of the classes:
    //eigensymm( betweenCovarA, withinCovarA, dvecA );
    eigensymm_obsolete( betweenCovarA, withinCovarA, dvecA );

     // Calculate numerical rank, "m" is eventually reduced.
     //
     // this is not used any more, because the dimension reduction
     // can lead to worse results, if the rank is previously known
     // exactly (and this applies here).
     /*
         k = rank( b, w, d, n );
         if( k < *m ) *m = k;
     */

    for( i = 0; i < betweenCovarA.dim( 0 ); ++i )
	{
        for( j = 0; j < m; ++j )
        {
            double t = 0.;

            for( k = 0; k < betweenCovarA.dim( 0 ); k++ )
			{
                t += transMatA( i, k ) * withinCovarA( k, j );
			}

            betweenCovarA( j, i ) = t;
        }
	}

    for( i = 0; i < m; i++ )
	{
        for( j = 0; j < betweenCovarA.dim( 0 ); ++j )
		{
            transMatA( i, j ) = betweenCovarA( i, j );
		}
	}
}






