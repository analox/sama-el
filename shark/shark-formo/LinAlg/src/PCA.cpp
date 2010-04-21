//===========================================================================
/*!
 *  \file PCA.cpp
 *
 *  \brief Methods for the "Principal Component Analysis" class.
 *
 *  This file offers methods for the \em Principal \em Component
 *  \em Analysis class, that is used to compress data for a
 *  better visualization and analysis.
 *
 *  \author  M. Kreutz
 *  \date    1998-10-15
 *
 *  \par Copyright (c) 1998-2000:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: PCA.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/11/19 15:45:39 $
 *
 *  \par Changes:
 *      $Log: PCA.cpp,v $
 *      Revision 2.1  2004/11/19 15:45:39  shark-admin
 *      eigensymm_obsolete(...) exchanged by eigensymm(...). (no sick internal modifications anymore of the matrix, which was called by reference)
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.7  2002/05/16 13:07:35  rudi
 *      doxygen comments added/modified.
 *
 *      Revision 1.6  2002/02/06 15:23:18  rudi
 *      Error in method "rtransform" fixed, new method for extraction of eigenvalues added.
 *
 *      Revision 1.5  2001/11/30 13:24:30  rudi
 *      Retransformation method totally reimplemented (didn't work before).
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


#include <cstdlib>
#include "Array/ArrayOp.h"
#include "LinAlg/linalg.h"
#include "LinAlg/PCA.h"


//===========================================================================
/*!
 *  \brief Constructs an "PCA" object with default values.
 *
 *  This constructor is defined in order to prevent the implicit definition
 *  by the compiler.
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
PCA::PCA( )
{
  modifiedM   = true;
  whiteningM  = false;
  removeMeanM = false;
  countM      = 0;
}


//===========================================================================
/*!
 *  \brief Destructs an "PCA" object.
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
PCA::~PCA()
{
}

//===========================================================================
/*!
 *  \brief Resets internal variables of a "PCA" object to initial values.
 *
 *  \return none.
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
void PCA::reset( )
{
    modifiedM = true;
    countM    = 0;
    meanM     = 0.;
    covarM    = 0.;
}


//===========================================================================
/*!
 *  \brief Updates the transformation matrix.
 *
 *  Based on the previous modified internal structures, the
 *  transformation matrix is updated.
 *  This method will be automatically called by #transform, #rtransform,
 *  #transMat and #eigVal, so you can be sure, that you are always
 *  using the most current transformation matrix.
 *
 *  \return none.
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
 *  \sa transform, rtransform, transMat, eigVal, train
 *
 */
void PCA::update( )
{
  // Only updating, when matrix is modified and
  // function "train" was called before: 
  if( modifiedM && meanM.nelem( ) > 0 && countM > 0 )
  {
    // The "real" covariance matrix:
    Array2D< double > cmatL( meanM.dim( 0 ), meanM.dim( 0 ) );
    // Allocate memory for eigenvalues-vector:
    dvecL.resize( meanM.dim( 0 ) );

    // Compute "real" covariance matrix:
    for( unsigned i = 0; i < meanM.dim( 0 ); ++i )
    {
      double t = meanM( i ) / countM;

            for( unsigned j = 0; j <= i; ++j )
			{
                cmatL( i, j ) = ( covarM( i, j ) - meanM( j ) * t ) / countM;
			}
        }

        //
        // compute transformation matrix
        //
		if( whiteningM )
		{
			rankDecomp( cmatL, transMatM, dvecL );
		}
		else
		{
		   eigensymm( cmatL, transMatM, dvecL );
		   //eigensymm_obsolete( cmatL, transMatM, dvecL );
		}

		modifiedM = false;
  }
}


//===========================================================================
/*!
 *  \brief If set to "true" a matrix will be created in method
 *  #update, that transforms a given sequence of data
 *  to white noise. 
 *
 *  Specifies, which kind of data transformation is performed
 *  at the next call of #transform. A value of \em true means,
 *  that a sequence of data values is transformed into white noise,
 *  i.e. the components of each transformed data vector is
 *  uncorrelated after the whitening. Additionally, the 
 *  variances of the transformed components are scaled to
 *  1.
 *  
 *  \param flagA \em true  - a whitening of data is performed,
 *               \em false - no whitening is performed.
 *  \return none
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
void PCA::setWhitening( bool flagA )
{
	whiteningM = flagA;
}


//===========================================================================
/*!
 *  setRemoveMean
 *
 *  \brief If set to "true" a matrix will be created in method
 *  #update, that will move the data in a way, that
 *  the mean value of the data lies in the zero-point.  
 *
 *  Specifies, which kind of data transformation is performed
 *  at the next call of #transform. A value of \em true means,
 *  that the mean value of the data lies in the zero-point
 *  of the coordinate system after the transformation.
 *  
 *  \param flagA \em true - the mean value of the data is moved,
 *               \em false - the mean value of the data is left untouched.
 *  \return none
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
void PCA::setRemoveMean( bool flagA )
{
	removeMeanM = flagA;
}




//===========================================================================
/*!
 *  \brief Calculates the estimates of the mean and the covariance
 *         matrix for the sample data "v".
 *
 *  The data for the estimation of the covariance and the mean matrix, that
 *  is used as base for the calculation of the transformation
 *  matrix is trained for the data in \em v. If \em v contains
 *  more than one vector, the training is performed for
 *  all data vectors. <br>
 *  Notice, that the transformation matrix itself is calculated only
 *  when calling method #update. This separation of the calculations
 *  is done due to efficiency. Guess, you are calling the train method
 *  at different times. Then the calculation of the transformation matrix
 *  will be performed only once at the end of the last training, instead
 *  of each time, the training method is called. 
 *
 *  \param  v Data vector(s) used for training.
 *  \return   none.
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em v is not one-
 *         or 2-dimensional
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
 *  \sa update
 *
 */
void PCA::train( const Array< double >& v )
{
  // Only for vector or 2-dimensional matrix:
  SIZE_CHECK( v.ndim( ) == 1 || v.ndim( ) == 2 )

  // 2-dimensional matrix:
  if ( v.ndim( ) == 2 )
  {
    // Calling "train" for each data vector:
    for ( unsigned i = 0; i < v.dim( 0 ); ++i )
    {
      train( v[ i ] );
    }
  } else {
  // Vector:
    if ( meanM.nelem( ) != v.dim( 0 ) )
    {
      if( countM > 0 )
      {
        throw 1;
      } else {
        meanM .resize( v.dim( 0 ) );
        meanM  = 0.;
        covarM.resize( v.dim( 0 ), v.dim( 0 ) );
        covarM = 0.;
      }
    }

    ++countM;

    for( unsigned i = 0; i < meanM.dim( 0 ); ++i )
    {
      double t = v( i );
      meanM( i ) += t;

      for( unsigned j = 0; j <= i; ++j )
      {
        covarM( i, j ) += v( j ) * t;
      }
    }

    modifiedM = true;
  }
}


//===========================================================================
/*!
 *  \brief Transformation of original data "v" to compressed data "w".
 *
 *  Using the original data in \em v and a previously calculated
 *  transformation matrix, the new compressed data \em w is produced.
 * 
 *  \param v       Original data.
 *  \param w       Compressed data.
 *  \param maxDimA Used to decrease the number of dimensions
 *                 for the compressed data; if set to zero
 *                 (= default value) the original number of dimensions
 *                 will be used .
 *  \return        none.
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
 *  \sa update
 *
 */
void PCA::transform
(
    const Array< double >& v,
    Array< double >& w,
    unsigned maxDimA
)
{
    unsigned                 i(0),         // counter
                             j(0),         // variables
                             v_dims,       // no. of dimensions of v.
                             v_dim1(0),    // x-dimension of v.
                             v_dim2(0),    // y-dimension of v.
                             mean_dim1(0); // x-dimension of mean value matrix.
    double                   t(0.);        // Stores single values for w.
    ArrayReference< double > wi;           // reference to a single 
                                           // transformed data vector


    v_dims = v.ndim();
    if (v_dims >= 1) v_dim1 = v.dim( 0 );
    if (v_dims >= 2) v_dim2 = v.dim( 1 );
    mean_dim1 = meanM.dim( 0 );
 
    // Check if transformation matrix exists:
    if ( meanM.nelem( ) == 0 || countM == 0 )
    {
        throw 1;
    }

    // Transform several data-vectors:
    if ( v_dims == 2 )
    {
        // Possibly set valid value for dimension reduction:
        if ( maxDimA == 0 || maxDimA > v_dim2 )
        {
            maxDimA = v_dim2;
        }

        // Allocate memory for transformed data-vectors: 
        w.resize( v_dim1, maxDimA );

        // Perform transformation for each single data vector:
        for ( i = 0; i < v_dim1; ++i )
        {
            wi.copyReference( w[ i ] );
            transform( v[ i ], static_cast< Array< double >& >( wi ), 
                       maxDimA );
        }
    }
    else
    {
    // Transform a single data-vector:

        // Update of transformation matrix for each single data-vector:
        update( );

        // Possibly set valid value for dimension reduction:
        if ( maxDimA == 0 || maxDimA > v_dim1 )
        {
            maxDimA = v_dim1;
        }

        // Allocate memory for transformed data-vector:  
        w.resize( maxDimA );

        if ( removeMeanM )
        {
	    // Move mean value of data to the zero point: 
            Array< double > vv( v - meanM / double( countM ) );

	    // Transform data vector by multiplying it with the
	    // transformation matrix:
            for ( i = maxDimA; i--; )
            {
                t = 0.;

                for ( j = mean_dim1; j--; )
                {
                    t += transMatM( j, i ) * vv( j );
                }
                w( i ) = t;
            }
        } 
	else
	{
	    // Transform data vector by multiplying it with the
	    // transformation matrix:
            for ( i = maxDimA; i--; )
            {
                t = 0.;

                for ( j = mean_dim1; j--; )
                {
                    t += transMatM( j, i ) * v( j );
                }
	        w( i ) = t;
            }
        }
    }
}


//===========================================================================
/*!
 *  \brief Retransformation of compressed data "v" to original data "w".
 *
 *  Using the compressed data in \em v and a previously calculated
 *  transformation matrix, the original data \em w is reconstructed.
 *  Notice, that there is a loss of information if a dimension reduction
 *  has taken place during the previous transformation.
 * 
 *  \param v       Compressed data.
 *  \param w       Original data.
 *  \param maxDimA Used to decrease the number of dimensions
 *                 for the original data; if set to zero
 *                 (= default value) the original number of dimensions
 *                 will be used .
 *  \return        none.
 *
 *  \author  R. Alberts
 *  \date    2001-11-07
 *
 *  \par Changes
 *      2003-07-18, ThB: <br>
 *      Dimension for resize of "w" and the limits for iteration
 *      were wrong. Bug fixed.
 *
 *      2002-01-30, ra: <br>
 *      Because of a misunderstanding, a retransformation was not
 *      possible for dimension reducted data before. This limitation
 *      was removed.
 *
 *  \par Status
 *      stable
 *
 *  \sa update
 *
 */
void PCA::rtransform
(
    const Array< double >& v,
    Array< double >& w,
    unsigned maxDimA
)
{
    unsigned                 i(0),         // counter
                             j(0),         // variables
                             v_dims,       // no. of dimensions of v.
                             v_dim1(0),    // x-dimension of v.
                             v_dim2(0),    // y-dimension of v.
                             mean_dim1(0); // x-dimension of mean value matrix.
    double                   t(0.);        // Stores single values for w.
    ArrayReference< double > wi;           // reference to a single 
                                           // retransformed data vector


    v_dims = v.ndim();
    if (v_dims >= 1) v_dim1 = v.dim( 0 );
    if (v_dims >= 2) v_dim2 = v.dim( 1 );
    mean_dim1 = meanM.dim( 0 );
 
    // Check if transformation matrix exists:
    if ( meanM.nelem( ) == 0 || countM == 0 )
    {
        throw 1;
    }

    // Retransform several data-vectors:
    if ( v_dims == 2 )
    {
        // Possibly set valid value for dimension reduction:
        if ( maxDimA == 0 )
        {
            maxDimA = v_dim2;
        }

        // Allocate memory for retransformed data-vectors: 
        w.resize( v_dim1, maxDimA );

        // Perform retransformation for each single data vector:
        for ( i = 0; i < v_dim1; ++i )
        {
            wi.copyReference( w[ i ] );
            rtransform( v[ i ], static_cast< Array< double >& >( wi ), 
                       maxDimA );
        }
    }
    else
    {
    // Retransform a single data-vector:

        // Update of transformation matrix for each single data-vector:
        update( );

        // Possibly set valid value for dimension reduction:
        if ( maxDimA == 0 )
        {
            maxDimA = v_dim1;
        }

        // Allocate memory for retransformed data-vector:  
        // w.resize( maxDimA );
	// ThB-CHANGE of prev. line:
        w.resize( mean_dim1 );

	// Retransform data vector by multiplying it with the
	// inverted (transposed because of orthogonality)
	// transformation matrix:
       
	// for ( i = maxDimA; i--; )
	// ThB: change of prev. line:
	for ( i = mean_dim1; i--; )
        {
            t = 0.;

           // for ( j = mean_dim1; j--; )
	   // ThB: change of prev. line:
           for ( j = maxDimA; j--; )
            {
	        if ( whiteningM )
		{
		    // For a retransformation after whitening
                    // the scaling of the eigenvectors must be
		    // undone additionally. "dvecL" contains
		    // the eigenvalues.
                    t += transMatM( i, j ) * v( j ) * dvecL( j );
		}
		else
		{
		    // "Normal" retransformation. 
                    t += transMatM( i, j ) * v( j );
		}
            }
	    if ( removeMeanM ) {
	        // Undo movement of mean value to zero point:
		w( i ) = t + meanM( i ) / double( countM );
            } else {
		w( i ) = t;
            }
        }
    }
}


//===========================================================================
/*!
 *  \brief Returns the transformation matrix.
 *
 *  The transformation matrix previously calculated is returned and
 *  possibly updated before.
 *
 *  \return The transformation matrix.
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
 *  \sa update
 *
 */
Array< double > PCA::transMat( )
{
	update( );
	return transMatM;
}



//===========================================================================
/*!
 *  \brief Returns the eigenvalues evaluated during the calculation
 *         of the transformation matrix.
 *
 *  The calculated eigenvalues are sometimes used for further processing.
 * 
 *  \return the eigenvalues.
 *
 *  \author  R. Alberts
 *  \date    2002-01-30
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa update
 *
 */
Array< double > PCA::eigVal( )
{
	update( );
	return dvecL;
}






