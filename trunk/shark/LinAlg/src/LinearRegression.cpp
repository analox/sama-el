//===========================================================================
/*!
 *  \file LinearRegression.cpp
 *
 *  \brief This file offers a class for creating a regression model that 
 *         will approximate data pairs (x,y) by finding a linear mapping for 
 *         them.
 *
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Copyright (c) 1995,1999:
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
 *      $RCSfile: LinearRegression.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:10 $
 *
 *  \par Changes:
 *      $Log: LinearRegression.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.12  2002/05/16 13:07:35  rudi
 *      doxygen comments added/modified.
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

#include "Array/ArrayOp.h"
#include "LinAlg/linalg.h"
#include "LinAlg/LinearRegression.h"


//===========================================================================
/*!
 *  \brief Creates a new Linear Regression object and resets
 *         internal values.
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
LinearRegression::LinearRegression( )
{
	//
	// the following initialization is not really necessary since
	// all comparision use sxM.nelem( ) instead of sxM.dim( 0 )
	//
	syM.resize( 0U );
	sxM.resize( 0U );

	reset( );
}


//===========================================================================
/*!
 *  \brief Destructs a Linear Regression object.
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
LinearRegression::~LinearRegression( )
{
}


//===========================================================================
/*!
 *  \brief Resets internal variables of the class to initial values.
 *
 *  \return none
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
void LinearRegression::reset( )
{
	modifiedM = true;
	countM    = 0;
	syM       = 0.;
	sxM       = 0.;
	rxxM      = 0.;
	rxyM      = 0.;
}


//===========================================================================
/*!
 *  \brief Calculates the current regression model.
 *
 *  Based on internal variables the current regression model
 *  is calculated if these internal variables have changed before
 *  due to adding/removing vector pair(s). This method is called
 *  automatically when you are calling the methods for returning
 *  the regression coefficients \f$A\f$ or \f$b\f$ or when
 *  you are calling the mehod for approximating an \f$y\f$ for a 
 *  given \f$x\f$ (methods #A, #b, #recall), so you can be sure
 *  that always the most current model is used. 
 *
 *  \return none
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
 *  \sa #train, #remove, #A, #b, #recall
 *
 */
void LinearRegression::update( )
{
    if( modifiedM && sxM.nelem( ) > 0 && countM > 0 )
	{
	    unsigned i, j;
		double   s, t;

        s = countM;

		Array  < double > dvecL( sxM.dim( 0 ) );
		Array  < double > mxL  ( sxM / s );
		Array  < double > myL  ( syM / s );
		Array2D< double > cxxL ( sxM.dim( 0 ), sxM.dim( 0 ) );
		Array2D< double > cxyL ( sxM.dim( 0 ), syM.dim( 0 ) );

        for( i = 0; i < sxM.dim( 0 ); ++i )
		{
            t = mxL( i );

            for( j = 0; j <= i; j++ )
			{
                cxxL( i, j ) = rxxM( i, j ) / s - t * mxL( j );
			}

            for( j = 0; j < syM.dim( 0 ); ++j )
			{
                cxyL( i, j ) = rxyM( i, j ) / s - t * myL( j );
			}
        }

        linearRegress( cxxL, cxyL, mxL, myL, amatM, bvecM, dvecL );

        modifiedM = false;
    }
}


//===========================================================================
/*!
 *  \brief Adds the data vector pair(s) "x"/"y" to the regression model.
 *
 *  Given the data vector \f$x\f$ and its corresponding \f$y\f$,
 *  the internal variables from which the model is calculated 
 *  are changed, i.e. the information given by \f$x\f$ and \f$y\f$
 *  is added to these variables. <br>
 *  You can also add several x/y-vector-pairs at once. <br>
 *  The model itself is not changed
 *  here, it will changed after calling method #update. <br>
 *  This separation is more efficient, e.g. when you are adding
 *  several data vectors at different time, because the model itself 
 *  is then calculated
 *  only once after all vectors are added instead of each time
 *  after a single vector pair is added.  
 *
 *  \param  x Vector(s) containing the first component of the data.
 *  \param  y Vector(s) containing the second component of the data.
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em x and \em y
 *         are not one- and 2-dimensional or that one of both
 *         is one-dimensional and the other one 2-dimensional
 *         or that the number of vectors in \em x is different
 *         to that in \em y
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
 *  \sa #update
 *
 */
void LinearRegression::train
(
	const Array< double >& x,
	const Array< double >& y
)
{
	SIZE_CHECK
	(
		( x.ndim( ) == 1 && y.ndim( ) == 1 ) ||
		( x.ndim( ) == 2 && y.ndim( ) == 2 )
	)

    if( x.ndim( ) == 2 )
	{
		SIZE_CHECK( x.dim( 0 ) == y.dim( 0 ) )

        for( unsigned i = 0; i < x.dim( 0 ); ++i )
		{
			train( x[ i ], y[ i ] );
		}
    }
	else
	{
		if( sxM.nelem( ) != x.dim( 0 ) || syM.nelem( ) != y.dim( 0 ) )
		{
			if( countM > 0 )
			{
				throw 1;
			}
			else
			{

				sxM  .resize( x.dim( 0 ) );
				syM  .resize( y.dim( 0 ) );
				rxxM .resize( x.dim( 0 ), x.dim( 0 ) );
				rxyM .resize( x.dim( 0 ), y.dim( 0 ) );
				amatM.resize( y.dim( 0 ), x.dim( 0 ) );
				bvecM.resize( y.dim( 0 ) );

				sxM  = 0.;
				syM  = 0.;
				rxxM = 0.;
				rxyM = 0.;
			}
		}

		//
		// train linear regression model
		//
		++countM;

		syM += y;

		for( unsigned i = 0, j; i < sxM.dim( 0 ); ++i )
		{
			double t  = x( i );
			sxM( i ) += t;

			for( j = 0; j <= i; ++j )
			{
				rxxM( i, j ) += t * x( j );
			}

			for( j = 0; j < syM.dim( 0 ); ++j )
			{
				rxyM( i, j ) += t * y( j );
			}
		}

		modifiedM = true;
	}
}



//===========================================================================
/*!
 *  \brief Removes the data vector pair(s) "x"/"y" from the regression model.
 *
 *  Given the data vector \f$x\f$ and its corresponding \f$y\f$,
 *  the internal variables from which the model is calculated 
 *  are changed, i.e. the information given by \f$x\f$ and \f$y\f$
 *  is removed from these variables. <br>
 *  You can also remove several x/y-vector-pairs at once. <br>
 *  The model itself is not changed
 *  here, it will changed after calling method #update. <br>
 *  This separation is more efficient, e.g. when you are removing
 *  several data vectors at different time, because the model itself 
 *  is then calculated
 *  only once after all vectors are removed instead of each time
 *  after a single vector pair is removed.  
 *
 *  \param  x Vector(s) containing the first component of the data.
 *  \param  y Vector(s) containing the second component of the data.
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em x and \em y
 *         are not one- and 2-dimensional or that one of both
 *         is one-dimensional and the other one 2-dimensional
 *         or that the number of vectors in \em x is different
 *         to that in \em y or that the model was not
 *         calculated yet or that it was not calculated for
 *         as many vectors as there are vectors in \em x and \em y
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
 *  \sa #update
 *
 */
void LinearRegression::remove
(
	const Array< double >& x,
	const Array< double >& y
)
{
	SIZE_CHECK
	(
		( x.ndim( ) == 1 && y.ndim( ) == 1 ) ||
		( x.ndim( ) == 2 && y.ndim( ) == 2 )
	)

    if( x.ndim( ) == 2 )
	{
		SIZE_CHECK( x.dim( 0 ) == y.dim( 0 ) )

        for( unsigned i = 0; i < x.dim( 0 ); ++i )
		{
			remove( x[ i ], y[ i ] );
		}
    }
	else
	{
		SIZE_CHECK
		(
			countM > 0 &&
			sxM.dim( 0 ) == x.dim( 0 ) &&
			syM.dim( 0 ) == y.dim( 0 )
		)

		//
		// remove vector from linear regression model
		//
		--countM;

		syM -= y;

		for( unsigned i = 0, j; i < sxM.dim( 0 ); ++i )
		{
			double t  = x( i );
			sxM( i ) -= t;

			for( j = 0; j <= i; ++j )
			{
				rxxM( i, j ) -= t * x( j );
			}

			for( j = 0; j < syM.dim( 0 ); ++j )
			{
				rxyM( i, j ) -= t * y( j );
			}
		}

		modifiedM = true;
	}
}


//===========================================================================
/*!
 *  \brief Given one/several data vector(s) "x" the corresponding 
 *         vector(s) "y" is/are calculated by using the current regression 
 *         model.
 *
 *  The current regression model \f$y = A \cdot x + b\f$ is used to 
 *  approximate \f$y\f$ for one/several given data vector(s) \f$x\f$.
 *
 *  \param  x Input vector(s) for which the corresponding y-vector(s)
 *            will be calculated.
 *  \param  y The corresponding (approximated) y-vector(s) for the
 *            input vector(s).
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em x is not
 *         one- or 2-dimensional or that the model was not calculated
 *         for any training vector yet 
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
 *  \sa #update
 *
 */
void LinearRegression::recall
(
	const Array< double >& x,
	Array< double >& y
)
{
	SIZE_CHECK( x.ndim( ) == 1 || x.ndim( ) == 2 )
    SIZE_CHECK( sxM.nelem( ) > 0 && countM > 0 )

    if( x.ndim( ) == 2 )
	{
		y.resize( x.dim( 0 ), syM.dim( 0 ) );

		for( unsigned i = 0; i < x.dim( 0 ); ++i )
		{
			ArrayReference< double > yi;
			yi.copyReference( y[ i ] );
			recall( x[ i ], static_cast< Array< double >& >( yi ) );
		}
    }
	else
	{
		update( );

		y.resize( syM.dim( 0 ) );

        for( unsigned i = 0; i < syM.dim( 0 ); ++i )
		{
            double t = bvecM( i );

            for( unsigned j = 0; j < sxM.dim( 0 ); ++j )
			{
                t += amatM( i, j ) * x( j );
			}

            y( i ) = t;
        }
	}
}


//===========================================================================
/*!
 *  \brief Temporarily removes the vector pair(s) "x"/"y" from the
 *         regression model and uses this modified model to approximate
 *         the vector(s) "y", when the vector(s) "x" is/are given.
 *
 *  This method can be used to check if some vector pairs \f$(x, y)\f$
 *  are necessary as information for the training of the regression model
 *  or not. <br>
 *  After the call of this method, the regression model is reset to the
 *  further state, where the \f$(x, y)\f$ pair(s) were included in the
 *  model.
 *
 *  \param x    The x-vector(s) that will be temporarily removed from the
 *              model. 
 *  \param y    The corresponding y-vector(s) that will be temporarily 
 *              removed from the model. 
 *  \param yout The y-vector(s) for the given vector(s) \em x, approximated
 *              by the modified model.
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em x and \em y
 *         are not one- and 2-dimensional or that one of both
 *         is one-dimensional and the other one 2-dimensional
 *         or that the number of vectors in \em x is different
 *         to that in \em y
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
void LinearRegression::leaveOneOut
(
	const Array< double >& x,
	const Array< double >& y,
	Array< double >&       yout
)
{
	SIZE_CHECK
	(
		( x.ndim( ) == 1 && y.ndim( ) == 1 ) ||
		( x.ndim( ) == 2 && y.ndim( ) == 2 )
	)

    SIZE_CHECK( sxM.nelem( ) > 0 && countM > 0 )

    if( x.ndim( ) == 2 )
	{
		SIZE_CHECK( x.dim( 0 ) == y.dim( 0 ) )

		yout.resize( y );

		for( unsigned i = 0; i < x.dim( 0 ); ++i )
		{
			ArrayReference< double > youti;
			youti.copyReference( yout[ i ] );
			leaveOneOut( x[ i ], y[ i ], static_cast< Array< double >& >( youti ) );
		}
    }
	else
	{
		remove( x, y );
		recall( x, yout );
		train( x, y );
	}
}


//===========================================================================
/*!
 *  \brief Returns the regression coefficient "A" of the regression model.
 *
 *  Given the current regression model \f$y = A \cdot x + b\f$, 
 *  this method returns \f$A\f$.
 *
 *  \return The regression coefficient \f$A\f$.
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
Array< double > LinearRegression::A( )
{
	update( );
	return amatM;
}


//===========================================================================
/*!
 *  \brief Returns the regression coefficient "b" of the regression model.
 *
 *  Given the current regression model \f$y = A \cdot x + b\f$, 
 *  this method returns \f$b\f$.
 *
 *  \return The regression coefficient \f$b\f$.
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
Array< double > LinearRegression::b( )
{
	update( );
	return bvecM;
}









