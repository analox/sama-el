//===========================================================================
/*!
 *  \file LinearClassifier.cpp
 *
 *  \brief This file offers a class with methods, that can perform a data 
 *         classification by using linear discriminant analysis (LDA) for 
 *         calculating a transformation of the data. 
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
 *      $RCSfile: LinearClassifier.cpp,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/07/09 12:14:18 $
 *
 *  \par Changes:
 *      $Log: LinearClassifier.cpp,v $
 *      Revision 2.1  2004/07/09 12:14:18  saviapbe
 *      A bag with the memory allocation was removed.
 *
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:07:35  rudi
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

#include "ZNconfig.h"
#include <cmath>
#include "Array/ArrayOp.h"
#include "LinAlg/linalg.h"
#include "LinAlg/LinearClassifier.h"


//! Initial value for the different number of classes.
static const unsigned initialClassNumberC = 50;

//! The number of classes that is added to the existing ones
//! when the class number of the current training vector
//! exceeds the number of classes already known.
static const unsigned classIncrC = 10;


//===========================================================================
/*!
 *  \brief Constructs a new Linear Classifier object with initial values.
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
LinearClassifier::LinearClassifier( ) : 
    classCountM( initialClassNumberC ),
    classAbsM  ( initialClassNumberC ),
    classRelM  ( initialClassNumberC ),
    classMeanM ( initialClassNumberC ),
    transMeanM ( initialClassNumberC )
{
    reset( );
}


//===========================================================================
/*!
 *  \brief Destructs an Linear Classifier object.
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
LinearClassifier::~LinearClassifier( )
{
}


//===========================================================================
/*!
 *  \brief Resets the current Linear Classifier object by assigning
 *         initial values to internal variables.
 *
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
void LinearClassifier::reset( )
{
    modifiedM    = true;
    dimM         = 0;
    outDimM      = 0;
    usedClassesM = 0;
    countM       = 0;
    classCountM  = 0U;
    meanM        = 0.;
    covarM       = 0.;
}


//===========================================================================
/*!
 *  \brief Vector "v" is used for training the classificator. 
 *
 *  Vector \em v that contains the characteristics of one individual
 *  is used for the training of the classification function.
 *  This original data with the class it is assigned to is given
 *  to modify the classificator. <br>
 *  In this method only internal variables from which the classificator
 *  is evaluated are calculated. The classificator
 *  itself is calculated by method #update. This separation of
 *  the calculation has the sense of increasing the calculation
 *  time when the training is called different times. Then, the
 *  classificator will only be calculated once at the end of the
 *  training instead of each time when the train method will be called.
 *
 *  \param v the training vector
 *  \param m the number of the class the individual belongs to
 *  \return none
 *  \throw check_exception the type of the exception will be "size mismatch"
 *         and indicates that \em v is not one-dimensional
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
void LinearClassifier::train
(
    const Array< double >& v,
    unsigned m
)
{
    // Array must be a vector:
    SIZE_CHECK( v.ndim( ) == 1 )

    unsigned c;


    // Dimension of training vector must be equal to the number
    // of attributes used for classification:	
    if ( dimM != v.dim( 0 ) )
    {
        if ( countM > 0 )
        {
	    // Not the first training vector, so this is definitely
	    // an error:
            throw 1;
        }
	else
	{
	    // This is the first training vector, so adapt
	    // internal data structures to its size:
            dimM   = v.dim( 0 );
            meanM .resize( dimM );
            covarM.resize( dimM, dimM );
            meanM  = 0.;
            covarM = 0.;
	}
    }   

    // If the class of the training vector exceeds the number
    // of known classes, space for new classes is created:
    if ( m >= classCountM.nelem( ) )
    {
        unsigned newSizeL = classCountM.nelem( ) + classIncrC;
        classCountM.resize( newSizeL, true );
        classAbsM  .resize( newSizeL, true );
        classRelM  .resize( newSizeL, true );
        classMeanM .resize( newSizeL, true );
        transMeanM .resize( newSizeL, true );

	 for( unsigned int i = classCountM.nelem( )-1;
             i >= classCountM.nelem( ) - classIncrC;
             i-- )
          classCountM( i ) = 0;
    }

  

    // Does the training vector belongs to a class
    // unknown before?
    if ( classCountM( m ) == 0 )
    {
        // New class, so assign a new save index
        // for this class and update internal variables: 
        c = usedClassesM++;

        classAbsM( c ) = m;
        classRelM( m ) = c;

        classMeanM( c ).resize( dimM );
        classMeanM( c ) = 0.;
    }
    else
    {
        // Class is already known, determine save index:
        c = classRelM( m );
    }

    // Update training vectors and class members counter:
    ++countM;
    ++classCountM( m );


    Array< double >& classMeanL = classMeanM( c );
    for ( unsigned i = 0; i < dimM; ++i )
    {
        double t = v( i );
        classMeanL( i ) += t;
        meanM     ( i ) += t;

        for ( unsigned j = 0; j <= i; ++j )
        {
            covarM( i, j ) += v( j ) * t;
        }
    }

    // The values on which the classificator is based was
    // modified, so the classificator needs an update:
    modifiedM = true;
}



//===========================================================================
/*!
 *  \brief The vectors in "vectorsA" are used for training the classificator. 
 *
 *  Similar to method #train(const Array< double >& v,unsigned m), but
 *  here not one vector, but several vectors are used for training. 
 *
 *  \param vectorsA Array containing the training vectors.
 *  \param classNumsA Vector with class numbers that correspond to
 *         the training vectors in \em vectorsA. classNumsA(i) is
 *         the number of the class, training vector vectorsA(i) belongs to.
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em vectorsA is not
 *         2-dimensional or \em classNumsA is not one-dimensional
 *         or that the number of class numbers in \em classNumsA
 *         is different from the number of training vectors in \em vectorsA
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
void LinearClassifier::train
(
	const Array< double >& vectorsA,
	const Array< unsigned >& classNumsA
)
{
	SIZE_CHECK
	(
		vectorsA  .ndim( ) == 2 &&
		classNumsA.ndim( ) == 1 &&
		classNumsA.dim( 0 ) == vectorsA.dim( 0 )
	)

    for( unsigned i = 0; i < vectorsA.dim( 0 ); ++i )
	{
		train( vectorsA[ i ], classNumsA( i ) );
	}
}


//===========================================================================
/*!
 *  \brief Classifies vector "v".
 * 
 *  Based on the characteristics stored in vector \em v, the trained 
 *  classificator will assign the individual with these characteristics 
 *  to a specific class.
 * 
 *  \param  v the vector that will be classified.
 *  \return the number of the class to which the classificator
 *          has assigned the vector.
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
unsigned LinearClassifier::classify
(
	const Array< double >& v
)
{
    double     d, t;
    double     dist = -1;
    unsigned   c, i, j, m = 0;

    if( usedClassesM == 1 )
	{
        m = classAbsM( 0 );
	}
    else if( dimM && usedClassesM >= 2)
    {
		update( );

		Array< double > dvecL( outDimM );

		//
		// Lc->d = Lc->TransMat * v
		//
		for( i = outDimM; i--; )
		{
			t = 0.;

			for( j = dimM; j--; )
			{
				t += transMatM( i, j ) * v( j );
			}

			dvecL( i ) = t;
		}


		distVecM.resize( usedClassesM );

		for( c = usedClassesM; c--; )
		{
			//
			// Abstand zu Klasse 'c'
			//
			d = sqrDistance( dvecL, transMeanM( c ) );

			//
			// Abstand kuerzer ?
			//
			if( dist < 0 || d < dist )
			{
				dist  = d;
				m     = classAbsM( c );
			}

			//
			// Abstandsquadrate speichern
			//
			distVecM( c ) = d;
		}
    }

    return m;
}


//===========================================================================
/*!
 *  \brief The vectors in "vectorsA" are classified.
 *
 *  Similar to method #classify(const Array< double >& v), but
 *  here not one vector, but several vectors are classified.
 *
 *  \param vectorsA Array containing the vectors, that will be
 *                  classified.
 *  \param classNumsA Vector with the resulting class numbers that
 *                    will be assigned to the single vectors
 *                    by the classificator.
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em vectorsA is
 *         not 2-dimensional
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
void LinearClassifier::classify
(
	const Array< double >& vectorsA,
	Array< unsigned >& classNumsA
)
{
	SIZE_CHECK( vectorsA.ndim( ) == 2 )

    if( dimM == 0 || countM == 0 )
	{
		throw 1;
	}

	classNumsA.resize( vectorsA.dim( 0 ) );

	for( unsigned i = 0; i < vectorsA.dim( 0 ); ++i )
	{
		classNumsA( i ) = classify( vectorsA[ i ] );
	}
}


//===========================================================================
/*!
 *  \brief Performs a "leaving-one-out" test with vector "v" taken
 *         out of the training set.
 *
 *  A \em leaving-one-out test is a special case of cross validation
 *  and used for a realistic evaluation of the classificator.
 *  For this evaluation, vector \em v is taken out of the training
 *  set and the classificator is updated. Then vector \em v is
 *  used as test vector for the new classificator and the result
 *  of this classification is returned. Last, the vector is 
 *  added to the training set again. To evaluate the error rate
 *  of a classificator you must perform the leaving-one-out test
 *  for all possible vectors in the training set, i.e. having
 *  \f$n\f$ vectors in the training set, the error \f$E\f$
 *  is calculated by:
 *
 *  \f$
 *      E = \frac{1}{n} \sum_{i=1}^n | y_i - C_{z(i)}(x_i)|
 *  \f$
 *
 *  Here \f$x_i\f$ is the vector, that is taken out of the training
 *  set, \f$y_i\f$ is the class \f$x_i\f$ belongs to and 
 *  \f$C_{z(i)}\f$ is the result of the classification of vector
 *  \f$x_i\f$ for the classificator based on the training
 *  set \f$z(i) = \{x_1, \dots, x_{i-1}, x_{i+1}, x_n\}\f$.
 *
 *  \param v vector, that is taken out of the training set
 *  \param m the class to which \em v belongs
 *  \return the class assigned to \em v by the classificator
 *          based on the training set without \em v
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em v is not
 *         one-dimensional or that its size is wrong or
 *         that \em m is no valid class number or was not
 *         used before 
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
unsigned LinearClassifier::leaveOneOut
(
	const Array< double >& v,
	unsigned m
)
{
	SIZE_CHECK
	(
		v.ndim( ) == 1 &&
		v.dim( 0 ) == dimM &&
		m < classCountM.nelem( ) &&
		classCountM( m ) > 0
	)

    unsigned  i, j, c, cn;
    double    t;

	//
	// relative Klassennummer holen
	//
	c = classRelM( m );

	//
	// Vektor aus dem Trainingssatz nehmen
	//
	--countM;
	--classCountM( m );

	Array< double >& classMeanL = classMeanM( c );
	for( i = 0; i < dimM; ++i )
	{
		t = v( i );
		classMeanL( i ) -= t;
		meanM     ( i ) -= t;

		for( j = 0; j <= i; ++j )
		{
			covarM( i, j ) -= v( j ) * t;
		}
	}

	modifiedM = true;
	cn = classify( v );

	//
	// Vektor wieder zum Trainingssatz hinzufuegen
	//
	++countM;
	++classCountM( m );

	for( i = 0; i < dimM; ++i )
	{
		t = v( i );
		classMeanL( i ) += t;
		meanM     ( i ) += t;

		for( j = 0; j <= i; ++j )
		{
			covarM( i, j ) += v( j ) * t;
		}
	}

	modifiedM = true;

	return cn;
}


//===========================================================================
/*!
 *  \brief Performs a "leaving-one-out" test for all
 *         vectors in "vectorsA".
 *
 *  Does the same as method #leaveOneOut(const Array< double >& v,unsigned m)
 *  but for all vectors in \em vectorsA.
 *
 *  \param vectorsA the vectors for which the test will be performed
 *  \param classNumsA the classes, the vectors in \em vectorsA
 *         belong to
 *  \param resultNumsA the results of the classification for the
 *         vectors in \em vectorsA
 *  \return none
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em vectorsA is
 *         is not 2-dimensional or that \em classNumsA is not
 *         one-dimensional or that the number of class numbers
 *         in \em classNumsA don't correspond to the number
 *         of vectors in \em vectorsA  
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
void LinearClassifier::leaveOneOut
(
	const Array< double >& vectorsA,
	const Array< unsigned >& classNumsA,
	Array< unsigned >& resultNumsA
)
{
	SIZE_CHECK
	(
		vectorsA.ndim( ) == 2 &&
		classNumsA.ndim( ) == 1 &&
		vectorsA.dim( 0 ) == classNumsA.dim( 0 )
	)

    if( dimM == 0 || countM == 0 )
	{
		throw 1;
	}

	resultNumsA.resize( vectorsA.dim( 0 ) );

	for( unsigned i = 0; i < vectorsA.dim( 0 ); ++i )
	{
		resultNumsA( i ) = leaveOneOut( vectorsA[ i ], classNumsA( i ) );
	}
}


//===========================================================================
/*!
 *  \brief Updates the classificator.
 *
 *  Here the current classificator is calculated by using the
 *  internal variables that were changed before. <br>
 *  This update is performed automatically by the #classify methods
 *  when a change of the internal variables has taken place
 *  before, so you can be sure, that always the most current
 *  classificator is used.
 *
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
 *  \sa classify, train
 *
 */
void LinearClassifier::update( )
{
    unsigned  i, j, c;
    double    cnt, t;

    if( modifiedM && dimM > 0 && usedClassesM >= 2 )
    {
		Array2D< double > betweenCovarL( dimM, dimM );
		Array2D< double > withinCovarL ( dimM, dimM );
		Array  < double > dvecL( dimM );

        //
        // Berechnung von B und W
        //
        for( i = 0; i < dimM; ++i )
		{
            for (j = 0; j <= i; j++)
			{
                betweenCovarL( i, j ) = 0.;
			}
		}

        for( c = usedClassesM; c--; )
		{
			//
			// Anzahl fuer absolute Klassennummer holen
			//
			cnt = double( classCountM( classAbsM( c ) ) );
			Array< double >& classMeanL = classMeanM( c );
			for( i = 0; i < dimM; ++i )
			{
				t = classMeanL( i ) / cnt;

				for( j = 0; j <= i; ++j )
				{
					betweenCovarL( i, j ) += classMeanL( j ) * t;
				}
			}
		}

        for( i = 0; i < dimM; ++i )
        {
            t = meanM( i ) / countM;

            for( j = 0; j <= i; ++j )
            {
                withinCovarL ( i, j )  = covarM( i, j ) - betweenCovarL( i, j );
                betweenCovarL( i, j ) -= meanM( j ) * t;
            }
		}

		// normalization necessary?
        for( i = 0; i < dimM; ++i )
		{
            for( j = 0; j <= i; ++j )
            {
                withinCovarL ( i, j ) /= countM;
                betweenCovarL( i, j ) /= countM;
            }
		}

		outDimM = dimM < usedClassesM - 1 ? dimM : usedClassesM - 1;

        discrimAnalysis
		(
			betweenCovarL,
			withinCovarL,
			transMatM,
			dvecL,
            outDimM
		);

        //
        // transformierte Klassenvektoren vorberechnen
        //
        for( c = usedClassesM; c--; )
		{
			Array< double >& classMeanL = classMeanM( c );
			Array< double >& transMeanL = transMeanM( c );
			transMeanL.resize( outDimM );

			//
			// Anzahl fuer absolute Klassennummer holen
			//
			cnt = double( classCountM( classAbsM( c ) ) );

			for( i = outDimM; i--; )
			{
				t = 0.0;

				for( j = dimM; j--; )
				{
					t += transMatM( i, j ) * classMeanL( j );
				}

				transMeanL( i ) = t / cnt;
			}
		}

        modifiedM = false;
    }
}



