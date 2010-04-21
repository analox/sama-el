//===========================================================================
/*!
 *  \file LinearClassifier.h
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
 *      $RCSfile: LinearClassifier.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: LinearClassifier.h,v $
 *      Revision 2.1  2004/06/16 15:22:20  saviapbe
 *
 *      Some bugs in the doxygen documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:20:31  rudi
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

#ifndef LINEARCLASSIFIER_H
#define LINEARCLASSIFIER_H

#include "ZNconfig.h"
#include <list>
#include "Array/Array2D.h"


//===========================================================================
/*!
 *  \brief With the methods in this class you can perform a data classification
 *         by using linear discriminant analysis (LDA) for calculating
 *         a transformation of the data.
 *
 *  If you have data that can be assigned to different classes, you
 *  need an efficient way to calculate a separation between the
 *  different classes, so you can easily classify new data. <br>
 *  Given two classes consider the separation as line between the
 *  two data sets, where all data points of set 1 are on one
 *  side of the line and all data points of set 2 are on the other
 *  side. This can be easily generalized for \f$m\f$ classes, where
 *  you have several separation lines. <br>
 *  LDA guarantees maximal separability of the classes by maximizing
 *  the ratio of between-class variance to within-class variance. <br>
 *  The location of the original data sets doesn't change, instead LDA
 *  tries to provide more class separability and draw a decision
 *  region between the given classes. <br>
 *  In this class you will find methods for training a classificator
 *  (a transformation matrix that will map data to the discriminant
 *  space with maximum class separability) using this classificator
 *  for classifying new data and testing the importance of data
 *  of the training set for the calculation of the classificator. <br>
 *  For more details about LDA, please refer to the reference of
 *  function discrimAnalysis, that is used as subroutine here. <br>
 *  The following example will show how the class and some of its
 *  methods are used: <br>
 *
 *  \example linearClassifier_test.cpp
 *
 *  You can find, compile and run this example in 
 *  "$SHARK_ROOTDIR/LinAlg/examples".
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 *  
 *  \sa discrimAnalysis
 *
 */
class LinearClassifier
{
  public:


      //! Constructs a new Linear Classifier object with initial values.
      LinearClassifier( );

      //! Destructs an Linear Classifier object.
      ~LinearClassifier( );

      //! Resets the current Linear Classifier object by assigning
      //! initial values to internal variables.
      void reset( );

      //! Updates the classificator.
      void update( );

      //! Vector "v" is used for training the classificator. 
      void train
      (
	  const Array< double >&,
	  unsigned
      );

      //! The vectors in "vectorsA" are used for training the classificator. 
      void train
      (
          const Array< double >&,
	  const Array< unsigned >&
      );

      //====================================================================== 
      /*!
       *  \brief The vectors in list "vectorsA" are used for training the 
       *         classificator. 
       *
       *  Similar to method #train(const Array< double >& v,unsigned m), but
       *  here not one vector, but several vectors, that are stored in
       *  a list, are used for training. The class numbers are also
       *  stored in a list. 
       *
       *  \param vectorsA List of training vectors.
       *  \param classNumsA List of class numbers that correspond to
       *         the training vectors in \em vectorsA. classNumsA(i) is
       *         the number of the class, training vector vectorsA(i) belongs 
       *         to.
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
       *  \sa update
       *
       */   
      void train
      (
          const std::list< Array< double > >& vectorsA,
	  const std::list< unsigned >&        classNumsA
      )
      {
          std::list< Array< double > >::const_iterator vecIterL;
	  std::list< unsigned >       ::const_iterator classIterL;
 
	  for ( vecIterL = vectorsA.begin( ),   
                classIterL = classNumsA.begin( );
		vecIterL != vectorsA.end( ) && classIterL != classNumsA.end( );
		++vecIterL, ++classIterL
	      )
	  {
	      train( *vecIterL, *classIterL );
	  }
      }

      //! Classifies vector "v".
      unsigned classify( const Array< double >& );

      //! The vectors in "vectorsA" are classified.
      void classify( const Array< double >&, Array< unsigned >& );

      //====================================================================== 
      /*!
       *  \brief The vectors in list "vectorsA" are classified.
       *
       *  Similar to method #classify(const Array< double >& v), but
       *  here not one vector, but several vectors, that are stored in
       *  a list, are classified. The resulting class numbers are also
       *  stored in a list. 
       *
       *  \param vectorsA List of vectors, that will be classified.
       *  \param classNumsA List of resulting class numbers.
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
       *  \sa update
       *
       */   
      void classify
      (
          const std::list< Array< double > >& vectorsA,
	  std::list< unsigned >&              classNumsA
      )
      {
	  std::list< Array< double > >::const_iterator vecIterL;

	  classNumsA.clear( );
	  for ( vecIterL  = vectorsA.begin( );
		vecIterL != vectorsA.end  ( );
		++vecIterL
	      )
	  {
	      classNumsA.push_back( classify( *vecIterL ) );
	  }
      }


      //! Performs a "leaving-one-out" test with vector "v" taken out of 
      //! the training set.
      unsigned leaveOneOut( const Array< double >&, unsigned );


      //! Performs a "leaving-one-out" test for all vectors in "vectorsA".
      void leaveOneOut
      (
          const Array< double >&,
	  const Array< unsigned >&,
	  Array< unsigned >&
      );

      //======================================================================
     /*!
      *  \brief Performs a "leaving-one-out" test for all
      *         vectors in the list "vectorsA".
      *
      *  Similar to method 
      *  #leaveOneOut(const Array< double >& v,unsigned m)
      *  but here all vectors in the list \em vectorsA are used.
      *  The predefined class numbers and the results of the
      *  classifications are also saved in lists.
      *
      *  \param vectorsA the list of vectors for which the test will be 
      *                  performed
      *  \param classNumsA the list of classes, the vectors in \em vectorsA
      *                    belong to
      *  \param resultNumsA the list of results of the classification for the
      *                     vectors in \em vectorsA
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
      void leaveOneOut
      (
          const std::list< Array< double > >& vectorsA,
	  const std::list< unsigned >&        classNumsA,
	  std::list< unsigned >&              resultNumsA
      )
      {
          std::list< Array< double > >::const_iterator vecIterL;
	  std::list< unsigned >       ::const_iterator classIterL;

	  resultNumsA.clear( );
	  for ( vecIterL = vectorsA.begin( ),   
                classIterL = classNumsA.begin( );
		vecIterL != vectorsA.end( ) && classIterL != classNumsA.end( );
			++vecIterL, ++classIterL
	      )
	  {
	      resultNumsA.push_back( leaveOneOut( *vecIterL, *classIterL ) );
	  }
      }



  protected:

    //! Is set to true, when the classificator needs an update
    //! (e.g. after training).
    bool modifiedM;

    //! The number of predefined attributes used for the classification
    //! of each individual. Because of the fact that training vectors
    //! enumerate the characteristics of individuals refering to the
    //! attributes, the dimension of each vector must be equal to
    //! the value of "dimM".
    unsigned dimM;

    //! When a training vector belongs to a class that is unknown
    //! then space not only for this but for several classes
    //! is added to the existing one to speed up the whole
    //! training process. To keep track of the space that is
    //! really used, "usedClassesM" saves the number of classes
    //! that are used.
    unsigned usedClassesM;

    //! The dimension of the discriminant space, which is equal
    //! to the number of used eigenvectors.
    unsigned outDimM;

    //! The number of vectors that were already used for the training.
    unsigned countM;

    //! The number of members for each class. 
    Array< unsigned > classCountM;

    //! Hashtable. The classes are saved here in the order they
    //! become known to the training method at run time. <br>
    //! classAbsM[ i ] will return the "absolute number" of the 
    //! the i-th introduced class. <br>
    //! This is the opposite hashtable to #classRelM.
    Array< unsigned > classAbsM;

    //! Hashtable. The classes are saved in the order they become
    //! known to the training method at run time. This hashtable
    //! maps the "absolute number" of the classes to the "relative
    //! number" under which they were saved. <br>
    //! classRelM[ i ] will return the relative number of class no. i, i.e.
    //! the index under which it is saved. This index
    //! then can be used to access data about this class
    //! which is saved in other structures under this index. <br>
    //! This is the opposite hashtable to #classAbsM.
    Array< unsigned > classRelM;

    //! The sense of this Array is unknown, because the only value
    //! that is ever assigned to it is zero.
    Array< Array< double > > classMeanM;

    //! The sense of this Array is unknown, because only 
    //! memory is allocated for it but no value is ever
    //! assigned to it.
    Array< Array< double > > transMeanM;

    //! If \f$n\f$ data vectors were already used for the training
    //! and each single vector has \f$d\f$ dimensions, then 
    //! \f$meanM_i = \sum_{j=1}^n {v_j}_i \mbox{,\ for\ } i = 1, \dots, d\f$
    Array< double > meanM;

    //! If \f$n\f$ data vectors were already used for the training
    //! and each single vector has \f$d\f$ dimensions, then this is
    //! a \f$d \times d\f$ matrix with 
    //! \f$covarM_{i, j} = \sum_{k=1}^n {v_k}_i \cdot {v_k}_j \mbox{,\ for\ } i, j = 1, \dots, d\f$
    Array2D< double > covarM;

    //! The transformation matrix, that will project the data
    //! to the discriminant space.
    Array2D< double > transMatM;

    //! Contains the euclidian distance from all classes
    //! to the last classified vector.
    Array< double > distVecM;

};

#endif  // LINEARCLASSIFIER_H




