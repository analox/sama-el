//===========================================================================
/*!
 *  \file LinearRegression.h
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
 *      <BR>
 *
 *  \par Project:
 *      LinAlg
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: LinearRegression.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2005/06/21 14:15:10 $
 *
 *  \par Changes:
 *      $Log: LinearRegression.h,v $
 *      Revision 2.2  2005/06/21 14:15:10  glasmtbl
 *
 *
 *      'private' members changed to 'protected' allowing for greater flexibility.
 *
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
 *      Revision 1.3  2002/05/16 13:20:31  rudi
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


#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H

#include "Array/Array2D.h"

//===========================================================================
/*!
 *  \brief Methods for creating a regression model that will approximate
 *         data pairs (x,y) by finding a linear mapping for them.
 *
 *  If you have data that consists of two components \f$x\f$ and \f$y\f$
 *  of different types where you know or suspect a correlation between 
 *  these two components (e.g. flow resistance and speed), you can use
 *  linear regression to describe this correlation and to summarize
 *  your data. <br>
 *  Linear regression will analyze the correlation of given input data
 *  and build a model for it, i.e. to create a linear mapping
 *  \f$y = A \cdot x + b\f$, that will approximate the data. <br>
 *  This class offers methods to build such a regression model
 *  for multidimensional data vectors \f$x\f$ and \f$y\f$. <br>
 *  You can add and remove data vector pairs \f$(x, y)\f$ to/from
 *  the model and test, how important a special vector pair is
 *  for the creation of a good approximation of the whole data. <br>
 *  If you are satisfied with your model, you can use it to approximate
 *  the data component \f$y\f$ for a given \f$x\f$. <br>
 *  The following example will show how the class and some of its
 *  methods are used: <br>
 *
 *  \example linearRegression_test.cpp
 *
 *  You can find, compile and run this example in 
 *  "$SHARK_ROOTDIR/LinAlg/examples".
 *  
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class LinearRegression
{
  public:

    //! Creates a new empty Linear Regression object and resets
    //! internal values.
    LinearRegression( );

    //! Destructs a Linear Regression object.
    ~LinearRegression( );

    //! Resets internal variables of the class to initial values.
    void reset( );

    //! Calculates the current regression model.
    void update ( );

    //! Adds the data vector pair(s) "x"/"y" to the regression model.
    void train
	(
		const Array< double >&,
		const Array< double >&
	);

    //! Removes the data vector pair(s) "x"/"y" from the regression model.
    void remove
	(
		const Array< double >&,
		const Array< double >&
	);

    //! Given one/several data vector(s) "x" the corresponding 
    //! vector(s) "y" is/are calculated by using the current regression model.
    void recall
	(
		const Array< double >&,
		Array< double >&
	);

    //! Temporarily removes the vector pair(s) "x"/"y" from the
    //! regression model and uses this modified model to approximate
    //! the vector(s) "y", when the vector(s) "x" is/are given.
    void leaveOneOut
	(
		const Array< double >&,
		const Array< double >&,
		Array< double >&
	);

    //! Returns the regression coefficient "A" of the regression model.
    Array< double > A( );

    //! Returns the regression coefficient "b" of the regression model.
    Array< double > b( );


#ifdef _WIN32
    // dummy !!!
    bool operator == ( const LinearRegression& ) const
    {
        return false;
    }

    // dummy !!!
    bool operator < ( const LinearRegression& ) const
    {
        return false;
    }
#endif


  protected:

        // Set to "true" if the internal variables sxM, syM, rxxM and rxyM
        // are changed due to adding/removing vector pair(s) x/y.
	bool modifiedM;

        // The no. of data vector pairs x/y used for training of the model.
	unsigned countM;

        // Given k data vectors x with dimension n and the corresponding
        // k vectors y with dimension m, this n-dimensional vector 
        // contains the sum of all x_l for l = 1,...,k.
	Array  < double > sxM;

        // Given k data vectors x with dimension n and the corresponding
        // k vectors y with dimension m, this n-dimensional vector 
        // contains the sum of all y_l for l = 1,...,k.
	Array  < double > syM;

        // Given k data vectors x with dimension n and the corresponding
        // k vectors y with dimension m, this n x n matrix contains the
        // sum of all s_l for l = 1,...,k with s_l = x_i * x_j
        // for i,j = 1,...,n.  
	Array2D< double > rxxM;

        // Given k data vectors x with dimension n and the corresponding
        // k vectors y with dimension m, this n x m matrix contains the
        // sum of all s_l for l = 1,...,k with s_l = x_i * y_j
        // for i = 1,...,n and j = 1,...,m.  
	Array2D< double > rxyM;

        // The transformation matrix "A" of the model y = Ax + b.
	Array2D< double > amatM;

        // The vector "b" of the model y = Ax + b.
	Array  < double > bvecM;


        friend class LocalLinearRegression;
};

#endif /* !__LINEARREGRESSION_H */


