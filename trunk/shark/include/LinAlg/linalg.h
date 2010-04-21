//===========================================================================
/*!
 *  \file linalg.h
 *
 *  \brief Some operations for matrices.
 *
 *
 *  \author  M. Kreutz
 *  \date    1995-01-01
 *
 *  \par Copyright (c) 1999-2001:
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
 *      $RCSfile: linalg.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/16 15:22:20 $
 *
 *  \par Changes:
 *      $Log: linalg.h,v $
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
 *      Revision 1.8  2002/08/13 14:59:37  thorsten
 *      variance for linear arrays
 *
 *      Revision 1.7  2002/08/01 10:18:22  arne
 *      g_inverse() now returns rank
 *
 *      Revision 1.6  2002/05/16 13:21:21  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.5  2002/02/06 14:27:37  rudi
 *      Error in method "transpose" removed, some doxygen comments added.
 *
 *      Revision 1.4  2001/11/30 14:10:14  rudi
 *      Merge of old "linalg.h" and "arraylinalg.h".
 *      doxygen comments added.
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

#ifndef LINALG_H
#define LINALG_H

#include "Array/Array2D.h"


// The commands in the following comments are used for the
// doxygen documentation program. doxygen is forced to include the
// source code of the named example programs into the documentation,
// that then can be linked directly.
// This "dirty trick" is necessary, because to put "\example [filename]"
// directly into method documentations has some unaesthetic side effects.
//
/*! \example rank_test.cpp */
/*! \example rank_decomp_test.cpp */
/*! \example svd_test.cpp */
/*! \example svdrank_test.cpp */
/*! \example svdsort_test.cpp */
/*! \example eigensymm_test.cpp */
/*! \example eigensymmJacobi_test.cpp */
/*! \example eigensymmJacobi2_test.cpp */
/*! \example eigensort_test.cpp */
/*! \example eigenerr_test.cpp */
/*! \example bfgs_test.cpp */
/*! \example bfgs2_test.cpp */
/*! \example g_inverse_matrix.cpp */
/*! \example detsymm_test.cpp */
/*! \example linalg_simple_test.cpp */
/*! \example covar_corrcoef_test.cpp */
/*! \example linearRegression_test.cpp */
/*! \example linearClassifier_test.cpp */


//! Sorts the eigenvalues in vector "dvecA" and the corresponding
//! eigenvectors in matrix "vmatA".
void eigensort
(
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);

//! Calculates the eigenvalues and the normalized eigenvectors of the
//! symmetric matrix "amatA" using the Jacobi method.
void eigensymmJacobi
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);

//! Calculates the eigenvalues and the normalized
//! eigenvectors of the symmetric matrix "amatA" using a modified
//! Jacobi method.
void eigensymmJacobi2
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);



//! Calculates the eigenvalues and the normalized
//! eigenvectors of a symmetric matrix "amatA" using the Givens
//! and Householder reduction, however, "amatA" is destroyed after application.
void eigensymm_obsolete
(
	const Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);

//! Used as frontend for
//! #eigensymm_obsolete(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA)
//! for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
//! and Householder reduction, however, "amatA" is destroyed after application.
void eigensymm_obsolete
(
    const Array< double >& A,
    Array< double >& G,
    Array< double >& l
);


//! Used as frontend for 
//! #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA) 
//! for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
//! and Householder reduction without corrupting "amatA" during application. Each time this frontend is 
//! called additional memory is allocated for intermediate results.

void eigensymm
(
	const Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);


//! Used as another frontend for 
//! #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA) for calculating 
//! the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
//! and Householder reduction without corrupting "A" during application. Each time this frontend is 
//! called additional memory is allocated for intermediate results.

void eigensymm
(
    const Array< double >& A,
    Array< double >& G,
    Array< double >& l
);


//! Calculates the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
//! and Householder reduction. This method works without corrupting "amatA" during application by demanding
//! another Array 'odvecA' as an algorithmic buffer instead of using the last row of 'amatA' to store 
//! intermediate algorithmic results.

void eigensymm
(
	const Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA,
	Array  < double >& odvecA
);


//! Used as frontend for 
//! #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA) 
//! for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the 
//! Givens and Householder reduction without corrupting "A" during application.

void eigensymm
(
    const Array< double >& A,
    Array< double >& G,
    Array< double >& l,
    Array< double >& od
);

//===========================================================================
/*!
 *  \brief Another frontend for 
 *  #eigensymm_obsolete(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA) for calculating
 *  the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
 *  and Householder reduction. 
 *
 *  Frontend for function 
 *  #eigensymm_obsolete(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA), 
 *  when using type \em ArrayReference instead of type \em Array2D.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom
 *                triangular matrix must contain values.
 *                Values below the diagonal will be destroyed
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column contains one
 *                eigenvector
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order
 *      \return none
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *     previously 'eigensymm', renamed by S. Wiegand 2003/10/01 
 *
 *  \par Status
 *      stable
 *
 */
inline void eigensymm_obsolete
(
    const ArrayReference< double > A,
	  ArrayReference< double > G,
          ArrayReference< double > l
)
{
    eigensymm_obsolete( static_cast< const Array< double >& >( A ),
                        static_cast< Array< double >& >( G ),
                        static_cast< Array< double >& >( l ) );
}

//===========================================================================
/*!
 *  \brief Another frontend for
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA) 
 *  for calculating  the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the Givens
 *  and Householder reduction without corrupting the input matrix "amatA" during application. Each time this 
 *  frontend is called additional memory is allocated for intermediate results.
 *
 *  Frontend for function 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA),
 *  when using type \em ArrayReference instead of type \em Array2D. Each time this frontend is 
 *  called additional memory is allocated for intermediate results.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom triangular matrix must contain values.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column contains one
 *                eigenvector
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order
 *      \return none
 *
 *  \author  M. Kreutz
 *  \date    1998
 *
 *  \par Changes
 *       none
 *
 *  \par Status
 *      stable
 *
 */
inline void eigensymm
(
    const ArrayReference< double > A,
	  ArrayReference< double > G,
          ArrayReference< double > l
)
{
    eigensymm( static_cast< const Array< double >& >( A ),
               static_cast< Array< double >& >( G ),
               static_cast< Array< double >& >( l ) );
}

//===========================================================================
/*!
 *  \brief Another frontend for
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA) 
 *  for calculating the eigenvalues and the normalized eigenvectors of a symmetric matrix "amatA" using the
 *  Givens and Householder reduction without corrupting the input matrix "amatA" during application.
 *
 *  Frontend for function 
 *  #eigensymm(const Array2D<double> &amatA,Array2D<double> &vmatA,Array<double> &dvecA,Array<double> &odvecA),
 *  when using type \em ArrayReference instead of type \em Array2D.
 *
 *      \param  A \f$ n \times n \f$ matrix, which must be symmetric, so
 *                only the bottom triangular matrix must contain values.
 *      \param	G \f$ n \times n \f$ matrix with the calculated normalized 
 *                eigenvectors, each column contains one
 *                eigenvector
 *	\param  l n-dimensional vector with the calculated 
 *                eigenvalues in descending order
 *	\param  od n-dimensional vector with the calculated 
 *              offdiagonal of the Householder transformation.
 *      \return none
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes
 *        none
 *
 *  \par Status
 *      stable
 *
 */
inline void eigensymm
(
    const ArrayReference< double > A,
	  ArrayReference< double > G,
          ArrayReference< double > l,
          ArrayReference< double > od
)
{
    eigensymm( static_cast< const Array< double >& >( A ),
               static_cast< Array< double >& >( G ),
               static_cast< Array< double >& >( l ),
               static_cast< Array< double >& >( od )  );
}

//! Calculates the relative error of one eigenvalue with no. "c".
double eigenerr
(
	const Array2D< double >& amatA,
	const Array2D< double >& vmatA,
	const Array  < double >& dvecA,
	unsigned c
);

//! Determines the numerical rank of the symmetric matrix "amatA".
unsigned rank
(
	const Array2D< double >& amatA,
	const Array2D< double >& vmatA,
	const Array  < double >& dvecA
);

//! Calculates the determinate of the symmetric matrix "amatA".
double detsymm
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);

//! Calculates the rank of the symmetric matrix "amatA", its eigenvalues and 
//! eigenvectors.
unsigned rankDecomp
(
	Array2D< double >& amatA,
	Array2D< double >& vmatA,
	Array  < double >& dvecA
);

//! Given "m" classes of data, the covariances between
//! the values within one class and the covariances
//! between all the classes, this method calculates
//! the transformation matrix that will project the
//! data in a way, that maximum separation of the 
//! different classes is given.
void discrimAnalysis
(
	Array2D< double >& betweenCovarA,
	Array2D< double >& withinCovarA,
	Array2D< double >& transMatA,
	Array  < double >& dvecA,
	unsigned& m
);


//! Given the correlations of the n-dimensional data vector "x" 
//! and the m-dimensional data vector "y" and also given
//! their mean values, this function summarizes the data by 
//! finding a linear mapping that will approximate the data. 
void linearRegress
(
	Array2D< double >& cxxMatA,
	Array2D< double >& cxyMatA,
	Array  < double >& mxVecA,
	Array  < double >& myVecA,
        Array2D< double >& amatA,
	Array  < double >& bvecA,
	Array  < double >& dvecA
);


//! Determines the numerical rank of a rectangular matrix "amatA",
//! when a singular value decomposition for "amatA" has taken place 
//! before.
unsigned svdrank
(
	const Array2D< double >& amatA,
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& wvecA
);


//! Calculates the generalized inverse matrix of input matrix "amatA".
unsigned g_inverse
(
	const Array2D< double >& amatA,
	Array2D< double >& bmatA
);


//! Determines the singular value decomposition of a rectangular
//! matrix "amatA".
void svd
(
	const Array2D< double >& amatA,
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& wvecA
);


//! Used as frontend for function #svd(const Array2D< double >& amatA, Array2D< double >& umatA, Array2D< double >& vmatA, Array  < double >& w),
//! calculates singular value decomposition for a square matrix "A".
void svd
(
    const Array< double >& A,
	  Array< double >& U,
	  Array< double >& V,
	  Array< double >& W
);


//! Sorts the singular values in vector "wvecA" by descending order.
void svdsort
(
	Array2D< double >& umatA,
	Array2D< double >& vmatA,
	Array  < double >& wvecA
);



//===========================================================================
/*!
 *  \brief Transpose of matrix "v".
 *
 *  Given the matrix \em v, the transposed matrix is created, i.e.
 *  the rows become columns and vice versa:
 *  \f[
 *      \left(
 *      \begin{array}{*{4}{c}}
 *          1  & 2  & 3  & 4\\
 *          5  & 6  & 7  & 8\\
 *          9  & 10 & 11 & 12\\
 *          13 & 14 & 15 & 16\\    
 *      \end{array}
 *      \right)
 *      \longrightarrow
 *      \left(
 *      \begin{array}{*{4}{c}}
 *          1 & 5 & 9  & 13\\
 *          2 & 6 & 10 & 14\\
 *          3 & 7 & 11 & 15\\
 *          4 & 8 & 12 & 16\\    
 *      \end{array}
 *      \right) 
 *  \f]
 *  The original matrix \em v will not be modified.
 *
 *      \param  v matrix, that will be transposed
 *      \return the transposed matrix
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em v is not
 *             2-dimensional
 *
 *  \example linalg_simple_test.cpp
 *
 *  Please follow the link to view the source code of the example.
 *  The example can be executed in the example directory
 *  of package LinAlg. 
 *
 *  \author  M. Kreutz
 *  \date    1995
 *
 *  \par Changes
 *      2002-01-23, ra: <br>
 *      Formerly worked only for square matrices - fixed.
 *
 *  \par Status
 *      stable
 *
 */
template < class T >
Array< T > transpose( const Array< T >& v )
{
    SIZE_CHECK( v.ndim( ) == 2 )
    Array< T > z;
    z.resize( v.dim( 1 ), v.dim( 0 ) );
    for( unsigned i = v.dim( 0 ); i--; )
        for( unsigned j = v.dim( 1 ); j--; )
	    z( j, i ) = v( i, j );
    return Array< T >( z, true );
}


//===========================================================================
/*!
 *  \brief Creates a new matrix with vector "v" lying on the
 *         diagonal.
 *
 *  Given the vector \em v, a new matrix is created. The diagonal
 *  of the new matrix adapt the values of \em v and all other
 *  other matrix positions are set to zero.<br>
 *  Example:
 *  
 *  \f[
 *      \left(
 *      \begin{array}{*{4}{c}}
 *          1 \\
 *          2 \\
 *          3 \\
 *          4 \\    
 *      \end{array}
 *      \right)
 *      \longrightarrow
 *      \left(
 *      \begin{array}{*{4}{c}}
 *          1 & 0 & 0 & 0\\
 *          0 & 2 & 0 & 0\\
 *          0 & 0 & 3 & 0\\
 *          0 & 0 & 0 & 4\\    
 *      \end{array}
 *      \right) 
 *  \f]
 *
 *      \param  v vector with values for the matrix diagonal
 *      \return the new matrix with \em v as diagonal
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em v is not
 *             one-dimensional
 *
 *  \example linalg_simple_test.cpp
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
template < class T >
Array< T > diagonal( const Array< T >& v )
{
    SIZE_CHECK( v.ndim( ) == 1 )
    Array< T > z( v.nelem( ), v.nelem( ) );
    z = 0;
    for( unsigned i = v.nelem( ); i--; )
        z( i, i ) = v.elem( i );
    return Array< T >( z, true );
}


//===========================================================================
/*!
 *  \brief Evaluates the sum of the values at the diagonal of
 *         matrix "v".
 *
 *  Example:
 *  \f[
 *      \left(
 *      \begin{array}{*{4}{c}}
 *          {\bf 1} & 5       & 9        & 13\\
 *          2       & {\bf 6} & 10       & 14\\
 *          3       & 7       & {\bf 11} & 15\\
 *          4       & 8       & 12       & {\bf 16}\\    
 *      \end{array}
 *      \right) 
 *      \longrightarrow 1 + 6 + 11 + 16 = 34
 *  \f]
 *
 *      \param  v square matrix
 *      \return the sum of the values at the diagonal of \em v
 *      \throw check_exception the type of the exception will be
 *             "size mismatch" and indicates that \em v is
 *             not a square matrix
 *
 *  \example  linalg_simple_test.cpp
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
template < class T >
T trace( const Array< T >& v )
{
    SIZE_CHECK( v.ndim( ) == 2 && v.dim( 0 ) == v.dim( 1 ) )

    T t( v( 0, 0 ) );
    for( unsigned i = 1; i < v.dim( 0 ); ++i )
        t += v( i, i );
    return t;
}


//! Calculates the mean vector of array "x".
Array< double > mean( const Array< double >& x );


//! Calculates the variance vector of array "x".
Array< double > variance( const Array< double >& x );


//! Calculates the angle between the vectors "x" and "y".
double angle( const Array< double >& x, const Array< double >& y );


//! Calculates the coefficient of correlation of the data 
//! vectors "x" and "y".
double corrcoef( const Array< double >& x, const Array< double >& y );


//! Calculates the coefficient of correlation matrix of the data 
//! vectors stored in matrix "x".
Array< double > corrcoef( const Array< double >& x );


//! Calculates the mean and variance values of matrix "x".
void meanvar
(
    const Array< double >& x,
    Array< double >&,
    Array< double >&
);

//! Calculates the mean and variance values of 1d-arrays p(x)
void meanvar
(
    const Array< double >& pxA,
    const Array< double >& xA,
	double &mA,
	double &vA,
	const int startA = -1,
	const int endA = -1
);


//! Calculates the covariance between the data vectors "x" and "y".
double covariance( const Array< double >& x, const Array< double >& y );


//! Calculates the covariance matrix of the data vectors stored in
//! matrix "x".
Array< double > covariance( const Array< double >& x );


//! Returns the generalized inverse matrix of input matrix
//! "A" by using singular value decomposition. Used as frontend
//! for metod #g_inverse when using type "Array" instead of
//! "Array2D".
Array< double > invert( const Array< double >& );

#endif  // LINALG_H
