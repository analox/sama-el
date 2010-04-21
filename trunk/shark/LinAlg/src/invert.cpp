//===========================================================================
/*!
 *  \file invert.cpp
 *
 *  \brief Determines the generalized inverse matrix of an input matrix
 *         by using singular value decomposition. Used as frontend
 *         for metod #g_inverse when using type "Array" instead of
 *         "Array2D".
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
 *      $RCSfile: invert.cpp,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:10 $
 *
 *  \par Changes:
 *      $Log: invert.cpp,v $
 *      Revision 2.0  2003/11/28 16:23:10  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:04  shark-admin
 *      INI Administration
 *
 *      Revision 1.4  2002/05/16 13:55:23  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.3  2001/11/30 13:25:45  rudi
 *      Doxygen comments added.
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
 *  \brief Returns the generalized inverse matrix of input matrix
 *         "A" by using singular value decomposition. Used as frontend
 *         for method #g_inverse when using type "Array" instead of
 *         "Array2D".
 *
 *  For a more exact description see documentation of method
 *  #g_inverse. 
 *  Here not only the usage of variable type "Array< double >"
 *  instead of "Array2D< double >" as storage for matrices
 *  is different, but also the resulting generalized inverse
 *  matrix is returned directly and not given back by assigning
 *  it to a second parameter.
 *
 *  \param  A The input matrix.
 *  \return   The generalized inverse matrix.
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em A is not a
 *         square matrix 
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
 *  \sa g_inverse.cpp, svd.cpp
 *
 */
Array< double > invert( const Array< double >& A )
{
    SIZE_CHECK( A.ndim( ) == 2 )

    Array< double > B( A.dim( 1 ), A.dim( 0 ) );
    Array2DReference< double > amatL( const_cast< Array< double >& >( A ) );
    Array2DReference< double > bmatL( B );

    g_inverse( amatL, bmatL );

    return Array< double >( B, true );
}
