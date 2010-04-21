//===========================================================================
/*!
 *  \file ArrayIndexOp.h
 *
 *  \brief Offers some operations for the index arrays defined in class
 *         ArrayIndex.
 *
 *  Here you will find operations for <br>
 *  <ul>
 *      <li>extracting indices of array elements, that fulfill
 *          a comparison to a given thresh value</li>
 *      <li>combining index arrays with AND and OR</li>
 *      <li>copying array elements indicated by an index array
 *          to another array</li>
 *      <li>testing index arrays for equality and inequality</li>
 *      <li>output of index arrays</li>
 *  </ul>
 *
 *  \author  D. Homberg
 *
 *  \par Copyright (c) 1995, 2003:
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
 *      Array
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: ArrayIndexOp.h,v $<BR>
 *      $Revision: 2.0 $<BR>
 *      $Date: 2003/11/28 16:23:21 $
 *
 *  \par Changes:
 *      $Log: ArrayIndexOp.h,v $
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.6  2002/05/16 13:18:35  rudi
 *      doxygen commands added/modified.
 *
 *
 *  This file is part of Array. This library is free software;
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

#ifndef _ARRAYINDEXOP_H

//! Set to "1" if this header file is read in. Used to avoid
//! multiple inclusion (and all the problems coming with it)
//! of this header file.
#define _ARRAYINDEXOP_H 1

#include "Array/Array.h"
#include "Array/ArrayIndex.h"


//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are greater than the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are greater than the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereGT(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr > thresh)
		{
			index.add(k);
		}
	}
	return index;
}



//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are greater than or equal to the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are greater than or equal to the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereGE(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr >= thresh)
		{
			index.add(k);
		}
	}
	return index;
}



//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are equal to the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are equal to the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereEQ(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr == thresh)
		{
			index.add(k);
		}
	}
	return index;
}



//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are not equal to the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are not equal to the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereNE(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr != thresh)
		{
			index.add(k);
		}
	}
	return index;
}


//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are less than the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are less than the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereLT(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr < thresh)
		{
			index.add(k);
		}
	}
	return index;
}


//===================================================================
/*! 
 *  \ingroup Extract
 *  \brief Adds the indices of all elements of array "vec"
 *         that are less than or equal to the given "thresh" to the
 *         index array "index" and returns this index array.
 *
 *  The indices that will be stored in \em index are related to
 *  the element vector Array::e of array \em vec.
 *
 *  \param index  stores the indices of all array elements,
 *                that are less than or equal to the \em thresh. Any
 *                existing elements are lost
 *  \param vec    the array with the elements, that will
 *                be checked
 *  \param thresh determines which elements of \em vec
 *                will be stored in \em index
 *  \return reference to the index array \em index
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T>
ArrayIndex& whereLE(ArrayIndex& index, Array<T>& vec, const T& thresh)
{
	unsigned k;
	index.setPos();
	T* vecPtr = vec.elemvec();
	
	for (k = 0; k < vec.nelem( ); k++, vecPtr++)
	{
		if ( *vecPtr <= thresh)
		{
			index.add(k);
		}
	}
	return index;
}


//===================================================================
/*! 
 *  \ingroup Math
 *  \brief Given two index arrays "a" and "b", where the indices are
 *         stored in ascending order, an index array with all indices 
 *         stored in "a" AND "b" is returned.
 *
 *  For sorting the index arrays, use function "sort" in ArraySort.h.
 *
 *  \param a the first index array, whose indices are compared to those
 *           of \em b
 *  \param b the second index array, whose indices are compared to those
 *           of \em a
 *  \return index array with all indices, that are stored in
 *          \em a AND \em b
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em a or \em b
 *         are not one-dimensional
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
inline ArrayIndex operator&&(ArrayIndex& a, ArrayIndex& b)
{
	unsigned i = a.getPos() > b.getPos() ? b.getPos() : a.getPos(), ia, ib;
	if (i == 0) i = 1;
	ArrayIndex v(i);

	for (ia = 0, ib = 0; ia < a.getPos() && ib < b.getPos(); )
	{
		if (a(ia) < b(ib))
		{
			ia++;
			continue;
		}
		if (a(ia) > b(ib))
		{
			ib++;
			continue;
		}

                // a(ia) == b(ib):
		v.add(a(ia));
		ia++; ib++;
	}
	return v;
}


//===================================================================
/*! 
 *  \ingroup Math
 *  \brief Given two index arrays "a" and "b", where the indices are
 *         stored in ascending order, an index array with all indices 
 *         stored in "a" OR "b" is returned.
 *
 *  If an index value is stored in \em a and \em b, then it is only 
 *  stored once in the returned index array. <br>
 *  For sorting the index arrays, use function "sort" in ArraySort.h.
 *
 *  \param a the first index array, whose indices are compared to those
 *           of \em b
 *  \param b the second index array, whose indices are compared to those
 *           of \em a
 *  \return index array with all indices, that are stored in
 *          \em a OR \em b
 *  \throw check_exception the type of the exception will be
 *         "size mismatch" and indicates that \em a or \em b
 *         are not one-dimensional
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
inline ArrayIndex operator||(ArrayIndex& a, ArrayIndex& b)
{
	unsigned ia, ib;
	ArrayIndex v(a.getPos() + b.getPos());

	for (ia = 0, ib = 0; ia < a.getPos() || ib < b.getPos(); )
	{
		if (ia < a.getPos() && ib < b.getPos())
		{
			if (a(ia) == b(ib))
			{
				v.add(a(ia));
				ia++; ib++;
			} else {
				if (a(ia) < b(ib))
				{
					v.add(a(ia));
					ia++;
				} else {
					v.add(b(ib));
					ib++;
				}
			}
		} else {
			if (ia < a.getPos()) {
				v.add(a(ia));
				ia++;
			} else {
				v.add(b(ib));
				ib++;
			}
		}
	}
	return v;
}


//===================================================================
/*! 
 *  \ingroup Copy
 *  \brief Copies all elements of "src", whose indices are stored in
 *         index array "index" to the array "dst".
 *
 *
 *  \param dst    contains copies of all elements of \em src,
 *                whose indices are stored in \em index.
 *                The copied elements from \em dst are casted
 *                to the element type of \em dst 
 *  \param src    array from which elements are copied to \em dst
 *  \param index  indicates which elements of \em src are stored
 *                in dst
 *  \return none
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
template <class T1, class T2>
void copy(Array<T1>& dst, const Array<T2>& src, const ArrayIndex& index)
{
	dst.resize(index.getPos());
	Array<T1>::iterator iterDst;
	ArrayIndex::iterator iter;
	for (iter = index.begin(), iterDst = dst.begin();
	     iter != index.end(); iter++, iterDst++)
	{
		*iterDst = (T1) src.elem(*iter);
	}
}


//===================================================================
/*! 
 *  \ingroup InOut
 *  \brief Writes the content of index array "a" to the output
 *         stream "os". 
 *
 *  First the type of the array \em a ("ArrayIndex<unsigned>") and
 *  the number of the indices that are stored in it are written
 *  to \em os. <br>
 *  Then all indices stored in \em a are written to the output stream,
 *  separated by a tabulator character.
 *
 *  \par Example
 *  \code
 *  #include "Array/ArrayIndex.h"
 *  #include "Array/ArrayIndexOp.h"
 *  
 *  void main()
 *  {
 *    ArrayIndex test( 5 );
 *  
 *    for ( unsigned i = 0; i < 5; i++ ) test.add( i );
 *  
 *    cout << test;
 *  }
 *  \endcode
 *
 *  The given source code will produce the output: <br>
 *
 *  \f$
 *  \mbox{\ }\\ \noindent
 *  \mbox{ArrayIndex}<\mbox{unsigned}>\mbox{(5)}\\
 *  \mbox{0\ \ \ \ 1\ \ \ \ 2\ \ \ \ 3\ \ \ \ 4}\\
 *  \f$
 *
 *  \param a  the index array, whose elements are written to \em os
 *  \param os the output stream to which the content of \em a is written
 *  \return reference to the output stream \em os
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
inline std::ostream& operator << ( std::ostream& os, const ArrayIndex& a )
{
    os << "ArrayIndex<unsigned>(" << a.getPos() << ")\n";

		for (ArrayIndex::iterator iter = a.begin();
		     iter != a.end(); iter++)
		{
			if( iter != a.begin() ) os << '\t';
			os << *iter;
		}

    return os << std::endl;
}

//===================================================================
/*! 
 *  \ingroup Compare
 *  \brief Returns "true" if the index arrays "v" and "w" are equal.
 *
 *  The index array \em v and \em w are equal, if they contain
 *  the same number of elements and each index stored in \em v
 *  is equal to the corresponding index stored in \em w.
 *
 *  \param v the first index array that is compared to \em w
 *  \param w the second index array that is compared to \em v
 *  \return "true", if both index arrays are equal, "false" otherwise
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
inline bool operator == ( const ArrayIndex& v, const ArrayIndex& w )
{
		if ( v.getPos() != w.getPos() ) return false;
		ArrayIndex::iterator iterV, iterW;
		for(iterV = v.begin(), iterW = w.begin(); iterV != v.end();
			  iterV++, iterW++)
		{
			if (*iterV != *iterW) return false;
		}
		return true;
}


//===================================================================
/*! 
 *  \ingroup Compare
 *  \brief Returns "true" if the index arrays "v" and "w" are not equal.
 *
 *  The index array \em v and \em w are not equal, if they contain
 *  different numbers of elements or at least one index stored in \em v
 *  is not equal to the corresponding index stored in \em w.
 *
 *  \param v the first index array that is compared to \em w
 *  \param w the second index array that is compared to \em v
 *  \return "true", if both index arrays are not equal, "false" otherwise
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */
inline bool operator != ( const ArrayIndex& v, const ArrayIndex& w )
{
	return !(v == w);
}

#endif /* _ARRAYINDEXOP_H */











