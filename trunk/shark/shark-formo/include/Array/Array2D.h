//===========================================================================
/*!
 *  \file Array2D.h
 *
 *  \brief Provides a more efficient handling of 2-dimensional arrays.
 *
 *  Using an additional pointer matrix that allows direct access
 *  to all array elements, especially the single array rows, the
 *  computation time for operations is decreased.
 *
 *  \author  M. Kreutz
 *  \date    2001-04-24
 *
 *  \par Copyright (c) 2001:
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
 *      $RCSfile: Array2D.h,v $<BR>
 *      $Revision: 2.3 $<BR>
 *      $Date: 2005/07/06 10:07:47 $
 *
 *  \par Changes:
 *      $Log: Array2D.h,v $
 *      Revision 2.3  2005/07/06 10:07:47  glasmtbl
 *      inserted some "this->" vor ANSI compatibility.
 *
 *      Revision 2.2  2005/05/23 10:02:18  christian_igel
 *      In a template definition, unqualified names will no longer find
 *      members of a dependent base (as specified by [temp.dep]/3 in the C++
 *      standard).
 *
 *      Hence we made the names dependent by prefixing them with this->.
 *
 *      Revision 2.1  2004/12/29 16:07:07  glasmtbl
 *      *** empty log message ***
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.3  2002/05/16 13:18:34  rudi
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

#ifndef ARRAY2D_H
#define ARRAY2D_H

// include this STL-header for template functions "min" and "max"
#include "Array/Array.h"

// forward declaration
template< class T >
class Array2DReference;


//===========================================================================
/*!
 *  \brief A class for handling the special case of 2-dimensional arrays.
 *         The class uses a pointer matrix to increase the efficiency.
 *
 *  Besides the normal internal structure of arrays as given in the
 *  template class #Array, a pointer matrix for direct access of
 *  each array element is administrated. Especially when accessing
 *  single array rows, this additional structure speeds up 
 *  operations.
 *
 *  \author  M. Kreutz
 *  \date    2001-04-24
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
template< class T >
class Array2D : public Array< T >
{
  public:

    //! Pointer to a single array element. Used for compatibility with 
    //! the other stdlib structures.
    typedef T* iterator;

    //! Pointer to a single array element. Used for compatibility with 
    //! the other stdlib structures.
    typedef T* const_iterator;




   //===================================================================
   /*! 
    *  \ingroup Create
    *  \brief Creates a new empty Array2D object.
    *
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D( ) : Array< T >( )
    {
	allocPtrArr( );
    }


   //===================================================================
   /*! 
    *  \ingroup Create
    *  \brief Creates a new 2-dimensional array with "i" rows and "j" 
    *         columns.
    *
    *  A 2-dimensional array is created and memory for the pointer
    *  matrix #ptrArrM is allocated.
    *
    *  \param i number of rows of the new array
    *  \param j number of columns of the new array
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D( unsigned i, unsigned j ) : Array< T >( i, j )
    {
	allocPtrArr( );
    }


   //=======================================================================
   /*! 
    *  \ingroup Create
    *  \brief Creates a new 2-dimensional array with the structure and 
    *         content of the constant array "v".
    *
    *  Memory for the pointer matrix #ptrArrM will be allocated.
    * 
    *  \param v array, whose content will be copied to the new array
    *  \return none
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" and indiactes that \em v is no 2-dimensional
    *         or non-dimensional array
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D( const Array< T >& v )
	  : Array< T >( v )
    {
		SIZE_CHECK( v.ndim( ) == 0 || v.ndim( ) == 2 )
		allocPtrArr( );
    }

   //=======================================================================
   /*! 
    *  \ingroup Create
    *  \brief Creates a new 2-dimensional array with the structure and 
    *         content of the constant 2-dimensional array "v".
    *
    *  Memory for the pointer matrix #ptrArrM will be allocated.
    * 
    *  \param v 2-dimensional array, whose content will be copied to the 
    *           new array
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D( const Array2D< T >& v )
	  : Array< T >( v )
    {
		allocPtrArr( );
    }


   //===================================================================
   /*! 
    *  \brief Destructs the current Array2D object and frees memory.
    *
    *  The pointer matrix #ptrArrM is destructed, too.
    *
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    ~Array2D( )
    {
		deallocPtrArr( );
    }

   //=======================================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns the element at row "i" and column "j" of the
    *         current array.
    *
    *  \param i row-position of the element that will be returned.
    *  \param j column-position of the element that will be returned.
    *  \return The i-th, j-th array element
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" if the current array is not 2-dimensional
    *         and "range check error" if \em i exceeds the size of the
    *         array's first dimension or \em j exceeds the size of the
    *         array's second dimension
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    T& operator ( ) ( unsigned i, unsigned j )
    {
	SIZE_CHECK( this->nd == 2 )
        RANGE_CHECK( i < this->d[ 0 ] )
        RANGE_CHECK( j < this->d[ 1 ] )
        return ptrArrM[ i ][ j ];
    }


   //=======================================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns the element at row "i" and column "j" of the
    *         current array as constant.
    *
    *  \param i row-position of the element that will be returned.
    *  \param j column-position of the element that will be returned.
    *  \return The i-th, j-th array element as constant
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" if the current array is not 2-dimensional
    *         and "range check error" if \em i exceeds the size of the
    *         array's first dimension or \em j exceeds the size of the
    *         array's second dimension
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    const T& operator ( ) ( unsigned i, unsigned j ) const
    {
	SIZE_CHECK( this->nd == 2 )
        RANGE_CHECK( i < this->d[ 0 ] )
        RANGE_CHECK( j < this->d[ 1 ] )
        return ptrArrM[ i ][ j ];
    }


   //=======================================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns row no. "i" of the current array.
    *
    *  \param i number of the row that will be returned
    *  \return pointer to the first element of the i-th row
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" if the current array is not 2-dimensional
    *         and "range check error" if \em i exceeds the array's
    *         first dimension
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    T* operator [ ] ( unsigned i )
	{
		SIZE_CHECK( this->nd == 2 )
                RANGE_CHECK( i < this->d[ 0 ] )
		return ptrArrM[ i ];
	}


   //=======================================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns row no. "i" of the current array as constant.
    *
    *  \param i number of the row that will be returned
    *  \return constant pointer to the first element of the i-th row
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" if the current array is not 2-dimensional
    *         and "range check error" if \em i exceeds the array's
    *         first dimension
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    const T* operator [ ] ( unsigned i ) const
	{
		SIZE_CHECK( this->nd == 2 )
                RANGE_CHECK( i < this->d[ 0 ] )
		return ptrArrM[ i ];
	}

   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the value "v" to all positions of the current
    *         2-dimensional array.
    *
    *  \param v the new value for all array elements
    *  \return the new 2-dimensional array with the new values
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D< T >& operator = ( const T& v )
    {
	Array< T >::operator = ( v );
        return *this;
    }


   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the values of array "v" to the current 2-dimensional
    *         array.
    *
    *  The size of the current 2-dimensional array is adopted to the size of 
    *  array \em v.
    *
    *  \param v array with the new values for the current 2-dimensional array
    *  \return the current 2-dimensional array with the new values
    *  \throw check_exception the type of the exception is "size mismatch"
    *         and indicates that \em v is no 2-dimensional array or
    *         that the current array is a static array reference whose
    *         size is different to that of \em v
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D< T >& operator = ( const Array< T >& v )
    {
		SIZE_CHECK( v.ndim( ) == 0 || v.ndim( ) == 2 )
		Array< T >::operator = ( v );
		reallocPtrArr( );
		return *this;
    }

   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the values of the 2-dimensional array "v" to the 
    *         current 2-dimensional array.
    *
    *  The size of the current 2-dimensional array is adopted to the size of 
    *  the 2-dimensional array \em v.
    *
    *  \param v 2-dimensional array with the new values for the current 
    *           2-dimensional array
    *  \return the current 2-dimensional array with the new values
    *  \throw check_exception the type of the exception is "size mismatch"
    *         and indicates that the current array is a static array 
    *         reference whose size is different to that of \em v
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2D< T >& operator = ( const Array2D< T >& v )
    {
		Array< T >::operator = ( v );
		reallocPtrArr( );
		return *this;
    }


   //===================================================================
   /*!
    *  \ingroup Information 
    *  \brief Returns the number of array rows.
    *
    *  If the array is 2-dimensional, then the number of rows
    *  is returned, otherwise zero is returned.
    *
    *  \return number of rows or zero
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
	unsigned rows( ) const
	{
		return (this->nd == 2) ? (this->d)[ 0 ] : 0;
	}


   //===================================================================
   /*! 
    *  \ingroup Information
    *  \brief Returns the number of array columns.
    *
    *  If the array is 2-dimensional, then the number of columns
    *  is returned, otherwise zero is returned.
    *
    *  \return number of columns or zero
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
	unsigned cols( ) const
	{
		return (this->nd) == 2 ? (this->d)[ 1 ] : 0;
	}



   //===================================================================
   /*! 
    *  \brief Returns a pointer to the internal pointer matrix.
    *
    *  Handle this method with care, because direct manipulation
    *  of the pointer matrix #ptrArrM can cause errors.
    *
    *  \return pointer to the pointer matrix
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    T** ptrArr( )
    {
        return ptrArrM;
    }


   //===================================================================
   /*! 
    *  \brief Returns a pointer to the internal pointer matrix as constant.
    *
    *  Handle this method with care, because direct manipulation
    *  of the pointer matrix #ptrArrM can cause errors.
    *
    *  \return constant pointer to the pointer matrix
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    T** const ptrArr( ) const
    {
        return ptrArrM;
    }

   //===================================================================
   /*! 
    *  \ingroup Copy
    *  \brief Returns an identical copy of this Array2D object.
    * 
    *  Creates a new 2-dimensional array that is identical to the current one
    *  and returns a pointer to this clone.
    *
    *  \return a copy of the current array
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    ArrayBase* clone( ) const { return new Array2D< T >( *this ); }


   //===================================================================
   /*! 
    *  \ingroup Create
    *  \brief Returns an empty 2-dimensional array with the same type 
    *         as this Array2D object.
    * 
    *  A new empty Array2D object with the same type as the current array is
    *  created and a pointer to this new object returned.
    *
    *  \return an empty copy of the current array
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    ArrayBase* empty( ) const { return new Array2D< T >( );       }


   //=======================================================================
   /*! 
    *  \brief Returns an iterator to the first array element.
    *
    *  A pointer to the first element of the element vector Array::e
    *  of the current array is returned as iterator. <br>
    *  This method is used for compatibility with the other
    *  stdlib structures.
    *
    *  \return Iterator to the first array element.
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    inline iterator begin( )
    {
        return this->e;
    }
         
   //=======================================================================
   /*! 
    *  \brief Returns a constant iterator to the first array element.
    *
    *  A constant pointer to the first element of the element vector
    *  Array::e of the current array is returned as iterator. <br>
    *  This method is used for compatibility with the other
    *  stdlib structures.
    *
    *  \return Constant iterator to the first array element.
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    inline const const_iterator begin( ) const
    {
        return this->e;
    }

   //=======================================================================
   /*! 
    *  \brief Returns an iterator to the last array element.
    *
    *  A constant pointer to the last element of the element vector
    *  Array::e of the current array is returned as iterator. <br>
    *  This method is used for compatibility with the other
    *  stdlib structures.
    *
    *  \return Iterator to the last array element.
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    inline iterator end( )
    {
        return this->eEnd;
    }
     
   //=======================================================================
   /*! 
    *  \brief Returns a constant iterator to the last array element.
    *
    *  A constant pointer to the last element of the element vector
    *  Array::e of the current array is returned as iterator. <br>
    *  This method is used for compatibility with the other
    *  stdlib structures.
    *
    *  \return Constant iterator to the last array element.
    * 
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */          
    inline const const_iterator end()   const
    {
        return this->eEnd;
    }

  protected:

        //! Additional pointer matrix to the single elements of 
        //! the current 2-dimensional array. <br>
        //! This pointer matrix is used for more efficiency, because
        //! you can use it to directly access e.g. single rows of the
        //! array.
        T** ptrArrM;

        //! Number of array rows.
	unsigned rowsM;

        //===================================================================
        /*! 
         *  \brief Creates a new pointer matrix and allocates memory
	 *         for it.
         *
         *  Based on the structure of the current template Array object, 
         *  memory for the pointer matrix is allocated. <br>
         *  If there is only an empty array, then also no memory for
         *  the pointer matrix is allocated. <br>
         *  If there is a non-empty array, then the corresponding
         *  memory for the pointer matrix is allocated and the 
         *  single pointers are set to the elements in the
         *  array element vector Array::e. <br>
         *  Be cautious! It is not checked whether the array
         *  is really 2-dimensional.
         *
         *  \return none
         *
         *  \author  M. Kreutz
         *  \date    2001-04-24
         *
         *  \par Changes
         *      none
         *
         *  \par Status
         *      stable
         *
         */
	void allocPtrArr( )
	{
		if( this->nd == 0 )
		{
			ptrArrM = new T*[ rowsM = 0 ];
		}
		else
		{
			ptrArrM = new T*[ rowsM = (this->d)[ 0 ] ];
			for( unsigned i = 0, j = 0; i < (this->d)[ 0 ]; ++i, j += (this->d)[ 1 ] )
			{
				ptrArrM[ i ] = this->e + j;
			}
		}
	}

        //===================================================================
        /*! 
         *  \brief Frees the memory for an existing pointer array.
         *
         *
         *  \return none
         *
         *  \author  M. Kreutz
         *  \date    2001-04-24
         *
         *  \par Changes
         *      none
         *
         *  \par Status
         *      stable
         *
         */
	void deallocPtrArr( )
	{
		delete[ ] ptrArrM;
	}

        //===================================================================
        /*! 
         *  \brief First frees the memory an existing pointer array to 
         *         allocate it new then.
         *
         *  \return none
         *
         *  \author  M. Kreutz
         *  \date    2001-04-24
         *
         *  \par Changes
         *      none
         *
         *  \par Status
         *      stable
         *
         */
	void reallocPtrArr( )
	{
		deallocPtrArr( );
		allocPtrArr( );
	}


   //===================================================================
   /*! 
    *  \brief Handles memory allocation for the other resize methods,
    *         including allocation of the pointer matrix. 
    *
    *  This function handles memory allocation of the
    *  respective template type and the pointer matrix in case of 
    *  resizing and is only for internal use.
    *
    *  \param dimA  one-dimensional vector of array-dimensions
    *  \param ndimA number of array-dimensions, must be zero or 2
    *  \param copyA flag which indicates whether existing elements of
    *               an array should be copied in case of resizing (as long
    *               as possible)
    *  \return none
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" and indicates that \em ndimA
    *         is not 2-dimensional and not non-dimensional 
    *         or that the current array
    *         is a static array reference and the sizes given
    *         by \em dimA and \em ndimA are not the same than
    *         the current array has
    *          
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    *  \sa Array::resize_i(unsigned*, unsigned, bool)
    *
    */
    void resize_i
	(
		unsigned* dimA,
		unsigned  ndimA,
		bool      copyA
	)
    {
		SIZE_CHECK( ndimA == 0 || ndimA == 2 )
		Array< T >::resize_i( dimA, ndimA, copyA );
		reallocPtrArr( );
    }
};



//===========================================================================
/*!
 *  \brief Class that implements references to objects of the
 *         Array2D template class.
 *         
 *  A reference is an Array2D object, that is identified
 *  by the flag "stat", which is set to "true".
 *
 *  \author  M. Kreutz
 *  \date    2001-04-24
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
template< class T >
class Array2DReference : public Array2D< T >
{
  public:

   //========================================================================
   /*!
    *  \ingroup Create
    *  \brief Creates a new empty Array2D Reference object.
    * 
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2DReference( )
	  : Array2D< T >( )
    {
        this->stat = true;
    }


   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the values of array "v" to the current array 2D 
    *         reference.
    *
    *  The sizes of the current array reference, including the pointer 
    *  matrix #ptrArrM are adopted to the size of array \em v 
    *
    *  \param v array with the new values for the current array reference
    *  \return the current array reference with the new values
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2DReference( Array< T >& v )
	  : Array2D< T >( )
    {
        this->d    = v.dimvec( );
        this->e    = v.elemvec( );
        this->nd   = v.ndim( );
        this->ne   = v.nelem( );
        this->stat = true;
        this->eEnd = this->e + this->ne;
	this->reallocPtrArr( );
    }


   //========================================================================
   /*!
    *  \brief Destroys the current array 2D reference.
    *
    *  It is important to define the destructor, because otherwise
    *  Visual C++ wouldn't call the destructors of the base classes.
    *
    *  \return none
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    ~Array2DReference( )
    {
    }

   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the value "v" to all positions of the current array
    *         2D reference.
    *
    *  \param v the new value for all array reference elements
    *  \return the array 2D reference with the new values
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2DReference< T > operator = ( const T& v )
    {
        Array2D< T >::operator = ( v );
        return *this;
    }

   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the values of array "v" to the current array 2D 
    *         reference.
    *
    *  The size of the current array reference is adopted to the size of 
    *  array \em v.
    *
    *  \param v array with the new values for the current array reference
    *  \return the current array 2D reference with the new values
    *  \throw check_exception the type of the exception is "size mismatch"
    *         and indicates that \em v is not a 2-dimensional or
    *         non-dimensional array or that 
    *         the current array is a static array 
    *         reference whose size is different to that of \em v
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2DReference< T > operator = ( const Array< T >& v )
    {
        Array2D< T >::operator = ( v );
        return *this;
    }


   //========================================================================
   /*!
    *  \ingroup Assign
    *  \brief Assigns the values of the 2-dimensional array "v" to the 
    *         current array 2D reference.
    *
    *  The size of the current array 2D reference is adopted to the size of 
    *  the 2-dimensional array \em v.
    *
    *  \param v 2-dimensional array with the new values for the current 
    *           array reference
    *  \return the current array 2D reference with the new values
    *  \throw check_exception the type of the exception is "size mismatch"
    *         and indicates that the current array is a static array 
    *         reference whose size is different to that of \em v
    *
    *  \author  M. Kreutz
    *  \date    2001-04-24
    *
    *  \par Changes
    *      none
    *
    *  \par Status
    *      stable
    *
    */
    Array2DReference< T > operator = ( const Array2D< T >& v )
    {
        Array2D< T >::operator = ( v );
        return *this;
    }
};

//===========================================================================

#endif /* !__ARRAY_H */













