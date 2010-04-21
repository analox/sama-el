//===========================================================================
/*!
 *  \file ArrayIndex.h
 *
 *  \brief Offers special index arrays (arrays with elements of type
 *         "unsigned") for usage with the programming language IDL.
 *
 *  IDL stands for Interactive Data Language and is a scientific
 *  computing environment that combines mathematics, advanced data
 *  visualization, scientific graphics, and a graphical user interface
 *  toolkit to analyze and visualize scientific data. <br>
 *  IDL often uses index vectors, so this class here was generated
 *  to allow the exchange of data between IDL and the array library and
 *  is used as base for some operations, familiar to IDL users. <br>
 *  The index arrays are used to address elements of other arrays 
 *  (the indices in index arrays are always related to the indices
 *  of the element vector Array::e of an array) and
 *  besides their type they also own another speciality, a pointer
 *  ArrayIndex::pos to the last occupied index array element (in the
 *  following named as "occupation pointer"). This is necessary,
 *  because when adding new elements, extra storage for further
 *  add operations is allocated to speed up these operations.  
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
 *
 *  \par Project:
 *      Array
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  \par File and Revision:
 *      $RCSfile: ArrayIndex.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/01 15:41:45 $
 *
 *  \par Changes:
 *      $Log: ArrayIndex.h,v $
 *      Revision 2.1  2004/06/01 15:41:45  saviapbe
 *
 *      The bugs in the doxygen's documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:21  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.10  2002/05/16 13:18:34  rudi
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


#ifndef _ARRAYINDEX_H

//! Set to "1" if this header file is read in. Used to avoid
//! multiple inclusion (and all the problems coming with it)
//! of this header file.
#define _ARRAYINDEX_H 1

#include "Array/Array.h"


//===========================================================================
/*!
 *  \brief This class defines index arrays for usage with the programming
 *         language IDL.
 *
 *  IDL stands for Interactive Data Language and is a scientific
 *  computing environment that combines mathematics, advanced data
 *  visualization, scientific graphics, and a graphical user interface
 *  toolkit to analyze and visualize scientific data. <br>
 *  IDL often uses index vectors, so this class here was generated
 *  to allow the exchange of data between IDL and the array library. <br>
 *  The index arrays (arrays of type unsigned) are used to address 
 *  elements of other arrays 
 *  (the indices in index arrays are always related to the indices
 *  of the element vector Array::e of an array) and
 *  besides their type they also own another speciality, a pointer
 *  ArrayIndex::pos to the last occupied index array element (in the
 *  following named as "occupation pointer"). This is necessary,
 *  because when adding new elements, extra storage for further
 *  add operations is allocated to speed up these operations. <br>
 *  This class is used as base for the functions in ArrayIndexOp.h.
 *
 *  \author  D. Homberg
 *  \date    ????
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class ArrayIndex : protected Array<unsigned>
{
	public:

                //! Pointer to a single element of an index array. 
                //! Used for compatibility with the other stdlib structures.
                //! This is a specialization of Array::iterator.
		typedef unsigned* iterator;

                //! Constant pointer to a single element of an index array. 
                //! Used for compatibility with the other stdlib structures.
                //! This is a specialization of Array::const_iterator.
		typedef unsigned* const_iterator;

	private:

	protected:

                //! Index of the last occupied index array position,
		//! in the following named as "occupation pointer".
		unsigned pos;
	
	public:

               //=============================================================
               /*! 
                *  \ingroup Create
                *  \brief Creates a new empty index array.
                *
                *  The occupation pointer #pos is set to zero.
                *
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
		ArrayIndex() : Array<unsigned>(), pos(0)
		{
			eEnd = e;
		}
		

               //=============================================================
               /*! 
                *  \ingroup Create
                *  \brief Creates a new empty index array with
                *         size "i".
                *
                *  The occupation pointer #pos is set to zero.
                *
		*  \param i of the new index array
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
		ArrayIndex(unsigned i) : Array<unsigned>(i), pos(0)
		{
			eEnd = e;
		}
		

               //=============================================================
               /*! 
                *  \ingroup Create
                *  \brief Creates a new index array with
                *         the content of index array "v".
                *
		*  \param v the index array, whose content will
		*           be used for the new created index array
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
		ArrayIndex(const ArrayIndex& v ) :
			Array<unsigned>(static_cast<const Array<unsigned>&>(v)) 
		{
			pos = v.pos;
			eEnd = e + pos;
		}

               //=============================================================
               /*! 
                *  \brief Destructs the current index array.
                *
                *  The destructor is only virtual!
                *
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
		virtual ~ArrayIndex() {};
    

               //=============================================================
               /*! 
                *  \ingroup Assign
                *  \brief Assigns the content of index array "v" to the
		*         current index array.
                *
                *  \param v the index array, whose content will be assigned
		*           to the calling index array.
                *  \return reference to the current index array with the
		*          new content
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
		ArrayIndex& operator = ( const ArrayIndex& v )
		{
			if (this == &v) return *this;
			Array<unsigned>::operator=(v);
			pos = v.pos;
			eEnd = e + pos;
			return *this;
		}
	
               //=============================================================
               /*! 
                *  \brief Returns an iterator to the first index array element.
                *
                *  A pointer to the first element of the element vector 
                *  Array::e
                *  of the current index array is returned as iterator. <br>
                *  This method is used for compatibility with the other
                *  stdlib structures.
                *
                *  \return Iterator to the first index array element.
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
		inline iterator begin()      
		{
			return e;
		}
		

               //=============================================================
               /*! 
                *  \brief Returns a constant iterator to the first index array 
                *         element.
                *
                *  A pointer to the first element of the element vector 
                *  Array::e
                *  of the current index array is returned as constant 
		*  iterator. <br>
                *  This method is used for compatibility with the other
                *  stdlib structures.
                *
                *  \return Constant iterator to the first index array
		*          element.
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
		inline const const_iterator begin() const
		{
			return e;
		}


               //=============================================================
               /*! 
                *  \brief Returns an iterator to the last index array 
                *         element.
                *
                *  A pointer to the last element of the element vector 
                *  Array::e
                *  of the current index array is returned as iterator. <br>
                *  This method is used for compatibility with the other
                *  stdlib structures.
                *
                *  \return Iterator to the last index array
		*          element.
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
		inline iterator       end()
		{
			return eEnd;
		}
		

               //=============================================================
               /*! 
                *  \brief Returns a constant iterator to the last index array 
                *         element.
                *
                *  A pointer to the last element of the element vector 
                *  Array::e
                *  of the current index array is returned as constant 
		*  iterator. <br>
                *  This method is used for compatibility with the other
                *  stdlib structures.
                *
                *  \return Constant iterator to the last index array
		*          element.
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
		inline const const_iterator end()   const
		{
			return eEnd;
		}

               //=============================================================
               /*! 
		*  \ingroup Information
                *  \brief Returns the current value of the occupation pointer 
		*         #pos.
                *
                *  \return the current value of the occupation pointer
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
		unsigned       getPos()       {return pos;}


               //=============================================================
               /*! 
		*  \ingroup Information
                *  \brief Returns the current value of occupation pointer 
		*         #pos as constant.
                *
                *  \return the current value of the occupation pointer as 
		*          constant
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
		const unsigned getPos() const {return pos;}
		

               //=============================================================
               /*! 
                *  \brief Sets occupation pointer #pos to the new value "i".
                *
                *  If #pos exceeds the number of elements for which
		*  memory is allocated, then the index array is
		*  resized and its content is kept if possible.
                *
                *  \param  i new occupation pointer value, by default set to 
		*            zero
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
		void setPos(const unsigned i = 0)
		{
			pos  = i;
			if (pos >= ne)
			{
				resize(pos, true);
			}
			eEnd = e + pos;
		}
		

               //=============================================================
               /*! 
                *  \brief Resets occupation pointer #pos to zero.
                *
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
		void reset()
		{
			setPos();
		}
		
               //=============================================================
               /*!
		*  \ingroup Resize 
                *  \brief Resizes the index array for method #add by adding
		*         at least memory for one element and also
		*         some extra memory depending on the number
		*         of elements stored in the index array.
                *
                *  Each resize operation consumes calculation time.
		*  If you're e.g. calling the add method many times,
		*  then a normal resize method would be called
		*  many times, too. To prevent this behaviour, this
		*  resize method always allocates some extra memory
		*  that is \f$frac{1}{5}\f$ of the number of elements stored
		*  in the index array. 
                *
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
                *  \sa add
                *
                */		
		void resize()
		{
			Array<unsigned>::resize((nelem() * 1.2 + 1), true);
			eEnd = e + pos;
		}
		

               //=============================================================
               /*! 
		*  \ingroup Resize
                *  \brief Resizes the calling index array to the size of the
		*         index array "ai".
                *
                *  If the occupation pointer value #pos of the calling index 
		*  array is greater than the occupation pointer value of 
		*  \em ai, it will be set to the occupation pointer value of 
		*  \em ai.
                *
                *  \return none
                *
                *  \author  D. Homberg
                *  \date    ????
                *
                *  \par Changes
                *      2002-03-14, ra: The occupation pointer "pos"
		*      was set wrong, when
		*      "ai" had less elements than the calling
		*      array. This error was fixed by adding
		*      the conditional. 
                *
                *  \par Status
                *      stable
                *
                *  \sa add
                *
                */		
		void resize(const ArrayIndex& ai)
		{
			Array<unsigned>::resize(ai.nelem());
                       
                        if ( ai.pos < pos ) pos = ai.pos;
          
			eEnd = e + pos;		
		}
		

               //=============================================================
               /*! 
		*  \ingroup Resize
                *  \brief Resizes the calling index array to the new size "i"
		*         and keeps the old content if "cp" is set to "true".
                *
		*  \param i the new size of the index array
		*  \param cp if set to "true", the old content of the
		*            index array is kept if possible, otherwise
		*            it is erased.
                *  \return none
                *
                *  \author  D. Homberg
                *  \date    ????
                *
                *  \par Changes
                *      2002-03-14, ra: Swapped the instructions from
		*      the then and else part of the conditional,
		*      because otherwise the old content of the
		*      index array is not kept if the copy-flag is
		*      set to "true". Also added the occupation pointer ("pos")
		*      reset in the new else part. 
                *
                *  \par Status
                *      stable
                *
                *  \sa add
                *
                */		
		void resize(unsigned i, bool cp = false)
		{
			unsigned _d[ 1 ];
			_d[ 0 ] = i;
			resize_i( _d, 1, cp );
			if (cp)
			{
                                if (pos >= ne) pos = ne;
				eEnd = e + pos;
			} else {
			        pos = 0;
				eEnd = e;
			}
		}
	
               //=============================================================
               /*! 
		*  \ingroup Add
                *  \brief Adds element "p" to the array.
                *
                *  The new value \em p is added at position #pos.
                *  If #pos exceeds the number of elements for which
                *  memory is already allocated, then the array is resized
		*  in a way, that not only memory for the new element
		*  is allocated, but also some extra memory for future 
		*  add operations, so the time consuming resize method
		*  has not to be called so often, which is important,
		*  because calling the add method is the only
		*  way to add new elements to the index array. <br>
		*  The size of extra memory depends on the number of
		*  elements that are already stored in the array.
		*  The more elements are stored, the greater the extra
		*  memory is.
                *
                *  \param  p new element that is added to the array
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
                *  \sa #resize()
                *
                */		
		void add(const unsigned& p)
		{
			if (pos >= nelem())
			{
			    ArrayIndex::resize();
			}
			e[pos] = p;
			pos++;
			eEnd++;
		}


   //=============================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns the i-th element of the array index.
    *
    *  \param  i index of the element that will be returned
    *  \return the i-th element
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" and indicates that the current index array
    *         is not one-dimensional or that \em i exceeds the
    *         size of the index array's dimension
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
    unsigned& operator ( ) ( unsigned i )
    {
			SIZE_CHECK( nd == 1 )
			RANGE_CHECK( i < d[ 0 ] )
			return e[ i ];
    }

   //=============================================================
   /*! 
    *  \ingroup Extract
    *  \brief Returns the i-th element of the array index as constant.
    *
    *  \param  i ndex of the element that will be returned
    *  \return the i-th element as constant
    *  \throw check_exception the type of the exception will be
    *         "size mismatch" and indicates that the current index array
    *         is not one-dimensional or that \em i exceeds the
    *         size of the index array's dimension
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
    const unsigned& operator ( ) ( unsigned i ) const
    {
			SIZE_CHECK( nd == 1 )
			RANGE_CHECK( i < d[ 0 ] )
			return e[ i ];
    }
};

#endif /* _ARRAYINDEX_H */










